#!/usr/bin/env python3
"""Run the built-in Poisson correctness check on valid boundary/grid layouts."""
import argparse
import itertools
import os
from pathlib import Path
import re
import signal
import subprocess

# Fixture, grid, process grid, pencil axis, distributed tridiagonal solve.
LAYOUTS = [
    ('even', (8, 10, 12), (1, 1), 1, False),
    ('odd', (9, 11, 13), (1, 2), 1, False),
    ('odd', (9, 11, 13), (2, 2), 2, False),
    ('odd', (9, 11, 13), (2, 3), 1, True),
    ('odd', (9, 11, 13), (3, 2), 2, True),
    ('odd', (9, 11, 13), (1, 2), 3, False),
    ('stretched', (9, 11, 9), (1, 1), 1, False),
    ('stretched', (9, 11, 9), (1, 2), 1, True),
    ('thin', (2, 5, 7), (2, 1), 2, False),
    ('thin', (5, 2, 7), (2, 1), 1, False),
    ('thin', (5, 7, 2), (1, 2), 1, False),
]
BCS = [(pair,)*3 for pair in ('PP', 'DD', 'NN', 'DN', 'ND')]
BCS += [('DD', 'NN', 'DN'), ('ND', 'DN', 'PP')]
INVERSE = dict(PP='PP', DD='NN', NN='DD', DN='ND', ND='DN')
BACKENDS = {'fft-fft': (True, True), 'gemm-fft': (False, True),
            'fft-gemm': (True, False), 'gemm-gemm': (False, False)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--backend', choices=['all', *BACKENDS], default='all')
    args = parser.parse_args()
    testdir = Path(__file__).resolve().parent
    rundir = testdir.parent.parent/'run'
    exe = rundir/'cans'
    output = rundir/'poisson'
    # Allow oversubscription with both Open MPI 4 (ORTE) and 5 (PRRTE).
    os.environ.setdefault('OMPI_MCA_rmaps_base_oversubscribe', '1')
    os.environ.setdefault('PRTE_MCA_rmaps_default_mapping_policy', ':oversubscribe')
    os.environ.setdefault('OMP_NUM_THREADS', '2')
    os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
    os.environ.setdefault('MKL_NUM_THREADS', '1')
    cases = list(itertools.product(LAYOUTS, BCS))
    cases += [(('thin', (9, 11, 1), (2, 1), 1, False), bc)
              for bc in (('PP', 'PP', 'PP'), ('ND', 'DN', 'PP'))]
    # EIGEN additionally supports Z/YZ diffusion and selectable FFT/GEMM axes.
    cases = [(layout, bc, 0, None, False) for layout, bc in cases]
    cases += [(('even', (8, 10, 12), (1, 1), 2, False), bc, mode, None, False)
              for mode, bc in itertools.product((1, 2, 3), (BCS[0], BCS[-1]))]
    cases += [(('odd', (9, 11, 13), (2, 3), 1, True), BCS[2], mode, None, False)
              for mode in (1, 3)]
    backends = BACKENDS if args.backend == 'all' else {args.backend: BACKENDS[args.backend]}
    count = 0
    for backend, fft in backends.items():
        selected = list(cases)
        if not all(fft):
            stretch = (0. if fft[0] else 1.5, 0. if fft[1] else 1.5, 1.5)
            selected += [(('odd', (9, 11, 13), (1, 1), 2, False), BCS[-1], mode, stretch, natural)
                         for mode, natural in ((2, False), (3, False), (3, True))]
            selected += [(('thin', grid, (1, 1), 2, False), BCS[0], 3, None, False)
                         for grid in ((1, 11, 13), (9, 1, 13))]
        for index, (layout, bc, mode, stretch, natural) in enumerate(selected, 1):
            fixture, grid, dims, axis, dtdma = layout
            name = f'{index:02d}-{fixture}-{"-".join(bc)}-{dims[0]}x{dims[1]}-p{axis}-d{int(dtdma)}'
            name += f'-m{mode}'
            case = output/backend/name
            (case/'data').mkdir(parents=True, exist_ok=True)
            text = (testdir/f'input-{fixture}.nml').read_text()
            updates = {
                'ng(1:3)': ', '.join(map(str, grid)),
                'dims(1:2)': f'{dims[0]}, {dims[1]}, ipencil_axis = {axis}',
                'is_poisson_dtdma': 'T' if dtdma else 'F',
                'impdiff_mode': str(mode),
                'is_poisson_fft': ', '.join('T' if value else 'F' for value in fft),
                'cbcscal(0:1,1:3,:)': ', '.join(repr(c) for pair in bc for c in pair),
                'cbcpre(0:1,1:3)': ', '.join(repr(c) for pair in bc for c in pair),
            }
            if stretch:
                updates['gtype(1:3)'] = '2, 3, 4, gr(1:3) = '+', '.join(map(str, stretch))
                updates['is_gridpoint_natural_channel(1:3)'] = ', '.join('T' if natural and value else 'F' for value in stretch)
            for component in range(1, 4):
                updates[f'cbcvel(0:1,1:3,{component})'] = ', '.join(repr(c) for pair in bc for c in INVERSE[pair])
            for key, value in updates.items():
                text = re.sub(r'^'+re.escape(key)+r'\s*=.*$', key+' = '+value, text, flags=re.M)
            (case/'input.nml').write_text(text)
            with (case/'run.log').open('w') as log:
                process = subprocess.Popen(['mpirun', '-n', str(dims[0]*dims[1]), str(exe)], cwd=case,
                                           stdout=log, stderr=subprocess.STDOUT, start_new_session=True)
                try:
                    status = process.wait(timeout=60)
                except subprocess.TimeoutExpired:
                    os.killpg(process.pid, signal.SIGKILL)
                    process.wait()
                    raise RuntimeError(f'MPI timeout: {case}/run.log') from None
            log = (case/'run.log').read_text()
            if status != 0 or '*** Fim ***' not in log or 'ERROR:' in log:
                raise RuntimeError(f'Failed {case}\n{log[-3000:]}')
            count += 1
            print(f'PASS {backend}/{name}', flush=True)
    print(f'PASS: {count} Poisson cases; logs in {output}')


if __name__ == '__main__':
    main()
