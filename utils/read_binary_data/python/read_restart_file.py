# -
#
# SPDX-FileCopyrightText: Pedro Costa and the CaNS contributors
# SPDX-License-Identifier: MIT
#
# -
#!/usr/bin/env python
def read_restart_file(filenamei):
    import os
    import numpy as np
    #
    # setting up some parameters
    #
    iprecision = 8            # precision of the real-valued data
    r0 = np.array([0.,0.,0.]) # domain origin
    precision  = 'float64'
    if(iprecision == 4): precision = 'float32'
    #
    # read geometry file
    #
    geofile  = "geometry.out"
    data = np.loadtxt(geofile, comments = "!", max_rows = 2)
    ng = data[0,:].astype('int')
    #
    # read grid
    #
    grids = []
    for axis, n, offset in zip(["x","y","z"],ng,r0):
        grid = np.fromfile("grid_"+axis+".bin",dtype=precision)
        grid = np.reshape(grid,(n,4),order='F')
        grids.append((offset + grid[:,2], offset + grid[:,3]))
    #
    # split centered and staggered grids
    #
    xp = grids[0][0] # centered grid
    xu = grids[0][1] # staggered grid
    yp = grids[1][0] # centered grid
    yv = grids[1][1] # staggered grid
    zp = grids[2][0] # centered grid
    zw = grids[2][1] # staggered grid
    #
    # read checkpoint binary file
    #
    disp = np.prod(ng)
    fldinfo = np.zeros([2])
    legacy_size = iprecision*(4*disp+2)
    def split_prefix(filename):
        root, ext = os.path.splitext(filename)
        if(ext == ".bin"):
            for suffix in ["_u","_v","_w","_p"]:
                if(root.endswith(suffix)):
                    return root[:-len(suffix)]
            return root
        return filename
    def read_one(filename):
        with open(filename,'rb') as f:
            fld = np.fromfile(f,dtype=precision,count=disp)
            info = np.fromfile(f,dtype=precision,count=2)
        return np.reshape(fld,(ng[0],ng[1],ng[2]),order='F'), info
    if(os.path.exists(filenamei) and os.path.getsize(filenamei) == legacy_size):
        #
        # u, v, w, p
        #
        data = np.zeros([ng[0],ng[1],ng[2],4])
        offset = 0
        with open(filenamei,'rb') as f:
            for q in range(4):
                f.seek(offset)
                fld = np.fromfile(f,dtype=precision,count=disp)
                data[:,:,:,q] = np.reshape(fld,(ng[0],ng[1],ng[2]),order='F')
                offset += iprecision*disp
            f.seek(offset)
            fldinfo[:] = np.fromfile(f,dtype=precision,count=2)
        u = data[:,:,:,0]
        v = data[:,:,:,1]
        w = data[:,:,:,2]
        p = data[:,:,:,3]
    else:
        prefix = split_prefix(filenamei)
        u, fldinfo = read_one(prefix+"_u.bin")
        v, _       = read_one(prefix+"_v.bin")
        w, _       = read_one(prefix+"_w.bin")
        p, _       = read_one(prefix+"_p.bin")
    time  =     fldinfo[0]
    istep = int(fldinfo[1])
    return u, v, w, p, time, istep

if __name__ == "__main__":
    filenamei = input("Name of the binary restart file written by CaNS [fld.bin]: ") or "fld.bin"
    u, v, w, p, time, istep = read_restart_file(filenamei)
