# -
#
# SPDX-FileCopyrightText: Pedro Costa and the CaNS contributors
# SPDX-License-Identifier: MIT
#
# -
#!/usr/bin/env python
def read_single_field_binary(data_dir,filenamei,iskip):
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
    geofile  = data_dir+"/geometry.out"
    geo = np.loadtxt(geofile, comments = "!", max_rows = 2)
    ng = geo[0,:].astype('int')
    #
    # read grid
    #
    grids = []
    for axis, n, offset in zip(["x","y","z"],ng,r0):
        grid = np.fromfile(os.path.join(data_dir,"grid_"+axis+".bin"),dtype=precision)
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
    # read binary file
    #
    iskip       = np.asarray(iskip,dtype=int)
    n           = ((ng[:]-1)//iskip[:] + 1).astype(int)
    fld         = np.fromfile(data_dir+"/"+filenamei,dtype=precision)
    if(fld.size != np.prod(n)):
        raise ValueError("expected {} values for ng={} and iskip={}, found {}".format(np.prod(n),ng,iskip,fld.size))
    data        = np.reshape(fld,(n[0],n[1],n[2]),order='F')
    #
    # reshape grid
    #
    xp = xp[0:ng[0]:iskip[0]]
    yp = yp[0:ng[1]:iskip[1]]
    zp = zp[0:ng[2]:iskip[2]]
    xu = xu[0:ng[0]:iskip[0]]
    yv = yv[0:ng[1]:iskip[1]]
    zw = zw[0:ng[2]:iskip[2]]
    return data, xp, yp, zp, xu, yv, zw

if __name__ == "__main__":
    import numpy as np
    filenamei   = input("Name of the binary file written by CaNS (e.g. vex_fld_0000000.bin)]: ")
    iskipx      = input("Data saved every (ix, iy, iz) points. Value of ix? [1]: ") or "1"
    iskipy      = input("Data saved every (ix, iy, iz) points. Value of iy? [1]: ") or "1"
    iskipz      = input("Data saved every (ix, iy, iz) points. Value of iz? [1]: ") or "1"
    iskip       = np.array([iskipx,iskipy,iskipz]).astype(int)
    data, xp, yp, zp, xu, yv, zw = read_single_field_binary("./",filenamei,iskip)
