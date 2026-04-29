#!/usr/bin/env python
#
# SPDX-FileCopyrightText: Pedro Costa and the CaNS contributors
# SPDX-License-Identifier: MIT
#
def read_single_field_hdf5(data_dir,filenamei,varname=""):
    import os
    import h5py
    import numpy as np
    #
    # setting up some parameters
    #
    r0 = np.array([0.,0.,0.]) # domain origin
    #
    # read field file
    #
    filepath = os.path.join(data_dir,filenamei)
    hf = h5py.File(filepath,'r')
    if(len(varname) == 0):
        fields = list(hf["fields"].keys())
        if(len(fields) != 1):
            raise ValueError("varname must be provided when the HDF5 file stores multiple fields")
        varname = fields[0]
    data = np.transpose(np.asarray(hf["fields/"+varname]))
    if("meta/lo" in hf):
        lo = np.asarray(hf["meta/lo"],dtype=int)
        hi = np.asarray(hf["meta/hi"],dtype=int)
        nskip = np.asarray(hf["meta/nskip"],dtype=int)
    else:
        ng_local = np.array(data.shape,dtype=int)
        lo = np.array([1,1,1],dtype=int)
        hi = ng_local.copy()
        nskip = np.array([1,1,1],dtype=int)
    hf.close()
    #
    # read grid
    #
    grids = []
    for axis, offset in zip(["x","y","z"],r0):
        with h5py.File(os.path.join(data_dir,"grid_"+axis+".h5"),"r") as hf:
            grids.append((offset + np.asarray(hf["rc"]), offset + np.asarray(hf["rf"])))
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
    # reshape grid
    #
    xp = xp[lo[0]-1:hi[0]:nskip[0]]
    yp = yp[lo[1]-1:hi[1]:nskip[1]]
    zp = zp[lo[2]-1:hi[2]:nskip[2]]
    xu = xu[lo[0]-1:hi[0]:nskip[0]]
    yv = yv[lo[1]-1:hi[1]:nskip[1]]
    zw = zw[lo[2]-1:hi[2]:nskip[2]]
    return data, xp, yp, zp, xu, yv, zw
if __name__ == "__main__":
    filenamei = input("Name of the HDF5 file written by CaNS (e.g. vex_fld_0000000.h5)]: ")
    varname = input("Name of the HDF5 field dataset []: ") or ""
    data, xp, yp, zp, xu, yv, zw = read_single_field_hdf5("./",filenamei,varname)
