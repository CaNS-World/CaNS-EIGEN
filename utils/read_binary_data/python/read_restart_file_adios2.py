#!/usr/bin/env python
#
# SPDX-FileCopyrightText: Pedro Costa and the CaNS contributors
# SPDX-License-Identifier: MIT
#
def read_restart_file_adios2(filenamei):
    import os
    os.environ.setdefault("ADIOS2_ALWAYS_USE_MPI","1")
    import adios2
    import numpy as np
    #
    # read checkpoint ADIOS2 file
    #
    FileReader = getattr(adios2,"FileReader",None)
    def open_bp(filename):
        return FileReader(filename) if FileReader else adios2.open(filename,"r")
    def read_bp5(filename,names):
        adios = adios2.ADIOS()
        io = adios.DeclareIO("reader")
        io.SetEngine("BP5")
        engine = io.Open(filename,adios2.Mode.Read)
        engine.BeginStep()
        available_variables = io.AvailableVariables()
        data = {}
        for name in names:
            var = io.InquireVariable(name)
            if(var is None):
                raise KeyError(name)
            shape = tuple(var.Shape())
            dtype = np.float32 if var.Type() == "float" else np.float64
            if("int" in var.Type()): dtype = np.int32
            arr = np.empty(shape,dtype=dtype)
            engine.Get(var,arr,adios2.Mode.Sync)
            data[name] = arr
        engine.EndStep()
        engine.Close()
        return available_variables,data
    def split_prefix(filename):
        root, ext = os.path.splitext(filename)
        if(ext == ".bp"):
            for suffix in ["_u","_v","_w","_p"]:
                if(root.endswith(suffix)):
                    return root[:-len(suffix)]
            return root
        return filename
    def read_one(filename,fieldname):
        if(FileReader):
            with open_bp(filename) as fh:
                data = np.transpose(np.asarray(fh.read(fieldname)))
                time = float(np.asarray(fh.read("time"))[0])
                istep = int(np.asarray(fh.read("istep"))[0])
        else:
            _, values = read_bp5(filename,[fieldname,"time","istep"])
            data = np.transpose(values[fieldname])
            time = float(values["time"][0])
            istep = int(values["istep"][0])
        return data, time, istep
    is_combined = False
    if(os.path.exists(filenamei)):
        if(FileReader):
            with open_bp(filenamei) as fh:
                available_variables = fh.available_variables()
        else:
            available_variables, _ = read_bp5(filenamei,[])
        is_combined = all(fieldname in available_variables for fieldname in ["u","v","w","p"])
    if(is_combined):
        if(FileReader):
            with open_bp(filenamei) as fh:
                u = np.transpose(np.asarray(fh.read("u")))
                v = np.transpose(np.asarray(fh.read("v")))
                w = np.transpose(np.asarray(fh.read("w")))
                p = np.transpose(np.asarray(fh.read("p")))
                time = float(np.asarray(fh.read("time"))[0])
                istep = int(np.asarray(fh.read("istep"))[0])
        else:
            _, values = read_bp5(filenamei,["u","v","w","p","time","istep"])
            u = np.transpose(values["u"])
            v = np.transpose(values["v"])
            w = np.transpose(values["w"])
            p = np.transpose(values["p"])
            time = float(values["time"][0])
            istep = int(values["istep"][0])
    else:
        prefix = split_prefix(filenamei)
        u, time, istep = read_one(prefix+"_u.bp","u")
        v, _   , _     = read_one(prefix+"_v.bp","v")
        w, _   , _     = read_one(prefix+"_w.bp","w")
        p, _   , _     = read_one(prefix+"_p.bp","p")
    return u, v, w, p, time, istep
if __name__ == "__main__":
    filenamei = input("Name of the ADIOS2 restart file written by CaNS [fld.bp]: ") or "fld.bp"
    u, v, w, p, time, istep = read_restart_file_adios2(filenamei)
