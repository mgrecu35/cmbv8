import glob
retFs=glob.glob("johnkV2/*HDF5")

retFs=sorted(retFs)

import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
for f in retFs:
    with nc.Dataset(f) as fh:
        print(fh['KuGMI'])
        cmb_sfc_precip=fh['KuGMI/nearSurfPrecipTotRate'][:,:]
        orbit=f.split(".")[-2]
        f_npy=glob.glob("npy_tmp_dir/%s/*npz"%orbit)
        f_npy=sorted(f_npy)
        print(f_npy)
        v=[]
        for f in f_npy:
            d=np.load(f)
            v.append(d['near_surface_onnx_precip_rate'])
            #break

    onnx_precip=np.concat(v,axis=0)
    nscans=cmb_sfc_precip.shape[0]
    onnx_precip=onnx_precip[:nscans,:]
    a=np.nonzero(cmb_sfc_precip>-0.01)
    b=np.nonzero(cmb_sfc_precip[a]<30)
    print(np.corrcoef(onnx_precip[a][b],cmb_sfc_precip[a][b]))
    plt.scatter(onnx_precip[a],cmb_sfc_precip[a])
