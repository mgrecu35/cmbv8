import glob
pFiles=glob.glob("paramFiles/*")
pFiles=sorted(pFiles)[16:18]
import os
for pf in pFiles[:]:
    cmd="./combAlg2.exe junk %s "%pf
    orbit=pf.split(".")[-2]
    print(cmd,orbit)
    #continue
    os.system(cmd)
    
    try:
        cfiles=glob.glob("npy_tmp_dir/*.npz")
        #orbit=cfiles[0].split(".")[-2]
        print(orbit)
        os.system("mkdir npy_tmp_dir/%s"%orbit)
        os.system("mv npy_tmp_dir/*%s*npz npy_tmp_dir/%s"%(orbit,orbit))
    except:
        continue
    #break
