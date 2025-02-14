#!/usr/bin/python
#makeparams
#Makes 2BCMB parameter files
#Set gmisource, dprsource, iteG, iteD

import os
import sys
import subprocess
import shutil
import string

def runCmd(str):
  import subprocess
  proc = subprocess.Popen(str, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
  (st, err) = proc.communicate()
  return st.strip()

def prevDay(date):
  #return previous day, yyyymmdd
  date = str(date)
  year1s = date[0:4]
  year = int(year1s)
  if(year%4 == 0): febdays = 29
  else: febdays = 28
  dlist = [31, febdays, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
  month = int(date[4:6])
  day = int(date[6:8])
  if(month < 1 or month > 12):
    print("Input date error - month = ", month)
    return -1
  if(day < 1 or day > dlist[month-1]):
    print("Input date error - day = ", day)
    return -1
  day -= 1
  if(day == 0):
    month -= 1
    if(month == 0):
      month = 12
      year -= 1
    day = dlist[month-1]
  return str(year * 10000 + month * 100 + day)

def nextDay(date):
  #return next day, yyyymmdd
  date = str(date)
  year1s = date[0:4]
  year = int(year1s)
  if(year%4 == 0): febdays = 29
  else: febdays = 28
  dlist = [31, febdays, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
  month = int(date[4:6])
  day = int(date[6:8])
  if(month < 1 or month > 12):
    print("Input date error - month = ", month)
    return -1
  if(day < 1 or day > dlist[month-1]):
    print("Input date error - day = ", day)
    return -1
  day += 1
  if(day > dlist[month-1]):
    month += 1
    if(month == 13):
      month = 1
      year += 1
    day = 1
  return  str(year * 10000 + month * 100 + day)

def filelist(date, dateType, path, subdir, filestring1, filestring2):
  year = date[0:4]
  month = date[4:6]
  day = date[6:8]
  if len(subdir) > 0: subdir = "/" + subdir
  if(dateType == 0): path2 = path + subdir
  elif(dateType == 1): path2 = path + "/" + date + subdir
  else: path2 = path + "/" + year + "/" + month + "/" + day + subdir
  path3 = path2 + "/" + filestring1 + "*" + filestring2
  print(path3)
  lsstr = runCmd("ls " + path3)
  #print lsstr
  return lsstr.split()


#Start makeparams

gmisource  = 1  # 0 - Production, 1 - ITE
dprsource  = 1  # 0 - Production, 1 - ITE, 2 - Special
iteG = "ITE753"
iteD = "ITE705"

nargs = len(sys.argv)
if(nargs < 3):
  print("")
  print("   Usage - makeparams dir date1 [date2]   Date format yyyymmdd")
  print("      or   makeparams dir date orbit")
  print("")
  print("   Make Level 2 GPM parameter files.")
  print("   Directory dir is created (in current directory if dir has no path).")
  print("   Parameter files are saved in dir.")
  print("   Set gmisource, dprsource, iteG, iteD in script.")  
  print("")
  sys.exit(0)

dir1 = sys.argv[1]
if dir1[0] != "/":
   cpath = os.getcwd()
   dir1 = cpath + "/" + dir1

date1 = sys.argv[2]
date2 = date1
orbit = 0
if(nargs > 3):
  if len(sys.argv[3]) == 8: date2 = sys.argv[3]
  else: orbit = int(sys.argv[3])

os.system("rm -rf " + dir1)
os.makedirs(dir1)
os.chdir(dir1)
dprsource=0
if dprsource == 0:
  #DPR-Production
  pathD = "/gpmdata"
  dateTypeD = 2  # 0 - no date, 1 - yyyymmdd, 2 - yyyy/mm/dd
  subdirD = "radar"
  string1DKu = "2A.GPM.Ku."
  string2DKu = "V07*.HDF5"
  string1DDPR = "2A.GPM.DPR*.V"
  string2DDPR = "V07*.HDF5"
  string1DKuE = "2A-ENV.GPM.Ku."
  string2DKuE = "V07*.HDF5"
if dprsource == 1:
  #DPR-ITE
  pathD = "/itedata/" + iteD
  #pathD = "/PANFS/ite/data/product/regular"
  dateTypeD = 2  # 0 - no date, 1 - yyyymmdd, 2 - yyyy/mm/dd
  subdirD = "radar"
  string1DKu = "2A.GPM.Ku"
  string2DKu = "." + iteD + ".HDF5"
  string1DDPR = "2A.GPM.DPR*.V9"
  string2DDPR = "." + iteD + ".HDF5"
  string1DKuE = "2A-ENV.GPM.Ku"
  string2DKuE = "." + iteD + ".HDF5"
gmisource=0
if gmisource == 0:
  #GMI-Production
  pathG = "/gpmdata"
  dateTypeG = 2  # 0 - no date, 1 - yyyymmdd, 2 - yyyy/mm/dd
  subdirG = "1C"
  string1G = "1C.GPM.GMI.XCAL2016"
  string2G = ".V07*.HDF5"
else:
  #GMI-ITE
  pathG = "/itedata/" + iteG
  dateTypeG = 2  # 0 - no date, 1 - yyyymmdd, 2 - yyyy/mm/dd
  subdirG = "1C"
  string1G = "1C.GPM.GMI."
  string2G = "." + iteG + ".HDF5"

while(1):
  if(date1 > date2): break
  print(date1)
  date1m = prevDay(date1)
  date1p = nextDay(date1)
  flKu = filelist(date1, dateTypeD, pathD, subdirD, string1DKu, string2DKu)
  if len(flKu) == 0:
    print("No Ku files")
    date1 = nextDay(date1)
    continue
  flgmi1 = filelist(date1, dateTypeG, pathG, subdirG, string1G, string2G)
  flgmi0 = filelist(date1m, dateTypeG, pathG, subdirG, string1G, string2G)
  flgmi2 = filelist(date1p, dateTypeG, pathG, subdirG, string1G, string2G)
  flDPR = filelist(date1, dateTypeD, pathD, subdirD, string1DDPR, string2DDPR)
  flKuE = filelist(date1, dateTypeD, pathD, subdirD, string1DKuE, string2DKuE)
  print((len(flDPR)))
  print((len(flKu)))
  print(flKu)
  print(flDPR)
  lastorbit = len(flKu) - 1

  i = 0
  for fKu in flKu:
    #print(fKu)
    fKu=fKu.decode("utf-8")
    #pos = fKu.find("-E") + 9
    #orbit1 = fKu[pos:pos+6]
    orbit1=fKu.split('.')[-3]
    if orbit > 0 and int(orbit1) != orbit:
      i += 1
      continue
    fgmi1 = flgmi1[i]
    if i == 0: fgmi0 = flgmi0[-1]
    else: fgmi0 = flgmi1[i-1]
    if i == lastorbit: fgmi2 = flgmi2[0]
    else: fgmi2 = flgmi1[i + 1]
    print((orbit,fKu,i))
    #print(flDPR[i])
    flDPRu = fKu.replace("2A.GPM.Ku","2A.GPM.DPR")
    flKuEu = fKu.replace("2A.GPM.Ku","2A-ENV.GPM.Ku")
    fDPR = flDPRu#[i]
    fKuE = flKuEu#[i]
    i += 1

    paramfile = date1 + "." + orbit1 + ".param0"
    fp = open(paramfile, "w")
    fp.write("ialg=1\n")
    fgmi0=fgmi0.decode('utf-8')
    fgmi1=fgmi1.decode('utf-8')
    fgmi2=fgmi2.decode('utf-8')
    fp.write("f1CGMI=" + fgmi0 + "\n")
    fp.write("f1CGMI=" + fgmi1 + "\n")
    fp.write("f1CGMI=" + fgmi2 + "\n")
    fp.write("f2AKu=" + fKu + "\n")
    fp.write("f2ADPR=" + fDPR + "\n")
    fp.write("f2KuENV=" + fKuE + "\n")
    fp.write("fCMB=cmb_output/2B.GPM.DPRGMI." + date1 +"." + orbit1 + ".HDF5" + "\n")
    fp.write("rseed=157983119\n")
    fp.write("rseed=835697248\n")
    fp.write("ifs=1\n")
    fp.write("iiad=0\n")
    fp.close()

  date1 = nextDay(date1)

print("")
exit(0)
#End

