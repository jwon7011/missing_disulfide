import os
import sys

#Usage: annotatePDB.py [input interXray file] [input DSB file] [outfile]
#
#Add DSB annotation to interXray file output from pdbDisulfide_interXray_gz.py
#Uses output directly from pdbDisulfide_interXray_gz.py and the appropriate annotations from DSB
#Output from DSB must be saved in tab-delimited format.

if len(sys.argv) < 4:
   print "Usage: annotatePDB.py [input interXray file] [input DSB file] [outfile]"
   exit()


DSBanno = {}     #{PDB&chain&Cys1-Cys2:count}

DSBfile = open(sys.argv[2],"r")
for line in DSBfile:
    tmp = line.split("\t")
    DSBanno[tmp[0]+tmp[5]+tmp[6]+"-"+tmp[10]]=line
DSBfile.close()
Xfile = open(sys.argv[1],"r")
out = open(sys.argv[3],"w")
for line in Xfile:
    tmp = line.split("\t")
    PDBs = tmp[3].replace("\'","").replace("[","").replace("]","").replace(" ","").split(",")
    POSs = tmp[4].replace("\'","").replace("[","").replace("]","").replace(" ","").split(",")
    found = 0
    pos = 0
    while found == 0 and pos < len(PDBs):
          try:
              out.write(line[:-1]+"\t"+DSBanno[PDBs[pos]+POSs[pos]])
              found = 1
          except KeyError:
              pos+=1
    if found == 0:
       out.write(line[:-1]+"\n")
Xfile.close()
out.close()

