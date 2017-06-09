"""Echo command line arguments


"""

__author__ = "Jason Wong"
__version__ = "$Revision: 1.0 $"
__date__ = "$Date: 2017/04/27 10:06:20 $"
__copyright__ = "Copyright (c) 2017 Jason WH Wong"
__license__ = "UNSW Faculty of Medicine"

import sys
import os
import getopt
import tokenize
import StringIO
import gzip
import copy

class pdb:
      prot = ""
      pos = {}
      revpos = {}
      start = 0
      end = 0
      chain = ""
      AA = {}

class uniprot:
      ID = ""
      SSBond = []
      pdbs = []

class exBedGraph:

    def __init__(self):
        self.outFile = None
        self.pdbDir = None
        self.mapFile = None
        self.Prots = {}
        self.pdbMap = {}


    # Reads the residue level map file from PDBSWS, such that for every PDB chain the corresponding Uniprot
    # start, end, and positions of cysteine residues are recorded.
    def readInputMap(self):
    
        print "Processing BigMap file..."+os.path.basename(self.bigmapFile)
        f = open(self.bigmapFile,'r')
        c = 0
        for line in f:
	  if len(line) > 5:
            if c%421977 == 0:
               print "."
            line = ' '.join(line.split())
            tmp = line.split(" ")
            tmpID = tmp[0].swapcase()+tmp[1]
	    if len(tmp) > 7:
             try:
                tmpPDB = self.pdbMap[tmpID]
                #if tmp[3] == 'CYS':
                if tmp[6] == 'C':                  # only store cysteine positions to save memory
                  try:
                      tmpPDB.pos[int(tmp[4])]=int(tmp[7])
                  except ValueError:
                      tmpPDB.pos[int(tmp[4][:-1])]=int(tmp[7])
                  try:
                      tmpPDB.revpos[int(tmp[7])]=int(tmp[4])
                  except ValueError:
                      tmpPDB.revpos[int(tmp[7])]=int(tmp[4][:-1])
                  try:
                      tmpPDB.AA[int(tmp[4])]=tmp[3]
                  except ValueError:
                      tmpPDB.AA[int(tmp[4][:-1])]=tmp[3]
                if tmpPDB.end < int(tmp[7]):
                  tmpPDB.end = int(tmp[7])
                self.pdbMap[tmpID] = tmpPDB
             except KeyError:
                try:
                  tmpPDB = pdb()
                  tmpPDB.prot = tmp[5]
                  #if tmp[3] == 'CYS':
                  if tmp[6] == 'C':
                     try:
                         tmpPDB.pos[int(tmp[4])]=int(tmp[7])
                     except ValueError:
                         tmpPDB.pos[int(tmp[4][:-1])]=int(tmp[7])
                     try:
                        tmpPDB.revpos[int(tmp[7])]=int(tmp[4])
                     except ValueError:
                        tmpPDB.revpos[int(tmp[7])]=int(tmp[4][:-1])
                  tmpPDB.start = int(tmp[7])
                  tmpPDB.end = 0
                  try:
                      tmpPDB.AA[int(tmp[4])]=tmp[3]
                  except ValueError:
                      tmpPDB.AA[int(tmp[4][:-1])]=tmp[3]
                  tmpPDB2 = pdb()
                  tmpPDB2.prot = tmpPDB.prot
                  tmpPDB2.pos = copy.deepcopy(tmpPDB.pos)
                  tmpPDB2.revpos = copy.deepcopy(tmpPDB.revpos)
                  tmpPDB2.start = tmpPDB.start
                  tmpPDB2.end = tmpPDB.end
                  tmpPDB2.AA = copy.deepcopy(tmpPDB.AA)
                  tmpPDB.pos.clear()
                  tmpPDB.revpos.clear()
                  self.pdbMap[tmpID] = tmpPDB2
                except:
                  a = None
             except IndexError:
                a = None
            c+=1
        print "PDBs: "+str(len(self.pdbMap))


    # Builds a list of all UniProt IDs in the database and for each of these store all disulfide bonds present and all PDB chains
    # associated with that Uniprot ID
    def processPDB(self):
       # Creates an object holder for each UniProt protein with a list of disulfide bonds and PDB chains
	print "Processing Input file..."+os.path.basename(self.mapFile)
        f = open(self.mapFile,'r')
        for line in f:
            line = line[:-1]
            tmp = line.split(" ")
            tmpProt = uniprot()
            tmpProt.ID = tmp[2]
            tmpProt2 = uniprot()
            tmpProt2.ID = tmpProt.ID
            tmpProt2.SSBond = tmpProt.SSBond[:]
            tmpProt2.pdbs = tmpProt.pdbs[:]
            self.Prots[tmpProt2.ID] = tmpProt2
	f.close()
	print "Uniprots: "+str(len(self.Prots))

        # Start populating the disulfide bond and PDB chains list for each Uniprot protein entry.
        print "Reading PDB..."+os.path.basename(self.pdbDir)
        dirList = os.listdir(self.pdbDir)
        for fname in dirList:
            if fname.find(".ent") != -1:
               print self.pdbDir+"//"+fname
               f = gzip.open(self.pdbDir+"//"+fname,'r')
               header = ""
               expt = ""
               starts = {}
               ends = {}
               hasSS = 0
               chains = []
               for line in f:
                   if line[0:6] == "HEADER":
                      line = line[:-1].strip()
                      header = line[-4:]
                      print "HEAD: "+header
                   elif line[0:6] =="EXPDTA":
                      line = ' '.join(line[:-1].split())
                      expt = line.split(" ")[1]
                      expt = expt.split(",")[0]
                      #print "EXPT: "+expt
                      if expt <> "X-RAY":
                         print "NOT XRAY!!!"
                         break
                   elif line[0:6] =="SSBOND":
                      hasSS = 1
                      line = ' '.join(line.split())
                      bond = line.split(" ")
                      pos1 = bond[4]
                      pos2 = bond[7]
                      #print line
		      try:
                      	tmpID1 = header+bond[3]
                      	tmpID2 = header+bond[6]
                      except:
                        #print "invalid SSBOND line"
			continue
		      tmpPDB1 = None
                      tmpPDB2 = None
                      try:
                          tmpPDB1 = self.pdbMap[tmpID1]
                      except:
                          #print "Uniprot C1 not found"
                          continue
                      try:
                          tmpPDB2 = self.pdbMap[tmpID2]
                      except KeyError, k:
                          #print "Uniprot C2 not found: "+str(k)
                          continue
                      SS = tmpPDB1.start
                      SE = tmpPDB1.end
                      try:
                          AAC1 = tmpPDB1.pos[int(bond[4])]
                      except KeyError:
                          #print "no C1"
                          continue
                      except ValueError:
                          try:
                              AAC1 = tmpPDB1.pos[int(bond[4][:-1])]
                          except KeyError:
                              #print "no C1"
                              continue
                      try:
                          AAC2 = tmpPDB2.pos[int(bond[7])]
                      except KeyError:
                          #print "no C2"
                          continue
                      except ValueError:
                          try:
                              AAC2 = tmpPDB2.pos[int(bond[7][:-1])]
                          except KeyError:
                              #print "no C2"
                              continue
                      if AAC2 < AAC1:
                         tmp = AAC1
                         AAC1 = AAC2
                         AAC2 = tmp
                      if tmpPDB1.prot == tmpPDB2.prot and bond[3] == bond[6]:     #Uniprot and Chains are the same
                         try:
                             curProt = self.Prots[tmpPDB1.prot]
                             if expt.find("X-RAY") != -1:
                                found = 0
                                for i in curProt.pdbs:
                                    if i == tmpID1:
                                       found = 1
                                if found == 0:
                                   curProt.pdbs.append(tmpID1)
                                found = 0
                                pos = 0
                                for i in curProt.SSBond:
                                    if i[0] == int(AAC1) and i[1] == int(AAC2):
                                      tmplist = i[2][:]
                                      tmplist.append(tmpID1)
                                      tmpposlist = i[3][:]
                                      tmpposlist.append(str(self.pdbMap[tmpID1].revpos[AAC1])+"-"+str(self.pdbMap[tmpID1].revpos[AAC2]))
                                      curProt.SSBond[pos] = ((i[0],i[1],tmplist[:],tmpposlist[:]))
                                      found = 1
                                    pos+=1
                                if found == 0:
                                    IDlist = []
                                    poslist = []
                                    IDlist.append(tmpID1)
                                    poslist.append(str(self.pdbMap[tmpID1].revpos[AAC1])+"-"+str(self.pdbMap[tmpID1].revpos[AAC2]))
                                    curProt.SSBond.append((int(AAC1),int(AAC2),IDlist[:],poslist[:]))
                             self.Prots[tmpPDB1.prot] = curProt
                         except KeyError:
                                print "Not found: "+tmpPDB1.prot
                   elif line[0:5] == "DBREF" and line.find("UNP") <> -1:
                        tmpc = line.split(" ")
                        #print line
                        chains.append(tmpc[3])
                   elif line[0:4] =="ATOM":
                        if hasSS == 0:
                           for k in chains:
                             try:
                               curProt = self.Prots[self.pdbMap[header+k].prot]
                               curProt.pdbs.append(header+k)
                             except KeyError:
                               print header+k
                        break
               f.close()


    # Outputs results and also work out whether each disulfide bond listed is missing in any other PDB chain and for what reason
    def outputRes(self):
        print "Output results file..."+os.path.basename(self.outFile)
        out = open(self.outFile,'w')
        for key in self.Prots.keys():
            curProt = self.Prots[key]
            uniqXray = []               # Bonds that are only in one structure because not covered by any other
            commonXray = []             # Bonds that appear more than once
            missXray = []               # Bonds that is in one structure but not another yet is covered
            for i in curProt.SSBond:
                for k in curProt.pdbs:
                    found = 0
                    for j in i[2]:
                        if k == j:
                           found =1
                    if found == 0:
                          curPDB = self.pdbMap[k]
                          if curPDB.start <= i[0] and curPDB.end >= i[1]:               # if PDB covers bond
                            found1 = 0
                            found2 = 0
                            for key,val in curPDB.pos.items():
                                if i[0] == val:
                                   found1 = 1
                                if i[1] == val:
                                   found2 = 1
                            out.write(curProt.ID+"\t"+str(i[0])+"\t"+str(i[1])+"\t"+str(i[2])+"\t"+str(i[3])+"\t")
                            try:
                                out.write(str(curPDB.revpos[i[0]])+"\t")
                            except KeyError:
                                out.write("X\t")
                            try:
                                out.write(str(curPDB.revpos[i[1]]))
                            except KeyError:
                                out.write("X")
                            if found1 ==1 and found2 ==1 and curPDB.AA[curPDB.revpos[i[0]]] == "CYS" and curPDB.AA[curPDB.revpos[i[1]]] == "CYS":
                                out.write("\tSSMissing"+"\t"+k)
                                found3 = 0
                                for x in i[2]:
                                    if x.find(k[:-1]) <> -1:
                                       found3 = 1
                                if found3 == 1:
                                   out.write("\t1\n")
                                else:
                                   out.write("\t0\n")
                            else:
                                out.write("\tSSMutant"+"\t"+k+"\t-\n")
                          else:
                            out.write(curProt.ID+"\t"+str(i[0])+"\t"+str(i[1])+"\t"+str(i[2])+"\t"+str(i[3])+"\t")
                            try:
                                out.write(str(curPDB.revpos[i[0]])+"\t")
                            except KeyError:
                                out.write("X\t")
                            try:
                                out.write(str(curPDB.revpos[i[1]]))
                            except KeyError:
                                out.write("X")
                            out.write("\tNonOverlap"+"\t"+k+"\t"+str(curPDB.start)+"-"+str(curPDB.end)+"\n")
        out.close()


    def ParseCommandLine(self, Arguments):
        try:
            (Options, Args) = getopt.getopt(Arguments, "hi:m:o:b:")
        except:
            print "Unknown option entered"
            print UsageInfo
            sys.exit(1)
        OptionsSeen = {}
        for (Option, Value) in Options:
            OptionsSeen[Option] = 1
            if Option == "-i":
                self.pdbDir = Value
            elif Option == "-m":
                self.mapFile = Value
            elif Option == "-b":
                self.bigmapFile = Value
            elif Option == "-o":
                self.outFile = Value
            else:
                print "** Unknown option:", Option, Value
        if not self.pdbDir:
            print "\n* Error: -i missing"
            print UsageInfo
            sys.exit(1)

UsageInfo = """
pdbDisulfide_interXray_gz.py - Checks for disulfide bonds that are in the PDB chain of a XRay structure but not in others.

Parameters:
 -i [FILENAME] input PDB directory (contains .ent.gz files directly from PDB)
 -m [FILENAME] Chain level map file from PDBSWS (http://www.bioinf.org.uk/pdbsws/pdbsws_chain.txt.gz)
 -b [FILENAME] Residue level map file from PDBSWS (http://www.bioinf.org.uk/pdbsws/pdbsws_res.txt.gz)
 -o [Output] Output file

Example:
  pdbDisulfide_interXray_gz -i ~/Z/DB/pdb/  -m ~/Z/DB/pdb/pdb_uniprot_chain_map.lst.2 -b ~/Z/DB/pdb_uniprot_map.lst.2 -o ~/Z/Data/disulfide/interXRay_all_out.txt
"""

def Main(exBed= None):
    global MAX_RESULTS_FILES_TO_PARSE
    if not exBed:
        exBed = exBedGraph()
        exBed.ParseCommandLine(sys.argv[1:])
        exBed.readInputMap()
        exBed.processPDB()
        exBed.outputRes()

if __name__ == "__main__":
    try:
        import psyco
        psyco.full()
    except:
        print ""
    #TestMain()
    Main()
