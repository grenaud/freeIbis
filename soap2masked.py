#!/usr/bin/env python

import sys,os
from optparse import OptionParser
from optparse import OptionGroup
import subprocess
import time
import re
import createTrainingSeqs_reader as Reader



if (len(sys.argv[0]) > 0)  and os.path.isdir(os.getcwd()+"/"+os.path.dirname(sys.argv[0])):
  root_bin = os.path.dirname(os.getcwd()+"/"+sys.argv[0])
elif (len(sys.argv[0]) > 0)  and os.path.isdir(os.path.dirname(sys.argv[0])):
  root_bin = os.path.dirname(sys.argv[0])
else:
  root_bin = os.getcwd()
f=open(root_bin+'/params.py')
c=compile(f.read(), root_bin+'/params.py', "exec")
eval(c)

#########################
### USER PARAMETER ...
#########################

percentageFraction=0.0;
parser = OptionParser(usage="usage: %prog [options]")
group = OptionGroup(parser, "General","General options")
group.add_option("-f", "--file", dest="filesoap", help="Soap output")
group.add_option("-r", "--ref", dest="reference", help="PhiX reference (fasta format)")
group.add_option("-p", "--percent", dest="percentage", help="Percentage of divergent bases to output base (default 10)",default=10,type="int")
group.add_option("-q", "--quality", dest="quality", help="Minimum quality for a mismatch to be considered (default 0)",default=0,type="int")
group.add_option("-o", "--outprefix", dest="outprefix", help="Prefix for output files")

parser.add_option_group(group)


(options, args) = parser.parse_args()

if (options.filesoap == None) or (not os.path.isfile(options.filesoap)):
  print "Need valid SOAP output\n"
  sys.exit()

if (options.reference == None) or (not os.path.isfile(options.reference)):
  print "Need valid reference sequence\n"
  sys.exit()


if (options.outprefix == None) :
  print "Need valid prefix for output files\n"
  sys.exit()


if (options.percentage >= 0) and (options.percentage <= 100):
  percentageFraction=float(options.percentage)/100.0;
else:
  print "Enter a percentage between 0 and 100\n"
  sys.exit()



## BEGIN reading phiX ref

phixfp = open(options.reference)

line = phixfp.readline();

if(not(line.startswith(">"))):
  print "Invalid reference fasta file"
  sys.exit()

#phiXrefgenome="";
#
#while 1:
#    line = phixfp.readline();
#    if(not(line)):
#        break;
#    if(line.startswith(">")):
#       print "Reference cannot be a multi-fasta file"
#       sys.exit()
#    phiXrefgenome= phiXrefgenome+line.rstrip();

phiXrefgenome = {};
phiXcoverage  = {};
phiXmmA       = {};
phiXmmC       = {};
phiXmmG       = {};
phiXmmT       = {};

for title,seq in Reader.read_fasta(options.reference):
  phiXrefgenome[title.split()[0]]  =seq;
  phiXcoverage[title.split()[0]]   =[0]*len(seq);
  phiXmmA[title.split()[0]]        =[0]*len(seq);
  phiXmmC[title.split()[0]]        =[0]*len(seq);
  phiXmmG[title.split()[0]]        =[0]*len(seq);
  phiXmmT[title.split()[0]]        =[0]*len(seq);


## END reading phiX ref




soapout = open(options.filesoap)

while 1:
    line = soapout.readline();
    if(not(line)):
        break;
    columns=line.split("\t");
    sequenceRead=columns[1];
    indexPhiX   =int(columns[8])-1;
    alignLength =int(columns[5]);
    mm          =int(columns[9]);
    controlName =columns[7];

    for i in range(indexPhiX,indexPhiX+alignLength):
      phiXcoverage[controlName][i]+=1;
    if(mm > 0):
      for mmFound in columns[10:]:
        #                    1      2     3        4
        mat=re.search("^([ACGT])->(\d+)([ACGT])-?(\d+)$", mmFound);
        if( mat ):
          if(int(mat.group(4))>=int(options.quality)):
            indexMMRead=int( mat.group(2) );

            if(phiXrefgenome[controlName][indexPhiX+indexMMRead].lower() != mat.group(1).lower()):
              print "Base pair on phiX does not match the one from the SOAP output in line "+line;
              sys.exit(1)

            if(sequenceRead[indexMMRead].lower() != mat.group(3).lower()):
              print "Base pair on read does not match the one from the SOAP output for line "+line;
              sys.exit(1)

            if(mat.group(3).lower()   == "a"):
              phiXmmA[controlName][indexPhiX+indexMMRead]+=1;
            elif(mat.group(3).lower() == "c"):
              phiXmmC[controlName][indexPhiX+indexMMRead]+=1;
            elif(mat.group(3).lower() == "g"):
              phiXmmG[controlName][indexPhiX+indexMMRead]+=1;
            elif(mat.group(3).lower() == "t"):
              phiXmmT[controlName][indexPhiX+indexMMRead]+=1;
            else:
              print "Wrong format for SOAP output, invalid bp for line = "+line;
              sys.exit(1)
        else:
           print "Wrong format for SOAP output for line = "+line;
           sys.exit(1)


    

setMaskedPositions = set([]);
try:
  fileHcoverage = open ( options.outprefix + ".covr", 'w' ) ;
except:
  print "Cannot write to file "+options.outprefix + ".covr";
  sys.exit(1)

try:
  fileHmasked   = open ( options.outprefix + ".mask", 'w' ) ;
except:
  print "Cannot write to file "+options.outprefix + ".mask";
  sys.exit(1)

for controlkeys in phiXrefgenome.keys():
  for i in range(len(phiXrefgenome[controlkeys])):
    fileHcoverage.write(controlkeys+"\t"+str(i)+"\t"+str(phiXcoverage[controlkeys][i])+"\t"+str(phiXmmA[controlkeys][i])+"\t"+str(phiXmmC[controlkeys][i])+"\t"+str(phiXmmG[controlkeys][i])+"\t"+str(phiXmmT[controlkeys][i])+"\n");
    if(phiXcoverage[controlkeys][i] > 0):
      if(float(phiXmmA[controlkeys][i])/float(phiXcoverage[controlkeys][i]) >= percentageFraction):
        setMaskedPositions.add(i);
      if(float(phiXmmC[controlkeys][i])/float(phiXcoverage[controlkeys][i]) >= percentageFraction):
        setMaskedPositions.add(i);
      if(float(phiXmmG[controlkeys][i])/float(phiXcoverage[controlkeys][i]) >= percentageFraction):
        setMaskedPositions.add(i);
      if(float(phiXmmT[controlkeys][i])/float(phiXcoverage[controlkeys][i]) >= percentageFraction):
        setMaskedPositions.add(i);

  setMaskedPositions=sorted(setMaskedPositions);

  for maskPos in setMaskedPositions:
    fileHmasked.write(controlkeys+"\t"+str(maskPos)+"\n");

fileHcoverage.close();
fileHmasked.close();
