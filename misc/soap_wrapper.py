#!/usr/bin/env python
# -*- coding: ASCII -*-

"""

:Author: Martin Kircher
:Contact: Martin.Kircher@eva.mpg.de
:Date: *22.06.2009

"""

import sys,os
from optparse import OptionParser
from optparse import OptionGroup
import gzip
import string
import time
import subprocess

max_length_soap=80
max_frags_created = 1
def_soap_path = '/mnt/solexa/bin/soap_1.11/soap'
def_temp = '/tmp/'

compSEQ = string.maketrans('CATG','GTAC')

def repair_mm(clist,offset):
  res = []
  for elem in clist:
    position = ''
    fields = elem.split("->")
    for char in fields[1]:
      if char.isdigit():
        position += char
      else:
        break
    res.append(fields[0]+"->"+str(int(position)+offset)+fields[1].lstrip(position))
  return res

def rev_compl(seq):
  global compSEQ
  nseq=seq.translate(compSEQ)
  nseq=list(nseq)
  nseq.reverse()
  return "".join(nseq)

def read_fasta(filename):
  if filename.endswith('.gz'):
    GZinfile = gzip.GzipFile(filename)
    infile = GZinfile.read().splitlines()
  else:
    infile = open(filename)
  seq=""
  header=""
  for line in infile:
    line = line.strip()
    if (line[0] == ">"):
      if (header != "") and (seq != ""):
        yield header,seq
      header = line[1:]
      seq = ""
    else:
      seq+=line
  if (header != ""):
    yield header,seq
  if filename.endswith('.gz'):
    GZinfile.close()
  else:
    infile.close()
  raise StopIteration

def read_fastq_file(filename):
  if filename.endswith('.gz'):
    GZinfile = gzip.GzipFile(filename)
    infile = GZinfile.read().splitlines()
  else:
    infile = open(filename)
  name,seq,qual = "","",""
  count = 0
  for line in infile:
    if line.startswith("@") and ((count == 0) or (count == 3)):
      count=1
      if len(seq) > 0:
        yield name,seq,qual
      name = line[1:].strip()
      seq = ""
      qual = ""
    elif (not line.startswith("+")) and (count == 1):
      seq += line.strip()
    elif line.startswith("+") and (count == 1):
      count=2
    elif (count >= 2):
      qual += line.strip()
      count = 3
    else:
      if line.startswith('>'):
        for header,seq in read_fasta(filename):
          yield header,seq,None
        raise StopIteration
      else:
        print "Unexpected line:",line.strip()
        raise StopIteration
  if len(seq) > 0:
    yield name,seq,qual
  if filename.endswith('.gz'):
    GZinfile.close()
  else:
    infile.close()
  raise StopIteration

def write_out_soap_fa(ofilestream,name,oseq):
  global max_length_soap
  global max_frags_created
  frags=len(oseq)//((len(oseq)//max_length_soap)+1)
  count = 1
  seq = oseq
  while len(seq) > 0:
    if len(seq) <= max_length_soap:
      cseq = seq
      seq = ""
    else:
      cseq = seq[:frags]
      seq = seq[frags:]
    ofilestream.write(">%s#%d-%d\n%s\n"%(name,count,len(oseq),cseq))
    if max_frags_created < count: max_frags_created=count
    count += 1

def eval_soap_output(infilename,outfilename):
  outfile = open(outfilename,'w')
  infile = open(infilename,'r')
  count = 0
  line = infile.readline()
  fields = line.split()
  cid = ""
  cnum = -1
  cdata = {}
  delta_bases = -1
  next_delta_bases = -1
  while line != "":
    nextline = infile.readline()
    nextfields = nextline.split()
    try:
      hfields = nextfields[0].split("#")
      nextcid = "#".join(hfields[:-1])
      hfields = hfields[-1].split("-")
      nextnum = int(hfields[0])
      next_delta_bases = int(hfields[1])
    except:
      nextcid = ""
      nextnum = -1
      next_delta_bases = -1

    if (nextcid != cid) and (cid != ""):
      cdata[cnum] = fields
      keys = cdata.keys()
      keys.sort()
      if len(cdata[keys[0]]) > 9:
        soap_line = cdata[keys[0]]
        soap_line[0] = cid
        fstart = int(cdata[keys[0]][8])-1
        lpos = int(cdata[keys[0]][8])
        lref = cdata[keys[0]][7]
        lstrand = cdata[keys[0]][6]
        if lstrand == "-":
          seq = rev_compl(cdata[keys[0]][1])
          soap_line=soap_line[:10] + repair_mm(soap_line[10:],delta_bases-len(seq))
        else:
          seq = cdata[keys[0]][1]
          lpos += len(cdata[keys[0]][1])
        consistent = True
        for key in keys[1:]:
          if (len(cdata[key]) <= 9) or (lref != cdata[key][7]) or (lstrand != cdata[key][6]):
            consistent = False
            break
          else:
            if (lstrand == '-') and (lpos-len(cdata[key][1]) == int(cdata[key][8])):
              lpos -= len(cdata[key][1])
              soap_line[8] = cdata[key][8]
              soap_line=soap_line[:10] + repair_mm(cdata[key][10:],delta_bases-len(seq)-len(cdata[key][1])) + soap_line[10:]
              seq += rev_compl(cdata[key][1])
            elif (lstrand == '+') and (lpos == int(cdata[key][8])):
              lpos += len(cdata[key][1])
              soap_line.extend(repair_mm(cdata[key][10:],len(seq)))
              seq += cdata[key][1]
            else:
              consistent = False
              #print lpos+1,lref,lstrand,'  vs  ',cdata[key][7],cdata[key][6],cdata[key][8],'  ',key
              break

        if consistent and (len(seq) == delta_bases):
          if lstrand == "-": seq=rev_compl(seq)
          soap_line[1]=seq
          soap_line[2]="h"*len(seq)
          soap_line[5]=str(len(seq))
          soap_line[9]=str(len(soap_line)-10)
          #if len(soap_line) <= 10:
          outfile.write("\t".join(soap_line)+"\n")
          #else:
            #outfile.write("\t".join(soap_line[:10])+'\t'+' '.join(soap_line[10:])+"\n")
          count += 1

      cdata = {}
      cid = nextcid
      delta_bases = next_delta_bases
    elif (cid == ""):
      hfields = nextfields[0].split("#")
      cid = "#".join(hfields[:-1])
      hfields = hfields[-1].split("-")
      cnum = int(hfields[0])
      delta_bases = int(hfields[1])
      cdata[cnum] = fields
    else:
      cdata[cnum] = fields
    cnum = nextnum
    line = nextline
    fields = nextfields
  infile.close()
  outfile.close()
  return count

#########################
### USER PARAMETER ...
#########################

parser = OptionParser(usage="usage: %prog [options]")
group = OptionGroup(parser, "SOAP","Parameters for SOAP mapping")
group.add_option("-s", "--soap", dest="soap", help="Path to SOAP binary (default "+def_soap_path+")",default=def_soap_path)
group.add_option("-p", "--params", dest="params", help="Additional SOAP params (do not use -g!,default \'-r 0\')",default='-r 0')
group.add_option("-c", "--cores", dest="cores", help="Number of CPU cores for SOAP (default 1)",default=1,type="int")
group.add_option("-r", "--reference", dest="reference", help="Reference genome file")
group.add_option("-m", "--mismatch", dest="mismatch", help="Number of mismatches allowed (default 5)",default=5,type="int")
parser.add_option_group(group)

group = OptionGroup(parser, "Files","Options for final/intermediate output")
group.add_option("-i", "--infile", dest="infile", help="Path for input file")
group.add_option("-o", "--outfile", dest="outfile", help="Path for output file")
group.add_option("--temp", dest="tmp", help="Path to temporary folder (default "+def_temp+")",default=def_temp)
group.add_option("--keep", dest="keep", help="Keep temporary files for debugging (default OFF)",action="store_true",default=False)
parser.add_option_group(group)

(options, args) = parser.parse_args()

################################
# EVALUATE READ STARTS AND ENDS
################################

if (options.soap == None) or not(os.path.isfile(options.soap)):
  print "Need path to SOAP binary."
  sys.exit()

if (options.reference == None) or not(os.path.isfile(options.reference)):
  print "Need path to reference Genome."
  sys.exit()

if (options.tmp == None) or not(os.path.isdir(options.tmp)):
  print "Need path to temporary folder."
  sys.exit()

if (options.outfile == None):
  print "Need name for output file."
  sys.exit()

if (options.infile == None):
  print "Need name for input file."
  sys.exit()
elif not os.path.isfile(options.infile):
  print "Input file not valid."
  sys.exit()

timestamp = str(time.time())
tstart = 0
tend = None

outfilename = options.tmp+"/"+timestamp+"_soap_input.txt"
print "Creating temporary file with sequences:",outfilename
outfile = open(outfilename,'w')
print "Reading",options.infile
for name,seqorg,qualorg in read_fastq_file(options.infile):
  if tend == None: tend=len(seqorg)
  write_out_soap_fa(outfile,name,seqorg)
outfile.close()

if tend == None:
  sys.exit()
seed_length = max(6,min(12,(min((tend-tstart)//(((tend-tstart)//max_length_soap)+1),max_length_soap)-3)//2))
#seed_length = 8

print ""
print "Calling SOAP..."
outfilename2 = options.tmp+"/"+timestamp+"_soap_output.txt"
cmdline = options.soap+" -s "+str(seed_length)+" -d "+options.reference+" -a "+outfilename+" -v "+str(options.mismatch)+" -g 0 -c 0 -p "+str(options.cores)+" -o "+outfilename2+" "+options.params
print cmdline
proc = subprocess.Popen(cmdline,shell=True)
proc.wait()

print ""
if (max_frags_created > 1) and (options.cores > 1):
   print "Sorting SOAP output by sequence ID"
   proc = subprocess.Popen("sort "+outfilename2+" > "+outfilename2+"_sort",shell=True)
   proc.wait()
   proc = subprocess.Popen("mv "+outfilename2+"_sort "+outfilename2,shell=True)
   proc.wait()
print "Iterating over SOAP output file"
count = eval_soap_output(outfilename2,options.outfile)
if not options.keep:
  print "Removing temporary SOAP output file",outfilename2
  proc = subprocess.Popen("rm -f "+outfilename2,shell=True)
  print "Removing temporary SOAP input file",outfilename
  proc = subprocess.Popen("rm -f "+outfilename,shell=True)
print ""
print count," sequences available in",options.outfile
