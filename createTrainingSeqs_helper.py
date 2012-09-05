#!/usr/bin/env python
# -*- coding: ASCII -*-

"""

:Author: Martin Kircher
:Contact: Martin.Kircher@eva.mpg.de
:Date: *23.06.2009

"""

import sys,os
from optparse import OptionParser
from optparse import OptionGroup
import string

import createTrainingSeqs_reader as Reader

max_length_soap=60
max_frags_created = 1

table = string.maketrans('.','N')

# CALCULATE IDENTITY
def seq_identity(seq1,seq2):
  length = min(30,len(seq1),len(seq2))
  seq1 = seq1[:length]
  seq2 = seq2[:length]
  dist = length
  for pos in range(length):
    if (seq1[pos] != seq2[pos]): dist-=1
  if length > 0:
    return dist/float(length)
  else:
    return 1.0

def write_out_soap_fa(ofilestream,name,seq,adapter=None):
  global max_length_soap

  if adapter != None:
    for pos in xrange(max(0,len(seq)-1)):
      if seq_identity(seq[pos:],adapter) >= 0.9:
        if ((len(seq)-pos) < len(adapter)):
          seq = seq[:pos]
          break
        else:
          seq = ""
          break

  frags=len(seq)//((len(seq)//max_length_soap)+1)
  outseqs = []
  while len(seq) > 0:
    if len(seq) <= max_length_soap:
      cseq = seq
      seq = ""
    else:
      cseq = seq[:frags]
      seq = seq[frags:]
    outseqs.append(cseq)

  for count,cseq in enumerate(outseqs):
    ofilestream.write(">%s#%d/%d\n%s\n"%(name,count+1,len(outseqs),cseq))


#########################
### USER PARAMETER ...
#########################

parser = OptionParser(usage="usage: %prog [options]")
parser.add_option("--expID", dest="expID", help="Name of experiment to use for sequence names",default="SoapTraining")

parser.add_option("--infile", dest="infile", help="Input filename")
parser.add_option("--ftype", dest="ftype", help="Filetype")

parser.add_option("--start",dest="start",help="First base to include")
parser.add_option("--end",dest="end",help="Last base to include")
parser.add_option("--srange",dest="srange",help="SRANGE string")

parser.add_option("--start2",dest="start2",help="Second first base to include")
parser.add_option("--end2",dest="end2",help="Second last base to include")
parser.add_option("--srange2",dest="srange2",help="SRANGE2 string")

parser.add_option("--numberN", dest="numberN", help="Maximum number of missing bases to be accepted",type="int")
parser.add_option("--adapter", dest="adapter", help="Adapter sequence",default='')
parser.add_option("--control_index", dest="control_index", help="Sequence of index identifying control reads used for training/test data set (default '' = no filter)",default="")

parser.add_option("--index_seq1", dest="index_seq1", help="Index read is stored attached to read1 input file",default=False,action="store_true")
parser.add_option("--index_seq2", dest="index_seq2", help="Index read is stored in a separate input file",default=False,action="store_true")

parser.add_option("--outfile", dest="outfile", help="Path for output file")
parser.add_option("--outfile2", dest="outfile2", help="Path for second output file")
(options, args) = parser.parse_args()


if (options.infile == None) or (options.ftype == None) or (options.start == None) or (options.end == None) or (options.numberN == None) or (options.adapter == None) or (options.outfile == None) or (options.srange == None):
  print "Not enough parameters."
  sys.exit()

options.adapter = options.adapter.upper()
if (len(options.adapter) == 0) or (options.adapter == 'NONE'): options.adapter = None

start = int(options.start)
start2 = None
end = int(options.end)
end2 = None

second_in_first = False
if (options.start2 != None) and (options.end2 != None) and (options.srange2 != None) and (options.outfile2 != None):
  second_in_first = True
  start2 = int(options.start2)
  end2 = int(options.end2)

outfile = open(options.outfile,'w')
if second_in_first:
  outfile_r2 = open(options.outfile2,'w')
else: # SHOULD NEVER BE NECCESSARY
  outfile_r2 = outfile

file_reader = None
if options.ftype == "qseq":
  file_reader = Reader.read_qseq_file
elif options.ftype == "bcl":
  file_reader = Reader.read_bcl_tile
elif options.ftype == "fastq":
  file_reader = Reader.read_fastq_file
elif options.ftype == "seq_prb":
  file_reader = Reader.read_seq_prb_file
else:
  print "Unexpected file type list."
  sys.exit()

ifile = None
ilength = len(options.control_index)
indexFilter = False
options.control_index = options.control_index.strip().upper()
if (len(options.control_index) > 1):
  indexFilter = True
  if (options.index_seq2):
    ifile = file_reader("_".join(options.infile.split("_")[:-3]+['2']+options.infile.split("_")[-2:]))
  elif (options.index_seq1):
    ifile = file_reader("_".join(options.infile.split("_")[:-3]+['1']+options.infile.split("_")[-2:]))
  #print ifile

cfailed = 0
ckept = 0
for name,seqorg,qualorg in file_reader(options.infile):
  if indexFilter:
    index= None
    if ifile != None:
      iname,iseqorg,iqualorg = ifile.next()
      index = iseqorg[-ilength:].upper()
    else:
      index=seqorg[end:end+ilength].upper()

    if index != options.control_index:
      cfailed += 1
      continue
    else:
      ckept += 1


  seq=seqorg[start:end]
  name = ":".join(name)
  if (seq.count("N") <= options.numberN):
    write_out_soap_fa(outfile,options.expID+":"+name+":"+options.srange,seq,options.adapter)

  if second_in_first:
    seq = seqorg[start2:end2]
    if (seq.count("N") <= options.numberN):
      write_out_soap_fa(outfile_r2,options.expID+":"+name+":"+options.srange2,seq,options.adapter)

outfile.close()
if second_in_first: outfile_r2.close()

if cfailed > 0:
  print "Excluded %d (kept %d) sequences because index did not match."%(cfailed,ckept)
