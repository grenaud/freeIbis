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
import gzip
import string
import time
import math
import subprocess

if (len(sys.argv[0]) > 0)  and os.path.isdir(os.getcwd()+"/"+os.path.dirname(sys.argv[0])):
  root_bin = os.path.dirname(os.getcwd()+"/"+sys.argv[0])
elif (len(sys.argv[0]) > 0)  and os.path.isdir(os.path.dirname(sys.argv[0])):
  root_bin = os.path.dirname(sys.argv[0])
else:
  root_bin = os.getcwd()
f=open(root_bin+'/params.py')
c=compile(f.read(), root_bin+'/params.py', "exec")
eval(c)

offset=64 #Solexa quality offset
table = string.maketrans('.','N')
compSEQ = string.maketrans('CATG','GTAC')

def entropy(seq):
  counts = []
  for x in 'ACGT':
    counts.append(seq.count(x))
  total = float(sum(counts))
  entropy = 0
  for elem in counts:
    if (total > 0) and (elem/total <> 0):
      entropy -= elem/total*math.log(elem/total,2)
  return entropy

def rev_compl(seq):
  nseq=seq.translate(compSEQ)
  nseq=list(nseq)
  nseq.reverse()
  return "".join(nseq)

def parse_rangestr(rangestr):
  res = []
  fields = rangestr.split(',')
  for elem in fields:
    if "-" in elem:
      se = elem.split('-')
      if len(se) == 2:
        try:
          start = int(se[0])
          end = int(se[1])
          res.extend(range(start,end+1))
        except: return None
      else: return None
    else:
      try:
        pos = int(elem)
        res.append(pos)
      except: return None
  if len(res) == 0: return None
  else:
    res = list(set(res))
    res.sort()
    return res

def read_fastq_file(filename):
  infile = open(filename)
  name,seq,qual = "","",""
  count = 0
  for line in infile:
    if line.startswith("@") and ((count == 0) or (count == 3)):
      count=1
      if len(seq) > 0:
        yield name,seq,qual
      name = tuple(line[1:].strip().split(":")[-4:])
      if len(name) != 4:
        print "Unexpected sequence ID:",line.strip()
        raise StopIteration
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
      print "Unexpected line:",line.strip()
      raise StopIteration
  if len(seq) > 0:
    yield name,seq,qual
  infile.close()
  raise StopIteration

def read_qseq_file(filename):
  global table
  infile = open(filename)
  flane = 2
  ftile = 3
  fx = 4
  fy = 5
  fseq = 8
  fqual = 9
  for line in infile:
    fields = line.split('\t')
    if len(fields) > fqual:
      yield (fields[flane],fields[ftile],fields[fx],fields[fy]),fields[fseq].translate(table),fields[fqual]
    else:
      print "Unexpected line",line.rstrip()
      raise StopIteration
  infile.close()
  raise StopIteration

def read_seq_prb_file(filename):
  global offset
  prb = open(filename[1],'r')
  base = open(filename[0],'r')
  cprb = " "
  cbase = " "
  while ((cprb <> []) and (cbase <> [])):
    cbase = base.readline().split()
    if (len(cbase) > 4):
      name=tuple(cbase[:4])
      cbase=cbase[4].translate(table)
      cprb = prb.readline().split("\t")
      prstr=''
      for x in cprb:
        prstr+=chr(max(map(int,x.split()))+offset)
      yield name,cbase,prstr
  prb.close()
  base.close()

#########################
### USER PARAMETER ...
#########################

parser = OptionParser(usage="usage: %prog [options]")
group = OptionGroup(parser, "Bustard","General options describing sequencing data")
group.add_option("-p", "--path", dest="path", help="Path to bustard folder (default .)", default=".")
group.add_option("-l", "--lanes", dest="lanes", help="Lanes, example: 1,3,4-8 (default 4)",default="4")
group.add_option("-t", "--tiles", dest="tiles", help="Lanes, example: 1-13,100 (default 1-120)",default="1-120")
parser.add_option_group(group)

group = OptionGroup(parser, "Filter/Paired End","Excluding sequences and bases for training")
group.add_option("--start",dest="start",help="First base to include/first base in each read (comma separate)")
group.add_option("--end",dest="end",help="Last base to include/last base in each read (comma separate)")
group.add_option("--adapter", dest="adapter", help="NON-FUNCTIONAL",default='')
group.add_option("-c","--lowestQuality", dest="qcutoff", help="Lowest accepted quality score (default 15)",default=15,type="int")
group.add_option("-e","--lowestEntropy", dest="ecutoff", help="Lowest accepted entropy value (default 0.85)",default=0.85,type="float")
group.add_option("--numberN", dest="numberN", help="Maximum number of missing bases to be accepted (default 3)",default=3,type="int")
group.add_option("--indexlength",dest="indexlength",help="Length of the index read (default 0/AUTO)",default=0,type="int")
group.add_option("--2nd_indexlength",dest="indexlength2", help="NON-FUNCTIONAL",type="int",default=0)
group.add_option("--control_index", dest="control_index", help="NON-FUNCTIONAL",default='')
parser.add_option_group(group)

group = OptionGroup(parser, "Files","Options for final/intermediate output")
group.add_option("-o", "--outfile", dest="outfile", help="Path for output file(s)")
group.add_option("--temp", dest="tmp", help="Path to temporary folder (default "+def_temp+")",default=def_temp)
parser.add_option_group(group)

(options, args) = parser.parse_args()

################################
# EVALUATE READ STARTS AND ENDS
################################

reads = 1
tstart = 1
tstart2 = None
tend = None
tend2 = None
if (options.start != None) and ((options.start).count(",") == 1):
  try:
    fields = (options.start).split(",")
    tstart = int(fields[0])
    tstart2 = int(fields[1])
    reads = 2
  except:
    print "Can't evaluate sequence starts:",options.start
    sys.exit()
elif (options.start != None) and ((options.start).count(",") == 0):
  try:
    tstart = int(options.start)
  except:
    print "Can't evaluate sequence start:",options.start
    sys.exit()
elif options.start != None:
  print "Can't evaluate sequence start:",options.start
  sys.exit()

if (options.end != None) and ((options.end).count(",") == 1):
  try:
    fields = (options.end).split(",")
    tend = int(fields[0])
    tend2 = int(fields[1])
    reads = 2
  except:
    print "Can't evaluate sequence ends:",options.end
    sys.exit()
elif (options.end != None) and ((options.end).count(",") == 0):
  try:
    tend = int(options.end)
  except:
    print "Can't evaluate sequence end:",options.end
    sys.exit()
elif options.end != None:
  print "Can't evaluate sequence end:",options.end
  sys.exit()

if tstart > 0: tstart -= 1
else: tstart=0
if tend <= 0: tend = None

if tstart2 > 0: tstart2 -= 1
else: tstart2=0
if tend2 <= 0: tend2 = None

if (reads == 2) and ((tend2 == None) or (tstart2 == None) or (tend == None)):
  print "Need both starts and ends for two reads."
  sys.exit()

print "Using bases",tstart+1,"to",tend
if (reads == 2):
  print "Using bases",tstart2+1,"to",tend2," of second read."

if (options.tmp == None) or not(os.path.isdir(options.tmp)):
  print "Need path to temporary folder."
  sys.exit()

if (options.path == None) or not(os.path.isdir(options.path)):
  print "Need path to Bustard folder."
  sys.exit()
options.path=options.path.rstrip("/")

if (options.outfile == None):
  print "Need name for training output file."
  sys.exit()

lanes = parse_rangestr(options.lanes)
if lanes == None:
  print "Need a valid range of lanes."
  sys.exit()
print "Using lanes:",lanes

tiles = parse_rangestr(options.tiles)
if tiles == None:
  print "Need a valid range of tiles."
  sys.exit()
print "Using tiles:",tiles

# SHOULD NEVER BE NECCESSARY, BUT YOU NEVER KNOW :)
if tstart2 == None: tstart2=tstart
if tend2 == None: tend2=tend

srange1=str(tstart)+":"+str(tend)
srange2=str(tstart2)+":"+str(tend2)

outfile = open(options.outfile,'w')
if reads == 2:
  outfile_r2 = open(options.outfile+"_r2",'w')

count_reads = 0
count_reads2 = 0
count_reads_filter = 0
count_reads2_filter = 0

print ""
for lane in lanes:
  lane = str(lane)
  files = []
  first_read = []
  second_read = []
  index_read2 = False

  for nr in tiles:
    nr = str(nr)
    while len(nr) < 4:
      nr="0"+nr
    if os.path.isfile(options.path+"/s_"+lane+"_1_"+nr+"_qseq.txt"):
      first_read.append(options.path+"/s_"+lane+"_1_"+nr+"_qseq.txt")
    if os.path.isfile(options.path+"/s_"+lane+"_3_"+nr+"_qseq.txt"):
      second_read.append(options.path+"/s_"+lane+"_3_"+nr+"_qseq.txt")
      index_read2 = True
    elif os.path.isfile(options.path+"/s_"+lane+"_2_"+nr+"_qseq.txt"):
      second_read.append(options.path+"/s_"+lane+"_2_"+nr+"_qseq.txt")

    if os.path.isfile(options.path+"/s_"+lane+"_"+nr+"_seq.txt") and os.path.isfile(options.path+"/s_"+lane+"_"+nr+"_prb.txt"):
      first_read.append((options.path+"/s_"+lane+"_"+nr+"_seq.txt",None)) #options.path+"/s_"+lane+"_"+nr+"_prb.txt"
    if os.path.isfile(options.path+"/s_"+lane+"_"+nr+"_seq.txt") and (not os.path.isfile(options.path+"/s_"+lane+"_"+nr+"_prb.txt")):
      first_read.append((options.path+"/s_"+lane+"_"+nr+"_seq.txt",None))

    if os.path.isfile(options.path+"/s_"+lane+"_"+nr+".fastq"):
      first_read.append(options.path+"/s_"+lane+"_"+nr+".fastq")
    if os.path.isfile(options.path+"/s_"+lane+"_1_"+nr+".fastq"):
      first_read.append(options.path+"/s_"+lane+"_1_"+nr+".fastq")
    if os.path.isfile(options.path+"/s_"+lane+"_3_"+nr+".fastq"):
      second_read.append(options.path+"/s_"+lane+"_3_"+nr+".fastq")
      index_read2 = True
    elif os.path.isfile(options.path+"/s_"+lane+"_2_"+nr+".fastq"):
      second_read.append(options.path+"/s_"+lane+"_2_"+nr+".fastq")

  start2ndread=None
  second_in_first = ((len(second_read) == 0) and (reads == 2))
  file_reader = None
  for filename in first_read:
    had_training_sequences = True
    if len(filename) != 2 and filename.endswith("qseq.txt"):
      file_reader = read_qseq_file
    elif len(filename) != 2 and filename.endswith("fastq.txt"):
      file_reader = read_fastq_file
    elif len(filename) == 2:
      file_reader = read_seq_prb_file
    else:
      print "Unexpected content of file list."
      sys.exit()

    print "Reading",filename
    for name,seqorg,qualorg in file_reader(filename):
      name = ":".join(name)
      if index_read2:
         for iname,iseqorg,iqualorg in file_reader(filename.replace("s_"+lane+"_1_","s_"+lane+"_2_")):
           start2ndread=len(seqorg)+len(iseqorg)
           break
      elif (start2ndread == None) and (tend == None):
        tend = len(seqorg)
        srange1 = str(tstart)+":"+str(tend)
      elif (start2ndread == None) and (tend != None):
        start2ndread = tend+options.indexlength
      if tend == None:
        tend = len(seqorg)
        srange1 = str(tstart)+":"+str(tend)

      if start2ndread == None: start2ndread=max(len(seqorg),tstart2)
      seq=list(seqorg[tstart:tend])
      qual=qualorg[tstart:tend]
      qual_int = map(lambda x:ord(x)-offset,list(qual))
      for pos,qscore in enumerate(qual_int):
        if qscore < options.qcutoff:
          seq[pos]="N"
      seq="".join(seq).upper()
      if (seq.count("N") <= options.numberN) and (entropy(seq) >= options.ecutoff):
        outfile.write("%s\t%d\t%d\t%s\n"%("\t".join(name),tstart,tend,seq))
        count_reads+=1
      else:
        count_reads_filter+=1

      if second_in_first:
        seq = list(seqorg[tstart2:tend2])
        qual=qualorg[tstart2:tend2]
        qual_int = map(lambda x:ord(x)-offset,list(qual))
        for pos,qscore in enumerate(qual_int):
          if qscore < options.qcutoff:
            seq[pos]="N"
        seq="".join(seq).upper()
        if (seq.count("N") <= options.numberN) and (entropy(seq) >= options.ecutoff):
          outfile_r2.write("%s\t%d\t%d\t%s\n"%("\t".join(name),tstart2,tend2,seq))
          count_reads2+=1
        else:
          count_reads2_filter+=1

  if start2ndread == None: start2ndread = tend
  for filename in second_read:
    if len(filename) != 2 and filename.endswith("qseq.txt"):
      file_reader = read_qseq_file
    elif len(filename) != 2 and filename.endswith("fastq.txt"):
      file_reader = read_fastq_file
    elif len(filename) == 2:
      file_reader = read_seq_prb_file
    else:
      print "Unexpected content of file list."
      sys.exit()

    print "Reading",filename
    for name,seqorg,qualorg in file_reader(filename):
      seq=list(seqorg[max(tstart2-start2ndread,0):max(tend2-start2ndread,0)])
      qual=qualorg[max(tstart2-start2ndread,0):max(tend2-start2ndread,0)]
      qual_int = map(lambda x:ord(x)-offset,list(qual))
      for pos,qscore in enumerate(qual_int):
        if qscore < options.qcutoff:
          seq[pos]="N"
      seq="".join(seq).upper()
      if (seq.count("N") <= options.numberN) and (entropy(seq) >= options.ecutoff):
        outfile_r2.write("%s\t%d\t%d\t%s\n"%("\t".join(name),tstart2,tend2,seq))
        count_reads2+=1
      else:
        count_reads2_filter+=1

print "Saved",count_reads,"sequences to",options.outfile,"(Reads not passing filter %d)"%count_reads_filter
outfile.close()
if reads == 2:
  outfile_r2.close()
  print "Saved",count_reads2,"sequences to",options.outfile+"_r2","(Reads not passing filter %d)"%count_reads2_filter
