# -*- coding: ASCII -*-

"""

:Author: Martin Kircher
:Contact: Martin.Kircher@eva.mpg.de
:Date: *27.01.2010

"""

## Documentation for createTrainingSeqs_reader.py
#
# This class contains various subroutines 
# to read different file formats:
#
# read_fasta
# read_fasq_file
# read_qseq_file (old Illumina basecalls)
# read_seq_prb_file (old Illumina format)
# read_bcl : generator to read a single file in bcl format
# get_bcl_cycles : call the read_bcl for each cycle
# read_bcl_tile : calls get_bcl_cycles for a given lane/tile
#
##

import sys,os
import string

offset=64 #Solexa quality offset
table = string.maketrans('.','N')
lendian = 1 ## use -1 for big endian systems
bases = 'ACGT'

def read_fasta(filename):
  infile = open(filename)
  seq=""
  header=""
  for line in infile:
    if (line[0] == ">"):
      line = line.strip()
      if (header != "") and (seq != ""):
        yield header,seq
      header = line[1:]
      seq = ""
    else:
      seq+=line.strip()
  if (header != ""):
    yield header,seq
  infile.close()
  raise StopIteration

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
  filename = filename.split(",")
  if filename[1] != "None": prb = open(filename[1],'r')
  base = open(filename[0],'r')
  cprb = " "
  cbase = " "
  while ((cprb <> []) and (cbase <> [])):
    cbase = base.readline().split()
    if (len(cbase) > 4):
      name=tuple(cbase[:4])
      cbase=cbase[4].translate(table)
      if filename[1] != "None": cprb = prb.readline().split("\t")
      prstr=''
      if filename[1] != "None":
        for x in cprb:
          prstr+=chr(max(map(int,x.split()))+offset)
      yield name,cbase,prstr
  if filename[1] != "None": prb.close()
  base.close()

def to_int(s):
  global lendian
  power = 1
  total = 0
  for elem in s[::lendian]:
    total += ord(elem)*power
    power *= 256
  return total

def read_bcl(filename,clusteridx=None):
  global bases
  infile = open(filename,'rb')
  try:
    nrclusters = to_int(infile.read(4))
    data = infile.read()
    if len(data) == nrclusters:
      iterclusters = xrange(nrclusters)
      if clusteridx != None:
        iterclusters = clusteridx
      for cluster in iterclusters:
        base = (ord(data[cluster]) & 3)
        quality = ord(data[cluster]) >> 2
        yield bases[base],quality
    else:
      print "File content does not reflect number of clusters.",nrclusters,len(data)
  except IOError:
    print "Error reading",filename
  infile.close()
  raise StopIteration

def get_bcl_cycles(root_path,lane,tile,start,end,clusteridx=None,qualities=True):
  #print root_path,lane,tile,start,end
  lanes = map(lambda x: "L%03d"%(x+1),range(8))
  res_bases = None
  res_quals = None
  if lane > 0 and lane <= len(lanes):
    cycles = map(lambda x: "C%d.1"%x,range(start,end+1))
    #print cycles
    if os.path.isdir(root_path+"/"+lanes[lane-1]):
      res_bases = []
      if qualities: res_quals = []
      for ind,cycle in enumerate(cycles):
        if os.path.isfile(root_path+"/"+lanes[lane-1]+"/"+cycle+"/s_%d_%d.bcl"%(lane,tile)):
          count = 0
          for cluster in read_bcl(root_path+"/"+lanes[lane-1]+"/"+cycle+"/s_%d_%d.bcl"%(lane,tile),clusteridx):
            if ind == 0:
              res_bases.append(cluster[0])
              if qualities: res_quals.append(chr(cluster[1]+64))
            else:
              res_bases[count]+=cluster[0]
              if qualities: res_quals[count]+=chr(cluster[1]+64)
            count+=1
        else:
          print "Error: Could not find",root_path+"/"+lanes[lane-1]+"/"+cycle
    else:
      print "Error: Could not find",root_path+"/"+lanes[lane-1]
      return None,None
  return res_bases,res_quals

def read_bcl_tile(filename,clusteridx=None,qualities=True):
  path,lane,tile,start,end = filename.split(';')
  #print path,int(lane),int(tile),int(start),int(end)
  reads,qualities = get_bcl_cycles(path,int(lane),int(tile),int(start),int(end),clusteridx,qualities)
  if reads != None:
    for ind,seq in enumerate(reads):
      yield (str(lane),str(tile),str(ind),"IDX"),seq,None
    raise StopIteration
