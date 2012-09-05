#!/usr/bin/env python
# -*- coding: ASCII -*-

"""

:Author: Martin Kircher
:Contact: Martin.Kircher@eva.mpg.de
:Date: *05.11.2008

"""
import sys,os
import string
import gzip
from optparse import OptionParser

table = string.maketrans('ATGC','TACG')

def compl_seq(seq):
  global table
  return seq.translate(table)

def revcompl_seq(seq):
  return compl_seq(rev_seq(seq))

def rev_seq(seq):
  res = list(seq)
  res.reverse()
  return "".join(res)

def read_sequence_file(filename):
  if filename.endswith('.gz'):
    GZinfile = gzip.GzipFile(filename)
    infile = GZinfile.read().splitlines()
  else:
    infile = open(filename)
  seqid = None
  seq = ""
  qual = -1
  fastq = None
  for line in infile:
    if (line.startswith("@") and (fastq != False) and ((qual == -1) or (len(qual)==len(seq)))) or (line.startswith(">") and (fastq != True)):
      # FIRST LINE DECIDES FILE FORMAT!
      if fastq == None:
        if line.startswith("@"): fastq = True
        else: fastq = False
      if len(seq) > 0:
        yield seqid,seq,qual
      seqid=line[1:].strip()
      seq = ""
      qual = None
    elif line.startswith("+") and (qual == None):
      qual = ""
    else:
      if (qual==None): seq+=line.strip()
      else: qual+=line.strip()
  if filename.endswith('.gz'):
    GZinfile.close()
  else:
    infile.close()
  if len(seq) > 0:
    yield seqid,seq,qual
  raise StopIteration

def parse_soap(filehandle,filter_length=-1,reverse_reads=False):
  def parse_edits(edit,rev=False,length=0,offset=0):
    fields = edit.split("->")
    if len(fields) == 2:
      res = [None,None,fields[0]]
      if rev:
        res = [None,None,compl_seq(fields[0])]
      pos = ""
      for x in range(len(fields[1])):
        if fields[1][x].isdigit():
          pos+=fields[1][x]
        else:
          break
      res[0] = offset+int(pos)+1
      if rev:
        res[0] = length-int(pos)
      res[1] = fields[1][x]
      if rev:
        res[1] = compl_seq(fields[1][x])
      return res
    else:
      print "Uups",edit
      return [1,None,None]
  while True:
    line = filehandle.readline()
    res = {}
    if not(line):
      raise StopIteration
    else:
      fields = line.split()
      if len(fields) >= 10:
        res['fragID']=fields[0]
        res['fragSeq']=fields[1]
        if reverse_reads: res['fragSeq']=revcompl_seq(fields[1])
        res['fragQual']=fields[2]
        res['nrHits']=int(fields[3])
        res['paired_ab']=fields[4]
        res['fragLen']=int(fields[5])
        res['refStrand']=fields[6]
        res['refName']=fields[7]
        res['refStart']=int(fields[8])
        res['nrEdits']=int(fields[9])
        res['edits']={}
        if res['refStrand'] == "+": # + STRAND
          if len(fields) > 10 and res['nrEdits'] < 100: # NO GAPS
            todo=fields[10:]
            edits = []
            for elem in range(len(todo)):
              if not reverse_reads: edit = parse_edits(todo[elem])
              else: edit = parse_edits(todo[elem],rev=True)
              res['edits'][edit[0]] = (edit[1],edit[2])
          elif (res['nrEdits'] < 200) and (len(fields) == 11): #INSERT
            length = res['nrEdits']-100
            res['edits']['insert']=(int(fields[10]),length)
          elif (res['nrEdits'] < 300) and (len(fields) == 11): #DELETION
            length = res['nrEdits']-200
            res['edits']['deletion']=(int(fields[10]),length)
          for ind,elem in enumerate(res['fragSeq']):
            if elem == 'N': res['edits'][ind+1]=('N','N')
        else: # - STRAND
          if len(fields) > 10 and res['nrEdits'] < 100: # NO GAPS
            todo=fields[10:]
            edits = []
            for elem in range(len(todo)):
              if not reverse_reads:  edit = parse_edits(todo[elem],rev=True,length=res['fragLen'])
              else: edit = parse_edits(todo[elem],length=res['fragLen'])
              res['edits'][edit[0]] = (edit[1],edit[2])
          elif (res['nrEdits'] < 200) and (len(fields) == 11): #INSERT
            length = res['nrEdits']-100
            res['edits']['insert']=(res['fragLen']-int(fields[10])-length,length)
          elif (res['nrEdits'] < 300) and (len(fields) == 11): #DELETION
            length = res['nrEdits']-200
            res['edits']['deletion']=(res['fragLen']-int(fields[10]),length)
          for ind,elem in enumerate(res['fragSeq'][::-1]):
            if elem == 'N': res['edits'][ind+1]=('N','N')
        if (filter_length < 0) or (filter_length <= res['fragLen']):
          yield res

parser = OptionParser()
parser.add_option("-1", "--bustard", dest="bustard", help="Bustard soap output file")
parser.add_option("-2", "--ibis", dest="ibis", help="Ibis soap output file")
parser.add_option("-i", "--input", dest="fasta", help="Input file with all reads (for counting)",default=None)
parser.add_option("-o", "--output", dest="output", help="Output file (default stdout)",default=None)
parser.add_option("-r", "--reverse", dest="reverse", help="Right align reads (default OFF)",default=False,action="store_true")
parser.add_option("-f", "--filter", dest="length_filter", help="Restrict analysis to reads longer than N nt (default OFF)",default=-1,type="int")
parser.add_option("-m", "--mismatch", dest="mismatch", help="Restrict analysis to reads with less or equal N mismatches (default 5)",default=5,type="int")
parser.add_option("-v", "--verbose", dest="verbose", help="Print out additional information (default OFF).",default=False,action="store_true")

(options, args) = parser.parse_args()

if (options.fasta == None or not(os.path.isfile(options.fasta))) or (options.bustard == None or not(os.path.isfile(options.bustard))) or (options.ibis == None or not(os.path.isfile(options.ibis))):
  print "Need valid input files."
  sys.exit()

outfile = sys.stdout
if options.output <> None:
  outfile = open(options.output,'w')
  if options.verbose: print "Open "+options.output+" for output..."

#Count raw reads
raw_count = 0
for seqid,seq,qual in read_sequence_file(options.fasta):
  raw_count += 1
outfile.write("Number of raw reads:\t%d\n"%raw_count)

#IDENTIFY READS MAPPED WITH BOTH APPROACHES:
infile = open(options.bustard)
bustard = set()
bustard_total = {}
for line in infile:
  fields = line.split()
  if len(fields) >= 10:
    edits = int(fields[9])+fields[1].count('N')
    if edits <= options.mismatch:
      if edits in bustard_total: bustard_total[edits]+=1
      else: bustard_total[edits]=1
      if fields[0][0].isdigit():
        bustard.add("_".join(fields[0].split('_')[1:]))
      else:
        bustard.add(fields[0])
infile.close()

infile = open(options.ibis)
ibis = set()
ibis_total = {}
for line in infile:
  fields = line.split()
  if len(fields) >= 10:
    edits = int(fields[9])+fields[1].count('N')
    if edits <= options.mismatch:
      if edits in ibis_total: ibis_total[edits]+=1
      else: ibis_total[edits]=1
      if fields[0][0].isdigit():
        ibis.add("_".join(fields[0].split('_')[1:]))
      else:
        ibis.add(fields[0])
infile.close()

to_sort = bustard_total.items()
to_sort.sort()
help = sum(bustard_total.values())
outfile.write("Bustard mapped:\t%d (%.2f%%)\n"%(help,(help/float(raw_count))*100))
for key,value in to_sort:
  outfile.write("%d\t%d\t%.2f%%\n"%(key,value,(value/float(raw_count))*100))

to_sort = ibis_total.items()
to_sort.sort()
help = sum(ibis_total.values())
outfile.write("Ibis mapped:\t%d (%.2f%%)\n"%(help,(help/float(raw_count))*100))
for key,value in to_sort:
  outfile.write("%d\t%d\t%.2f%%\n"%(key,value,(value/float(raw_count))*100))

intersection = bustard.intersection(ibis)
outfile.write("Number of mapped reads in intersection:\t%d\n"%len(intersection))

for filename in [options.bustard,options.ibis]:
  outfile.write("\n\n%s\n"%filename)
  infile = open(filename,'r')
  if options.verbose: print "Open "+filename+" as input..."

  errors = {}
  all_pairs = set()
  for elem1 in ['A','T','G','C']:
    for elem2 in ['A','T','G','C']:
      if elem1 != elem2:
        all_pairs.add((elem1,elem2))
  all_pairs.add(('N','N'))

  counts = {}
  for record in parse_soap(infile,options.length_filter,options.reverse):
    if record["fragID"][0].isdigit():
      record["fragID"]="_".join(record["fragID"].split('_')[1:])
    if record["fragID"] in intersection:
      if 'edits' in record:
        for pos,bpair in record['edits'].iteritems():
          if (pos in errors) and (bpair in errors[pos]):
            errors[pos][bpair]+=1
          elif (pos in errors) and (bpair not in errors[pos]):
            errors[pos][bpair]=1
          else:
            errors[pos]={}
            errors[pos][bpair]=1
      for elem in range(record['fragLen']):
        if elem+1 in counts:
          counts[elem+1]+=1
        else:
          counts[elem+1]=1

  all_pairs = list(all_pairs)
  all_pairs.sort()
  outfile.write('#Pos\t'+'\t'.join(map(lambda (x,y):x+"->"+y,all_pairs))+'\tTotal\n')
  ave_error = 0.0
  for cycle in range(1,len(counts)+1):
    outstr = str(cycle)
    for pair in all_pairs:
      if (cycle in errors) and (pair in errors[cycle]):
        outstr += "\t"+str(errors[cycle][pair]/float(counts[cycle]))
      else:
        outstr += "\t0"
    total = sum(errors[cycle].values())/float(counts[cycle])
    ave_error += total
    outfile.write(outstr+"\t"+str(total)+"\n")
  if len(counts) > 0: ave_error=ave_error/float(len(counts))
  outfile.write("\nAverage error rate:\t%.2f%%\n"%(ave_error*100))
  if options.verbose: print "Read",count,'results from',filename
  infile.close()

if options.output <> None:
  outfile.close()
