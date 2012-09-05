#!/usr/bin/env python
# -*- coding: ASCII -*-

"""

:Author: Martin Kircher
:Contact: Martin.Kircher@eva.mpg.de
:Date: *23.02.2009

"""

import sys,os
import shelve
from optparse import OptionParser
from optparse import OptionGroup
import gzip
import time
import subprocess
import math

if (len(sys.argv[0]) > 0)  and os.path.isdir(os.getcwd()+"/"+os.path.dirname(sys.argv[0])):
  root_bin = os.path.dirname(os.getcwd()+"/"+sys.argv[0])
elif (len(sys.argv[0]) > 0)  and os.path.isdir(os.path.dirname(sys.argv[0])):
  root_bin = os.path.dirname(sys.argv[0])
else:
  root_bin = os.getcwd()
f=open(root_bin+'/params.py')
c=compile(f.read(), root_bin+'/params.py', "exec")
eval(c)

sub_script = def_ibis_path+'trainingSeqs2SVMlight_helper.py'
nr_id_fields_seqs = 6

def mean_sd(vector):
  try:
    nr = float(len(vector))
    mean = sum(vector)/nr
    sd = math.sqrt(1.0/(nr-1.0)*sum(map(lambda x: (x-mean)*(x-mean),vector)))
    return mean,sd
  except:
    return None,None

#####################
# USER PARAMETER ...
#####################

parser = OptionParser(usage="usage: %prog [options]")
parser.add_option("-i", "--infile", dest="infile", help="Sequence input file(s) to use for training (max 2, use comma as separator)")
parser.add_option("--NoTestDataSet", dest="NoTestDataSet", help="Don't separate data for a test data set",action="store_true",default=False)

group = OptionGroup(parser, "Parameter","Parameters for creating the training data")
group.add_option("--indexlength", dest="indexlength", help="Length of Index read (default 0)",type="int",default=0)
group.add_option("-r", "--readlength", dest="readlength", help="Total read length (incl R1,Index,R2; default None)",type="int")
group.add_option("-u", "--use", dest="use", help="Use every Nth sequence to create training data (default 1)",type="int",default=1)
group.add_option("-s", "--skip", dest="skip", help="Use every Nth tile to create training data (default 1)",type="int",default=1)
parser.add_option_group(group)
(options, args) = parser.parse_args()

infiles = []
if (options.infile == None):
  print "Need input file with sequences for training."
  sys.exit()
elif (options.infile.count(",") == 1):
   fields = options.infile.split(",")
   if not(os.path.isfile(fields[0])) or not(os.path.isfile(fields[1])):
     print "Need valid input files with sequences for training."
     sys.exit()
   else:
     infiles.append(fields[0])
     infiles.append(fields[1])
elif os.path.isfile(options.infile):
  infiles.append(options.infile)
else:
  print "Need (maximum two) valid input files with sequences for training."
  sys.exit()

readlength = options.readlength
lstart = 0
ranges=[]
for filename in list(infiles): ## WE MAY CHANGE INFILES IN CASE OF INDEX READS - MAKE SURE WE DO NOT INFLUENCE THE FOR LOOP
  try:
    infile = open(filename)
    fields = map(int,infile.readline().split()[nr_id_fields_seqs-2:nr_id_fields_seqs])
    if (lstart > 0) and (options.indexlength > 0): #FIX INDEX WHEN SEEING SECOND TRAINING DATA SET FILE
      ranges.append((lstart,ranges[-1][1],min(lstart+options.indexlength,fields[0])-lstart+ranges[-1][1],min(lstart+options.indexlength,fields[0])))
      lstart = min(lstart+options.indexlength,fields[0])
      options.indexlength=0
      infiles.insert(0,infiles[0])
    if options.readlength == None:
      readlength = max(fields[1]+options.indexlength,readlength)
    ranges.append((fields[0],fields[0],fields[1],fields[1]))
    lstart=fields[1]
    infile.close()
  except None:
    print "Unexpected file format found when extracting sequence starts and ends."
    sys.exit()

## SINGLE READ RUN, HAVE TO FIX INDEX IN ANOTHER WAY
if (lstart > 0) and (options.indexlength > 0):
  ranges.append((lstart,ranges[-1][1],min(lstart+options.indexlength,readlength)-lstart+ranges[-1][1],min(lstart+options.indexlength,readlength)))
  lstart = min(lstart+options.indexlength,readlength)
  options.indexlength=0
  infiles.insert(0,infiles[0])

## MAKE SURE LAST ELEMENT IN RANGE IS READLENGTH
if (ranges[-1][3] != readlength):
  ranges[-1]=(ranges[-1][0],ranges[-1][1],ranges[-1][2],readlength)

#print ranges
#sys.exit()
have_test = False
reads = {}
for ind,infilename in enumerate(infiles):
  print "Reading",infilename
  crange = ranges[ind]
  infile = open(infilename)
  cur_lane_tile = None
  line_stack = {}
  count = 0
  count_tiles = 0
  reads[ind] = ([],[])
  for line in infile:
    if (count % options.use == 0):
      fields = line.split()
      if len(fields) == nr_id_fields_seqs+1:
        tile = fields[1]
        while len(tile) < 4:
          tile="0"+tile
        if (cur_lane_tile == None):
          cur_lane_tile = fields[0]+"_"+tile
        elif (cur_lane_tile != fields[0]+"_"+tile):
          if (options.NoTestDataSet or (count_tiles % 2 == 0)) and (count_tiles % options.skip == 0):
            #print "Reading",cur_lane_tile,"for training"
            #eval_lines(line_stack,cur_lane_tile,train_core_files,crange)
            reads[ind][0].append(cur_lane_tile)
          elif (count_tiles % options.skip == 0) and len(line_stack) > 0:
              #print "Reading",cur_lane_tile,"as test data set"
              reads[ind][1].append(cur_lane_tile)
              #eval_lines(line_stack,cur_lane_tile,test_core_files,crange)
              have_test = True
          line_stack = {}
          cur_lane_tile = fields[0]+"_"+tile
          count_tiles+=1
        line_stack[tuple(fields[:(nr_id_fields_seqs-2)])]=fields[nr_id_fields_seqs]
      else:
        print "Unexpected line in training input:",line
    count+=1
  if len(line_stack) > 0:
    if (options.NoTestDataSet or (count_tiles % 2 == 0)) and (count_tiles % options.skip == 0):
      #print "Reading",cur_lane_tile,"for training"
      #eval_lines(line_stack,cur_lane_tile,train_core_files,crange)
      reads[ind][0].append(cur_lane_tile)
    elif (count_tiles % options.skip == 0):
      #"Reading",cur_lane_tile,"as test data set"
      reads[ind][1].append(cur_lane_tile)
      #eval_lines(line_stack,cur_lane_tile,test_core_files,crange)
      have_test = True
  infile.close()

for read,value in reads.iteritems():
  print "Read",read+1
  print "Training:"
  print value[0]
  print "Test:"
  print value[1]
  print