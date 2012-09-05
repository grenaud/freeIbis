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
import random

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
group = OptionGroup(parser, "Paths","Location of folders and files")
group.add_option("-i", "--infile", dest="infile", help="Sequence input file(s) to use for training (max 2, use comma as separator)")
group.add_option("-o", "--outfile", dest="outfiles", help="Path for output files (default .)",default=".")
group.add_option("-p", "--path", dest="path", help="Path to Firecrest/IPAR/Intensities folder (default .)", default=".")
group.add_option("-e", "--epsilon", dest="epsilon", help="Termination criterion", type='float',default=0.01)
group.add_option("--training", dest="SVMlight", help="SVM light multiclass training program (default "+def_svm_train+")",default=def_svm_train)
group.add_option("--prediction", dest="SVMlightTest", help="SVM light multiclass program used for estimating test error (default "+def_svm_train+")",default=def_svm_test)
group.add_option("--temp", dest="tmp", help="Path to temporary folder (default "+def_temp+")",default=def_temp)
group.add_option("--keep", dest="keep", help="Keep temporary prediction files in separate subfolder",default=False,action="store_true")
group.add_option("--mock", dest="mock", help="Do a mock run for testing",default=False,action="store_true")
group.add_option("--verbose", dest="verbose", help="Print all status messages",default=False,action="store_true")
parser.add_option_group(group)

group = OptionGroup(parser, "Parameter","Parameters for creating the training data")
group.add_option("--indexlength", dest="indexlength", help="Length of Index read (default 0)",type="int",default=0)
group.add_option("-r", "--readlength", dest="readlength", help="Total read length (incl R1,Index,R2; default None)",type="int")
group.add_option("-u", "--use", dest="use", help="Use every Nth sequence to create training data (default 1)",type="int",default=1)
group.add_option("-c", "--cores", dest="cores", help="Number of CPU processes to be used (default 1)",type="int",default=1)
parser.add_option_group(group)
(options, args) = parser.parse_args()

jobs = []
train_names = []
test_names = []
test_temp_names = []
remove_files = []
outfilenames = []
timestamp = str(time.time())
external_jobs = []
free_ids = None

def handle_jobs(cjob):
  global options
  global jobs
  if options.mock:
    print cjob
    return None
  if len(jobs) < options.cores:
    jobs.append(subprocess.Popen(cjob,shell=True))
  else:
    njobs = []
    for elem in jobs:
      if not elem.poll():
        njobs.append(elem)
    jobs = njobs
    if len(jobs) >= options.cores:
      proc = jobs.pop(0)
      proc.wait()
    jobs.append(subprocess.Popen(cjob,shell=True))

def wait_jobs():
  global jobs
  while len(jobs) > 0:
    proc = jobs.pop(0)
    proc.wait()

def eval_lines(trainlines,lane_tile,files_dict,crange):
  global external_jobs, free_ids
  global options, timestamp, remove_files
  global sub_script

  if free_ids == None: 
    free_ids = range(options.cores)
    for elem in free_ids:
      remove_files.append(options.tmp+"/"+timestamp+"_lines_"+str(elem)+".shelve")

  if options.mock:
    cur_jobid = 0
    name = options.tmp+"/"+timestamp+"_lines_"+str(cur_jobid)+".shelve"
    print "Calling",sub_script+" --lines='"+name+"' --outfilesindex='"+files_dict[cur_jobid][0]+"' --lane_tile='"+lane_tile+"' --crange='"+str(crange)+"' -p "+options.path
  else:
    if len(free_ids) > 0:
      cur_jobid = free_ids.pop(0)
      name = options.tmp+"/"+timestamp+"_lines_"+str(cur_jobid)+".shelve"
      outfile = shelve.open(name,'n')
      for elem,value in trainlines.iteritems():
        outfile[str(elem)]=value
      outfile.close()
      cjob = sub_script+" --lines='"+name+"' --outfilesindex='"+files_dict[cur_jobid][0]+"' --lane_tile='"+lane_tile+"' --crange='"+str(crange)+"' -p "+options.path
      external_jobs.append((cur_jobid,subprocess.Popen(cjob,shell=True)))
    else:
      njobs = []
      for cur_jobid,proc in external_jobs:
        if not proc.poll():
          njobs.append((cur_jobid,proc))
        else:
          free_ids.append(cur_jobid)
      external_jobs = njobs
      if  len(free_ids) == 0:
        cur_jobid,proc = external_jobs.pop(0)
        proc.wait()
        free_ids.append(cur_jobid)

      cur_jobid = free_ids.pop(0)
      name = options.tmp+"/"+timestamp+"_lines_"+str(cur_jobid)+".shelve"
      outfile = shelve.open(name,'n')
      for elem,value in trainlines.iteritems():
        outfile[str(elem)]=value
      outfile.close()
      cjob = sub_script+" --lines='"+name+"' --outfilesindex='"+files_dict[cur_jobid][0]+"' --lane_tile='"+lane_tile+"' --crange='"+str(crange)+"' -p "+options.path
      external_jobs.append((cur_jobid,subprocess.Popen(cjob,shell=True)))

def wait_external_jobs():
  global external_jobs,free_ids
  while len(external_jobs) > 0:
    cid,proc = external_jobs.pop(0)
    proc.wait()
  free_ids = None

if options.mock:
  options.verbose=True

if (options.tmp == None) or not(os.path.isdir(options.tmp)):
  print "Need path to temporary folder."
  sys.exit()
options.tmp=options.tmp.rstrip("/")

if (options.path == None) or not(os.path.isdir(options.path)):
  print "Need path to Firecrest folder."
  sys.exit()
options.path=options.path.rstrip("/")

if (options.outfiles == None) or not(os.path.isdir(options.outfiles)):
  print "Need path for output files."
  sys.exit()
options.outfiles=options.outfiles.rstrip("/")

if not(os.path.isfile(options.SVMlight)):
  print "Can't access training routine",options.SVMlight
  sys.exit()

if not(os.path.isfile(options.SVMlightTest)):
  print "Can't access prediction routine",options.SVMlightTest
  sys.exit()

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

missing = set()
assignments = {}
reassignments = {}
if len(ranges) > 0: print "Have the following sequence parts:"
lstart = 0
for fragStart,seqStart,seqEnd,fragEnd in ranges:
  print "%d-%d -> trained on %d-%d"%(fragStart+1,fragEnd,seqStart+1,seqEnd)
  for elem in range(fragStart,fragEnd):
    if elem == fragStart:
      assignments[elem]=(elem,'B')
    elif elem == fragEnd-1:
      assignments[elem]=(elem,'E')
    else:
      assignments[elem]=(elem,'M')
  if (lstart < fragStart):
    missing.add(lstart)
    reassignments[lstart] = (fragStart,'B')
    assignments[lstart] = (fragStart,'B')
    reassignments[fragStart] = (fragStart+1,'M')
    assignments[fragStart] = (fragStart+1,'M')
    lstart += 1
  while (lstart < fragStart):
    missing.add(lstart)
    reassignments[lstart] = (fragStart+1,'M')
    assignments[lstart] = (fragStart+1,'M')
    lstart += 1
  lstart = fragEnd
if len(ranges) > 0: print ""
if lstart < readlength:
  missing.add(readlength-1)
  reassignments[readlength-1] = (fragEnd,'E')
  assignments[readlength-1] = (fragEnd,'E')
  reassignments[lstart-1] = (fragEnd-1,'M')
  assignments[lstart-1] = (fragEnd-1,'M')
  for elem in range(lstart,readlength-1):
    missing.add(elem)
    reassignments[elem] = (fragEnd-1,'M')
    assignments[elem] = (fragEnd-1,'M')

print "Readlength:",readlength
if len(reassignments) > 0:
  print ""
  print "For %d bases (%s) no traning data is assigned, will do the following (re)assignments:"%(len(missing),",".join(map(lambda x:str(x+1),missing)))
  for pos,replace in reassignments.iteritems():
    print "  Pos:",pos+1,"Using base:",replace[0]+1,"Type:",replace[1]
  print ""

for elem in range(readlength):
  # TRAINING DATA SET
  name = options.tmp+"/"+timestamp+"_SVMlight_train_cycle_"+str(elem+1)+".txt"
  train_names.append(name)
  remove_files.append(name)
  # TEST DATA SET AND PREDICTION FILES
  name = options.tmp+"/"+timestamp+"_SVMlight_test_cycle_"+str(elem+1)+".txt"
  test_names.append(name)
  remove_files.append(name)
  name = options.tmp+"/"+timestamp+"_SVMlight_test_res_cycle_"+str(elem+1)+".txt"
  test_temp_names.append(name)
  remove_files.append(name)
  # MODEL OUTPUT FILE
  outfilenames.append(options.outfiles+"/SVMlight_model_cycle_"+str(elem+1)+".txt")

train_core_files = {}
for core in range(options.cores):
  # TRAINING DATA SET
  name = options.tmp+"/"+timestamp+"_train_core_"+str(core)+".lst"
  train_core_files[core]=[name]
  if not options.mock: outfile = open(name,'w')
  remove_files.append(name)
  for elem in range(readlength):
    name = options.tmp+"/"+timestamp+"_SVMlight_train_cycle_"+str(elem+1)+"_"+str(core)+".txt"
    remove_files.append(name)
    if not options.mock: outfile.write(name+'\n')
    train_core_files[core].append(name)
  if not options.mock: outfile.close()

if len(ranges)!=len(infiles):
  print "Error: Something off with selecting parts of the complete read."
  sys.exit()

for ind,infilename in enumerate(infiles):
  print "Reading",infilename
  crange = ranges[ind]
  infile = open(infilename)
  cur_lane_tile = None
  line_stack = {}
  count = 0
  count_tiles = 0
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
          if options.verbose: print "Reading",cur_lane_tile,"for training"
          eval_lines(line_stack,cur_lane_tile,train_core_files,crange)
          line_stack = {}
          cur_lane_tile = fields[0]+"_"+tile
          count_tiles+=1
        line_stack[tuple(fields[:(nr_id_fields_seqs-2)])]=fields[nr_id_fields_seqs]
      else:
        print "Unexpected line in training input:",line
    count+=1
  if len(line_stack) > 0:
    if options.verbose: print "Reading",cur_lane_tile,"for training"
    eval_lines(line_stack,cur_lane_tile,train_core_files,crange)
  infile.close()
wait_external_jobs()

print "Merging independent training files..."
for core,values in train_core_files.iteritems():
  values=values[1:]
  for ind,trainingfile in enumerate(train_names):
    if not options.mock and os.path.isfile(values[ind]):
      infile=open(values[ind],'r')
      outfile=open(trainingfile,'a')
      for line in infile:
        outfile.write(line)
      outfile.close()
      infile.close()
  #print "Last index",ind

print "Training models..."
for ind, file in enumerate(train_names):
  if ind not in missing:
    if options.verbose: print "Calling",options.SVMlight+" -c 1 -v 0 "+file+" "+outfilenames[ind]
    handle_jobs(options.SVMlight+" -c 1 -v 0 "+file+" "+outfilenames[ind])
wait_jobs()
print ""

iteration = 0
SVMfinished = set()
scaling_logit = []
for ind,orgfilename in enumerate(train_names):
  scaling_logit.append([])

while True:
  iteration += 1
  total_miss = 0.0
  total_count = 0
  print "Creating predictions for training data..."
  for ind,orgfilename in enumerate(train_names):
    if ind not in missing and (ind not in SVMfinished):
      if options.verbose: print "Calling",options.SVMlightTest+" -v 0 "+orgfilename+" "+outfilenames[ind]+" "+test_temp_names[ind]
      handle_jobs(options.SVMlightTest+" -v 0 "+orgfilename+" "+outfilenames[ind]+" "+test_temp_names[ind])
  wait_jobs()
  print ""

  for ind,orgfilename in enumerate(train_names):
    if (ind not in missing) and (not options.mock) and (ind not in SVMfinished):
      infile = open(orgfilename)
      classes = []
      distancies = [([],[],[],[]),([],[],[],[]),([],[],[],[]),([],[],[],[])]
      for line in infile:
        cclass = int(line.split()[0])
        classes.append((cclass,line))
      infile.close()

      infile = open(test_temp_names[ind])
      message = True
      for line in infile:
        fields = line.split()
        if len(fields) == 5:
          cclass2 = int(fields[0])
          if cclass2 > 0:
            distancies[cclass2-1][0].append(float(fields[1]))
            distancies[cclass2-1][1].append(float(fields[2]))
            distancies[cclass2-1][2].append(float(fields[3]))
            distancies[cclass2-1][3].append(float(fields[4]))
          else:
            print "Problems when extracting error."
        elif message:
          print "Unexpected line in",test_temp_names[ind],". May be class is missing in training data..."
          message = False
      infile.close()
      A_mean,A_sd = mean_sd(distancies[0][0])
      C_mean,C_sd = mean_sd(distancies[1][1])
      G_mean,G_sd = mean_sd(distancies[2][2])
      T_mean,T_sd = mean_sd(distancies[3][3])
      cutoff_called = [0,0,0,0]
      cutoff_called[0]=A_mean-1.5*A_sd
      cutoff_called[1]=C_mean-1.5*C_sd
      cutoff_called[2]=G_mean-1.5*G_sd
      cutoff_called[3]=T_mean-1.5*T_sd

      miss = 0
      infile = open(test_temp_names[ind])
      outfile = open(orgfilename,'w')
      #hadclass = [0,0,0,0,0]
      for cclass,pline in classes:
        fields = infile.readline().split()
        if len(fields) == 5:
          cclass2 = int(fields[0])
          if cclass2 > 0:
            if (cclass != cclass2): miss += 1
            if (cclass != cclass2) and (cutoff_called[cclass2-1]*random.random() <= float(fields[cclass2])):
              #hadclass[cclass2]+=1
              outfile.write(fields[0]+' '+' '.join(pline.split(' ')[1:]))
            else:
              outfile.write(pline)
          else:
            print "Problems when extracting error."
        elif message:
          print "Unexpected line in",test_temp_names[ind],". May be class is missing in training data..."
          message = False
      infile.close()
      ##TRY TO NOT LOSS A CLASS IN THIS PROCESS, IF NECESSARY GET TRAINING DATA FROM NEIGHBORING CYLCE
      #for cclass,ccount in enumerate(hadclass):
        #if (cclass > 0) and (ccount/float(sum(hadclass)) < 0.10) and (ind > 0) and (ind < len(train_names)-1):
          #goal = round((0.1*sum(hadclass)-ccount)/0.9)
          #infile = open(train_names[ind-1])
          #for line in infile:
            #cclass2 = int(line.split()[0])
            #if cclass2 == cclass:
              #outfile.write(line)
              #hadclass[cclass2]+=1
              #goal-=1
            #if goal == 0: break
          #infile.close()
      outfile.close()

      rate = 1.0
      if len(classes) > 0:
        rate = miss/float(len(classes))
        print "Reclassification rate for %d. model in iteration %2d: %.4f%%"%(ind+1,iteration,rate*100.0)
        del classes

      if (rate < options.epsilon) and (iteration > 1):
        SVMfinished.add(ind)
      else:
        total_miss += rate
        total_count += 1
        if options.verbose: print "Calling",options.SVMlight+" -c 1 -v 0 "+orgfilename+" "+outfilenames[ind]
        handle_jobs(options.SVMlight+" -c 1 -v 0 "+orgfilename+" "+outfilenames[ind])

      helper = []
      for cclass in range(4):
        A_mean,A_sd = mean_sd(distancies[cclass][0])
        C_mean,C_sd = mean_sd(distancies[cclass][1])
        G_mean,G_sd = mean_sd(distancies[cclass][2])
        T_mean,T_sd = mean_sd(distancies[cclass][3])
        ra,rc,rg,rt = 0.01,0.01,0.01,0.01
        if cclass == 0: ra = 0.97
        elif cclass == 1: rc = 0.97
        elif cclass == 2: rg = 0.97
        else: rt = 0.97
        if A_mean == None or C_mean == None or G_mean == None or T_mean == None or (len(distancies[cclass][0]) < 2) or (len(distancies[cclass][1]) < 2) or (len(distancies[cclass][2]) < 2)or (len(distancies[cclass][3]) < 2):
          print "Estimating parameters of the normal distributions failed (use a bigger test data set!). Falling back to arbitrary defaults."
          A_mean,A_sd = -125.0,100.0
          C_mean,C_sd = -125.0,100.0
          G_mean,G_sd = -125.0,100.0
          T_mean,T_sd = -125.0,100.0
          ra,rc,rg,rt = 0.02,0.02,0.02,0.02
          if cclass == 0:
            A_mean = 200
            ra = 0.94
          elif cclass == 1:
            C_mean = 200
            rc = 0.94
          elif cclass == 2:
            G_mean = 200
            rg = 0.94
          else:
            T_mean = 200
            rt = 0.94
        helper.append(A_mean)
        helper.append(A_sd)
        helper.append(ra)
        helper.append(C_mean)
        helper.append(C_sd)
        helper.append(rc)
        helper.append(G_mean)
        helper.append(G_sd)
        helper.append(rg)
        helper.append(T_mean)
        helper.append(T_sd)
        helper.append(rt)
      helper.append(rate)
      scaling_logit[ind]=helper

  wait_jobs()
  #print iteration,total_miss,total_count,total_miss/(total_count+iteration-1.0)
  print len(SVMfinished),"models are finished after this iteration."
  if (total_miss/(total_count+iteration*round(math.log(max(iteration-5,1),3))+iteration-1.0) < options.epsilon) and (iteration > 1): break

if not(options.keep):
  print "Removing temporary files..."
  for elem in remove_files:
    if os.path.isfile(elem):
      handle_jobs("rm "+elem)
else:
  if (not os.path.isdir(options.outfiles+"/Training/")) and (not options.mock):
    os.makedirs(options.outfiles+"/Training/")
  for elem in remove_files:
    if os.path.isfile(elem):
      handle_jobs("mv "+elem+" "+options.outfiles+"/Training/")

order = assignments.items()
order.sort()

if not options.mock:
  outfile = open(options.outfiles+"/SVMlight_models.index",'w')
  for cycle,(imodel,ctype) in order:
    outfile.write(str(cycle+1)+"\t"+outfilenames[imodel]+"\t"+ctype+"\t"+"\t".join(map(str,scaling_logit[imodel]))+"\n")
  outfile.close()
print "Index file available as",options.outfiles+"/SVMlight_models.index"
wait_jobs()
