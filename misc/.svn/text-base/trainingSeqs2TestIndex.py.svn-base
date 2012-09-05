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
group = OptionGroup(parser, "Paths","Location of folders and files")
group.add_option("-i", "--infile", dest="infile", help="Sequence input file(s) to use for training (max 2, use comma as separator)")
group.add_option("-o", "--outfile", dest="outfiles", help="Path for output files (default .)",default=".")
group.add_option("-p", "--path", dest="path", help="Path to Firecrest/IPAR/Intensities folder (default .)", default=".")
group.add_option("--training", dest="SVMlight", help="SVM light multiclass training program (default "+def_svm_train+")",default=def_svm_train)
group.add_option("--prediction", dest="SVMlightTest", help="SVM light multiclass program used for estimating test error (default "+def_svm_train+")",default=def_svm_test)
group.add_option("--temp", dest="tmp", help="Path to temporary folder (default "+def_temp+")",default=def_temp)
group.add_option("--keep", dest="keep", help="Keep temporary prediction files in separate subfolder",default=False,action="store_true")
group.add_option("--NoTestDataSet", dest="NoTestDataSet", help="Don't separate data for a test data set",action="store_true",default=False)
group.add_option("--mock", dest="mock", help="Do a mock run for testing",default=False,action="store_true")
group.add_option("--verbose", dest="verbose", help="Print all status messages",default=False,action="store_true")
parser.add_option_group(group)

group = OptionGroup(parser, "Parameter","Parameters for creating the training data")
group.add_option("--indexlength", dest="indexlength", help="Length of Index read (default 0)",type="int",default=0)
group.add_option("-r", "--readlength", dest="readlength", help="Total read length (incl R1,Index,R2; default None)",type="int")
group.add_option("-u", "--use", dest="use", help="Use every Nth sequence to create training data (default 1)",type="int",default=1)
group.add_option("-s", "--skip", dest="skip", help="Use every Nth tile to create training data (default 1)",type="int",default=1)
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

def wait_jobs():
  global jobs
  while len(jobs) > 0:
    proc = jobs.pop(0)
    proc.wait()

def handle_jobs(cjob):
  global options
  global jobs
  if options.mock:
    print cjob
    return None
  if len(jobs) < options.cores:
    jobs.append(subprocess.Popen(cjob,shell=True))
  else:
    iteration = 0
    while len(jobs) >= options.cores:
      njobs = []
      for elem in jobs:
        if None == elem.poll():
          njobs.append(elem)
      jobs = njobs
      iteration += 1
      if len(jobs) >= options.cores:
        #DEAD LOCK?
        if iteration >= 20:
          wait_jobs()
        else: time.sleep(30)
    jobs.append(subprocess.Popen(cjob,shell=True))

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
options.path=options.path.rstrip("/")+"/"

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
  #name = options.tmp+"/"+timestamp+"_SVMlight_train_cycle_"+str(elem+1)+".txt"
  #train_names.append(name)
  #remove_files.append(name)
  # TEST DATA SET AND PREDICTION FILES
  name = options.tmp+"/"+timestamp+"_SVMlight_test_cycle_"+str(elem+1)+".txt"
  test_names.append(name)
  remove_files.append(name)
  name = options.tmp+"/"+timestamp+"_SVMlight_test_res_cycle_"+str(elem+1)+".txt"
  test_temp_names.append(name)
  remove_files.append(name)
  # MODEL OUTPUT FILE
  outfilenames.append(options.outfiles+"/SVMlight_model_cycle_"+str(elem+1)+".txt")

#train_core_files = {}
test_core_files = {}
for core in range(options.cores):
  # TRAINING DATA SET
  #name = options.tmp+"/"+timestamp+"_train_core_"+str(core)+".lst"
  #train_core_files[core]=[name]
  #if not options.mock: outfile = open(name,'w')
  #remove_files.append(name)
  #for elem in range(readlength):
    #name = options.tmp+"/"+timestamp+"_SVMlight_train_cycle_"+str(elem+1)+"_"+str(core)+".txt"
    #remove_files.append(name)
    #if not options.mock: outfile.write(name+'\n')
    #train_core_files[core].append(name)
  #if not options.mock: outfile.close()
  # TEST DATA SET
  name = options.tmp+"/"+timestamp+"_test_core_"+str(core)+".lst"
  test_core_files[core]=[name]
  if not options.mock: outfile = open(name,'w')
  remove_files.append(name)
  for elem in range(readlength):
    name = options.tmp+"/"+timestamp+"_SVMlight_test_cycle_"+str(elem+1)+"_"+str(core)+".txt"
    remove_files.append(name)
    if not options.mock: outfile.write(name+'\n')
    test_core_files[core].append(name)
  if not options.mock: outfile.close()

if len(ranges)!=len(infiles):
  print "Error: Something off with selecting parts of the complete read."
  sys.exit()

#print ranges
#sys.exit()
have_test = False
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
          if (options.NoTestDataSet or (count_tiles % 2 == 0)) and (count_tiles % options.skip == 0):
            #if options.verbose: print "Reading",cur_lane_tile,"for training"
            #eval_lines(line_stack,cur_lane_tile,train_core_files,crange)
            pass
          elif (count_tiles % options.skip == 0) and len(line_stack) > 0:
            if options.verbose: print "Reading",cur_lane_tile,"as test data set"
            eval_lines(line_stack,cur_lane_tile,test_core_files,crange)
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
      #if options.verbose: print "Reading",cur_lane_tile,"for training"
      #eval_lines(line_stack,cur_lane_tile,train_core_files,crange)
      pass
    elif (count_tiles % options.skip == 0):
      if options.verbose: print "Reading",cur_lane_tile,"as test data set"
      eval_lines(line_stack,cur_lane_tile,test_core_files,crange)
      have_test = True
  infile.close()
wait_external_jobs()

#print "Merging independent training files..."
#for core,values in train_core_files.iteritems():
  #values=values[1:] #REMOVE FILE WITH LIST
  #for ind,trainingfile in enumerate(train_names):
    #if not options.mock and os.path.isfile(values[ind]):
      #infile=open(values[ind],'r')
      #outfile=open(trainingfile,'a')
      #for line in infile:
        #outfile.write(line)
      #outfile.close()
      #infile.close()
  ##print "Last index",ind

print "Merging independent test files..."
for core,values in test_core_files.iteritems():
  values=values[1:]
  for ind,testfile in enumerate(test_names):
    if not options.mock and os.path.isfile(values[ind]):
      infile=open(values[ind],'r')
      outfile=open(testfile,'a')
      for line in infile:
        outfile.write(line)
      outfile.close()
      infile.close()
  #print "Last index",ind

#print "Training models..."
#for ind, file in enumerate(train_names):
  #if ind not in missing:
    #if options.verbose: print "Calling",options.SVMlight+" -c 1 -v 0 "+file+" "+outfilenames[ind]
    #handle_jobs(options.SVMlight+" -c 1 -v 0 "+file+" "+outfilenames[ind])
#wait_jobs()
#print ""

#if options.mock: print "Number of training file names:",len(train_names)
#names = train_names
names = None
if have_test:
  names = test_names
  if options.mock: print "Number of testing file names:",len(test_names)
if names == None:
  print "Missing test data set. This script is not supposed to work here in this situation!"
  sys.exit()

print "Creating test data..."
for ind,orgfilename in enumerate(names):
  if ind not in missing:
    if options.verbose: print "Calling",options.SVMlightTest+" -v 0 "+orgfilename+" "+outfilenames[ind]+" "+test_temp_names[ind]
    handle_jobs(options.SVMlightTest+" -v 0 "+orgfilename+" "+outfilenames[ind]+" "+test_temp_names[ind])
wait_jobs()
print ""

scaling_logit = []
for ind,orgfilename in enumerate(names):
  if (ind not in missing) and (not options.mock):
    infile = open(orgfilename)
    classes = []
    distancies = {}
    cor_class = {}
    for line in infile:
      cclass = int(line.split()[0])
      classes.append(cclass)
      if cclass not in distancies: 
        distancies[cclass]=([],[],[],[])
        distancies[cclass]=([],[],[],[])
      if cclass not in cor_class:
        cor_class[cclass]=[1,1,1,1] # finite sample estimate - next one is in class
    infile.close()

    miss = 0
    infile = open(test_temp_names[ind])
    for cclass in classes:
      fields = infile.readline().split()
      if len(fields) == 5:
        if cclass > 0:
          cclass2 = int(fields[0])
          if cclass != cclass2: miss += 1
          cor_class[cclass][cclass2-1] += 1
          distancies[cclass][0].append(float(fields[1]))
          distancies[cclass][1].append(float(fields[2]))
          distancies[cclass][2].append(float(fields[3]))
          distancies[cclass][3].append(float(fields[4]))
        else:
          print "Problems when extracting error."
      else:
        print "Unexpected line in",test_temp_names[ind],". May be class is missing in training data..."
    infile.close()

    if len(classes) > 0:
      rate = miss/float(len(classes))
      print "Misclassification rate for %d. model: %.4f%%"%(ind+1,rate*100.0)
      del classes

    helper = []
    order = distancies.keys()
    order.sort()
    for cclass in order:
      #print "Determined the following parameters for the normal distribution in class "+str(cclass)
      A_mean,A_sd = mean_sd(distancies[cclass][0])
      C_mean,C_sd = mean_sd(distancies[cclass][1])
      G_mean,G_sd = mean_sd(distancies[cclass][2])
      T_mean,T_sd = mean_sd(distancies[cclass][3])
      total_obs = float(sum(cor_class[cclass]))
      ra = cor_class[cclass][0]/total_obs
      rc = cor_class[cclass][1]/total_obs
      rg = cor_class[cclass][2]/total_obs
      rt = cor_class[cclass][3]/total_obs
      #print "A:\t"+str(A_mean)+"\t"+str(A_sd)+"\t"+str(ra)
      #print "C:\t"+str(C_mean)+"\t"+str(C_sd)+"\t"+str(rc)
      #print "G:\t"+str(G_mean)+"\t"+str(G_sd)+"\t"+str(rg)
      #print "T:\t"+str(T_mean)+"\t"+str(T_sd)+"\t"+str(rt)
      #print "Estimates based on:\t"+str(len(distancies[cclass][0]))+" observations"
      if A_mean == None or C_mean == None or G_mean == None or T_mean == None or (len(distancies[cclass][0]) < 2) or (len(distancies[cclass][1]) < 2) or (len(distancies[cclass][2]) < 2)or (len(distancies[cclass][3]) < 2):
        print "Estimating parameters of the normal distributions failed (use a bigger test data set!). Falling back to arbitrary defaults."
        A_mean,A_sd = -125.0,100.0
        C_mean,C_sd = -125.0,100.0
        G_mean,G_sd = -125.0,100.0
        T_mean,T_sd = -125.0,100.0
        ra,rc,rg,rt = 0.02,0.02,0.02,0.02
        if cclass == 1:
          A_mean = 200
          ra = 0.94
        elif cclass == 2:
          C_mean = 200
          rc = 0.94
        elif cclass == 3:
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
    scaling_logit.append(helper)
  else:
    scaling_logit.append([])

order = assignments.items()
order.sort()

if not options.mock:
  outfile = open(options.outfiles+"/SVMlight_models.index",'w')
  for cycle,(imodel,ctype) in order:
    outfile.write(str(cycle+1)+"\t"+outfilenames[imodel]+"\t"+ctype+"\t"+"\t".join(map(str,scaling_logit[imodel]))+"\n")
  outfile.close()
print "Index file available as",options.outfiles+"/SVMlight_models.index"
wait_jobs()

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
