#!/usr/bin/env python
# -*- coding: ASCII -*-

"""

:Author: Martin Kircher
:Contact: Martin.Kircher@eva.mpg.de
:Date: *16.04.2010

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
maximum_nr_clusters = 250000 
training_seq_limit = 20*maximum_nr_clusters 

def timeString():
    return str(time.strftime(" on %d %b %Y %H:%M:%S", time.localtime())+" secs = "+str(time.time()));


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
group.add_option("--prediction", dest="SVMlightTest", help="SVM light multiclass program used for estimating test error (default "+def_svm_test+")",default=def_svm_test)
group.add_option("--temp", dest="tmp", help="Path to temporary folder (default "+def_temp+")",default=def_temp)
group.add_option("--keep", dest="keep", help="Keep temporary prediction files in separate subfolder",default=False,action="store_true")
group.add_option("--NoTestDataSet", dest="NoTestDataSet", help="Don't separate data for a test data set",action="store_true",default=False)
group.add_option("--mock", dest="mock", help="Do a mock run for testing",default=False,action="store_true")
group.add_option("--verbose", dest="verbose", help="Print all status messages",default=False,action="store_true")
group.add_option("--quick", dest="quick", help="Use linear regression for prediction, quicker but less accurate",action="store_true",default=False)
parser.add_option_group(group)

group = OptionGroup(parser, "Parameter","Parameters for creating the training data")
group.add_option("--indexlength", dest="indexlength", help="Length of Index read (default 0)",type="int",default=0)
group.add_option("--2nd_indexlength", dest="indexlength2", help="Length of second index read (default 0)",type="int",default=0)
group.add_option("-t", "--coordianteType", dest="coordtype", help="Type of cluster coordinate transformation: round = Pipeline <= 1.4: round(x), floor = Pipeline 1.5.*: floor(abs(x)), shift_floor (default) = Pipeline 1.6.*: round(x*10+1000)",choices=["round","floor","shift_round"],default="shift_round")
group.add_option("-r", "--readlength", dest="readlength", help="Total read length (incl R1,Index,R2; default None)",type="int")
group.add_option("-u", "--use", dest="use", help="Use every Nth sequence to create training data (default 1)",type="int",default=1)
group.add_option("-s", "--skip", dest="skip", help="Use every Nth tile to create training data (default 1)",type="int",default=1)
group.add_option("-c", "--cores", dest="cores", help="Number of CPU processes to be used (default 1)",type="int",default=1)
parser.add_option_group(group)

group = OptionGroup(parser, "Recalibration","Parameters for the recalibration of quality scores")
group.add_option("--recalibration", dest="recalibration", help="Trigger quality score recalibration using control reads ",action="store_true",default=False)
group.add_option("--plotqual", dest="plotqual", help="Plot quality score recalibration lines (requires R)",action="store_true",default=False)

group.add_option("--trainingrecal",   dest="SVMlightRECAL", help="SVM light binary training program for recalibration (default "+def_svm_trainbinary+")",default=def_svm_trainbinary)
group.add_option("--predictionrecal", dest="SVMlightTestRECAL", help="SVM light binary program used for recalibration (default "+def_svm_testbinary+")",default=def_svm_testbinary)
group.add_option("--cgroup", dest="cgroup", help="Group consecutive this many cycles together for quality score recalibration default:1 (useful for low amounts of control sequences ex: miseq)",type="int",default=1)

parser.add_option_group(group)
(options, args) = parser.parse_args()

jobs = []
train_names = []
train_temp_names = []

test_names = []
test_temp_names = []

recal_names = []
recal_names1 = []
recal_names2 = []
recal_names3 = []
recal_names4 = []

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

def launchSimpleJob(cjob):
  while(True):
    if(options.mock):
      break;
    myjob=subprocess.Popen(cjob,shell=True);
    if (myjob.wait() == 0):
      break;
    else:
      print "WARNING: command "+myjob+" failed, restarting";


def handle_jobs(cjob,timeToSleep=30):
  global options
  global jobs
  if options.mock:
    print cjob
    return None
  if len(jobs) < options.cores:
    #jobs.append(subprocess.Popen(cjob,shell=True));
    jobcreated=subprocess.Popen(cjob,shell=True);
    jobs.append(jobcreated);
    return jobcreated;
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
        else: 
          time.sleep(timeToSleep)
    #jobs.append(subprocess.Popen(cjob,shell=True))
    jobcreated=subprocess.Popen(cjob,shell=True);
    jobs.append(jobcreated);
    return jobcreated;


def handleListOfjobs(alljobs,timeToSleep=30):
  unresolved=alljobs[:]; #copy
  numberOfIterations=0;
  while(len(unresolved) != 0):
    numberOfIterations+=1;
    arrayOfJobs=[]; #array of 2-upple (cmd string,Popen object)

    #launch jobs
    for jobToAdd in unresolved:
      arrayOfJobs.append([jobToAdd,handle_jobs(jobToAdd,timeToSleep)]);

    #wait for jobs
    if options.mock:
      break;
    while(1):
      allFinished=True;
      for toverify in arrayOfJobs:
        allFinished = allFinished and (toverify[1].poll() != None); #if one fails, they will all be false
      if(allFinished):
        break;
      else:
        time.sleep(4);

    #all are done, check return code
    unresolved=[]; 
    for proc in arrayOfJobs:
      if(proc[1].returncode != 0): #wrong code
        print "WARNING: process "+proc[0]+" failed, will be relaunched";
        unresolved.append(proc[0]);# we will loop on the failed jobs

##
# eval lines()
# called using the data accumulated from the training_seqs file and calls trainingSeqs2SVMlight_helper.py
# trainlines is a dictionary tuple - > sequence i.e : ('4', '2', '1055', '1211'): 'AACAATTTTAATTGCAGGGGCTTCGGCCCCTTACTTGAGGATAAATTATGT' 
#                                                         (lane,tile,x,y) - > sequence from training_seqs 
#
##
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
    print "Calling",sub_script+" -t "+options.coordtype+" --lines='"+name+"' --outfilesindex='"+files_dict[cur_jobid][0]+"' --lane_tile='"+lane_tile+"' --crange='"+str(crange)+"' -p "+options.path
  else:
    if len(free_ids) > 0:
      cur_jobid = free_ids.pop(0)
      name = options.tmp+"/"+timestamp+"_lines_"+str(cur_jobid)+".shelve"
      outfile = shelve.open(name,'n')
      for elem,value in trainlines.iteritems():
        outfile[str(elem)]=value
      outfile.close()
      cjob = sub_script+" -t "+options.coordtype+" --lines='"+name+"' --outfilesindex='"+files_dict[cur_jobid][0]+"' --lane_tile='"+lane_tile+"' --crange='"+str(crange)+"' -p "+options.path
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
      cjob = sub_script+" -t "+options.coordtype+" --lines='"+name+"' --outfilesindex='"+files_dict[cur_jobid][0]+"' --lane_tile='"+lane_tile+"' --crange='"+str(crange)+"' -p "+options.path
      external_jobs.append((cur_jobid,subprocess.Popen(cjob,shell=True)))

def wait_external_jobs():
  global external_jobs,free_ids
  while len(external_jobs) > 0:
    cid,proc = external_jobs.pop(0)
    proc.wait()
  free_ids = None

if options.mock:
  options.verbose=True


##
#
# Begin Parsing user parameters
#
##

print "Starting training module "+timeString();

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

if not(os.path.isfile(options.SVMlight.split()[0])):
  print "Can't access training routine",options.SVMlight
  sys.exit()

if not(os.path.isfile(options.SVMlightTest.split()[0])):
  print "Can't access prediction routine",options.SVMlightTest
  sys.exit()

infiles = [] #array to training_seqs,training_seqs_r2
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

##
#
# end parsing user parameters
#
##


# FIXES infiles for indices


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
      infiles.insert(1,"INDEX")
    if options.readlength == None:
      readlength = max(fields[1]+options.indexlength+options.indexlength2,readlength)
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
  #options.indexlength=0
  infiles.insert(1,"INDEX")

# SECOND INDEX READ
if (options.indexlength2 > 0):
  ranges.append((lstart,ranges[-1][1],min(lstart+options.indexlength2,readlength)-lstart+ranges[-1][1],min(lstart+options.indexlength2,readlength)))
  lstart = min(lstart+options.indexlength2,readlength+options.indexlength)
  #options.indexlength2=0
  infiles.append("INDEX")


## MAKE SURE LAST ELEMENT IN RANGE IS READLENGTH
if (ranges[-1][3] != readlength):
  ranges[-1]=(ranges[-1][0],ranges[-1][1],ranges[-1][2],readlength)

missing = set()
assignments = {}
reassignments = {}
#print ranges
#print missing, assignments, reassignments
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
    reassignments[lstart] = (seqStart,'B')
    assignments[lstart] = (seqStart,'B')
    reassignments[fragStart] = (seqStart+1,'M')
    assignments[fragStart] = (seqStart+1,'M')
    lstart += 1
  while (lstart < fragStart):
    missing.add(lstart)
    reassignments[lstart] = (seqStart+1,'M')
    assignments[lstart] = (seqStart+1,'M')
    lstart += 1
  lstart = seqEnd
  if ((seqEnd-seqStart) == (fragEnd-fragStart)): lstart = fragEnd

if len(ranges) > 0: print ""
if lstart < readlength:
  missing.add(readlength-1)
  reassignments[readlength-1] = (lstart-1,'E')
  assignments[readlength-1] = (lstart-1,'E')
  reassignments[lstart-1] = (lstart-2,'M')
  assignments[lstart-1] = (lstart-2,'M')
  for elem in range(lstart,readlength-1):
    missing.add(elem)
    reassignments[elem] = (lstart-2,'M')
    assignments[elem] = (lstart-2,'M')

print "Readlength:",readlength
if len(reassignments) > 0:
  print ""
  print "For %d bases (%s) no training data is assigned, will do the following (re)assignments:"%(len(missing),",".join(map(lambda x:str(x+1),missing)))
  for pos,replace in reassignments.iteritems():
    print "  Pos:",pos+1,"Using base:",replace[0]+1,"Type:",replace[1]
  print ""

if(options.recalibration):
  myrecalname2cycle = {};
  cycle2set = {};
  

#for binning recalibration data
denominator= int(math.floor(float(readlength)/float(options.cgroup)));
smallset= int(math.floor( float(readlength)/float(denominator)));
largeset= int(math.ceil(  float(readlength)/float(denominator)));
nofsmallset=(denominator-(readlength % denominator));
noflargeset=(readlength%denominator);

indexbin   =1;
numberinbin=0;
    
for elem in range(readlength):
  # TRAINING DATA SET
  name = options.tmp+"/"+timestamp+"_SVMlight_train_cycle_"+str(elem+1)+".txt"
  train_names.append(name)
  remove_files.append(name)


  if(options.recalibration):
  
        
    cycle2set[elem+1]=indexbin;       
    name = options.tmp+"/"+timestamp+"_SVMlight_train_res_cycle_"+str(elem+1)+".txt";
    train_temp_names.append(name);
    remove_files.append(name);

    name = "SVMlight_recal_cycle_"+str( indexbin )+"_1";

    if(name in myrecalname2cycle):
      myrecalname2cycle[name].append("SVMlight_recal_cycle_"+str( (elem+1) )+"_1");
    else:
      myrecalname2cycle[name]= ["SVMlight_recal_cycle_"+str( (elem+1) )+"_1"];

    name = "SVMlight_recal_cycle_"+str( indexbin )+"_2";
    if(name in myrecalname2cycle):
      myrecalname2cycle[name].append("SVMlight_recal_cycle_"+str( (elem+1) )+"_2");
    else:
      myrecalname2cycle[name]= ["SVMlight_recal_cycle_"+str( (elem+1) )+"_2"];

    name = "SVMlight_recal_cycle_"+str( indexbin )+"_3";
    if(name in myrecalname2cycle):
      myrecalname2cycle[name].append("SVMlight_recal_cycle_"+str( (elem+1) )+"_3");
    else:
      myrecalname2cycle[name]= ["SVMlight_recal_cycle_"+str( (elem+1) )+"_3"];

    name = "SVMlight_recal_cycle_"+str( indexbin )+"_4";
    if(name in myrecalname2cycle):
      myrecalname2cycle[name].append("SVMlight_recal_cycle_"+str( (elem+1) )+"_4");
    else:
      myrecalname2cycle[name]= ["SVMlight_recal_cycle_"+str( (elem+1) )+"_4"];

    if(options.plotqual): #if we plot, we need the original .est file
      remove_files.append(options.tmp+"/"+timestamp+"_SVMlight_recal_cycle_"+str( (elem+1) )+"_1.estp");
      remove_files.append(options.tmp+"/"+timestamp+"_SVMlight_recal_cycle_"+str( (elem+1) )+"_2.estp");
      remove_files.append(options.tmp+"/"+timestamp+"_SVMlight_recal_cycle_"+str( (elem+1) )+"_3.estp");
      remove_files.append(options.tmp+"/"+timestamp+"_SVMlight_recal_cycle_"+str( (elem+1) )+"_4.estp");

    if( numberinbin == 0):
      name = "SVMlight_recal_cycle_"+str(indexbin)+"_1";
      recal_names1.append(name);
      recal_names.append(name);
      remove_files.append(options.tmp+"/"+timestamp+"_"+name+".par");
      remove_files.append(options.tmp+"/"+timestamp+"_"+name+".svm");
      remove_files.append(options.tmp+"/"+timestamp+"_"+name+".est");
      remove_files.append(options.tmp+"/"+timestamp+"_"+name+".out");
      remove_files.append(options.tmp+"/"+timestamp+"_"+name+".txt");

      name = "SVMlight_recal_cycle_"+str(indexbin)+"_2";
      recal_names2.append(name);
      recal_names.append(name);
 
      remove_files.append(options.tmp+"/"+timestamp+"_"+name+".est");
      remove_files.append(options.tmp+"/"+timestamp+"_"+name+".out");
      remove_files.append(options.tmp+"/"+timestamp+"_"+name+".txt");

      name = "SVMlight_recal_cycle_"+str(indexbin)+"_3";
      recal_names3.append(name);
      recal_names.append(name);
  
      remove_files.append(options.tmp+"/"+timestamp+"_"+name+".est");
      remove_files.append(options.tmp+"/"+timestamp+"_"+name+".out");
      remove_files.append(options.tmp+"/"+timestamp+"_"+name+".txt");

      name = "SVMlight_recal_cycle_"+str(indexbin)+"_4";
      recal_names4.append(name);
      recal_names.append(name);

      remove_files.append(options.tmp+"/"+timestamp+"_"+name+".est");
      remove_files.append(options.tmp+"/"+timestamp+"_"+name+".out");
      remove_files.append(options.tmp+"/"+timestamp+"_"+name+".txt");

    numberinbin+=1;
    if(indexbin <= nofsmallset):
        if( numberinbin == smallset):
            indexbin+=1;
            numberinbin=0;
    else:
        if( numberinbin == largeset):
            indexbin+=1;
            numberinbin=0;
    
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
test_core_files = {}
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

have_test = False
index_cycles = {}

# This reads the different training_seqs and training_seqs_r2
# each time it encounters a new lane_tile number
# even tiles are picked for training and odd ones for testing
# 
for ind,infilename in enumerate(infiles):
  if infilename != "INDEX":
    print "Reading",infilename
    crange = ranges[ind]
    infile = open(infilename)
    cur_lane_tile = None
    line_stack = {}
    count = 0
    count_tiles = 0
    for line in infile:
      if count > training_seq_limit: break 
      if (count % options.use == 0):
        fields = line.split()
        if len(fields) == nr_id_fields_seqs+1:
          tile = fields[1]
          while len(tile) < 4:
            tile="0"+tile
          if (cur_lane_tile == None):
            cur_lane_tile = fields[0]+"_"+tile
          elif (cur_lane_tile != fields[0]+"_"+tile): #upon finding a new lane_tile
            if (options.NoTestDataSet or (count_tiles % 2 == 0)) and (count_tiles % options.skip == 0):
              if options.verbose: print "Reading",cur_lane_tile,"for training"
              eval_lines(line_stack,cur_lane_tile,train_core_files,crange)
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
      if len(line_stack) > maximum_nr_clusters:
        if (options.NoTestDataSet or (count_tiles % 2 == 0)) and (count_tiles % options.skip == 0):
          if options.verbose: print "Reading",cur_lane_tile,"for training"
          eval_lines(line_stack,cur_lane_tile,train_core_files,crange)
        elif (count_tiles % options.skip == 0) and len(line_stack) > 0:
          if options.verbose: print "Reading",cur_lane_tile,"as test data set"
          eval_lines(line_stack,cur_lane_tile,test_core_files,crange)
          have_test = True
        line_stack = {}
        count_tiles+=1
      count+=1
    if len(line_stack) > 0:
      if (options.NoTestDataSet or (count_tiles % 2 == 0)) and (count_tiles % options.skip == 0):
        if options.verbose: print "Reading",cur_lane_tile,"for training"
        eval_lines(line_stack,cur_lane_tile,train_core_files,crange)
      elif (count_tiles % options.skip == 0):
        if options.verbose: print "Reading",cur_lane_tile,"as test data set"
        eval_lines(line_stack,cur_lane_tile,test_core_files,crange)
        have_test = True
    infile.close()
  else:
    crange = ranges[ind]
    #print crange
    for key,value in zip(range(crange[0],crange[3]),range(crange[1],crange[2])):
      index_cycles[key] = value
wait_external_jobs()

#print ranges
#print index_cycles,len(index_cycles)

print "Merging independent training files..."+timeString();
for core,values in train_core_files.iteritems():
  values=values[1:] #REMOVE FILE WITH LIST
  for ind,trainingfile in enumerate(train_names):
    if not options.mock and (ind in index_cycles) and os.path.isfile(values[index_cycles[ind]]): # COPYING DATA FOR INDEX CYCLES / TRIMMING BACK
      stripL4 = assignments[ind][1] == 'E' and assignments[index_cycles[ind]][1] == 'M'
      infile=open(values[index_cycles[ind]],'r')
      outfile=open(trainingfile,'a')
      for line in infile:
        if stripL4: outfile.write(" ".join(line.split()[:-4])+"\n")
        else: outfile.write(line)
      outfile.close()
      infile.close()
  for ind,trainingfile in enumerate(train_names):
    if not options.mock and os.path.isfile(values[ind]):
      infile=open(values[ind],'r')
      outfile=open(trainingfile,'a')
      for line in infile:
        outfile.write(line)
      outfile.close()
      infile.close()
      #SET INPUT FILE EMPTY TO SAVE DISK SPACE, WILL DELETE IT LATER...
      infile=open(values[ind],'w')
      infile.close()

print "Merging independent test files..."+timeString();
for core,values in test_core_files.iteritems():
  values=values[1:]
  for ind,testfile in enumerate(test_names):
    if not options.mock and (ind in index_cycles) and os.path.isfile(values[index_cycles[ind]]):  # COPYING DATA FOR INDEX CYCLES / TRIMMING BACK
      stripL4 = assignments[ind][1] == 'E' and assignments[index_cycles[ind]][1] == 'M'
      infile=open(values[index_cycles[ind]],'r')
      outfile=open(testfile,'a')
      for line in infile:
        if stripL4: outfile.write(" ".join(line.split()[:-4])+"\n")
        else: outfile.write(line)
      outfile.close()
      infile.close()
  for ind,testfile in enumerate(test_names):
    if not options.mock and os.path.isfile(values[ind]):
      infile=open(values[ind],'r')
      outfile=open(testfile,'a')
      for line in infile:
        outfile.write(line)
      outfile.close()
      infile.close()
      #SET INPUT FILE EMPTY TO SAVE DISK SPACE, WILL DELETE IT LATER...
      infile=open(values[ind],'w')
      infile.close()

print "Training models..."+timeString();
arrayOfJobsToSend=[];
for ind, file in enumerate(train_names):
  if ind not in missing:
    if options.verbose: print "Calling",options.SVMlight+" "+file+" "+outfilenames[ind]
    #handle_jobs(options.SVMlight+" "+file+" "+outfilenames[ind])
    arrayOfJobsToSend.append(options.SVMlight+" "+file+" "+outfilenames[ind]);
#wait_jobs()
handleListOfjobs(arrayOfJobsToSend);
print ""

if options.mock: print "Number of training file names:",len(train_names)
names = train_names
if have_test:
  names = test_names
  if options.mock: print "Number of testing file names:",len(test_names)

print "Creating test data..."+timeString();
arrayOfJobsToSend=[];
for ind,orgfilename in enumerate(names):
  if ind not in missing:
    if options.verbose: print "Calling",options.SVMlightTest+" "+orgfilename+" "+outfilenames[ind]+" > "+test_temp_names[ind]
    #handle_jobs(options.SVMlightTest+" "+orgfilename+" "+outfilenames[ind]+" > "+test_temp_names[ind])
    arrayOfJobsToSend.append(options.SVMlightTest+" "+orgfilename+" "+outfilenames[ind]+" > "+test_temp_names[ind]);
#wait_jobs()
handleListOfjobs(arrayOfJobsToSend);
print "done "+timeString();




#########################
#  QUALITY SCORE RECAL ##
#########################

if(options.recalibration):

  arrayOfJobsToSend=[];
  for ind,file in enumerate(train_names):
    if ind not in missing:
      if options.verbose: 
        print "Calling",options.SVMlightTest+" "+file+" "+outfilenames[ind]+" > "+train_temp_names[ind]
      #handle_jobs(options.SVMlightTest+" "+file+" "+outfilenames[ind]+" > "+train_temp_names[ind])
      arrayOfJobsToSend.append(options.SVMlightTest+" "+file+" "+outfilenames[ind]+" > "+train_temp_names[ind]);
  #wait_jobs()
  handleListOfjobs(arrayOfJobsToSend);


  print "Creating recalibration data..."+timeString();
  #if(not options.mock):
  for ind in range(readlength):
    if options.quick: 
      if options.verbose: 
        print "Calling "+("paste "+train_temp_names[ind]+" "+ train_names[ind]+" |awk '{if($1==1){print log( (1/(1+exp(-$3))+1/(1+exp(-$4))+1/(1+exp(-$5)))/(1/(1+exp(-$2))+1/(1+exp(-$3))+1/(1+exp(-$4))+1/(1+exp(-$5))) )\" \"($1==$6) }}' >> "+options.tmp+"/"+timestamp+"_"+recal_names1[cycle2set[(ind+1)]-1]+".txt");
        print "Calling "+("paste "+train_temp_names[ind]+" "+ train_names[ind]+" |awk '{if($1==2){print log( (1/(1+exp(-$2))+1/(1+exp(-$4))+1/(1+exp(-$5)))/(1/(1+exp(-$2))+1/(1+exp(-$3))+1/(1+exp(-$4))+1/(1+exp(-$5))) )\" \"($1==$6) }}' >> "+options.tmp+"/"+timestamp+"_"+recal_names2[cycle2set[(ind+1)]-1]+".txt");
        print "Calling "+("paste "+train_temp_names[ind]+" "+ train_names[ind]+" |awk '{if($1==3){print log( (1/(1+exp(-$2))+1/(1+exp(-$3))+1/(1+exp(-$5)))/(1/(1+exp(-$2))+1/(1+exp(-$3))+1/(1+exp(-$4))+1/(1+exp(-$5))) )\" \"($1==$6) }}' >> "+options.tmp+"/"+timestamp+"_"+recal_names3[cycle2set[(ind+1)]-1]+".txt");
        print "Calling "+("paste "+train_temp_names[ind]+" "+ train_names[ind]+" |awk '{if($1==4){print log( (1/(1+exp(-$2))+1/(1+exp(-$3))+1/(1+exp(-$4)))/(1/(1+exp(-$2))+1/(1+exp(-$3))+1/(1+exp(-$4))+1/(1+exp(-$5))) )\" \"($1==$6) }}' >> "+options.tmp+"/"+timestamp+"_"+recal_names4[cycle2set[(ind+1)]-1]+".txt");
      if(not options.mock):
        arrayOfJobsToSend=[];
        arrayOfJobsToSend.append("paste "+train_temp_names[ind]+" "+ train_names[ind]+" |awk '{if($1==1){print log( (1/(1+exp(-$3))+1/(1+exp(-$4))+1/(1+exp(-$5)))/(1/(1+exp(-$2))+1/(1+exp(-$3))+1/(1+exp(-$4))+1/(1+exp(-$5))) )\" \"($1==$6) }}' >> "+options.tmp+"/"+timestamp+"_"+recal_names1[cycle2set[(ind+1)]-1]+".txt");
        arrayOfJobsToSend.append("paste "+train_temp_names[ind]+" "+ train_names[ind]+" |awk '{if($1==2){print log( (1/(1+exp(-$2))+1/(1+exp(-$4))+1/(1+exp(-$5)))/(1/(1+exp(-$2))+1/(1+exp(-$3))+1/(1+exp(-$4))+1/(1+exp(-$5))) )\" \"($1==$6) }}' >> "+options.tmp+"/"+timestamp+"_"+recal_names2[cycle2set[(ind+1)]-1]+".txt");
        arrayOfJobsToSend.append("paste "+train_temp_names[ind]+" "+ train_names[ind]+" |awk '{if($1==3){print log( (1/(1+exp(-$2))+1/(1+exp(-$3))+1/(1+exp(-$5)))/(1/(1+exp(-$2))+1/(1+exp(-$3))+1/(1+exp(-$4))+1/(1+exp(-$5))) )\" \"($1==$6) }}' >> "+options.tmp+"/"+timestamp+"_"+recal_names3[cycle2set[(ind+1)]-1]+".txt");
        arrayOfJobsToSend.append("paste "+train_temp_names[ind]+" "+ train_names[ind]+" |awk '{if($1==4){print log( (1/(1+exp(-$2))+1/(1+exp(-$3))+1/(1+exp(-$4)))/(1/(1+exp(-$2))+1/(1+exp(-$3))+1/(1+exp(-$4))+1/(1+exp(-$5))) )\" \"($1==$6) }}' >> "+options.tmp+"/"+timestamp+"_"+recal_names4[cycle2set[(ind+1)]-1]+".txt");
        handleListOfjobs(arrayOfJobsToSend,8);

    else:
      if options.verbose:        
        print "Calling "+("paste "+train_temp_names[ind]+" "+ train_names[ind]+" |awk '{if($1==1){print ($1==$6)\" 1:\"$2\" 2:\"$3\" 3:\"$4\" 4:\"$5 }}' >> "+options.tmp+"/"+timestamp+"_"+recal_names1[cycle2set[(ind+1)]-1]+".txt");
        print "Calling "+("paste "+train_temp_names[ind]+" "+ train_names[ind]+" |awk '{if($1==2){print ($1==$6)\" 1:\"$2\" 2:\"$3\" 3:\"$4\" 4:\"$5 }}' >> "+options.tmp+"/"+timestamp+"_"+recal_names2[cycle2set[(ind+1)]-1]+".txt");
        print "Calling "+("paste "+train_temp_names[ind]+" "+ train_names[ind]+" |awk '{if($1==3){print ($1==$6)\" 1:\"$2\" 2:\"$3\" 3:\"$4\" 4:\"$5 }}' >> "+options.tmp+"/"+timestamp+"_"+recal_names3[cycle2set[(ind+1)]-1]+".txt");
        print "Calling "+("paste "+train_temp_names[ind]+" "+ train_names[ind]+" |awk '{if($1==4){print ($1==$6)\" 1:\"$2\" 2:\"$3\" 3:\"$4\" 4:\"$5 }}' >> "+options.tmp+"/"+timestamp+"_"+recal_names4[cycle2set[(ind+1)]-1]+".txt");

      if(not options.mock):
        arrayOfJobsToSend=[];
        arrayOfJobsToSend.append("paste "+train_temp_names[ind]+" "+ train_names[ind]+" |awk '{if($1==1){print ($1==$6)\" 1:\"$2\" 2:\"$3\" 3:\"$4\" 4:\"$5 }}' >> "+options.tmp+"/"+timestamp+"_"+recal_names1[cycle2set[(ind+1)]-1]+".txt");
        arrayOfJobsToSend.append("paste "+train_temp_names[ind]+" "+ train_names[ind]+" |awk '{if($1==2){print ($1==$6)\" 1:\"$2\" 2:\"$3\" 3:\"$4\" 4:\"$5 }}' >> "+options.tmp+"/"+timestamp+"_"+recal_names2[cycle2set[(ind+1)]-1]+".txt");
        arrayOfJobsToSend.append("paste "+train_temp_names[ind]+" "+ train_names[ind]+" |awk '{if($1==3){print ($1==$6)\" 1:\"$2\" 2:\"$3\" 3:\"$4\" 4:\"$5 }}' >> "+options.tmp+"/"+timestamp+"_"+recal_names3[cycle2set[(ind+1)]-1]+".txt");
        arrayOfJobsToSend.append("paste "+train_temp_names[ind]+" "+ train_names[ind]+" |awk '{if($1==4){print ($1==$6)\" 1:\"$2\" 2:\"$3\" 3:\"$4\" 4:\"$5 }}' >> "+options.tmp+"/"+timestamp+"_"+recal_names4[cycle2set[(ind+1)]-1]+".txt");
        handleListOfjobs(arrayOfJobsToSend,8);

  #handleListOfjobs(arrayOfJobsToSend,8);

  #arrayOfJobsToSend=[];
  if(not(options.NoTestDataSet)):
    for ind in range(readlength):
      if options.quick: 
        if options.verbose: 
          print "Calling "+("paste "+test_temp_names[ind]+" "+ test_names[ind]+" |awk '{if($1==1){print log( (1/(1+exp(-$3))+1/(1+exp(-$4))+1/(1+exp(-$5)))/(1/(1+exp(-$2))+1/(1+exp(-$3))+1/(1+exp(-$4))+1/(1+exp(-$5))) )\" \"($1==$6) }}' >> "+options.tmp+"/"+timestamp+"_"+recal_names1[cycle2set[(ind+1)]-1]+".txt");
          print "Calling "+("paste "+test_temp_names[ind]+" "+ test_names[ind]+" |awk '{if($1==2){print log( (1/(1+exp(-$2))+1/(1+exp(-$4))+1/(1+exp(-$5)))/(1/(1+exp(-$2))+1/(1+exp(-$3))+1/(1+exp(-$4))+1/(1+exp(-$5))) )\" \"($1==$6) }}' >> "+options.tmp+"/"+timestamp+"_"+recal_names2[cycle2set[(ind+1)]-1]+".txt");
          print "Calling "+("paste "+test_temp_names[ind]+" "+ test_names[ind]+" |awk '{if($1==3){print log( (1/(1+exp(-$2))+1/(1+exp(-$3))+1/(1+exp(-$5)))/(1/(1+exp(-$2))+1/(1+exp(-$3))+1/(1+exp(-$4))+1/(1+exp(-$5))) )\" \"($1==$6) }}' >> "+options.tmp+"/"+timestamp+"_"+recal_names3[cycle2set[(ind+1)]-1]+".txt");
          print "Calling "+("paste "+test_temp_names[ind]+" "+ test_names[ind]+" |awk '{if($1==4){print log( (1/(1+exp(-$2))+1/(1+exp(-$3))+1/(1+exp(-$4)))/(1/(1+exp(-$2))+1/(1+exp(-$3))+1/(1+exp(-$4))+1/(1+exp(-$5))) )\" \"($1==$6) }}' >> "+options.tmp+"/"+timestamp+"_"+recal_names4[cycle2set[(ind+1)]-1]+".txt");

        if(not options.mock):
          arrayOfJobsToSend=[];
          arrayOfJobsToSend.append("paste "+test_temp_names[ind]+" "+ test_names[ind]+" |awk '{if($1==1){print log( (1/(1+exp(-$3))+1/(1+exp(-$4))+1/(1+exp(-$5)))/(1/(1+exp(-$2))+1/(1+exp(-$3))+1/(1+exp(-$4))+1/(1+exp(-$5))) )\" \"($1==$6) }}' >> "+options.tmp+"/"+timestamp+"_"+recal_names1[cycle2set[(ind+1)]-1]+".txt");
          arrayOfJobsToSend.append("paste "+test_temp_names[ind]+" "+ test_names[ind]+" |awk '{if($1==2){print log( (1/(1+exp(-$2))+1/(1+exp(-$4))+1/(1+exp(-$5)))/(1/(1+exp(-$2))+1/(1+exp(-$3))+1/(1+exp(-$4))+1/(1+exp(-$5))) )\" \"($1==$6) }}' >> "+options.tmp+"/"+timestamp+"_"+recal_names2[cycle2set[(ind+1)]-1]+".txt");
          arrayOfJobsToSend.append("paste "+test_temp_names[ind]+" "+ test_names[ind]+" |awk '{if($1==3){print log( (1/(1+exp(-$2))+1/(1+exp(-$3))+1/(1+exp(-$5)))/(1/(1+exp(-$2))+1/(1+exp(-$3))+1/(1+exp(-$4))+1/(1+exp(-$5))) )\" \"($1==$6) }}' >> "+options.tmp+"/"+timestamp+"_"+recal_names3[cycle2set[(ind+1)]-1]+".txt");
          arrayOfJobsToSend.append("paste "+test_temp_names[ind]+" "+ test_names[ind]+" |awk '{if($1==4){print log( (1/(1+exp(-$2))+1/(1+exp(-$3))+1/(1+exp(-$4)))/(1/(1+exp(-$2))+1/(1+exp(-$3))+1/(1+exp(-$4))+1/(1+exp(-$5))) )\" \"($1==$6) }}' >> "+options.tmp+"/"+timestamp+"_"+recal_names4[cycle2set[(ind+1)]-1]+".txt");
          handleListOfjobs(arrayOfJobsToSend,8);

      else:
        if options.verbose: 
          print "Calling "+("paste "+test_temp_names[ind]+" "+ test_names[ind]+" |awk '{if($1==1){print ($1==$6)\" 1:\"$2\" 2:\"$3\" 3:\"$4\" 4:\"$5 }}' >> "+options.tmp+"/"+timestamp+"_"+recal_names1[cycle2set[(ind+1)]-1]+".txt");
          print "Calling "+("paste "+test_temp_names[ind]+" "+ test_names[ind]+" |awk '{if($1==2){print ($1==$6)\" 1:\"$2\" 2:\"$3\" 3:\"$4\" 4:\"$5 }}' >> "+options.tmp+"/"+timestamp+"_"+recal_names2[cycle2set[(ind+1)]-1]+".txt");
          print "Calling "+("paste "+test_temp_names[ind]+" "+ test_names[ind]+" |awk '{if($1==3){print ($1==$6)\" 1:\"$2\" 2:\"$3\" 3:\"$4\" 4:\"$5 }}' >> "+options.tmp+"/"+timestamp+"_"+recal_names3[cycle2set[(ind+1)]-1]+".txt");
          print "Calling "+("paste "+test_temp_names[ind]+" "+ test_names[ind]+" |awk '{if($1==4){print ($1==$6)\" 1:\"$2\" 2:\"$3\" 3:\"$4\" 4:\"$5 }}' >> "+options.tmp+"/"+timestamp+"_"+recal_names4[cycle2set[(ind+1)]-1]+".txt");

        if(not options.mock):
          arrayOfJobsToSend=[];
          arrayOfJobsToSend.append("paste "+test_temp_names[ind]+" "+ test_names[ind]+" |awk '{if($1==1){print ($1==$6)\" 1:\"$2\" 2:\"$3\" 3:\"$4\" 4:\"$5 }}' >> "+options.tmp+"/"+timestamp+"_"+recal_names1[cycle2set[(ind+1)]-1]+".txt");
          arrayOfJobsToSend.append("paste "+test_temp_names[ind]+" "+ test_names[ind]+" |awk '{if($1==2){print ($1==$6)\" 1:\"$2\" 2:\"$3\" 3:\"$4\" 4:\"$5 }}' >> "+options.tmp+"/"+timestamp+"_"+recal_names2[cycle2set[(ind+1)]-1]+".txt");
          arrayOfJobsToSend.append("paste "+test_temp_names[ind]+" "+ test_names[ind]+" |awk '{if($1==3){print ($1==$6)\" 1:\"$2\" 2:\"$3\" 3:\"$4\" 4:\"$5 }}' >> "+options.tmp+"/"+timestamp+"_"+recal_names3[cycle2set[(ind+1)]-1]+".txt");
          arrayOfJobsToSend.append("paste "+test_temp_names[ind]+" "+ test_names[ind]+" |awk '{if($1==4){print ($1==$6)\" 1:\"$2\" 2:\"$3\" 3:\"$4\" 4:\"$5 }}' >> "+options.tmp+"/"+timestamp+"_"+recal_names4[cycle2set[(ind+1)]-1]+".txt");
          handleListOfjobs(arrayOfJobsToSend,8);

  #handleListOfjobs(arrayOfJobsToSend,8);


  if not(options.quick): 
    print "Training recalibration models ..."+timeString();
    arrayOfJobsToSend=[];
    for ind,file in enumerate(recal_names):
      if options.verbose: 
        print "Calling",options.SVMlightRECAL+" "+options.tmp+"/"+timestamp+"_"+file+".txt "+options.tmp+"/"+timestamp+"_"+file+".svm";
        #handle_jobs(options.SVMlightRECAL+" "+options.tmp+"/"+timestamp+"_"+file+".txt "+options.tmp+"/"+file+".svm");
      arrayOfJobsToSend.append(options.SVMlightRECAL+" "+options.tmp+"/"+timestamp+"_"+file+".txt "+options.tmp+"/"+timestamp+"_"+file+".svm");
    handleListOfjobs(arrayOfJobsToSend,15);
    #wait_jobs()
    print "done"

  print "Testing observed error rate ..."+timeString();
  arrayOfJobsToSend=[];
  for ind,file in enumerate(recal_names):
    if options.quick: 
      if options.verbose: 
        print "Calling","cat "+options.tmp+"/"+timestamp+"_"+file+".txt  |sort -rn > "+options.tmp+"/"+timestamp+"_"+file+".out";
      arrayOfJobsToSend.append("cat "+options.tmp+"/"+timestamp+"_"+file+".txt  |sort -rn > "+options.tmp+"/"+timestamp+"_"+file+".out");

    else:
      if options.verbose: 
        print "Calling",options.SVMlightTestRECAL+" "+options.tmp+"/"+timestamp+"_"+file+".txt "+options.tmp+"/"+timestamp+"_"+file+".svm |sort -n > "+options.tmp+"/"+timestamp+"_"+file+".out";
      arrayOfJobsToSend.append(options.SVMlightTestRECAL+" "+options.tmp+"/"+timestamp+"_"+file+".txt "+options.tmp+"/"+timestamp+"_"+file+".svm  |sort -n > "+options.tmp+"/"+timestamp+"_"+file+".out");

  handleListOfjobs(arrayOfJobsToSend,15);
  #wait_jobs()
  print "done"

  print "Estimating error rate ..."+timeString();
  arrayOfJobsToSend=[];
  for ind,file in enumerate(recal_names):
    if options.verbose: 
      print "Calling",def_estimateError+"  "+options.tmp+"/"+timestamp+"_"+file+".out > "+options.tmp+"/"+timestamp+"_"+file+".est";
    #handle_jobs(def_estimateError+"  "+options.tmp+"/"+timestamp+"_"+file+".out > "+options.tmp+"/"+timestamp+"_"+file+".est");
    arrayOfJobsToSend.append(def_estimateError+"  "+options.tmp+"/"+timestamp+"_"+file+".out > "+options.tmp+"/"+timestamp+"_"+file+".est");
    #wait_jobs()
  handleListOfjobs(arrayOfJobsToSend,10);
  print "done"

  print "Linear regression ... "+timeString();
  arrayOfJobsToSend=[];
  for ind,file in enumerate(recal_names):
    if options.quick: 
      if options.verbose: 
        print "Calling",def_linearRegression+"  "+options.tmp+"/"+timestamp+"_"+file+".est 1 > "+options.tmp+"/"+timestamp+"_"+file+".par";
      arrayOfJobsToSend.append(def_linearRegression+"  "+options.tmp+"/"+timestamp+"_"+file+".est 1 > "+options.tmp+"/"+timestamp+"_"+file+".par");
    else:
      if options.verbose: 
        print "Calling",def_linearRegression+"  "+options.tmp+"/"+timestamp+"_"+file+".est 0 > "+options.tmp+"/"+timestamp+"_"+file+".par";
      arrayOfJobsToSend.append(def_linearRegression+"  "+options.tmp+"/"+timestamp+"_"+file+".est 0 > "+options.tmp+"/"+timestamp+"_"+file+".par");

  handleListOfjobs(arrayOfJobsToSend,5);
  print "done"+timeString();

  arrayOfJobsToSend=[];
  for ind,file in enumerate(recal_names):
    for filetomoveitto in myrecalname2cycle[file]:
      if options.verbose: 
        print "Calling","cp -f "+  options.tmp+"/"+timestamp+"_"+file+".par "+  options.outfiles+"/"+filetomoveitto+".par"
      arrayOfJobsToSend.append("cp -f "+  options.tmp+"/"+timestamp+"_"+file+".par "+  options.outfiles+"/"+filetomoveitto+".par");
      if(options.plotqual): #if we plot, we need the original .est file
        if options.verbose: 
          print "Calling","cp -f "+options.tmp+"/"+timestamp+"_"+file+".est "+options.tmp+"/"+timestamp+"_"+filetomoveitto+".estp";
        arrayOfJobsToSend.append("cp -f "+options.tmp+"/"+timestamp+"_"+file+".est "+options.tmp+"/"+timestamp+"_"+filetomoveitto+".estp");
      if options.verbose: 
        print "Calling","cp -f "+  options.tmp+"/"+timestamp+"_"+file+".svm "+  options.outfiles+"/"+filetomoveitto+".svm"
      arrayOfJobsToSend.append("cp -f "+  options.tmp+"/"+timestamp+"_"+file+".svm "+  options.outfiles+"/"+filetomoveitto+".svm");
  handleListOfjobs(arrayOfJobsToSend,5);




#########################
#   END RECALIBRATION  ##
#########################

if(options.recalibration):
  print "Evaluating training procedure"+timeString();

  order = assignments.items()
  order.sort()

#
#  for ind,orgfilename in enumerate(names):
#    if (ind not in missing) and (not options.mock):
#      infile = open(orgfilename)
#      classes = []
#      for line in infile:
#        cclass = int(line.split()[0])
#        classes.append(cclass)
#      infile.close()
#
#      miss = 0
#      infile = open(test_temp_names[ind])
#      for cclass in classes:
#        fields = infile.readline().split()
#        if len(fields) == 5:
#          if cclass > 0:
#            cclass2 = int(fields[0])
#            if cclass != cclass2: miss += 1
#          else:
#            print "Problems when extracting error."
#        else:
#          print "Unexpected line in",test_temp_names[ind],". May be class is missing in training data..."
#      infile.close()
#
#      if len(classes) > 0:
#        rate = miss/float(len(classes))
#        print "Misclassification rate for %d. model: %.4f%%"%(ind+1,rate*100.0)
#        del classes
  arrayErrorRates   =[];
  arrayMisclassRates=[];

  for ind,orgfilename in enumerate(names):
    if (ind not in missing) and (not options.mock):
      misclassrate=[[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]];

      awkjob="paste "+str(test_temp_names[ind])+" "+str(orgfilename)+"  |awk '{print $1==$6}' |sort -n |uniq -c";
      myjobawk=subprocess.Popen(awkjob,shell=True,stdout=subprocess.PIPE);
      if (myjobawk.wait() == 0):
        awkoutput = myjobawk.communicate()[0];
        arrayout=awkoutput.split();
        rate=float(arrayout[0])/(float(arrayout[0])+float(arrayout[2]));
        arrayErrorRates.append(rate);
        if(arrayout[1] == '0' and arrayout[3] == '1' ):
          print "Misclassification rate for %d. model: %.4f%%"%(ind+1,rate*100.0);
        else:
          print "WARNING: command "+awkjob+" failed";
      else:
        print "WARNING: command "+awkjob+" failed";

      for mynuc in range(1,5):
        awkjob="paste "+str(test_temp_names[ind])+" "+str(orgfilename)+"  |awk '{if($6=="+str(mynuc)+"){print $1}}' |sort -n |uniq -c";
        myjobawk=subprocess.Popen(awkjob,shell=True,stdout=subprocess.PIPE);
        if (myjobawk.wait() == 0):
          awkoutput = myjobawk.communicate()[0];
          arrayout=awkoutput.split();
          
          nuctoconsider=0;
          totalobsbases=0.0;
          while(nuctoconsider<len(arrayout)):
            totalobsbases+=float(arrayout[nuctoconsider])
            nuctoconsider+=2;
          nuctoconsider=1;
          while(nuctoconsider<=len(arrayout)):
            misclassrate[mynuc-1][int(arrayout[nuctoconsider])-1]=float(arrayout[ int(nuctoconsider)-1 ])/totalobsbases;
            nuctoconsider+=2;
          #if(arrayout[1] == '1' and arrayout[3] == '2' and arrayout[5] == '3' and arrayout[7] == '4' ):
          #  misclassrate[mynuc-1][0]=float(arrayout[0])/(float(arrayout[0])+float(arrayout[2])+float(arrayout[4])+float(arrayout[6]));
          #  misclassrate[mynuc-1][1]=float(arrayout[2])/(float(arrayout[0])+float(arrayout[2])+float(arrayout[4])+float(arrayout[6]));
          #  misclassrate[mynuc-1][2]=float(arrayout[4])/(float(arrayout[0])+float(arrayout[2])+float(arrayout[4])+float(arrayout[6]));
          #  misclassrate[mynuc-1][3]=float(arrayout[6])/(float(arrayout[0])+float(arrayout[2])+float(arrayout[4])+float(arrayout[6]));
          #else:
          #  print "WARNING: command "+awkjob+" failed";
        else:
          print "WARNING: command "+awkjob+" failed";

      arrayMisclassRates.append(misclassrate);
        
  if not options.mock:
    outfile = open(options.outfiles+"/SVMlight_models.index",'w')
    for cycle,(imodel,ctype) in order:
      outfile.write(str(cycle+1)+"\t"+outfilenames[imodel]+"\t"+ctype+"\t");#.join([str(0.0)]*49)+"\n") #4 x 4 x 3 +1 
      for mynuc1 in range(1,5):
        for mynuc2 in range(1,5):
          outfile.write("\t".join([str(0.0)]*2)+"\t"  + str(arrayMisclassRates[cycle][mynuc1-1][mynuc2-1])+"\t");
      outfile.write(str(arrayErrorRates[cycle])+"\n");           
    outfile.close()

else:

  print "Evaluating training procedure"+timeString();

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
          cor_class[cclass]=[1,1,1,1] # finite sample estimate - next one is in class, overstates the counters by one, the different in estimated error rate is negligible
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
        #print "1:\t"+str(cor_class[cclass][0])+" observations"
        #print "2:\t"+str(cor_class[cclass][1])+" observations"
        #print "3:\t"+str(cor_class[cclass][2])+" observations"
        #print "4:\t"+str(cor_class[cclass][3])+" observations"
        #print "Estimates based on:\t"+str(total_obs)+" observations"
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

  print "done"+timeString();

  if not options.mock:
    outfile = open(options.outfiles+"/SVMlight_models.index",'w')
    for cycle,(imodel,ctype) in order:
      outfile.write(str(cycle+1)+"\t"+outfilenames[imodel]+"\t"+ctype+"\t"+"\t".join(map(str,scaling_logit[imodel]))+"\n")
    outfile.close()

if(options.recalibration):
  if(options.plotqual):
    arrayOfJobsToSend=[];
    arrayOfJobsToSend.append("R CMD BATCH --no-save --no-restore '--args "+options.outfiles+"/ "+options.tmp+"/"+timestamp+" "+options.outfiles+"/' "+str(plot_recal)+" /dev/stdout");
    
    handleListOfjobs(arrayOfJobsToSend,15);
#sys.exit(1);
   
print "Index file available as",options.outfiles+"/SVMlight_models.index"+timeString();

if not(options.keep):
  print "Removing temporary files..."+timeString();
  for elem in remove_files:
    if os.path.isfile(elem):
      handle_jobs("rm "+elem)
else:
  arrayOfJobsToSend=[];
  if (not os.path.isdir(options.outfiles+"/Training/")) and (not options.mock):
    os.makedirs(options.outfiles+"/Training/")
  for elem in remove_files:
    if os.path.isfile(elem):
      #handle_jobs("mv "+elem+" "+options.outfiles+"/Training/")
      arrayOfJobsToSend.append("mv "+elem+" "+options.outfiles+"/Training/")
    #else:
    #  print "warning: file "+elem+" not found";
#wait_jobs()
  handleListOfjobs(arrayOfJobsToSend);


