#!/usr/bin/env python
# -*- coding: ASCII -*-

"""

:Author: Martin Kircher
:Contact: Martin.Kircher@eva.mpg.de
:Date: *01.02.2010

"""

import sys,os
from optparse import OptionParser
from optparse import OptionGroup
import subprocess
import time

if (len(sys.argv[0]) > 0)  and os.path.isdir(os.getcwd()+"/"+os.path.dirname(sys.argv[0])):
  root_bin = os.path.dirname(os.getcwd()+"/"+sys.argv[0])
elif (len(sys.argv[0]) > 0)  and os.path.isdir(os.path.dirname(sys.argv[0])):
  root_bin = os.path.dirname(sys.argv[0])
else:
  root_bin = os.getcwd()
f=open(root_bin+'/params.py')
c=compile(f.read(), root_bin+'/params.py', "exec")
eval(c)

prog_seqs = def_extract_dataset
prog_training = def_training

prog_training_int = def_svm_train
prog_training_int_quick = def_log_train

prog_training_int2 = def_svm_test
prog_prediction_py = def_prediction_py
prog_prediction = def_svm_prediction


if not os.path.isfile(prog_seqs): 
  print "Cannot find the following program : ",prog_seqs
  sys.exit()
if not os.path.isfile(prog_training): 
  print "Cannot find the following program : ",prog_training
  sys.exit()
if not os.path.isfile(plot_recal): 
  print "Cannot find the following program : ",plot_recal
  sys.exit()
if not os.path.isfile(prog_training_int.split()[0]): 
  print "Cannot find the following program : ",prog_training_int.split()[0]
  sys.exit()
if not os.path.isfile(prog_training_int_quick.split()[0]): 
  print "Cannot find the following program : ",prog_training_int_quick.split()[0]
  sys.exit()
if not os.path.isfile(prog_training_int2.split()[0]): 
  print "Cannot find the following program : ",prog_training_int2.split()[0]
  sys.exit()
if not os.path.isfile(prog_prediction.split()[0]): 
  print "Cannot find the following program : ",prog_prediction.split()[0]
  sys.exit()
if not os.path.isfile(prog_prediction_py.split()[0]): 
  print "Cannot find the following program : ",prog_prediction.split()[0]
  sys.exit()
if not os.path.isfile( def_estimateError_bin ): 
  print "Cannot find the following program : ",def_estimateError_bin;
  sys.exit()
if not os.path.isfile( def_linearRegression_bin ):
  print "Cannot find the following program : ",def_linearRegression_bin;
  sys.exit()
if not os.path.isfile( def_svm_trainbinary_bin ):
  print "Cannot find the following program : ",def_svm_trainbinary_bin
  sys.exit()
if not os.path.isfile( def_svm_testbinary_bin ):
  print "Cannot find the following program : ",def_svm_testbinary_bin
  sys.exit()

#codecmd, outputcmd = subprocess.getstatusoutput()
jobcreated=subprocess.Popen("gzip -h > /dev/null",shell=True);
jobcreated.wait();
if(jobcreated.returncode != 0):
  print "Cannot find the following program : gzip"
  print "It must be installed on your system prior to running freeIbis";
  sys.exit()

  jobcreated=subprocess.Popen("bgzip -h > /dev/null",shell=True);
jobcreated.wait();
if(jobcreated.returncode != 0):
  print "Cannot find the following program : bgzip"
  print "Please download from http://samtools.sourceforge.net/tabix.shtml"
  sys.exit()





#########################
### USER PARAMETER ...
#########################

parser = OptionParser(usage="usage: %prog [options]")
group = OptionGroup(parser, "General","General options")
group.add_option("-c", "--cores", dest="cores", help="Number of CPU cores to be used (default 1)",default=1,type="int")
group.add_option("--corespred", dest="corespred", help="Number of CPU cores for prediction only (default :use same as -c)",default=None,type="int")

group.add_option("-e", "--expID", dest="expID", help="Name of experiment to use for sequence names")
group.add_option("-o", "--outpath", dest="outpath", help="Path for output files")
group.add_option("-b", "--bustard_path", dest="bustard", help="Path to bustard folder")
group.add_option("--NoFinishCheck", dest="NoFinishCheck", help="Do not check for s_[1-8]_finished.txt files in Bustard folder.",action="store_true",default=False)
group.add_option("--coordianteType", dest="coordtype", help="Type of cluster coordinate transformation: round = Pipeline <= 1.4: round(x), floor = Pipeline 1.5.*: floor(abs(x)), shift_round (default) = Pipeline 1.6.*: round(x*10+1000)",choices=["round","floor","shift_round"],default="shift_round")
group.add_option("--temp", dest="tmp", help="Path to temporary folder (default "+def_temp+")",default=def_temp)
group.add_option("--NoTestDataSet", dest="NoTestDataSet", help="Don't separate data for a test data set",action="store_true",default=False)
group.add_option("--mock", dest="mock", help="Do a mock run for testing",default=False,action="store_true")
group.add_option("--quick", dest="quick", help="Use logistic regression for prediction instead of SVMs, quicker but less accurate",action="store_true",default=False)
parser.add_option_group(group)

group = OptionGroup(parser, "Filter/Paired End","Excluding parts of sequences/Handle paired end runs")
group.add_option("--start",dest="start",help="First base to include/first base in each read (comma separated, not including index)")
group.add_option("--end",dest="end",help="Last base to include/last base in each read (comma separated, not including index)")
group.add_option("--numberN", dest="numberN", help="Maximum number of missing bases to be accepted (default 3)",default=3,type="int")
group.add_option("--indexlength", dest="indexlength", help="Length of index read (default 0)",type="int",default=0)
group.add_option("--2nd_indexlength", dest="indexlength2", help="Length of a second index read following the reverse read (default 0)",type="int",default=0)
group.add_option("--control_index", dest="control_index", help="Sequence of index (only first index!) identifying control reads used for training/test data set (default '' = no filter)",default='')
parser.add_option_group(group)

group = OptionGroup(parser, "Training","Parameters for training")
group.add_option("-l", "--lanes", dest="lanes", help="Lanes, example: 1,3,4-8 (default 4)",default="4")
group.add_option("-t", "--tiles", dest="tiles", help="Tiles, example: 1-13,100 (default 1-120)",default="1-120")
group.add_option("-a", "--align", dest="align", help="Path to SOAP/Bowtie binary",default=def_align_path)
group.add_option("-r", "--reference", dest="reference", help="Reference genome file",default=def_align_ref)
group.add_option("-m", "--mismatch", dest="mismatch", help="Number of mismatches allowed (default 5/bowtie max 3)",default=5,type="int")
group.add_option("--maskMM", dest="maskMM", help="Mask mismatches from mapping output in training data",action="store_true",default=False)
group.add_option("--nomask", dest="nomask", help="Do not mask variable positions on the control reads",action="store_true",default=False);
parser.add_option_group(group);

group = OptionGroup(parser, "Prediction","Parameters for prediction")
group.add_option("--NoPrediction", dest="NoPrediction", help="Don't do prediction, just do training",action="store_true",default=False)
group.add_option("--prediction_lanes", dest="prediction_lanes", help="Lanes, example: 1,3,4-8 (default all)",default=None)
group.add_option("--prediction_tiles", dest="prediction_tiles", help="Tiles, example: 1-13,100 (default all)",default=None)
group.add_option("--fastq", dest="fastq", help="Create output in fastq format",action="store_true",default=False)
group.add_option("--out4Q", dest="out4Q", help="Create output in 4Q format",action="store_true",default=False)
group.add_option("--predictionprog", dest="predprogram", help="Use this program for prediction",default=prog_prediction)
parser.add_option_group(group)

group = OptionGroup(parser, "Quality Scores","Parameters for quality scores")
group.add_option("--recalibration", dest="recalibration", help="Trigger quality score recalibration using control reads",action="store_true",default=False)
group.add_option("--plotqual", dest="plotqual", help="Plot quality score recalibration lines (requires R)",action="store_true",default=False)

group.add_option("--maxqual", dest="maxqual", help="Set maximum allowed quality score (PHRED scale default="+str(def_MAXQUALSCORE)+")",default=def_MAXQUALSCORE,type="int")
group.add_option("--cgroup", dest="cgroup", help="Group consecutive cycles together for quality score recalibration (useful for low amounts of control sequences ex: miseq)",type="int")

parser.add_option_group(group)


(options, args) = parser.parse_args()

jobs = []


def timeString():
    return str(time.strftime(" %d %b %Y %H:%M:%S", time.localtime())+" secs = "+str(time.time()));

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

if (options.outpath == None):
  print "Need path for output files"
  sys.exit()
elif os.path.isdir(options.outpath):
  pass
else:
  try: os.makedirs(options.outpath)
  except: pass
  if not os.path.isdir(options.outpath):
    print "Could not create",options.outpath
    sys.exit()
options.outpath=options.outpath.rstrip("/")
options.outpath+="/"

if not os.path.isdir(options.tmp):
  print "Need valid path for temporary files"
  sys.exit()
options.tmp=options.tmp.rstrip("/")
options.tmp+="/"

lanes = parse_rangestr(options.lanes)
if lanes == None:
  print "Need a valid range of lanes for training."
  sys.exit()
print "Using training lanes:",lanes
lanes = set(lanes)

if (options.bustard == None) or (not os.path.isdir(options.bustard)):
  print "Need valid path to bustard folder"
  sys.exit()
else:
  options.bustard=options.bustard.rstrip("/")
  options.bustard+="/"
  finished = True
  for lane in lanes:
    if not os.path.isfile(options.bustard+"s_"+str(lane)+"_finished.txt"):
      finished = False
      break
  if not finished and not options.NoFinishCheck:
    print "Need path to finished bustard folder. The given folder misses s_[1-8]_finished.txt files. If starting from RTA or unpacked SRF archive, please run with --NoFinishCheck."
    sys.exit()

if (options.expID == None):
  fields = options.bustard.split("/")
  if len(fields) >= 5:
    options.expID="_".join(fields[-5].split("_")[1:])
  if (options.expID == None):
    print "Extracting experiment name from bustard path failed."
    sys.exit()
  else:
    print "Extracted experiment name:",options.expID

if ( options.end   == None) :
  print "Plesae specify the end cycles of the reads using --end "
  sys.exit()

if ( options.start == None) :
  print "Plesae specify the start cycles of the reads using --start "
  sys.exit()
  


if ((options.end != None) and (options.start == None)):
  #no longer supported ! left for legacy 
  ends = []
  starts = [1]
  for elem in (options.end).split(","):
    try:
      ends.append(int(elem))
      starts.append(int(elem)+1)
    except: pass
  options.ends = ",".join(map(str,ends))
  if (options.indexlength != 0) and (len(starts) > 0):
    starts[1]+= options.indexlength
  else:
    print "Error fixing start/end parameters. Define completely."
    sys.exit()
  options.start = ",".join(map(str,starts[:-1]))
  print options.start,options.ends

  

if options.out4Q and options.fastq:
  print "Cannot specify both out4Q and fastq"
  sys.exit()
  

IsPairedEnd = (options.start != None) and ("," in options.start)

if (((options.indexlength != 0) or (options.indexlength2 != 0)) and (options.end == None)):
  print "If index flag is active, the end of the read has to be defined!"
  sys.exit()

lengthForForward   = int( (options.end).split(",")[0]) - int( (options.start).split(",")[0] )+1;
if(IsPairedEnd):
  lengthForReverse = int( (options.end).split(",")[1]) - int( (options.start).split(",")[1] ) +1 ;
else:
  lengthForReverse = 0;


options.control_index = options.control_index.strip()
if (options.control_index != '') and (len(options.control_index) != options.indexlength):
  print "Defined control index sequence does not match defined length of index read."
  sys.exit()

if (options.corespred == None) :
  options.corespred=options.cores;

if (options.predprogram != None) :
  prog_prediction=options.predprogram;


print "Extracting sequences for training..."+timeString();
params = ""
params += " --cores="+str(options.cores)
params += " --expID="+options.expID
params += " --path="+options.bustard
params += " --lanes="+options.lanes
params += " --tiles="+options.tiles
if options.start != None: params += " --start="+options.start
if options.end != None: params += " --end="+options.end
params += " --index="+str(options.indexlength)
params += " --2nd_index="+str(options.indexlength2)
params += " --control_index="+options.control_index
params += " -a "+options.align
params += " --reference="+options.reference
params += " --numberN="+str(options.numberN)
params += " --mismatch="+str(options.mismatch)
if options.maskMM: params += " --maskMM"
if options.nomask: params += " --nomask"

params += " --temp="+options.tmp
params += " --outfile="+options.outpath+"training_seqs"

timeBeforeTrainingSets= time.time();



##
modeldir=options.outpath+"Models/"
timeAfterSetsBefTraining= time.time(); #will get redefined if the index files does not exist

if not os.path.isfile(modeldir+"SVMlight_models.index"):
  if not os.path.isfile(options.outpath+"training_seqs") and (not IsPairedEnd or not os.path.isfile(options.outpath+"training_seqs_r2")):
    print "Calling",prog_seqs+params
    if not options.mock:
      proc = subprocess.Popen(prog_seqs+params,shell=True)
      proc.wait()
  else:
    print "Starting from earlier training data set."

  if (not(IsPairedEnd) and os.path.isfile(options.outpath+"training_seqs")) or (IsPairedEnd and os.path.isfile(options.outpath+"training_seqs") and os.path.isfile(options.outpath+"training_seqs_r2")) or options.mock:
    if options.cores > 1:
      print "Sorting training sequence output files."
      if not options.mock:
        if IsPairedEnd:
          proc2 = subprocess.Popen("sort -T %s -k 1,2 %straining_seqs_r2 > %straining_seqs_r2_sort"%(options.tmp,options.outpath,options.outpath),shell=True)
        proc = subprocess.Popen("sort -T %s -k 1,2 %straining_seqs > %straining_seqs_sort"%(options.tmp,options.outpath,options.outpath),shell=True)
        proc.wait()
        proc = subprocess.Popen("mv %straining_seqs_sort %straining_seqs"%(options.outpath,options.outpath),shell=True)
        proc.wait()
        if IsPairedEnd:
          proc2.wait()
          proc2 = subprocess.Popen("mv %straining_seqs_r2_sort %straining_seqs_r2"%(options.outpath,options.outpath),shell=True)
          proc2.wait()
    print "Training models based on training data..."

    if not os.path.isdir(modeldir):
      try: os.mkdir(modeldir)
      except:
        print "Could not create ",modeldir
        sys.exit()
    params = ""
    if IsPairedEnd: params += " --infile="+options.outpath+"training_seqs,"+options.outpath+"training_seqs_r2"
    else: params += " --infile="+options.outpath+"training_seqs"
    params += " --outfile="+modeldir
    firecrest_folder = os.path.abspath(options.bustard+'..').rstrip('/')+"/"
    params += " -t "+options.coordtype
    params += " --path="+firecrest_folder
    if options.quick:
      params += " --training='%s'"%prog_training_int_quick
    else:
      params += " --training='%s'"%prog_training_int
    params += " --prediction='%s'"%prog_training_int2
    params += " --temp="+options.tmp
    params += " --cores="+str(options.cores)
    params += " --indexlength="+str(options.indexlength)
    params += " --2nd_indexlength="+str(options.indexlength2)
    #params += " --keep" #used for debugging
    if options.recalibration:
      params += " --recalibration"
      if options.plotqual:
        params += " --plotqual"

      if options.cgroup:
        params += " --cgroup="+str(options.cgroup)

    if options.quick:
      params += " --quick"
      

    timeAfterSetsBefTraining= time.time();
    if options.NoTestDataSet:
      params += " --NoTestDataSet"
    print "on "+timeString()+" calling ",prog_training+params
    if not options.mock:
      proc = subprocess.Popen(prog_training+params,shell=True)
      proc.wait()
  else:
    print "Creating training data failed. Check parameters."
    sys.exit()

timeAfterTrainingBeforePred= time.time();



if ((not options.NoPrediction) and os.path.isfile(modeldir+"SVMlight_models.index")) or ((options.mock) and (not options.NoPrediction)):
  firecrest_folder = os.path.abspath(options.bustard+'..').rstrip('/')+"/"
  params = ""
  params += " --outfile="+options.outpath
  params += " --path="+firecrest_folder
  params += " --expID="+options.expID
  params += " --cores="+str(options.corespred)
  params += " --coordianteType="+options.coordtype
  params += " --prediction="+prog_prediction
  params += " --maxqual="+str(options.maxqual);
  params += " --forwardl="+str(lengthForForward);
  params += " --reversel="+str(lengthForReverse);
  params += " --index1l="+str(options.indexlength);
  params += " --index2l="+str(options.indexlength2);

  if options.prediction_lanes != None:
    params += " --lanes="+options.prediction_lanes
  if options.prediction_tiles != None:
    params += " --tiles="+options.prediction_tiles

  if options.recalibration:
    params += " --recalibration";
  if options.quick:
    params += " --quick"

  if options.out4Q:
    params += " --out4Q"
  if options.fastq:
    params += " --fastq"

  if options.mock:
    params += " --mock"

  print "on "+timeString()+" calling ",prog_prediction_py+params

  if not options.mock:
    proc = subprocess.Popen(prog_prediction_py+params,shell=True)
    proc.wait()
elif not os.path.isfile(modeldir+"SVMlight_models.index"):
  print "Training models failed. Check parameters."

timeAfterPred= time.time();
#if not options.mock:
#  print "Time to create training sets : "+str(timeAfterSetsBefTraining-timeBeforeTrainingSets)+" s";
#  print "Time to train                : "+str(timeAfterTrainingBeforePred-timeAfterSetsBefTraining)+" s";
#  print "Time to basecall             : "+str(timeAfterPred-timeAfterTrainingBeforePred)+" s";

