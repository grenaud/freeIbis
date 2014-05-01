#!/usr/bin/env python

"""

:Author: Gabriel
:Date: 2014

"""

import sys,os
from optparse import OptionParser
from optparse import OptionGroup
import string
import time
import subprocess

import createTrainingSeqs_reader as Reader


if (len(sys.argv[0]) > 0)  and os.path.isdir(os.getcwd()+"/"+os.path.dirname(sys.argv[0])):
  root_bin = os.path.dirname(os.getcwd()+"/"+sys.argv[0])
elif (len(sys.argv[0]) > 0)  and os.path.isdir(os.path.dirname(sys.argv[0])):
  root_bin = os.path.dirname(sys.argv[0])
else:
  root_bin = os.getcwd()
f=open(root_bin+'/params.py')
c=compile(f.read(), root_bin+'/params.py', "exec")
eval(c)
sub_script = def_ibis_path+'createTrainingSeqs_helper.py'

max_length_soap=60
max_frags_created = 1 # DETERMINED FROM FIRST SEQUENCE OF FIRST AND SECOND READ
jobs = []



def timeString():
    return str(time.strftime(" on %d %b %Y %H:%M:%S", time.localtime())+" secs = "+str(time.time()));

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

compSEQ = string.maketrans('CATG','GTAC')

def rev_compl(seq):
  nseq=seq.translate(compSEQ)
  nseq=list(nseq)
  nseq.reverse()
  return "".join(nseq)




remove_files_r1 = []
remove_files_r2 = []
timestamp = str(time.time())
external_jobs = []
free_ids = None




def handle_jobs(cjob,timeToSleep=30):
  global options
  global jobs
  if options.mock:
    print cjob
    return None
  if len(jobs) < options.cores:
    #jobs.append(subprocess.Popen(cjob,shell=True));
    print "launching "+str(cjob);
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

def wait_external_jobs():
  global external_jobs,free_ids
  global remove_files_r1,remove_files_r2
  global outfile,outfile_r2

  while len(external_jobs) > 0:
    cur_jobid,proc,ltype = external_jobs.pop(0)
    proc.wait()
    free_ids.append(cur_jobid)
    #COPY OUTPUT
    if ltype == 1 or ltype == 3:
      infile = open(remove_files_r1[cur_jobid],'r')
      outfile.write(infile.read())
      infile.close()
    if ltype == 2 or ltype == 3:
      infile = open(remove_files_r2[cur_jobid],'r')
      outfile_r2.write(infile.read())
      infile.close()

  #CLEAN ID LIST AND TEMP FILES
  for elem in remove_files_r1:
    if os.path.isfile(elem): os.remove(elem)
  for elem in remove_files_r2:
    if os.path.isfile(elem): os.remove(elem)

#########################
### USER PARAMETER ...
#########################

parser = OptionParser(usage="usage: %prog [options]")
group = OptionGroup(parser, "Bustard","General options describing sequencing data")
group.add_option("-p", "--path", dest="path", help="Path to bustard folder (default .)", default=".")
group.add_option("-e", "--expID", dest="expID", help="Name of experiment to use for sequence names",default="Training")
group.add_option("-l", "--lanes", dest="lanes", help="Lanes, example: 1,3,4-8 (default 4)",default="4")
group.add_option("-t", "--tiles", dest="tiles", help="Lanes, example: 1-13,100 (default 1-120)",default="1-120")
parser.add_option_group(group)

group = OptionGroup(parser, "Filter/Paired End","Excluding sequences and bases for training")
group.add_option("--start",dest="start",help="First base to include/first base in each read (comma separated)")
group.add_option("--end",dest="end",help="Last base to include/last base in each read (comma separated)")
group.add_option("--adapter", dest="adapter", help="Adapter sequence (default '')",default='')
group.add_option("--numberN", dest="numberN", help="not used",default=3,type="int")
group.add_option("--indexlength",dest="indexlength",help="Length of the index read (default 0/AUTO)",default=0,type="int")
group.add_option("--2nd_indexlength", dest="indexlength2", help="Length of a second index read following the reverse read (default 0)",type="int",default=0)
group.add_option("--control_index", dest="control_index", help="Sequence of index (only first index!) identifying control reads used for training/test data set (default '' = no filter)",default='')
parser.add_option_group(group)

group = OptionGroup(parser, "SOAP","Parameters for mapping")
group.add_option("-a", "--soap", dest="soap", help="not used",default=def_soap_path)
group.add_option("-c", "--cores", dest="cores", help="Number of CPU cores for SOAP (default 1)",default=1,type="int")
group.add_option("-r", "--reference", dest="reference", help="Reference genome file",default=def_soap_ref)
group.add_option("-m", "--mismatch", dest="mismatch", help="Number of mismatches allowed (default 5)",default=5,type="int")
group.add_option("--maskMM", dest="maskMM", help="not used",action="store_true",default=False)
group.add_option("--nomask", dest="nomask", help="Do not mask variable positions on the control reads",action="store_true",default=False)

parser.add_option_group(group)

group = OptionGroup(parser, "Files","Options for final/intermediate output")
group.add_option("-o", "--outfile", dest="outfile", help="Path for output file(s)")
group.add_option("--temp", dest="tmp", help="Path to temporary folder (default "+def_temp+")",default=def_temp)
group.add_option("--keep", dest="keep", help="Keep temporary files for debugging (default OFF)",action="store_true",default=False)
group.add_option("--mock", dest="mock", help="Do a mock run for testing",default=False,action="store_true")
parser.add_option_group(group)

(options, args) = parser.parse_args()

options.adapter = options.adapter.upper()
if len(options.adapter) == 0: options.adapter = None

##################################
#  BEGIN EVALUTATING PARAMETERS  #
##################################

##################################
# determine read starts and ends #
##################################

reads = 1        #number of reads
tstart = 1       # first cycle first read
tstart2 = None   # first cycle second read
tend = None      # last cycle first read
tend2 = None     # last cycle second read
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

print "Using bases",tstart+1,"to",tend,"for SOAP mapping."
if (reads == 2):
  print "Using bases",tstart2+1,"to",tend2,"for SOAP mapping of second read."

if (options.expID == None):
  print "Need an name for the experiment."
  sys.exit()

if (options.soap == None) or not(os.path.isfile(options.soap)):
  print "Need path to SOAP binary."
  sys.exit()

if (options.reference == None) or not(os.path.isfile(options.reference)):
  print "Need path to reference Genome."
  sys.exit()

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



##################################
#  END   EVALUTATING PARAMETERS  #
##################################


if True:
  #find positions
  # 2413  time ./bcl2phix -n 150 -f 1 -l 1 -t 1101 -i 151 -s CGATTCG /mnt/solexa/140417_M00518_0234_000000000-A78LB_RS_MAGE_1_2/Data/Intensities/BaseCalls/ /mnt/solexa/Genomes/phiX/control/whole_genome.fa        |wc -l 
  arrayOfJobsToSend=[];
  filesRead1=[];
  filesRead2=[];

  for lane in lanes:
    for tile in tiles:
      cmd = def_bcl2phix+" -n "+str(tend-tstart)+" -f "+str(tstart+1)+"  -l "+str(lane)+" -t "+str(tile)+" -d "+str(options.mismatch);
      if( len(options.control_index) != 0 ):
        cmd+= " -i "+str(tend-tstart+1);
        cmd+= " -s "+str(options.control_index);
      cmd+= " "+str(options.path)+"  "+str(options.reference) + " >  "+options.tmp+"/"+timestamp+"_"+str(lane)+"_"+str(tile)+"_phix.txt";
      #print cmd;
      filesRead1.append(options.tmp+"/"+timestamp+"_"+str(lane)+"_"+str(tile)+"_phix.txt");

      arrayOfJobsToSend.append(cmd);
      if reads == 2:
        cmd = def_bcl2phix+" -n "+str(tend2-tstart2)+" -f "+str(tstart2+1)+" -l "+str(lane)+" -t "+str(tile)+" -d "+str(options.mismatch);
        if( len(options.control_index) != 0 ):
          cmd+= " -i "+str(tend-tstart+1);
          cmd+= " -s "+str(options.control_index);
        cmd+= " "+str(options.path)+"  "+str(options.reference)  + " >  "+options.tmp+"/"+timestamp+"_"+str(lane)+"_"+str(tile)+"_phix_r2.txt";
        #print cmd;
        filesRead2.append(options.tmp+"/"+timestamp+"_"+str(lane)+"_"+str(tile)+"_phix_r2.txt");
        arrayOfJobsToSend.append(cmd);
  
  handleListOfjobs(arrayOfJobsToSend);
 
  #combine
  arrayOfJobsToSend=[];

  cmd = def_bcl2phixcombine+" -f "+str(float(def_percentageForMask)/float(100.0))+" -c "+str(options.outfile)+".covr -o "+str(options.outfile)+".usort ";
  if( options.nomask):
    cmd+= " -n ";
  else:
    cmd+=" -m "+str(options.outfile)+".mask";
  cmd += " "+str(options.reference) +" "+(" ".join(filesRead1));


  arrayOfJobsToSend.append(cmd);

  if reads == 2:
    cmd = def_bcl2phixcombine+" -f "+str(float(def_percentageForMask)/float(100.0))+" -c "+str(options.outfile)+"_r2.covr -o "+str(options.outfile)+"_r2.usort ";

    if( options.nomask):
      cmd+= " -n ";
    else:
      cmd+=" -m "+str(options.outfile)+"_r2.mask";
    cmd += " "+str(options.reference) +" "+(" ".join(filesRead2));

    arrayOfJobsToSend.append(cmd);

  handleListOfjobs(arrayOfJobsToSend);

#sorting
  arrayOfJobsToSend=[];

  cmd = "sort -T "+str(options.tmp)+" "+str(options.outfile)+".usort > "+str(options.outfile);
  arrayOfJobsToSend.append(cmd);

  if reads == 2:
    cmd = "sort -T "+str(options.tmp)+" "+str(options.outfile)+"_r2.usort > "+str(options.outfile)+"_r2 ";  
    arrayOfJobsToSend.append(cmd);
  
  handleListOfjobs(arrayOfJobsToSend);

  handle_jobs("rm -f "+str(options.outfile)+".usort");
  handle_jobs("rm -f "+str(options.outfile)+"_r2.usort");


  #remove temp file
  if not options.keep:  
    print "Removing temporary files..."+timeString();

    for elem in filesRead1:
      if os.path.isfile(elem):
        handle_jobs("rm -f "+elem)

    if reads == 2:
      for elem in filesRead2:
        if os.path.isfile(elem):
          handle_jobs("rm -f "+elem)




     


print "Done creating training sets "+timeString();
