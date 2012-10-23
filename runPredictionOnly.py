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
from struct import * #for pack()

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


#Putting the PG tag in the bam file
commitversionibis="unknown";
bamtag="@PG\tID:freeIbis";
if os.path.exists(root_bin+"/.git/logs/HEAD"): #read the last github that was done if it exists
  fhGit = open ( root_bin+"/.git/logs/HEAD" );
  while 1:
      line=fhGit.readline();
      line=line.rstrip('\n');
      if(not(line) ):
          break;
      commitversionibis=line.split(" ")[1];
      
  fhGit.close();
else:
  print "No version number available";
if(commitversionibis != "unknown"):
  bamtag=bamtag+"\tVN:"+str(commitversionibis);
bamtag=bamtag+"\n";


#####################
# USER PARAMETER ...
#####################

parser = OptionParser(usage="usage: %prog [options]")
group = OptionGroup(parser, "Paths","Location of folders and files")
group.add_option("-o", "--outfile", dest="outfiles", help="Path for output files")
group.add_option("-p", "--path", dest="path", help="Path to Firecrest/IPAR/Intensities folder")
group.add_option("-e", "--expID", dest="expID", help="Name of experiment to use for sequence names")
group.add_option("-l", "--lanes", dest="lanes", help="Restrict to subset of lanes (default all)",default=None)
group.add_option("-s", "--tiles", dest="tiles", help="Restrict to subset of tiles (default all)",default=None)
group.add_option("-c", "--cores", dest="cores", help="Number of CPU processes to be used (default 1)",type="int",default=1)
group.add_option("-t", "--coordianteType", dest="coordtype", help="Type of cluster coordinate transformation: round = Pipeline <= 1.3: round(x), floor = Pipeline 1.5.*: floor(abs(x)), shift_floor (default) = Pipeline 1.6.*: round(x*10+1000)",choices=["round","floor","shift_round"],default="shift_round")
group.add_option("--prediction", dest="predictor", help="Prediction program used (default "+def_svm_prediction+")",default=def_svm_prediction)
group.add_option("--nooverwrite", dest="nooverwrite", help="Do not overwrite existing files",default=False,action="store_true")
group.add_option("--keeptiles", dest="keeptiles", help="Keep per tile files",default=False,action="store_true")
group.add_option("--mock", dest="mock", help="Do a mock run for testing",default=False,action="store_true")
group.add_option("--recalibration", dest="recalibration", help="Trigger quality score recalibration using control reads ",action="store_true",default=False)
group.add_option("--quick", dest="quick", help="Use linear regression for prediction, quicker but less accurate",action="store_true",default=False)
group.add_option("--verbose", dest="verbose", help="Produce a more verbose output",action="store_true",default=False)
group.add_option("--maxqual", dest="maxqual", help="Set maximum allowed quality score (PHRED scale)",default=def_MAXQUALSCORE,type="int")
parser.add_option_group(group)
group = OptionGroup(parser, "Output","Output format (default:BAM)")
group.add_option("--bcl", dest="outbcl", help="Produce BCL output",action="store_true",default=False)
group.add_option("--out4Q", dest="out4Q", help="Produce output in 4Q format",action="store_true",default=False)
group.add_option("--fastq", dest="fastq", help="Produce output in fastq format",action="store_true",default=False)

parser.add_option_group(group);
group = OptionGroup(parser, "Cycles","Details regarding the cycles");
group.add_option("--forwardl", dest="forwardl", help="Number of cycles for forward read",type="int",default=None);
group.add_option("--reversel", dest="reversel", help="Number of cycles for reverse read",type="int",default=None);
group.add_option("--index1l",  dest="index1l",  help="Number of cycles for first index", type="int",default=None);
group.add_option("--index2l",  dest="index2l",  help="Number of cycles for second index",type="int",default=None);


parser.add_option_group(group)


(options, args) = parser.parse_args()

if( (options.outbcl) ):
  print "BCL format is no longer supported";
  sys.exit(1)


if( (options.outbcl   and options.out4Q)  or 
    (options.outbcl   and options.fastq)  or 
    (options.out4Q    and options.fastq)  ):
  print "Please specify only one output format";
  sys.exit(1)

if( (options.outbcl or options.out4Q or options.fastq) ):
  options.outbam=False;
else:
  options.outbam=True;

if( options.forwardl == None):
  print "Please specify the number of cycles for forward read";
  sys.exit(1)

if( options.reversel == None):
  print "Please specify the number of cycles for reverse read";
  sys.exit(1)

if( options.index1l == None):
  print "Please specify the number of cycles for first index";
  sys.exit(1)

if( options.index2l == None):
  print "Please specify the number of cycles for first index";
  sys.exit(1)
  


def compareTileNumber(x, y):
  return int(x.split("_")[2])-int(y.split("_")[2]);


def timeString():
    return str(time.strftime(" on %d %b %Y %H:%M:%S", time.localtime())+" secs = "+str(time.time()));

def total_free_memory():
  proc = subprocess.Popen('top -n 1 | grep "^[S|M]"',shell=True,stdout=subprocess.PIPE)
  stdout_list = (proc.communicate())[0].split()
  helper = []
  for elem in stdout_list:
    if elem[0].isdigit() and elem[-1]=="k":
      helper.append(elem[:-1])
  if len(helper) == 8:
    return int(helper[2])+int(helper[7])
  else:
    return None

def total_memory():
  proc = subprocess.Popen('top -n 1 | grep "^[S|M]"',shell=True,stdout=subprocess.PIPE)
  stdout_list = (proc.communicate())[0].split()
  helper = []
  for elem in stdout_list:
    if elem[0].isdigit() and elem[-1]=="k":
      helper.append(elem[:-1])
  if len(helper) == 8:
    return int(helper[0])
  else:
    return None

def memory_pid(cid):
  proc = subprocess.Popen('ps -p %s -o rss'%cid,shell=True,stdout=subprocess.PIPE)
  try:
    return int((proc.communicate())[0].splitlines()[1])
  except:
    return None

jobs = []
job_count = 0
max_instance = 500000
printed_mem_warning = False

def check_memory():
  global jobs, max_instance, printed_mem_warning
  memory = []
  count_running = 0
  for elem in jobs:
    if None == elem.poll():
      count_running += 1
      mem = memory_pid(str(elem.pid))
      if mem != None and mem > max_instance:
        max_instance = mem
        #print max_instance

  mem = total_memory()
  free_mem = total_free_memory()
  if (mem == None) or ((max_instance*(count_running+1) < mem) and (max_instance <  free_mem)):
    return True
  else:
    if not printed_mem_warning:
      print "Not enough memory... Waiting for jobs to finish. Consider running with fewer cores next time."
      printed_mem_warning = True
    return False

def wait_jobs():
  global jobs
  while len(jobs) > 0:
    proc = jobs.pop(0)
    proc.wait()

def handle_jobs(cjob,params=None):
  global options, jobs, job_count

  if options.mock:
    if params == None:
      print cjob
    else:
      print " ".join(params)
    return None

  if params == None:
    if (len(jobs) < options.cores):
      #jobs.append(subprocess.Popen(cjob,shell=True))
      jobcreated=subprocess.Popen(cjob,shell=True);
      jobs.append(jobcreated);
      return jobcreated;
    else:
      iteration = 0
      while (len(jobs) >= options.cores):
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
      #jobs.append(subprocess.Popen(cjob,shell=True))
      jobcreated=subprocess.Popen(cjob,shell=True);
      jobs.append(jobcreated);
      return jobcreated;
  else:
    if (len(jobs) < options.cores) and check_memory():
      #jobcreated=subprocess.Popen(params,executable=cjob);
      #print "JOBS "+cjob+" PARAM "+str(" ".join(params));
      jobcreated=subprocess.Popen(str(" ".join(params)),shell=True);
      jobs.append(jobcreated);
      job_count += 1
      if job_count < options.cores:
        time.sleep(15) # GIVE PROCESS SOME TIME TO DEVELOP ITS MEMORY FOOTPRINT
      return jobcreated;
    else:
      iteration = 0
      while (len(jobs) >= options.cores) or (not check_memory()):
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

      #jobcreated=subprocess.Popen(params,executable=cjob);

      jobcreated=subprocess.Popen(str(" ".join(params)),shell=True);
      jobs.append(jobcreated);
      job_count += 1
      return jobcreated;



def handleListOfjobs(alljobs,params=False):
  global options

  unresolved=alljobs[:]; #copy
  numberOfIterations=0;
  while(len(unresolved) != 0):
    numberOfIterations+=1;
    arrayOfJobs=[];

    #launch jobs
    for jobToAdd in unresolved:
      if(options.verbose):
        print "launching "+str(jobToAdd);
      if(params):
        arrayOfJobs.append([jobToAdd,handle_jobs(jobToAdd[0],jobToAdd[1])]);
      else:
        arrayOfJobs.append([jobToAdd,handle_jobs(jobToAdd)]);

    if options.mock:
      return None;

    #wait for jobs
    while(1):
      allFinished=True;
      for toverify in arrayOfJobs:
        allFinished = allFinished and (toverify[1].poll() != None);

      if(allFinished):
        break;
      else:
        time.sleep(4);

    #all are done, check return code
    unresolved=[]; 
    for proc in arrayOfJobs:
      if(proc[1].returncode != 0): #wrong code
        print "WARNING: process "+str(proc[0])+" failed, will be relaunched";
        unresolved.append(proc[0]);
  
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

print "Starting prediction module "+timeString();

if options.lanes == None:
  lanes = range(1,9)
else:
  lanes = parse_rangestr(options.lanes)
str_lanes = set(map(str,lanes))

if options.tiles == None:
  tiles = None
else:
  tiles = set(map(str,parse_rangestr(options.tiles)))

if (options.path == None) or not(os.path.isdir(options.path)):
  print "Need path to Firecrest folder."
  sys.exit()
firecrest_folder=options.path.rstrip("/")+"/"

if (options.expID == None):
  fields = firecrest_folder.split("/")
  if len(fields) >= 4:
    options.expID="_".join(fields[-4].split("_")[1:])
  if (options.expID == None):
    print "Extracting experiment name from bustard path failed."
    sys.exit()
  else:
    print "Extracted experiment name:",options.expID

if (options.outfiles == None) or not(os.path.isdir(options.outfiles)):
  print "Need path for output files."
  sys.exit()
svmpath=options.outfiles.rstrip("/")+"/"
modeldir=svmpath+"Models/"

if not(os.path.isfile(options.predictor.split()[0])):
  print "Can't access prediction routine",options.predictor
  sys.exit()
prog_prediction=options.predictor


arrayOfFileInfix={};
for lane in lanes:
  arrayOfFileInfix[str(lane)]=[];

if os.path.isfile(modeldir+"SVMlight_models.index"):
  act_lanes = set()
  act_tiles = set()
  print "Predicting bases for each tile..."+timeString();

  #COMMON PARAMETERS
  params = [prog_prediction]
  params.extend(["-v","0","-t"])
  if options.coordtype == "round": 
    params.append("1")
  elif options.coordtype == "floor": 
    params.append("2")
  else: 
    params.append("3")
  if options.out4Q:
    params.extend(["-f","1"])
  if options.outbam:
    params.extend(["-f","2"])

  params.extend(["-e",options.expID]);
  params.extend(["-m",str(options.maxqual)]);
  params.extend(["-n",str(def_MINQUALSCORE)]);

  if(options.recalibration):
    params.extend(["-r"]);
  if(options.quick):
    params.extend(["-q"]);

  #print list(params);

  firecrest_files = os.listdir(firecrest_folder)
  firecrest_files.sort()
  found_IPAR_Firecrest = False

  arrayOfJobsToSend=[];


  for filename in firecrest_files: ## IPAR AND FIRECREST FILES...
    fields = filename.split("_")
    if (len(fields) > 2) and (fields[1] not in str_lanes): continue
    if (len(fields) > 2) and (tiles != None) and (fields[2].lstrip('0') not in tiles): continue
   

    if (filename.endswith("_int.txt.gz") or filename.endswith("_int.txt.p.gz")):
      found_IPAR_Firecrest = True
      myparams = list(params)
      currentLane=fields[1];
      currentTile=fields[2];
      myparams.append("-1 "+str(options.forwardl));
      myparams.append("-2 "+str(options.reversel));
      myparams.append("-3 "+str(options.index1l));
      myparams.append("-4 "+str(options.index2l));

      if filename.endswith("_int.txt.p.gz"):
        fields = filename.split("_")

        if os.path.isfile(firecrest_folder+"Firecrest/L00"+fields[1]+"/s_"+fields[1]+"_"+fields[2]+"_idx.txt"):
          myparams.extend(["-i",firecrest_folder+"Firecrest/L00"+fields[1]+"/s_"+fields[1]+"_"+fields[2]+"_idx.txt"])
          act_lanes.add(fields[1])
          act_tiles.add(fields[2])
     
        elif os.path.isfile(firecrest_folder+"s_"+fields[1]+"_"+fields[2]+"_pos.txt"):
          myparams.extend(["-i",firecrest_folder+"s_"+fields[1]+"_"+fields[2]+"_pos.txt"])
          act_lanes.add(fields[1])
          act_tiles.add(fields[2])
  
        else:
          print "Could not find index file for tile:",filename,"Skipping..."
          continue


      if(options.outbcl):
        myparams.append("-b");
        myparams.append("-l "+str(currentLane));
        myparams.append("-k "+str(currentTile));

        myparams.append(firecrest_folder+filename)
        myparams.append(modeldir+"SVMlight_models.index")
        myparams.append(svmpath);
        arrayOfJobsToSend.append([prog_prediction,myparams]);

      else:
        myparams.append(firecrest_folder+filename);
        myparams.append(modeldir+"SVMlight_models.index");

        if options.out4Q:
          houtname = svmpath+"_".join(filename.split("_")[:3])+".4Q";
        elif options.fastq:
          houtname = svmpath+"_".join(filename.split("_")[:3])+".fastq";
        else:
          houtname = svmpath+"_".join(filename.split("_")[:3])+".ubam";

        arrayOfFileInfix[fields[1]].append("_".join(filename.split("_")[:3]));

        if ( options.out4Q or options.fastq):
          if not (options.nooverwrite and os.path.isfile(houtname+".gz")):
            myparams.append(" /dev/stdout | gzip -c > "+houtname+".gz");
        else:
          if not (options.nooverwrite and os.path.isfile(houtname)):
            myparams.append(" "+houtname);

        arrayOfJobsToSend.append([prog_prediction,myparams]);



  if not found_IPAR_Firecrest: ## CIF FILES...
    for lane in lanes:

      if os.path.isdir(firecrest_folder+"L00%d/"%lane):

        cycles_folders = filter(lambda x: x.startswith('C') and x.endswith('.1') and os.path.isdir(firecrest_folder+"L00%d/%s"%(lane,x)), os.listdir(firecrest_folder+"L00%d/"%lane))

        if len(cycles_folders) > 0:
          cif_files = filter(lambda x: x.endswith('.cif'),os.listdir(firecrest_folder+"L00%d/%s/"%(lane,cycles_folders[0])))

          for cif in cif_files:
            fields = cif.split('.')[0].split("_")
            
            if (len(fields) > 2):
              myparams = list(params)
              tile = int(fields[2].lstrip('0'));
            
              if (tiles != None) and (str(tile) not in tiles): continue
              myparams.append("-1 "+str(options.forwardl));
              myparams.append("-2 "+str(options.reversel));
              myparams.append("-3 "+str(options.index1l));
              myparams.append("-4 "+str(options.index2l));

              if(options.outbcl):
                myparams.append("-b");
                myparams.append("-l "+str(lane));
                myparams.append("-k "+str(fields[2]));

              if os.path.exists(firecrest_folder+"/s_%d_%04d_pos.txt"%(lane,tile)):
                myparams.extend(["-i", firecrest_folder+"/s_%d_%04d_pos.txt"%(lane,tile), firecrest_folder+"L00%d/"%lane ])
                act_lanes.add("%s"%lane)
                act_tiles.add("%s"%tile)
              elif os.path.exists(firecrest_folder+"/L00%d/s_%d_%d.locs"%(lane,lane,tile)):
                myparams.extend(["-i", firecrest_folder+"/L00%d/s_%d_%d.locs"%(lane,lane,tile), firecrest_folder+"L00%d/"%lane ])
                act_lanes.add("%s"%lane)
                act_tiles.add("%s"%tile)
              elif os.path.exists(firecrest_folder+"/L00%d/s_%d_%d.clocs"%(lane,lane,tile)):
                myparams.extend(["-i", firecrest_folder+"/L00%d/s_%d_%d.clocs"%(lane,lane,tile), firecrest_folder+"L00%d/"%lane ])
                act_lanes.add("%s"%lane)
                act_tiles.add("%s"%tile)
              else:
                continue


              if(options.outbcl):
                myparams.append(modeldir+"SVMlight_models.index")
                myparams.append(svmpath);
                arrayOfJobsToSend.append([prog_prediction,myparams]);               
              else:
                myparams.append(modeldir+"SVMlight_models.index")
                if options.out4Q:
                  houtname = svmpath+"_".join(fields[:3])+".4Q"
                elif options.fastq:
                  houtname = svmpath+"_".join(fields[:3])+".fastq"
                else:
                  houtname = svmpath+"_".join(fields[:3])+".ubam"
                  
                arrayOfFileInfix[fields[1]].append("_".join(fields[:3]));

                #myparams.append(houtname)
                if ( options.out4Q or options.fastq):
                  if not (options.nooverwrite and os.path.isfile(houtname+".gz")):
                    myparams.append(" /dev/stdout | gzip -c > "+houtname+".gz");
                else:
                  if not (options.nooverwrite and os.path.isfile(houtname)):
                    myparams.append("  "+houtname);

                arrayOfJobsToSend.append([prog_prediction,myparams]);
                #handle_jobs(prog_prediction,myparams)
        else: 
          continue

  handleListOfjobs(arrayOfJobsToSend,True);   

  if(options.outbcl):
    print "Per lane raw sequence files created."+timeString();
    sys.exit(0);

  print "Prediction finished, merging files "+timeString();


  if not options.mock:
    for lane in lanes:
      #outfile = open(svmpath+"s_"+str(lane)+"_finished.txt",'w');       
      #outfile.close();
      if(options.outbam):
        if not(options.mock):
          fileHandleWrite = open ( svmpath+"s_"+str(lane)+"_sequence.ubam", 'wb' ) ;
          fileHandleWrite.write("BAM");
          fileHandleWrite.write('\1');        
          #for elem in range(8):
          #  fileHandleWrite.write('\0');


          fileHandleWrite.write(pack('i',int( len(bamtag) )));#writing length of bam tag
          fileHandleWrite.write(bamtag);
          fileHandleWrite.write(pack('i',int(0)));
          fileHandleWrite.close();
      else:
          outfile = open(svmpath+"s_"+str(lane)+"_sequence.txt.gz",'w');       
          outfile.close();

  #merging files for each tile
  for lane in lanes:
    if(options.outbam):
      cmdCombine="cat "+svmpath+"s_"+str(lane)+"_sequence.ubam ";
      
    for index in sorted(arrayOfFileInfix[str(lane)],cmp=compareTileNumber):  
      filetomove=svmpath+index+".";

      if options.out4Q:       
        filetomove=filetomove+"4Q.gz";
      elif options.fastq:
        filetomove=filetomove+"fastq.gz";
      else:
        filetomove=filetomove+"ubam";

      if( os.path.isfile(filetomove) or options.mock):
        if(options.outbam):
          cmdCombine+=" "+filetomove+" ";
          continue;
        while(True):
          cmdcat="cat "+filetomove+"  >> "+svmpath+"s_"+str(lane)+"_sequence.txt.gz";
          print "Calling cat "+filetomove+"  >> "+svmpath+"s_"+str(lane)+"_sequence.txt.gz";
          if(options.mock):
            break;
          myjobcat=subprocess.Popen(cmdcat,shell=True);
          if (myjobcat.wait() == 0):
            break;
          else:
            print "WARNING: command "+str(myjobcat)+" failed, restarting";
      else:
        print "ERROR: file "+str(filetomove)+" not found";
        
    if(options.outbam):
      cmdCombine+=" | bgzip -c  > "+svmpath+"s_"+str(lane)+"_sequence.bam ";
      handle_jobs(cmdCombine);


  wait_jobs();

  if not options.keeptiles:
    for lane in lanes:
      if(options.outbam):
        handle_jobs("rm "+svmpath+"s_"+str(lane)+"_sequence.ubam");

      for index in arrayOfFileInfix[str(lane)]:  
        filetomove=svmpath+index+".";
        if options.out4Q:       
          filetomove=filetomove+"4Q.gz";
        elif options.fastq:
          filetomove=filetomove+"fastq.gz";
        else:
          filetomove=filetomove+"ubam";

        handle_jobs("rm "+filetomove);

    wait_jobs()

  if not options.mock:
    for lane in lanes:
      outfile = open(svmpath+"s_"+str(lane)+"_finished.txt",'w');       
      outfile.close();
  print "Per lane raw sequence files created."+timeString();
else:
  print "Could not access training data. Check parameters."
  sys.exit()
