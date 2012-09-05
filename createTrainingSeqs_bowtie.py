#!/usr/bin/env python
# -*- coding: ASCII -*-

"""

:Author: Martin Kircher
:Contact: Martin.Kircher@eva.mpg.de
:Date: *07.06.2010

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

max_length_bowtie=40
max_frags_created = 1 # DETERMINED FROM FIRST SEQUENCE OF FIRST AND SECOND READ


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

def eval_bowtie_output(infilename,outfilename,delta_bases,options):
  outfile = open(outfilename,'w')
  infile = open(infilename,'r')
  count = 0
  line = infile.readline()
  fields = line.rstrip('\n').split('\t')
  cid = ""
  cnum = -1
  ctotal = 1
  cdata = {}
  while line != "":
    if (cid == ""):
      hfields = fields[0].split("#")
      cid = "#".join(hfields[:-1])
      cnum,ctotal = map(int,hfields[-1].split('/'))
      cdata = {}
      cdata[cnum] = fields
    else:
      cdata[cnum] = fields

    nextline = infile.readline()
    nextfields = nextline.rstrip('\n').split('\t')
    try:
      hfields = nextfields[0].split("#")
      nextcid = "#".join(hfields[:-1])
      nextnum,nexttotal = map(int,hfields[-1].split('/'))
    except:
      nextcid = ""
      nextnum = -1
      nexttotal = 1

    if (nextcid != cid) and (cid != ""):
      keys = cdata.keys()
      cdata[cnum] = fields
      keys.sort()
      if (len(keys) == ctotal) and len(cdata[keys[0]]) > 7:
        lref = cdata[keys[0]][2]
        lstrand = cdata[keys[0]][1]
        cmmatch = cdata[keys[0]][7].split(',')
        if len(cmmatch) == 1 and cmmatch[0] == '': cmmatch = []
        total_mmatch = len(cmmatch)
        if lstrand == "-":
          hseq = list(rev_compl(cdata[keys[0]][4]))
          lpos = int(cdata[keys[0]][3])
        else:
          hseq = list(cdata[keys[0]][4])
          lpos = int(cdata[keys[0]][3])+len(cdata[keys[0]][4])
        for elem in cmmatch:
          mpos = int(elem.split(':')[0])
          if options.maskMM:
            hseq[mpos]='N'
          else:
            hseq[mpos]=elem.split(':')[1][0]
        seq = ''.join(hseq)
        consistent = True
        for key in keys[1:]:
          if (len(cdata[key]) <= 7) or (lref != cdata[key][2]) or (lstrand != cdata[key][1]) or ((lstrand == '+') and (lpos != int(cdata[key][3]))) or ((lstrand == '-') and (lpos != int(cdata[key][3])+len(cdata[key][4]))):
            #print 'FAIL',lref,cdata[key][2],lstrand,cdata[key][1],lpos,cdata[key][3],len(cdata[key][4])
            consistent = False
            break
          else:
            #print 'PASS',lref,cdata[key][2],lstrand,cdata[key][1],lpos,cdata[key][3],len(cdata[key][4])
            cmmatch = cdata[key][7].split(',')
            if len(cmmatch) == 1 and cmmatch[0] == '': cmmatch = []
            total_mmatch += len(cmmatch)
            if (lstrand == '-'):
              hseq = list(rev_compl(cdata[key][4]))
              lpos -= len(cdata[key][4])
            else:
              hseq = list(cdata[key][4])
              lpos += len(cdata[key][4])
            for elem in cmmatch:
              mpos = int(elem.split(':')[0])
              if options.maskMM:
                hseq[mpos]='N'
              else:
                hseq[mpos]=elem.split(':')[1][0]
            seq += ''.join(hseq)
        if consistent and (len(seq) == delta_bases):
          outfile.write("\t".join(cid.split(":")[1:])+"\t"+seq+"\n")
          count += 1
      cdata = {}
      cid = nextcid

    cnum = nextnum
    ctotal = nexttotal
    line = nextline
    fields = nextfields
  infile.close()
  outfile.close()
  return count



remove_files_r1 = []
remove_files_r2 = []
timestamp = str(time.time())
external_jobs = []
free_ids = None

def eval_file(filename,filetype,start,end,srange,start2=None,end2=None,srange2=None,isread2=False,index_seq=None):
  global external_jobs, free_ids
  global options, timestamp, remove_files_r1, remove_files_r2
  global sub_script
  global outfile,outfile_r2
  if free_ids == None:
    free_ids = range(options.cores)
    for elem in free_ids:
      remove_files_r1.append(options.tmp+"/"+timestamp+"_trainSeqs_"+str(elem)+"_r1.out")
      remove_files_r2.append(options.tmp+"/"+timestamp+"_trainSeqs_"+str(elem)+"_r2.out")

  if options.mock:
    cur_jobid = 0
    outfilename = remove_files_r1[cur_jobid]
    if isread2:
      outfilename = remove_files_r2[cur_jobid]

    param_str = "--infile='%s' --ftype=%s --expID=%s --numberN=%d --adapter='%s' --start=%d --end=%d --srange='%s' --outfile=%s --control_index=%s"%(filename,filetype,options.expID,options.numberN,options.adapter,start,end,srange,outfilename,options.control_index)
    if index_seq == "1":
      param_str += " --index_seq1"
    elif index_seq == "2":
      param_str += " --index_seq2"
    if start2 != None:
      param_str += " --start2=%d --end2=%d --srange2='%s' --outfile2=%s"%(start2,end2,srange2,remove_files_r2[cur_jobid])

    print "Calling",sub_script+" "+param_str
  else:
    # CHECK FOR FINISHED JOBS...
    njobs = []
    for cur_jobid,proc,ltype in external_jobs:
      if (proc.poll() == None):
        njobs.append((cur_jobid,proc,ltype))
      else:
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
    external_jobs = njobs

    # HAVE FREE JOB ID?...
    if (len(free_ids) == 0): # NO?
      cur_jobid,proc,ltype = external_jobs.pop()
      proc.wait()
      #COPY OUTPUT
      if ltype == 1 or ltype == 3:
        infile = open(remove_files_r1[cur_jobid],'r')
        outfile.write(infile.read())
        infile.close()
      if ltype == 2 or ltype == 3:
        infile = open(remove_files_r2[cur_jobid],'r')
        outfile_r2.write(infile.read())
        infile.close()
    else:
      cur_jobid = free_ids.pop(0)

    outfilename = remove_files_r1[cur_jobid]
    ctype = 1
    if isread2:
      outfilename = remove_files_r2[cur_jobid]
      ctype = 2

    param_str = "--infile='%s' --ftype=%s --expID=%s --numberN=%d --adapter='%s' --start=%d --end=%d --srange='%s' --outfile=%s --control_index=%s"%(filename,filetype,options.expID,options.numberN,options.adapter,start,end,srange,outfilename,options.control_index)
    if index_seq == "1":
      param_str += " --index_seq1"
    elif index_seq == "2":
      param_str += " --index_seq2"
    if start2 != None:
      ctype = 3
      param_str += " --start2=%d --end2=%d --srange2='%s' --outfile2=%s"%(start2,end2,srange2,remove_files_r2[cur_jobid])

    cjob = sub_script+" "+param_str
    external_jobs.append((cur_jobid,subprocess.Popen(cjob,shell=True),ctype))

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
group.add_option("-e", "--expID", dest="expID", help="Name of experiment to use for sequence names",default="BowtieTraining")
group.add_option("-l", "--lanes", dest="lanes", help="Lanes, example: 1,3,4-8 (default 4)",default="4")
group.add_option("-t", "--tiles", dest="tiles", help="Lanes, example: 1-13,100 (default 1-120)",default="1-120")
parser.add_option_group(group)

group = OptionGroup(parser, "Filter/Paired End","Excluding sequences and bases for training")
group.add_option("--start",dest="start",help="First base to include/first base in each read (comma separate)")
group.add_option("--end",dest="end",help="Last base to include/last base in each read (comma separate)")
group.add_option("--adapter", dest="adapter", help="Adapter sequence (default '')",default='')
group.add_option("--numberN", dest="numberN", help="Maximum number of missing bases to be accepted (default 3)",default=3,type="int")
group.add_option("--indexlength",dest="indexlength",help="Length of the index read (default 0/AUTO)",default=0,type="int")
group.add_option("--2nd_indexlength", dest="indexlength2", help="Length of a second index read following the reverse read (default 0)",type="int",default=0)
group.add_option("--control_index", dest="control_index", help="Sequence of index (only first index!) identifying control reads used for training/test data set (default '' = no filter)",default='')
parser.add_option_group(group)

group = OptionGroup(parser, "Bowtie","Parameters for Bowtie mapping")
group.add_option("-a", "--bowtie", dest="bowtie", help="Path to Bowtie binary (default "+def_bowtie_path+")",default=def_bowtie_path)
group.add_option("-c", "--cores", dest="cores", help="Number of CPU cores for Bowtie (default 1)",default=1,type="int")
group.add_option("-r", "--reference", dest="reference", help="Index prefix for reference genome file",default=def_bowtie_ref)
group.add_option("-m", "--mismatch", dest="mismatch", help="Number of mismatches allowed (default 3)",default=3,type="int")
group.add_option("--maxlength", dest="maxlength", help="Maximum sequence length to be submitted to bowtie (default %d)"%max_length_bowtie,default=max_length_bowtie,type="int")
group.add_option("--maskMM", dest="maskMM", help="Mask mismatches from Bowtie output in training data",action="store_true",default=False)
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

if options.mismatch > 3: options.mismatch = 3
max_length_bowtie = options.maxlength

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

print "Using bases",tstart+1,"to",tend,"for Bowtie mapping."
if (reads == 2):
  print "Using bases",tstart2+1,"to",tend2,"for Bowtie mapping of second read."

if (options.expID == None):
  print "Need an name for the experiment."
  sys.exit()

if (options.bowtie == None) or not(os.path.isfile(options.bowtie)):
  print "Need path to Bowtie binary."
  sys.exit()

if (options.reference == None) or not(os.path.isfile(options.reference+'.1.ebwt')):
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

print ""
outfilename = options.tmp+"/"+timestamp+"_bowtie_input.txt"
print "Creating temporary file with sequences:",outfilename
outfile = open(outfilename,'w')
if reads == 2:
  outfilename_r2 = options.tmp+"/"+timestamp+"_bowtie_input_2.txt"
  print "Creating temporary file with sequences:",outfilename_r2
  outfile_r2 = open(outfilename_r2,'w')
else: # SHOULD NEVER BE NECCESSARY
  outfile_r2 = outfile

had_training_sequences = False
for lane in lanes:
  slane = str(lane)
  files = []
  first_read = []
  second_read = []
  index_read2 = False
  ctype = None
  min_bcl,max_bcl = 1,None
  for tile in tiles:
    nr = str(tile)
    while len(nr) < 4:
      nr="0"+nr
    if os.path.isfile(options.path+"/s_"+slane+"_1_"+nr+"_qseq.txt") and (ctype == None or ctype == "qseq"):
      first_read.append(options.path+"/s_"+slane+"_1_"+nr+"_qseq.txt")
      ctype = "qseq"
    if os.path.isfile(options.path+"/s_"+slane+"_3_"+nr+"_qseq.txt") and (ctype == None or ctype == "qseq"):
      second_read.append(options.path+"/s_"+slane+"_3_"+nr+"_qseq.txt")
      index_read2 = True
      ctype = "qseq"
    elif os.path.isfile(options.path+"/s_"+slane+"_2_"+nr+"_qseq.txt") and (ctype == None or ctype == "qseq"):
      second_read.append(options.path+"/s_"+slane+"_2_"+nr+"_qseq.txt")
      ctype = "qseq"

    if os.path.isfile(options.path+"/L%03d/C1.1/s_%d_%d.bcl"%(lane,lane,tile)) and (ctype == None or ctype == "bcl"):
      if max_bcl == None:
        max_bcl = max(map(lambda x: int(x[1:-2]),filter(lambda x: x.startswith("C") and x.endswith(".1") and os.path.isdir(options.path+"/L%03d/"%lane+x), os.listdir(options.path+"/L%03d/"%lane))))
      first_read.append(options.path+";%d;%d;%d;%d"%(lane,tile,min_bcl,max_bcl))
      ctype = "bcl"

    if os.path.isfile(options.path+"/s_"+slane+"_"+nr+"_seq.txt") and os.path.isfile(options.path+"/s_"+slane+"_"+nr+"_prb.txt") and (ctype == None or ctype == "seq"):
      first_read.append(options.path+"/s_"+slane+"_"+nr+"_seq.txt,None") #options.path+"/s_"+slane+"_"+nr+"_prb.txt"
      ctype = "seq"
    if os.path.isfile(options.path+"/s_"+slane+"_"+nr+"_seq.txt") and (not os.path.isfile(options.path+"/s_"+slane+"_"+nr+"_prb.txt")) and (ctype == None or ctype == "seq"):
      first_read.append(options.path+"/s_"+slane+"_"+nr+"_seq.txt,None")
      ctype = "seq"

    if os.path.isfile(options.path+"/s_"+slane+"_"+nr+".fastq") and (ctype == None or ctype == "fastq"):
      first_read.append(options.path+"/s_"+slane+"_"+nr+".fastq")
      ctype = "fastq"
    if os.path.isfile(options.path+"/s_"+slane+"_1_"+nr+".fastq") and (ctype == None or ctype == "fastq"):
      first_read.append(options.path+"/s_"+lane+"_1_"+nr+".fastq")
      ctype = "fastq"
    if os.path.isfile(options.path+"/s_"+slane+"_3_"+nr+".fastq") and (ctype == None or ctype == "fastq"):
      second_read.append(options.path+"/s_"+slane+"_3_"+nr+".fastq")
      index_read2 = True
      ctype = "fastq"
    elif os.path.isfile(options.path+"/s_"+slane+"_2_"+nr+".fastq") and (ctype == None or ctype == "fastq"):
      second_read.append(options.path+"/s_"+slane+"_2_"+nr+".fastq")
      ctype = "fastq"

  start2ndread=None
  second_in_first = ((len(second_read) == 0) and (reads == 2))
  file_reader = None
  file_reader_type = "none"
  # DO SOME AUTODETECTION OF PARAMETERS (FIRST READ FILES)
  if len(first_read) > 0:
    had_training_sequences = True
    filename = first_read[0]
    if filename.endswith("qseq.txt"):
      file_reader = Reader.read_qseq_file
      file_reader_type = "qseq"
    elif filename.endswith("fastq.txt"):
      file_reader = Reader.read_fastq_file
      file_reader_type = "fastq"
    elif ";" in filename:
      file_reader = Reader.read_bcl_tile
      file_reader_type = "bcl"
    elif "," in filename:
      file_reader = Reader.read_seq_prb_file
      file_reader_type = "seq_prb"
    else:
      print "Unexpected content of file list."
      sys.exit()
    for name,seqorg,qualorg in file_reader(filename):
      name = ":".join(name)
      if index_read2:
         for iname,iseqorg,iqualorg in file_reader(filename.replace("s_"+slane+"_1_","s_"+slane+"_2_")):
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

      seq=seqorg[tstart:tend]
      if len(seq) != tend-tstart:
        print "Length ranges defined do not match length of sequences."
        sys.exit()
      frags=len(seq)//((len(seq)//max_length_bowtie)+1)
      outseqs = []
      while len(seq) > 0:
        if len(seq) <= max_length_bowtie:
          cseq = seq
          seq = ""
        else:
          cseq = seq[:frags]
          seq = seq[frags:]
        outseqs.append(cseq)
      if max_frags_created < len(outseqs): max_frags_created=len(outseqs)

      if second_in_first:
        seq = seqorg[tstart2:tend2]
        if len(seq) != tend2-tstart2:
          print "Length ranges defined do not match length of sequences."
          sys.exit()
        frags=len(seq)//((len(seq)//max_length_bowtie)+1)
        outseqs = []
        while len(seq) > 0:
          if len(seq) <= max_length_bowtie:
            cseq = seq
            seq = ""
          else:
            cseq = seq[:frags]
            seq = seq[frags:]
          outseqs.append(cseq)
        if max_frags_created < len(outseqs): max_frags_created=len(outseqs)

      break

  # DO SOME AUTODETECTION OF PARAMETERS (SECOND READ FILES)
  if start2ndread == None: start2ndread = tend
  if len(second_read) > 0:
    filename = second_read[0]
    if filename.endswith("qseq.txt"):
      file_reader = Reader.read_qseq_file
      file_reader_type = "qseq"
    elif filename.endswith("fastq.txt"):
      file_reader = Reader.read_fastq_file
      file_reader_type = "fastq"
    elif "," in filename:
      file_reader = Reader.read_seq_prb_file
      file_reader_type = "seq_prb"
    else:
      print "Unexpected content of file list."
      sys.exit()

    for name,seqorg,qualorg in file_reader(filename):
      if tend2 == None:
        tend2 = start2ndread+len(seqorg)
        tstart2 = start2ndread
        srange2=str(tstart2)+":"+str(tend2)
      seq=seqorg[max(tstart2-start2ndread,0):max(tend2-start2ndread,0)]
      if len(seq) != max(tend2-start2ndread,0)-max(tstart2-start2ndread,0):
        print "Length ranges defined do not match length of sequences."
        sys.exit()
      frags=len(seq)//((len(seq)//max_length_bowtie)+1)
      outseqs = []
      while len(seq) > 0:
        if len(seq) <= max_length_bowtie:
          cseq = seq
          seq = ""
        else:
          cseq = seq[:frags]
          seq = seq[frags:]
        outseqs.append(cseq)
      if max_frags_created < len(outseqs): max_frags_created=len(outseqs)

      break

  # PROCESS FIRST READ FILES EXTERNAL
  for filename in first_read:
    print "Reading",filename
    if second_in_first:
      eval_file(filename,file_reader_type,tstart,tend,srange1,start2=tstart2,end2=tend2,srange2=srange2)
    else:
      if index_read2: eval_file(filename,file_reader_type,tstart,tend,srange1,index_seq="2")
      else: eval_file(filename,file_reader_type,tstart,tend,srange1)

  # PROCESS SECOND READ FILES EXTERNAL
  for filename in second_read:
    print "Reading",filename
    if index_read2:
      eval_file(filename,file_reader_type,max(tstart2-start2ndread,0),max(tend2-start2ndread,0),srange2,isread2=True,index_seq="2")
    else:
      eval_file(filename,file_reader_type,max(tstart2-start2ndread,0),max(tend2-start2ndread,0),srange2,isread2=True,index_seq="1")

  wait_external_jobs()

outfile.close()
if reads == 2: outfile_r2.close()

if not had_training_sequences:
  print "Did not read a single input file. Check path to Bustard folder:", options.path
  sys.exit()

if options.mock: sys.exit()

print ""
print "Calling Bowtie..."
outfilename2 = options.tmp+"/"+timestamp+"_bowtie_output.txt"
cmdline = options.bowtie+" -v "+str(options.mismatch)+" -p "+str(options.cores)+" -f "+options.reference+" "+outfilename+" > "+outfilename2
print cmdline
proc = subprocess.Popen(cmdline,shell=True)
proc.wait()

if reads == 2:
  outfilename2_r2 = options.tmp+"/"+timestamp+"_bowtie_output_2.txt"
  cmdline = options.bowtie+" -v "+str(options.mismatch)+" -p "+str(options.cores)+" -f "+options.reference+" "+outfilename_r2+" > "+outfilename2_r2
  print cmdline
  proc = subprocess.Popen(cmdline,shell=True)
  proc.wait()

if not options.keep:
  print ""
  print "Removing temporary Bowtie input file:",outfilename
  proc = subprocess.Popen("rm -f "+outfilename,shell=True)
  proc.wait()
  if reads == 2:
    print "Removing temporary Bowtie input file:",outfilename_r2
    proc = subprocess.Popen("rm -f "+outfilename_r2,shell=True)
    proc.wait()

delta_bases = tend-tstart
if (max_frags_created > 1) and (options.cores > 1):
   print "Sorting Bowtie output by sequence ID"
   proc = subprocess.Popen("sort -T %s %s > %s_sort"%(options.tmp,outfilename2,outfilename2),shell=True)
   proc.wait()
   proc = subprocess.Popen("mv %s_sort %s"%(outfilename2,outfilename2),shell=True)
   proc.wait()
print "Iterating over Bowtie output file"
count = eval_bowtie_output(outfilename2,options.outfile,delta_bases,options)
if not options.keep:
  print "Removing temporary Bowtie output file",outfilename2
  proc = subprocess.Popen("rm -f "+outfilename2,shell=True)
print ""
print count,"training sequences available in",options.outfile

if reads == 2:
  print ""
  delta_bases = tend2-tstart2
  if (max_frags_created > 1) and (options.cores > 1):
    print "Sorting Bowtie output by sequence ID for second read"
    proc = subprocess.Popen("sort -T %s %s > %s_sort"%(options.tmp,outfilename2_r2,outfilename2_r2),shell=True)
    proc.wait()
    proc = subprocess.Popen("mv %s_sort %s"%(outfilename2_r2,outfilename2_r2),shell=True)
    proc.wait()
  print "Iterating over Bowtie output file"
  count = eval_bowtie_output(outfilename2_r2,options.outfile+"_r2",delta_bases,options)
  if not options.keep:
    print "Removing temporary Bowtie output file",outfilename2_r2
    proc = subprocess.Popen("rm -f "+outfilename2_r2,shell=True)
  print ""
  print count,"training sequences available in",options.outfile+"_r2"

