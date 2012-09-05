#!/usr/bin/env python
# -*- coding: ASCII -*-

"""

:Author: Martin Kircher
:Contact: Martin.Kircher@eva.mpg.de
:Date: *23.02.2009

"""

import sys,os
from optparse import OptionParser
import gzip,shelve
import math
import Numeric as N
import array

#######################################
# SVM LIGHT FEATURE / CLASS ENCODING
#######################################

#CLASSES - #1: A #2: C #3: G #4: T
bases2class = {"A":"1","C":"2","G":"3","T":"4"}

#FEATURES - #1: SA #2: SC #3: SG #4: ST
#NOT FOR FIRST CYCLE IN EACH READ - #5: SBA #6: SBC #7: SBG #8: SBT
#NOT FOR LAST CYCLE IN EACH READ - #9: SAA #10: SAC #11: SAG #12: SAT

nr_id_fields = 4

lendian = 1 ## use -1 for big endian systems
nbases = 4 ## Number of bases/pictures
blocksize = 2
maxsignedintblock = None
maxintblock = None

def coord_round(x): return "%d"%(round(float(x)))
def coord_floor(x): return "%d"%(math.floor(abs(float(x))))
def coord_shift_round(x): return "%d"%(round(float(x)*10+1000))
coordconv = coord_shift_round

def to_int(s,signed=False):
  global lendian
  global maxintblock, maxsignedintblock, blocksize
  power = 1
  total = 0
  for elem in s[::lendian]:
    total += ord(elem)*power
    power *= 256
  if signed:
    if maxsignedintblock == None:
      maxsignedintblock=2**(blocksize*8-1)
      maxintblock=2**(blocksize*8)
    if total > maxsignedintblock:
      return total-maxintblock
    elif total == maxsignedintblock:
      return None
    else:
      return total
  else:
    return total

def read_cif(filename,clusteridx=None):
  global blocksize,nbases
  infile = open(filename,'rb')
  brange = range(nbases)
  try:
    CIFstr = infile.read(3)
    if CIFstr == "CIF":
      version = infile.read(1)
      #print "CIF version",ord(version)
      blocksize = ord(infile.read(1))
      #print "Block size",blocksize
      firstcycle = to_int(infile.read(2))
      nrcycles = to_int(infile.read(2))
      nrclusters = to_int(infile.read(4))
      data = infile.read()
      if len(data) == nbases * nrclusters * blocksize:
        #print "Data read fits with estimate."

        # EITHER ALL, OR ONLY CLUSTER SPECIFIED ARE EXTRACTED
        iterclusters = xrange(nrclusters)
        if clusteridx != None:
          iterclusters = clusteridx
        # DISTANCE BETWEEN THE DIFFERENT BASES
        jump = blocksize*nrclusters
        for cluster in iterclusters:
          pos = cluster*blocksize
          res = map(lambda x:to_int(data[x*jump+pos:x*jump+pos+blocksize],True),brange)
          if None in res:
            yield [0]*nbases;
          else:
            yield res
      else:
        print "File content does not reflect number of cycles and channels.",nbases * nrclusters,len(data)
    else:
      print "No valid CIF file"
  except IOError:
    print "Error reading",filename
  infile.close()
  raise StopIteration

def read_position_txt(filename):
  global coordconv
  infile = open(filename)
  for line in infile:
    fields = map(float,line.split())
    yield coordconv(fields[0]),coordconv(fields[1])
  infile.close()
  raise StopIteration

def read_locs(filename):
  global coordconv
  infile = open(filename,'rb')
  infile.read(8) # First 8 Byte are unused
  clusters = to_int(infile.read(4))
  binvalues = array.array('f')
  binvalues.read(infile, clusters * 2)
  data = N.array(binvalues, typecode=N.Float)
  data = N.reshape(data, ( clusters ,2 ))
  for x,y in data:
    yield coordconv(x),coordconv(y)
  infile.close()
  raise StopIteration

def read_clocs(filename):
  EXPECTED_CLOCS_VERSION = 1
  BLOCK_SIZE = 25
  IMAGE_WIDTH = 2048
  BLOCKS_PER_LINE = (IMAGE_WIDTH + BLOCK_SIZE - 1) / BLOCK_SIZE
  totalBlocks = 0
  currentBlock = 0
  currentBlockUnreadClusters = 0

  infile = open(filename,'rb')
  clocsVersion = ord(infile.read(1))
  totalBlocks = to_int(infile.read(4))
  currentBlockUnreadClusters = ord(infile.read(1))
  currentBlock+=1

  while (currentBlock < totalBlocks or ( currentBlock == totalBlocks and currentBlockUnreadClusters > 0)):
     while (currentBlockUnreadClusters == 0 and currentBlock < totalBlocks):
        currentBlockUnreadClusters = ord(infile.read(1))
        currentBlock += 1
     dx = ord(infile.read(1))
     dy = ord(infile.read(1))
     x = 10 * BLOCK_SIZE * ((currentBlock - 1) % BLOCKS_PER_LINE) + dx + 1000;
     y = 10 * BLOCK_SIZE * ((currentBlock - 1) / BLOCKS_PER_LINE) + dy + 1000;
     yield x,y
     currentBlockUnreadClusters -= 1
  infile.close()
  raise StopIteration


def get_cif_cycle(root_path,lane,tile,cycle,clusteridx=None):
  lanestr = "L%03d"%lane
  cyclestr = "C%d.1"%cycle
  if os.path.isdir(root_path+"/"+lanestr) and os.path.exists(root_path+"/"+lanestr+"/"+cyclestr+"/s_%d_%d.cif"%(lane,tile)):
    helper = []
    count = 0
    for cluster in read_cif(root_path+"/"+lanestr+"/"+cyclestr+"/s_%d_%d.cif"%(lane,tile),clusteridx):
      helper.append(cluster)
      count+=1
    return helper
  else:
    print "Error: Could not find",root_path+"/"+lanestr+"/"+cyclestr+"/s_%d_%d.cif"%(lane,tile)
  return None


def get_cif_cycles(root_path,lane,tile,start,end,clusteridx=None):
  lanes = map(lambda x: "L%03d"%(x+1),range(8))
  res = None
  lcyclecount = None
  if lane > 0 and lane <= len(lanes):
    cycles = map(lambda x: "C%d.1"%x,range(start,end+1))
    if os.path.isdir(root_path+"/"+lanes[lane-1]):
      res = []
      for cycle in cycles:
        if os.path.exists(root_path+"/"+lanes[lane-1]+"/"+cycle+"/s_%d_%d.cif"%(lane,tile)):
          helper = []
          count = 0
          for cluster in read_cif(root_path+"/"+lanes[lane-1]+"/"+cycle+"/s_%d_%d.cif"%(lane,tile),clusteridx):
            helper.append(cluster)
            count+=1
          if lcyclecount == None: lcyclecount=len(helper)
          if len(helper) != lcyclecount:
            print "Number of extracted cycles inconsistent!"
            return None
          res.append(helper)
        else:
          print "Error: Could not find",root_path+"/"+lanes[lane-1]+"/"+cycle
    else:
      print "Error: Could not find",root_path+"/"+lanes[lane-1]
      return None
  return res


parser = OptionParser(usage="usage: %prog [options]")
parser.add_option("-p", "--path", dest="path", help="Path to Firecrest folder (default .)", default=".")
parser.add_option("-t", "--coordianteType", dest="coordtype", help="Type of cluster coordinate transformation: round = Pipeline <= 1.4: round(x), floor = Pipeline 1.5.*: floor(abs(x)), shift_floor (default) = Pipeline 1.6.*: round(x*10+1000)",choices=["round","floor","shift_round"],default="shift_round")
parser.add_option("--lines",dest="lines", help="Shelve with lines to be processed", default=None)
parser.add_option("--outfilesindex",dest="outfilesindex",help="File listing output files",default=None)
parser.add_option("--lane_tile",dest="lane_tile",help="Lane/Tile string",default=None)
parser.add_option("--crange",dest="crange",help="crange string",default=None)
(options, args) = parser.parse_args()

if options.coordtype == "round":  coordconv = coord_round
elif options.coordtype == "floor": coordconv = coord_floor
else: coordconv = coord_shift_round


def conv_range(start,fields):
  res = ""
  for elem in fields:
    res += str(start)+":"+elem+" "
    start += 1
  return res

def eval_lines(trainlines,lane_tile,filestreams,crange):
  global coordconv
  lane = lane_tile.split("_")[0]
  tile = lane_tile.split("_")[1].lstrip("0")
  if os.path.exists(options.path+"s_"+lane_tile+"_int.txt.gz"): ## FIRECREST
    filename = options.path+"s_"+lane_tile+"_int.txt.gz"
    print "Reading",filename,"for",lane_tile,"with",len(trainlines),"sequences"
    infile = gzip.GzipFile(filename)

    have_coords = True
    for ckey,seq in trainlines.iteritems():
      if ckey.endswith("'IDX')"): have_coords = False
      break
    lcounter = 0
    for line in infile.read().splitlines():
      fields = line.split('\t')
      if have_coords: ckey = str(tuple(fields[:nr_id_fields]))
      else:
        ckey = str(tuple(lane,tile,str(lcounter),'IDX'))
        lcounter += 1
      if ckey in trainlines:
        seq = trainlines[ckey][:(crange[-1]-crange[0])]
        lseq = len(seq)
        if (lseq <= (len(fields)-nr_id_fields)):
          signal = fields[(crange[1]+nr_id_fields):(crange[2]+nr_id_fields)]
          for ind,base in enumerate(seq):
            notNull = False
            for cint in signal[ind].split():
              if int(float(cint)) != 0:
                notNull = True
                break
            if not notNull: continue

            if (ind == 0) and (base != "N") and (ind+1 < lseq) and (seq[ind+1] != "N"):
              filestreams[ind+crange[0]].write(bases2class[base]+" "+conv_range(1,signal[ind].split())+
                                conv_range(5,signal[ind+1].split())+"\n")
            elif (ind+1 == lseq) and (base != "N") and (ind > 0) and (seq[ind-1] != "N"):
              filestreams[ind+crange[0]].write(bases2class[base]+" "+conv_range(1,signal[ind].split())+
                                conv_range(5,signal[ind-1].split())+"\n")
            elif (base != "N") and (ind > 0) and (seq[ind-1] != "N") and (ind+1 < lseq) and (seq[ind+1] != "N"):
                filestreams[ind+crange[0]].write(bases2class[base]+" "+conv_range(1,signal[ind].split())+
                                conv_range(5,signal[ind-1].split())+conv_range(9,signal[ind+1].split())+"\n")
        else:
          print "Error: Sequences longer than intensities."
          sys.exit()
    infile.close()
  elif os.path.exists(options.path+"s_"+lane_tile+"_int.txt.p.gz"):  ## IPAR
    print "Reading intensity file for",lane_tile,"with",len(trainlines),"sequences"
    have_coords = True
    lcounter = 0
    for ckey,seq in trainlines.iteritems():
      if ckey.endswith("'IDX')"): have_coords = False
      break
    if have_coords:
      if os.path.exists(options.path+"s_"+lane_tile+"_pos.txt"):
        filename = options.path+"s_"+lane_tile+"_pos.txt"
        infile = read_position_txt(filename)
        pos = 0
        positions = []
        lane,tile = map(str,map(int,lane_tile.split("_")))
        for x,y in infile:
          pos += 1
          ckey = str((lane,tile,x,y))
          if ckey in trainlines:
            positions.append((pos,ckey))
        nrclusters = pos
      else:
        filename = options.path+"Firecrest/"+"L00"+lane+"/s_"+lane_tile+"_idx.txt"
        infile = open(filename)
        pos = 0
        positions = []
        for line in infile:
          fields = line.split()
          pos += 1
          ckey = str(tuple(fields[:nr_id_fields]))
          if ckey in trainlines:
            positions.append((pos,ckey))
        nrclusters = pos
        infile.close()
    else:
      positions = []
      for ckey,seq in trainlines.iteritems():
        positions.append((int(ckey.split(',')[2].strip().strip("\'")),ckey))


    filename = options.path+"s_"+lane_tile+"_int.txt.p.gz"
    #print "Reading",filename,"for",lane_tile,"with",len(trainlines),"sequences"
    infile = gzip.GzipFile(filename)
    filecontent = infile.read().splitlines()
    infile.close()
    maxpos = len(filecontent)
    #print "Lines",maxpos
    for pos,ckey in positions:
      signal = []
      seq = trainlines[ckey][:(crange[-1]-crange[0])]
      lseq = len(seq)
      for cycle in range(crange[1],crange[2]):
        cpos = pos+((nrclusters+1)*cycle)
        #print pos,cpos,cycle+1
        if cpos < maxpos:
          signal.append(filecontent[cpos])
        else:
          print "Error: Sequences longer than intensities."
          sys.exit()
      for ind,base in enumerate(seq):
        notNull = False
        for cint in signal[ind].split():
          if int(float(cint)) != 0:
            notNull = True
            break
        if not notNull: continue

        if (ind == 0) and (base != "N") and (ind+1 < lseq) and (seq[ind+1] != "N"):
          filestreams[ind+crange[0]].write(bases2class[base]+" "+conv_range(1,signal[ind].split())+
                            conv_range(5,signal[ind+1].split())+"\n")
        elif (ind+1 == lseq) and (base != "N") and (ind > 0) and (seq[ind-1] != "N"):
          filestreams[ind+crange[0]].write(bases2class[base]+" "+conv_range(1,signal[ind].split())+
                            conv_range(5,signal[ind-1].split())+"\n")
        elif (base != "N") and (ind > 0) and (seq[ind-1] != "N") and (ind+1 < lseq) and (seq[ind+1] != "N"):
          filestreams[ind+crange[0]].write(bases2class[base]+" "+conv_range(1,signal[ind].split())+
                            conv_range(5,signal[ind-1].split())+conv_range(9,signal[ind+1].split())+"\n")
    del filecontent
    del positions
  elif os.path.exists(options.path+"L00"+lane+"/C%d.1/s_"%crange[-1]+lane+"_"+tile+".cif"):  ## RTA CIF FILES
    print "Reading intensity files for",lane_tile,"with",len(trainlines),"sequences"
    have_coords = True
    lcounter = 0
    #print trainlines.keys()[:10]
    for ckey,seq in trainlines.iteritems():
      if ckey.endswith("'IDX')"): have_coords = False
      break
    if have_coords:
      if os.path.exists(options.path+"s_"+lane_tile+"_pos.txt"):
        filename = options.path+"s_"+lane_tile+"_pos.txt"
        infile = read_position_txt(filename)
      elif os.path.exists(options.path+"/L00%s/s_%s_%s.locs"%(lane,lane,tile.lstrip('0'))):
        filename = options.path+"/L00%s/s_%s_%s.locs"%(lane,lane,tile.lstrip('0'))
        infile = read_locs(filename)
      elif os.path.exists(options.path+"/L00%s/s_%s_%s.locs"%(lane,lane,tile.lstrip('0'))):
        filename = options.path+"/L00%s/s_%s_%s.clocs"%(lane,lane,tile.lstrip('0'))
        infile = read_clocs(filename)
      else:
        print "Error finding cluster position file for",lane_tile
        return None
      positions = []
      pos = 0
      #IDENTIFY CLUSTERS TO BE USED FOR TRAINING
      for x,y in infile:
        ckey = str((lane,tile,x,y))
        if ckey in trainlines:
          positions.append((pos,ckey))
        pos += 1
      nrclusters = pos
      #print len(trainlines),len(positions),nrclusters
    else:
      positions = []
      for ckey,seq in trainlines.iteritems():
        positions.append((int(ckey.split(',')[2].strip().strip("\'")),ckey))
      nrclusters = len(positions)

    lane,tile = int(lane),int(tile)
    intensities = [[]]
    for posind,(pos,ckey) in enumerate(positions):
      seq = trainlines[ckey][:(crange[-1]-crange[0])] ## Training reference sequence of a certain cluster
      intensities[0].append(seq)
    ccycle,ind = crange[1]+1,0
    hpositions = map(lambda (x,y):x,positions)
    while ccycle <= crange[2]:
      intensities.append(get_cif_cycle(options.path,lane,tile,ccycle,hpositions))
      if intensities[-1] == None:
        print "Error getting intensities from CIF files (lane %d, tile %d, cycle %d, #positions %d)."%(lane,tile,ccycle,len(positions))
        break

      if len(intensities) > 4: intensities.pop(1)
      if len(intensities) > 2:
        for posind,(pos,ckey) in enumerate(positions):
          seq = intensities[0][posind]
          lseq = len(seq)

          isNull = False
          for ipos in range(1,len(intensities)):
            for cint in intensities[ipos][posind]:
              if int(float(cint)) == 0: 
                isNull = True
                break
          if isNull: continue

          if (ind == 1) and (seq[0] != "N") and (seq[1] != "N"): # FIRST CYCLE IN READ
            filestreams[crange[0]].write(bases2class[seq[0]]+" "+conv_range(1,map(str,intensities[1][posind]))+
                              conv_range(5,map(str,intensities[2][posind]))+"\n")
          if len(intensities) == 4:
            if (ind+1 == lseq) and (seq[ind] != "N") and (ind > 0) and (seq[ind-1] != "N"): # LAST CYCLE IN READ
              filestreams[ind+crange[0]].write(bases2class[seq[ind]]+" "+conv_range(1,map(str,intensities[3][posind]))+
                                conv_range(5,map(str,intensities[2][posind]))+"\n")
            if (seq[ind] != "N") and (ind > 1) and (seq[ind-2] != "N") and (seq[ind-1] != "N"): # ANY OTHER CYCLE
              filestreams[ind-1+crange[0]].write(bases2class[seq[ind-1]]+" "+conv_range(1,map(str,intensities[2][posind]))+
                                conv_range(5,map(str,intensities[1][posind]))+conv_range(9,map(str,intensities[3][posind]))+"\n")
      ccycle += 1
      ind += 1
    del intensities
    del positions

  else:
    print "Can't find Firecrest/IPAR/CIF files:"
    print options.path+"s_"+lane_tile+"_int.txt.gz",os.path.exists(options.path+"/s_"+lane_tile+"_int.txt.gz")
    print options.path+"s_"+lane_tile+"_int.txt.p.gz",os.path.exists(options.path+"/s_"+lane_tile+"_int.txt.p.gz")
    print options.path+"s_"+lane_tile+"_pos.txt",os.path.exists(options.path+"/s_"+lane_tile+"_pos.txt")
    print options.path+"L00"+lane+"/C%d.1/s_"%crange[-1]+lane+"_"+tile+".cif",os.path.exists(options.path+"L00"+lane+"/C%d.1/s_"%crange[-1]+lane+"_"+tile+".cif")
    return None

if (options.lines != None) and (options.outfilesindex != None) and (options.lane_tile != None) and (options.crange != None) and os.path.exists(options.outfilesindex) and os.path.exists(options.lines):
  lines = shelve.open(options.lines,'r')
  infile = open(options.outfilesindex,'r')
  filestreams=[]
  for line in infile:
    filestreams.append(open(line.strip(),'a'))
  infile.close()
  eval_lines(lines,options.lane_tile,filestreams,eval(options.crange))
  for elem in filestreams:
    elem.close()
else:
  print "Incomplete arguments."
