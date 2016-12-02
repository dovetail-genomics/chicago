#! /usr/bin/env python
import getopt
import sys
import random
import fnmatch
import os
from ntpath import basename

class Unbuffered(object):
   def __init__(self, stream):
       self.stream = stream
   def write(self, data):
       self.stream.write(data)
       self.stream.flush()
   def __getattr__(self, attr):
       return getattr(self.stream, attr)

sys.stdout = Unbuffered(sys.stdout)

def str2bool(v):
  return v.lower() in ("yes", "true", "t", "1")

minFragLen=150
maxFragLen=40000
maxLBrownEst = 1.5e6
binsize = 20000
removeB2B=True
removeAdjacent=True
rmapfile = ""
baitmapfile = ""
outfile = ""
designDir= ""


def usage():
  print "Usage: python makeProxOEFile.py [--minFragLen=%d] [--maxFragLen=%d] [--maxLBrownEst=%d] [--binsize=%d] [--removeb2b=True] [--removeAdjacent=True]\n\t[--rmapfile=designDir/*.rmap]\n\t[--baitmapfile=designDir/*.baitmap]\n\t[--outfile=designDir/<rmapfileName>.poe]\n\t[--designDir=.]\n\nIf designDir is provided and contains a single <baitmapfile>.baitmap and <rmapfile>.rmap, these will be used unless explicitly specified.\nLikewise, the output file will be saved as designDir/proxOEfile.poe unless explicitly specified." \
  % (minFragLen, maxFragLen, maxLBrownEst, binsize)


try:
  opts, args = getopt.getopt(sys.argv[1:], 'm:x:l:b:rja:b:f:t:o:d:', \
['minFragLen=', 'maxFragLen=', 'maxLBrownEst=', 'binsize=', \
'removeb2b=', 'removeAdjacent=', 'rmapfile=', 'baitmapfile=', 'outfile=', 'designDir='])
except getopt.GetoptError:
  usage()
  sys.exit(120)
   
for opt, arg in opts: 
  if opt in ('--minFragLen', '-m'):
    minFragLen = long(arg)
  elif opt in ('--maxFragLen', '-x'):
    maxFragLen = long(arg)
  elif opt in ('--maxLBrownEst', '-l'):
    maxLBrownEst = long(arg)
  elif opt in ('--binsize', '-b'):
    binsize = long(arg)
  elif opt == '--removeb2b':
    removeB2B = str2bool(arg)
  elif opt == '--removeAdjacent':
    removeAdjacent = str2bool(arg)
  elif opt == '-b':
    removeB2B = True  
  elif opt == '-j':
    removeAdjacent = True
  elif opt in ('--rmapfile', '-f'):
    rmapfile = arg
  elif opt in ('--baitmapfile', '-t'):
    baitmapfile = arg
  elif opt in ('--outfile', '-o'):
    outfile = arg
  elif opt in ('--designDir', '-d'):
    designDir = arg


if designDir != "":
  if os.path.isdir(designDir):
    print "\nUsing designDir %s" % designDir;
  else:
    print "\nError: designDir does not exist.\n";
    usage()
    sys.exit(1)
else:
  designDir = "."

if baitmapfile == "":
  files = os.listdir(designDir)
  whichFiles = []
  for file in files:
    if fnmatch.fnmatch(file, '*.baitmap'):
        whichFiles.append(file)
  if len(whichFiles)==1:
    baitmapfile=os.path.join(designDir, whichFiles[0])
    print "Located baitmapfile %s in %s" % (whichFiles[0], designDir)
  else:
    print "\nError: could not unambiguously locate baitmapfile in designDir.\n"
    usage()
    sys.exit(1)

if rmapfile == "":
  files = os.listdir(designDir)
  whichFiles = []
  for file in files:
    if fnmatch.fnmatch(file, '*.rmap'):
        whichFiles.append(file)
  if len(whichFiles)==1:
    rmapfile=os.path.join(designDir, whichFiles[0])
    print "Located rmapfile %s in %s" % (whichFiles[0], designDir)
  else:
    print "\nError: could not unambiguously locate baitmapfile in designDir.\n"
    usage()
    sys.exit(1)


if outfile == "":
  rmapfileName = os.path.splitext(basename(rmapfile))[0]
  outfile = os.path.join(designDir, rmapfileName + ".poe")
  print "Output fill be saved as %s\n" % outfile


print "Using options:\n\tminFragLen=%d, maxFragLen=%d, maxLBrownEst=%d, binsize=%d, removeb2b=%r, removeAdjacent=%r\n\trmapfile=%s\n\tbaitmapfile=%s\n\toutfile=%s\n" \
% (minFragLen, maxFragLen, maxLBrownEst, binsize, removeB2B, removeAdjacent, rmapfile, baitmapfile, outfile)

a = open(rmapfile)
print "Reading rmap...."
chr = []
st = []
end = []
id = []
for line in a:
  line = line.strip()
  l = line.split("\t")
  chr.append(l[0])
  st.append(int(l[1]))
  end.append(int(l[2]))
  id.append(int(l[3]))
a.close()

b = open(baitmapfile)
print "Reading baitmap..."
bid = []
for line in b:
  line = line.strip()
  l = line.split("\t")
  bid.append(int(l[3]))
b.close()

bid = set(bid)

print "Sorting rmap..."

oldchr = chr
oldst = st
chr = [x for (x,y) in sorted(zip(oldchr, oldst))]
st = [y for (x,y) in sorted(zip(oldchr, oldst))]
end = [z for (x,y,z) in sorted(zip(oldchr, oldst, end))]
id = [z for (x,y,z) in sorted(zip(oldchr, oldst, id))]
del oldchr
del oldst

of = open(outfile, "wt")
of.write("#\tminFragLen=%d\tmaxFragLen=%d\tmaxLBrownEst=%d\tbinsize=%d\tremoveb2b=%r\tremoveAdjacent=%r\trmapfile=%s\tbaitmapfile=%s\n" % \
(minFragLen, maxFragLen, maxLBrownEst, binsize, removeB2B, removeAdjacent, rmapfile, baitmapfile))

print "Looping through baits..."

for i in xrange(len(st)):
  if not id[i] in bid:
    continue
    
  for j in xrange(i-1,0,-1):
   if chr[j] != chr[i]:
    break
   if removeB2B:
     if id[j] in bid:
       continue 
   if removeAdjacent:
     if j==i-1:
       continue
   if (end[j]-st[j])<minFragLen:
       continue
   if (end[j]-st[j])>maxFragLen:
       continue       
   d = st[i]+(end[i]-st[i])/2-(st[j]+(end[j]-st[j])/2)
   if d>=maxLBrownEst:
    break
   of.write("%d\t%d\t%d\n" % (id[i], id[j], d))
  
  for j in xrange(i+1,len(st),1):
   if chr[j] != chr[i]:
    break
   if removeB2B:
     if id[j] in bid:
       continue
   if removeAdjacent:
     if j==i+1:
       continue 
   if (end[j]-st[j])<minFragLen:
       continue
   if (end[j]-st[j])>maxFragLen:
       continue       
   
   d = st[j]+(end[j]-st[j])/2-(st[i]+(end[i]-st[i])/2)
   if d>=maxLBrownEst:
    break
   of.write("%d\t%d\t%d\n" % (id[i], id[j], d))
   
  if int(random.uniform(0,100))==1: 
    print "%d " % i

of.close()

