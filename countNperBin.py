import getopt
import sys
import random
import cPickle as p

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
rmapfile = "/bi/home/spivakov/g/CHIC/Digest_Human_HindIII.bed"
baitmapfile = "/bi/home/spivakov/g/CHIC/Digest_Human_HindIII_baits.bed"
outfile = "/bi/home/spivakov/g/CHIC/NperBin_Baits_out.txt"
picklefile = "/bi/home/spivakov/g/CHIC/NperBin_Baits_out.pickle"

def usage():
  print "Usage: countNperBin.py [--minFragLen=%d] [--maxFragLen=%d] [--maxLBrownEst=%d] [--binsize=%d] [--remove2b=%r] [--removeAdjacent=%r]\n\t[--rmapfile=%s]\n\t[--baitmapfile=%s]\n\t[--outfile=%s]\n\t[--picklefile=%s]\n" \
  % (minFragLen, maxFragLen, maxLBrownEst, binsize, removeB2B, removeAdjacent, rmapfile, baitmapfile, outfile, picklefile)


try:
  opts, args = getopt.getopt(sys.argv[1:], 'm:x:l:b:rja:b:f:t:o:p:', \
['minFragLen=', 'maxFragLen=', 'maxLBrownEst=', 'binsize=', \
'removeb2b=', 'removeAdjacent=', 'rmapfile=', 'baitmapfile=', 'outfile=', 'picklefile='])
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
  elif opt in ('--picklefile', '-p'):
    picklefile = arg


print "Using options:\n\tminFragLen=%d, maxFragLen=%d, maxLBrownEst=%d, binsize=%d, removeb2b=%r, removeAdjacent=%r\n\trmapfile=%s\n\tbaitmapfile=%s\n\toutfile=%s\n\tpicklefile=%s\n" \
% (minFragLen, maxFragLen, maxLBrownEst, binsize, removeB2B, removeAdjacent, rmapfile, baitmapfile, outfile, picklefile)

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

print "Looping through baits..."

n={}
for i in range(len(st)):
  if not id[i] in bid:
    continue
    
  n[id[i]] = [0]*int(maxLBrownEst/binsize)
  
  for j in range(i-1,0,-1):
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
   n[id[i]][d/binsize] += 1
  
  for j in range(i+1,len(st),1):
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
   n[id[i]][d/binsize] += 1
   
  if int(random.uniform(0,100))==1: 
   print "%d " % i,

print "\nWriting out text file..."

of = open(outfile, "wt")
of.write("#\tminFragLen=%d\tmaxFragLen=%d\tmaxLBrownEst=%d\tbinsize=%d\tremoveb2b=%r\tremoveAdjacent=%r\trmapfile=%s\tbaitmapfile=%s\n" % \
(minFragLen, maxFragLen, maxLBrownEst, binsize, removeB2B, removeAdjacent, rmapfile, baitmapfile))
for k in sorted(n.keys()):
 of.write("%d\t" % k)
 for i in range(len(n[k])): 
  of.write("%d" % n[k][i])
  if i!=len(n[k])-1:
   of.write("\t")
 of.write("\n")
of.close()

if picklefile!=None:
  print "Writing out pickle..."
  pf = open(picklefile, "wb")
  p.dump(n, pf)
  pf.close()
