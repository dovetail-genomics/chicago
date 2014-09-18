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

maxLBrownEst = 1.5e6
binsize = 20000
rmapfile = "/bi/home/spivakov/g/CHIC/Digest_Human_HindIII.bed"
baitmapfile = "/bi/home/spivakov/g/CHIC/Digest_Human_HindIII_baits.bed"
outfile = "/bi/home/spivakov/g/CHIC/NBaitsPerBin_out.txt"
picklefile = "/bi/home/spivakov/g/CHIC/NBaitsPerBin_out.pickle"
removeAdjacent = True

def usage():
  print "Usage: countNbaitsPerBin.py [--maxLBrownEst=%d] [--binsize=%d] [--removeAdjacent]\n\t[--rmapfile=%s]\n\t[--baitmapfile=%s]\n\t[--outfile=%s]\n\t[--picklefile=%s]\n" \
  % (maxLBrownEst, binsize, removeAdjacent, rmapfile, baitmapfile, outfile, picklefile)


try:
  opts, args = getopt.getopt(sys.argv[1:], 'l:b:jr:b:o:p:', \
['maxLBrownEst=', 'binsize=', 'removeAdjacent=', 'rmapfile=', 'baitmapfile=', 'outfile=', 'picklefile='])
except getopt.GetoptError:
  usage()
  sys.exit(120)
   
for opt, arg in opts: 
  if opt in ('--maxLBrownEst', '-l'):
    maxLBrownEst = long(arg)
  elif opt in ('--binsize', '-b'):
    binsize = long(arg)
  elif opt in ('--rmapfile', '-f'):
    rmapfile = arg
  elif opt in ('--baitmapfile', '-t'):
    baitmapfile = arg
  elif opt in ('--outfile', '-o'):
    outfile = arg
  elif opt in ('--picklefile', '-p'):
    picklefile = arg
  elif opt == '--removeAdjacent':
    removeAdjacent = str2bool(arg)
  elif opt == '-j':
    removeAdjacent = True


print "Using options:\nmaxLBrownEst=%d, binsize=%d removeAdjacent=%r\n\trmapfile=%s\n\tbaitmapfile=%s\n\toutfile=%s\n\tpicklefile=%s\n" \
% (maxLBrownEst, binsize, removeAdjacent, rmapfile, baitmapfile, outfile, picklefile)

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

print "Looping through other ends..."

n={}
for i in range(len(st)):
    
  n[id[i]] = [0]*int(maxLBrownEst/binsize)
  
  if removeAdjacent:
    iSt=i-2
  else:
    iSt=i-1

  for j in range(iSt,0,-1):
   if chr[j] != chr[i]:
    break
   d = st[i]+(end[i]-st[i])/2-(st[j]+(end[j]-st[j])/2)
   if d>=maxLBrownEst:
    break
   if id[j] in bid:
    n[id[i]][d/binsize] += 1
  
  if removeAdjacent:
    iSt=i+2
  else:
    iSt=i+1

  for j in range(iSt,len(st),1):
   if chr[j] != chr[i]:
    break
   d = st[j]+(end[j]-st[j])/2-(st[i]+(end[i]-st[i])/2)
   if d>=maxLBrownEst:
    break
   if id[j] in bid:
    n[id[i]][d/binsize] += 1
   
  if int(random.uniform(0,1000))==1: 
   print "%d " % i,

print "\nWriting out text file..."

of = open(outfile, "wt")
of.write("#\tmaxLBrownEst=%d\tbinsize=%d\trmapfile=%s\n" % (maxLBrownEst, binsize, rmapfile))
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
