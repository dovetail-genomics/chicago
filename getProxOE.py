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

minSize=150
maxSize=40000
maxL = 1.5e6
bin = 20000
removeB2B=True
removeAdjacent=True
rmapfile = "/bi/home/spivakov/g/CHIC/Digest_Human_HindIII.bed"
bfile = "/bi/home/spivakov/g/CHIC/Digest_Human_HindIII_baits.bed"
outfile = "/bi/home/spivakov/g/CHIC/proxOE_out.txt"
#picklefile = "/bi/home/spivakov/g/CHIC/proxOE_out.pickle"

def usage():
  print "Usage: getProxOE.py [--minsize=%d] [--maxsize=%d] [--maxl=%d] [--binsize=%d] [--removeb2b=%r] [--removeAdjacent=%r]\n\t[--rmapfile=%s]\n\t[--bfile=%s]\n\t[--outfile=%s]\n" \
  % (minSize, maxSize, maxL, bin, removeB2B, removeAdjacent, rmapfile, bfile, outfile)


try:
  opts, args = getopt.getopt(sys.argv[1:], 'm:x:l:b:rja:b:f:t:o:', \
['minsize=', 'maxsize=', 'maxl=', 'binsize=', \
'removeb2b=', 'removeAdjacent=', 'rmapfile=', 'baitmapfile=', 'outfile='])
except getopt.GetoptError:
  usage()
  sys.exit(120)
   
for opt, arg in opts: 
  if opt in ('--minsize', '-m'):
    minSize = long(arg)
  elif opt in ('--maxsize', '-x'):
    maxSize = long(arg)
  elif opt in ('--maxl', '-l'):
    maxL = long(arg)
  elif opt in ('--binsize', '-b'):
    bin = long(arg)
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
    bfile = arg
  elif opt in ('--outfile', '-o'):
    outfile = arg

print "Using options:\n\tminsize=%d, maxsize=%d, maxl=%d, bin=%d, removeb2b=%r, removeAdjacent=%r\n\trmapfile=%s\n\tbfile=%s\n\toutfile=%s\n" \
% (minSize, maxSize, maxL, bin, removeB2B, removeAdjacent, rmapfile, bfile, outfile)

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

b = open(bfile)
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
of.write("#\tminsize=%d\tmaxsize=%d\tmaxl=%d\tbinsize=%d\tremoveb2b=%r\tremoveAdjacent=%r\trmapfile=%s\tbaitmapfile=%s\n" % \
(minSize, maxSize, maxL, bin, removeB2B, removeAdjacent, rmapfile, bfile))

print "Looping through baits..."

#n={}
for i in range(len(st)):
  if not id[i] in bid:
    continue
    
  #n[id[i]] = [0]*int(maxL/bin)
  
  for j in range(i-1,0,-1):
   if chr[j] != chr[i]:
    break
   if removeB2B:
     if id[j] in bid:
       continue 
   if removeAdjacent:
     if j==i-1:
       continue
   if (end[j]-st[j])<minSize:
       continue
   if (end[j]-st[j])>maxSize:
       continue       
   d = st[i]+(end[i]-st[i])/2-(st[j]+(end[j]-st[j])/2)
   if d>=maxL:
    break
   of.write("%d\t%d\t%d\n" % (id[i], id[j], d))
   #n[id[i]][d/bin] += 1
  
  for j in range(i+1,len(st),1):
   if chr[j] != chr[i]:
    break
   if removeB2B:
     if id[j] in bid:
       continue
   if removeAdjacent:
     if j==i+1:
       continue 
   if (end[j]-st[j])<minSize:
       continue
   if (end[j]-st[j])>maxSize:
       continue       
   
   d = st[j]+(end[j]-st[j])/2-(st[i]+(end[i]-st[i])/2)
   if d>=maxL:
    break
   of.write("%d\t%d\t%d\n" % (id[i], id[j], d))
   #n[id[i]][d/bin] += 1
   
  if int(random.uniform(0,100))==1: 
    print "%d " % i

#print "\nWriting out text file..."

#of = open(outfile, "wt")
#of.write("#\tminsize=%d\tmaxsize=%d\tmaxl=%d\tbinsize=%d\tremoveb2b=%r\tremoveAdjacent=%r\trmapfile=%s\tbaitmapfile=%s\n" % \
#(minSize, maxSize, maxL, bin, removeB2B, removeAdjacent, rmapfile, bfile))
#for k in sorted(n.keys()):
# of.write("%d\t" % k)
# for i in range(len(n[k])): 
#  of.write("%d" % n[k][i])
#  if i!=len(n[k])-1:
#   of.write("\t")
# of.write("\n")
#of.close()

#if picklefile!=None:
#  print "Writing out pickle..."
#  pf = open(picklefile, "wb")
#  p.dump(n, pf)
#  pf.close()
