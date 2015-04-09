library(data.table)
library(matrixStats)
library(cluster)
library(argparser)
library(Hmisc)

cat("\n")
args = commandArgs(trailingOnly=T)
spec = matrix(c("<names-file>", "Full path to a tab-separated file with sample names (1st column) and full paths to input Rda files (2nd column)", 
                "<output-prefix>", "Output file names prefix (can contain path to folders)", 
                "<digest-map>", "Full path to digest map file", 
                "<bait-map>", "Full path to bait map ID file"),  byrow=T, ncol=2)
p = arg.parser("Generate a peak matrix and a sample tree from CHiCAGO output Rda files.", name="Rscript makePeakMatrix.R")
p = add.argument(p, arg=spec[,1], help=spec[,2])
p = add.argument(p, arg="--cutoff", help = "Score cutoff to use", default = 5, type = "numeric")
p = add.argument(p, arg="--subset", help = "Number of interactions to randomly subset for clustering", default = 100000)
p = add.argument(p, arg="--maxdist", help = "Max distance from bait to include into the peak matrix and clustering", default = NA, type="numeric")
p = add.argument(p, arg="--scorecol", help = "Column name for the scores", default = "score")
p = add.argument(p, arg="--clustmethod", help = "The clustering method to use (average/ward.D2/complete)", default = "average")
p = add.argument(p, arg="--usetrans", help = "Also include trans-interactions", flag = T)
p = add.argument(p, arg="--twopass", help="Should the list of significant interactions be generated first (saves memory, but requires two passes of reading the input", flag=T)
opts = parse.args(p, args)

namesfile = opts[["<names-file>"]]
prefix = opts[["<output-prefix>"]]
cutoff = opts[["cutoff"]]
sampsize = opts[["subset"]]
rmapfile = opts[["<digest-map>"]]
baitmapfile = opts[["<bait-map>"]]
maxdist = opts[["maxdist"]]
scorecol = opts[["scorecol"]]
clMethod = opts[["clustmethod"]]
useTrans = opts[["usetrans"]]
twoPass = opts[["twopass"]]

bm = fread(baitmapfile)
setnames(bm, c("V1", "V2", "V3", "V4", "V5"), c("baitChr", "baitStart", "baitEnd", "baitID", "baitName"))

input = read.table(namesfile,stringsAsFactors = F)
names(input) = c("name", "file")

sel = data.table(baitID=numeric(0), otherEndID=numeric(0))
if(twoPass){
  cat("Creating a list of significant interactions in at least one sample...\n")
  for(i in 1:nrow(input)){
    f = input[i,"file"]
    
    cat("Loading", f, "...\n")
    load(f)
  
    cat("Filtering and adding to the list...\n")  
    setDT(x)
    
    x = x[scorecol>=cutoff]  
    
    if(!useTrans){
      x = x[is.na(distSign)==F]
    }
    
    if (!is.na(maxdist)){
      x = x[abs(distSign)<=maxdist]
    }
    
    sel = rbindlist(list(sel, x[, c("baitID", "otherEndID"), with=F]))
    setkey(sel, baitID, otherEndID)
    sel = unique(sel)
  }
  
  cat("Reloading samples...\n")
}

data = vector("list")
for (i in 1:nrow(input)){
  f = input[i,"file"]    
  
  cat("Loading", f, "...\n")
  load(f)
  
  cat("Filtering and processing...\n")

  setDT(x)
  
  if(twoPass){
    setkey(x, baitID, otherEndID)  
    x = x[sel]    
  }
  else{
    x = x[scorecol>=cutoff]  

    if(!useTrans){
      x = x[is.na(distSign)==F]
    }
    
    if (!is.na(maxdist)){
      x = x[abs(distSign)<=maxdist]
    }    
  }
  
  name = input[i, "name"]
  data[[name]] = x[, c("baitID", "otherEndID", scorecol), with=F]
  setnames(data[[name]], scorecol, name)
  setkey(data[[name]], baitID, otherEndID)
}

cat("Merging...\n")

z = Reduce(function(x,y) merge(x,y,all=TRUE), data)

if(!twoPass){
  cat("Retaining only interactions exceeding score cutoff", cutoff, "in at least one sample...\n")
  scoreCols = names(z)[3:ncol(z)]
  z[, rmax:=eval(parse(text=paste0("pmax(", paste0(scoreCols, collapse=","),", na.rm=T)")))]
  z = z[rmax>=cutoff]
}

cat("Replacing NAs and negative scores with zeros...\n")
for (i in 3:ncol(z)){
  set(z, which(is.na(z[[i]]) | z[[i]]<0), i, 0)
}

cat("Adding bait coordinates and annotation...\n")
setkey(z, baitID)

#baitNames = gsub("\\-\\d{3}", "", bm$baitName) # remove transcript IDs
#baitNamesList = strsplit(baitNames, ",")
#baitNamesDedup = sapply(baitNamesList, function(x)paste(x[!duplicated(x)], collapse="," )) # put back together without duplicates
#bm$baitName = baitNamesDedup

setkey(bm, baitID)
z = merge(z, bm)

cat("Adding otherEnd coordinates...\n")
setkey(z, otherEndID)
rm = fread(rmapfile)
setnames(rm, c("V1", "V2", "V3", "V4"), c("oeChr", "oeStart", "oeEnd", "otherEndID"))
setkey(rm, otherEndID)
z = merge(z, rm)

cat("Annotating bait2bait interations...\n")
#setkey(z, otherEndID) - still remains one from the above code
setnames(bm, c("baitID", "baitName"), c("otherEndID", "otherEndName")) # sic!
setkey(bm, otherEndID)
z = merge(z, bm[, c("otherEndID", "otherEndName"), with=F], all.x=T, all.y=F)
set(z, which(is.na(z$otherEndName)), "otherEndName", ".")

cat("Adding distances and sorting...\n")

z[, dist := NA_real_]
z[baitChr==oeChr, dist := ceiling(oeStart+(oeEnd-oeStart)/2-(baitStart+(baitEnd-baitStart)/2))]

setnames(z, "otherEndID", "oeID")
setnames(z, "otherEndName", "oeName")
setDF(z)
z = z[, c("baitChr", "baitStart", "baitEnd", "baitID", "baitName", "oeChr", "oeStart", "oeEnd", "oeID", "oeName", "dist", input[,"name"])]

z = z[order(z$baitChr, z$baitStart, z$oeChr, z$oeStart), ]

rdaname = paste0(prefix, ".Rda")
cat(paste0("Saving the result image as ", rdaname, "...\n"))
save(z, file=rdaname)

txtname = paste0(prefix, ".txt")
cat(paste0("Writing out the result as ", txtname, "...\n"))
write.table(z, txtname, quote = F, sep = "\t", col.names = T, row.names=F)

if (nrow(input)>2){
  cat("Clustering samples based on", sampsize,  "random interactions...\n")
  
#   cat("Using continuous signals...\n") 
  zsamp = z[sample(1:nrow(z), sampsize),]
  d = dist(t(zsamp[,12:ncol(zsamp)]))
  h = hclust(d, method=clMethod)
  
  pdfname = paste0(prefix, "_tree.pdf")
  cat(paste0("Saving the sample dendrogram as ", pdfname, "...\n"))
  pdf(pdfname, width=20, height=10)
  plot(h)
  dev.off()
  
  #cat("Using binary signals ( cutoff =", cutoff, ")...\n")
  #zsbin = apply(zsamp[, 12:ncol(zsamp)],1,function(x){x[x<cutoff]=0; x[x>=cutoff];x})
  #db = daisy(t(zsbin), metric="gower")
  #hb = hclust(db, method=clMethod)
  #pdfname = paste0(prefix, "_tree_binary.pdf")
  #cat(paste0("Saving the sample dendrogram as ", pdfname, "...\n"))
  #pdf(pdfname, width=20, height=10)
  #plot(hb)
  #dev.off()
}else{
  cat("Clustering not performed as n<=2\n")
}
cat("Done!\n")

