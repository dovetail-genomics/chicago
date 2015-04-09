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
p = add.argument(p, arg="--chunks", help="The number of chunks to split the merging process to save memory at the expense of time", default=1)
opts = parse.args(p, args)


# currdir = getwd()
# setwd(opts[["<data-root>"]])
# files = list.files(opts[["<data-root>"]], opts[["<file-pattern>"]], recursive=T)

namesfile = opts[["<names-file>"]]
prefix = opts[["<output-prefix>"]]
cutoff = opts[["cutoff"]]
sampsize = opts[["subset"]]
rmapfile = opts[["<digest-map>"]]
baitmapfile = opts[["<bait-map>"]]
maxdist = opts[["maxdist"]]
scorecol = opts[["scorecol"]]
clMethod = opts[["clustmethod"]]
nChunks = opts[["chunks"]] 

bm = fread(baitmapfile)
setnames(bm, c("V1", "V2", "V3", "V4", "V5"), c("baitChr", "baitStart", "baitEnd", "baitID", "baitName"))
baits = unique(bm$baitID)
baitsCut = cut2(baits, g=nChunks)

input = read.table(namesfile,stringsAsFactors = F)
names(input) = c("name", "file")

zlist = vector("list")
for (i in 1:nChunks){
  cat("Processing baits ", baitsCut[i], "\n...")
  cBaits = baits[as.numeric(baitsCut)==i]
  data = vector("list")
  for (i in 1:nrow(input)){
    f = input[i,"file"]
    cat("Loading", f, "...\n")
    load(f)
    
    cat("Processing",  f, "...")
    x = x[!is.na(x$distSign),]
    
    if (!is.na(maxdist)){
      x = x[abs(x$distSign)<=maxdist,]
      if (nChunks>1){
        setDT(x)
        setkey(x, baitID)
        x = x[J(cBaits)]
        setDF(x)
      }
      
    }
    
    name = input[i, "name"]
    data[[name]] = x[, c("baitID", "otherEndID", scorecol)]
    
    if (scorecol!="score"){
      names(data[[name]])[names(data[[name]])==scorecol] = "score"
    }
    cat("done\n")    
  }
  
  cat("Converting to data.table and indexing...\n")
  for (d in names(data)){
    setDT(data[[d]])
    setkey(data[[d]], baitID, otherEndID)
    setnames(data[[d]], "score", d)
  }
  
  cat("Merging...\n")
  z = Reduce(function(x,y) merge(x,y,all=TRUE), data)
  
  cat("Retaining only interactions exceeding score cutoff", cutoff, "in at least one sample...\n")
  scoreCols = names(z)[3:ncol(z)]
  z[, rmax:=eval(parse(text=paste0("pmax(", paste0(scoreCols, collapse=","),", na.rm=T)")))]
  z = z[rmax>=cutoff]
  
  cat("Replacing NAs and negative scores with zeros...\n")
  for (i in 3:ncol(z)){
    set(z, which(is.na(z[[i]]) | z[[i]]<0), i, 0)
  }
  
  cat("Adding bait coordinates and annotation...\n")
  setkey(z, baitID)
  
  baitNames = gsub("\\-\\d{3}", "", bm$baitName) # remove transcript IDs
  baitNamesList = strsplit(baitNames, ",")
  baitNamesDedup = sapply(baitNamesList, function(x)paste(x[!duplicated(x)], collapse="," )) # put back together without duplicates
  bm$baitName = baitNamesDedup
  
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
  
  cat("Reordering and sorting...\n")
  
  z[, dist := NA_real_]
  z[baitChr==oeChr, dist := ceiling(oeStart+(oeEnd-oeStart)/2-(baitStart+(baitEnd-baitStart)/2))]
  
  setnames(z, "otherEndID", "oeID")
  setnames(z, "otherEndName", "oeName")
  setDF(z)
  z = z[, c("baitChr", "baitStart", "baitEnd", "baitID", "baitName", "oeChr", "oeStart", "oeEnd", "oeID", "oeName", "dist", input[,"name"])]
  
  z = z[order(z$baitChr, z$baitStart, z$oeChr, z$oeStart), ]
  zlist[[i]] = z
}

z = rbindlist(zlist)

rdaname = paste0(prefix, ".Rda")
cat(paste0("Saving the result image as ", rdaname, "...\n"))
save(z, file=rdaname)

txtname = paste0(prefix, ".txt")
cat(paste0("Writing out the result as ", txtname, "...\n"))
write.table(z, txtname, quote = F, sep = "\t", col.names = T, row.names=F)

if (nrow(input)>2){
  cat("Clustering samples based on", sampsize,  "random interactions...\n")
  
  
  cat("Using continuous signals...\n") 
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

