library(argparser)
cat("\n")
args = commandArgs(trailingOnly=T)
spec = matrix(c("<names-file>", "Full path to a tab-separated file with sample names (1st column) and full paths to input Rda files (2nd column)", 
                "<output-prefix>", "Output file names prefix (can contain path to folders)", 
                "<digest-map>", "Full path to digest map file", 
                "<bait-map>", "Full path to bait map ID file"),  byrow=T, ncol=2)
p = arg.parser("Generate a peak matrix and a sample tree from CHiCAGO output Rda files.", name="Rscript makePeakMatrix.R")
p = add.argument(p, arg=spec[,1], help=spec[,2])
p = add.argument(p, arg="--cutoff", help = "Score cutoff to use", default = 12, type = "numeric")
p = add.argument(p, arg="--subset", help = "Number of interactions to randomly subset for clustering", default =10000, type="numeric")
opts = parse.args(p, args)


# currdir = getwd()
# setwd(opts[["<data-root>"]])
# files = list.files(opts[["<data-root>"]], opts[["<file-pattern>"]], recursive=T)

namesfile = opts[["<names-file>"]]
prefix = opts[["<output-prefix>"]]
cutoff = opts[["cutoff"]]
sampsize = opts[["sample"]]
rmapfile = opts[["<digest-map>"]]
baitmapfile = opts[["<bait-map>"]]

input = read.table(namesfile,stringsAsFactors = F)
names(input) = c("name", "file")

data = vector("list")
for (i in 1:nrow(input)){
     f = input[i,"file"]
     cat("Loading", f, "...\n")
     load(f)

     cat("Processing",  f, "...")
     x = x[!is.na(x$distSign),]
     
#      t1 = gsub("step2_data_processed_separately/res_\\S+_chicago2/data/", "", f)
#      t2 = gsub("_step2.RDa", "", t1)
  
     name = input[i, "name"]
     data[[name]] = x[, c("baitID", "otherEndID", "score")]
     cat("done\n")    
}

# setwd(currdir)

#save(data, file="data.Rda")

library(data.table)
library(matrixStats)
#cat("Loading data...\n")
#load("data.Rda")

cat("Converting to data.table and indexing...\n")
for (d in names(data)){
  data[[d]] = data.table(data[[d]])
  setkey(data[[d]], baitID, otherEndID)
  setnames(data[[d]], "score", d)
}

cat("Merging...\n")
z = Reduce(function(x,y) merge(x,y,all=TRUE), data)

cat("Retaining only interactions exceeding score cutoff", cutoff, "in at least one sample...\n")
z = as.data.frame(z)[rowMaxs(as.data.frame(z)[,3:ncol(z)], na.rm = T)>=cutoff ,]

cat("Replacing NAs and negative scores with zeros...\n")
for (i in 3:ncol(z)){
   z[is.na(z[,i]) | z[,i]<0, i] = 0
}


cat("Adding bait coordinates and annotation...\n")
z = data.table(z)
setkey(z, baitID)
bm = fread(baitmapfile)
setnames(bm, c("V1", "V2", "V3", "V4", "V5"), c("baitChr", "baitStart", "baitEnd", "baitID", "baitName"))

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
z$otherEndName[is.na(z$otherEndName)] = "."

cat("Reordering and sorting...\n")

setnames(z, "otherEndID", "oeID")
setnames(z, "otherEndName", "oeName")
z = as.data.frame(z)
z = z[, c("baitChr", "baitStart", "baitEnd", "baitID", "baitName", "oeChr", "oeStart", "oeEnd", "oeID", "oeName", input[,"name"])]

z = z[order(z$baitChr, z$baitStart, z$oeChr, z$oeStart), ]

rdaname = paste0(prefix, ".Rda")
cat(paste0("Saving the result image as ", rdaname, "...\n"))
save(z, file=rdaname)

txtname = paste0(prefix, ".txt")
cat(paste0("Writing out the result as ", txtname, "...\n"))
write.table(z, txtname, quote = F, sep = "\t", col.names = T, row.names=F)

cat("Clustering samples based on", sampsize,  "random interactions...\n")

zsamp = z[sample(1:nrow(z), sampsize),]
d = dist(t(zsamp[,3:ncol(zsamp)]))
h = hclust(d)

pdfname = paste0(prefix, "_tree.pdf")
cat(paste0("Saving the sample dendrogram as ", pdfname, "...\n"))

pdf(pdfname)
plot(h)
dev.off()

cat("Done!\n")


