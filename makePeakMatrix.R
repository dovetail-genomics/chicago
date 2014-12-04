args = commandArgs(trailingOnly=T)
if(!length(args)){
  stop("Usage: Rscript pullData.R <data-root-folder> <image-file-pattern>\nSee the code of this script for an example.\n")
  
}

#setwd("~/g/SNP_Project//CHiCAGOv2//DataRelease_141014//EachReplicate_141014/")
#files = list.files("~/g/SNP_Project//CHiCAGOv2//DataRelease_141014//EachReplicate_141014/", ".+(Macro|Mega|CD4_Naive|Erythro|Mono)\\S+RDa$", recursive=T)

setwd(args[1])
files = list.files(args[1], args[2], recursive=T)

cat("Input files to be used:\n")
print(files)

data = vector("list")
for (f in files){
     cat("loading", f, "...\n")
     load(f)

     cat("processing",  f, "...")
     x = x[!is.na(x$distSign),]
     t1 = gsub("step2_data_processed_separately/res_\\S+_chicago2/data/", "", f)
     t2 = gsub("_step2.RDa", "", t1)
     data[[t2]] = x[, c("baitID", "otherEndID", "score")]
     cat("done\n")    
}

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

cat("Filtering...\n")
z = as.data.frame(z)[rowMaxs(as.data.frame(z)[,3:ncol(z)], na.rm = T)>=12 ,]

cat("Replacing NAs with zeros...\n")
for (i in 3:ncol(z)){
   z[is.na(z[,i]), i] = 0
}
cat("Saving the result as scoreMatrix.Rda...\n")
save(z, file="scoreMatrix.Rda")

cat("Writing out the result as scoreMatrix.txt")
write.table(z, "scoreMatrix.txt", quote = F, sep = "\t", col.names = T, row.names=F)

cat("Clustering samples...\n")

zsamp = z[sample(1:nrow(z), 10000),]
d = dist(t(zsamp[,3:ncol(zsamp)]))
h = hclust(d)

pdf("samplesClustered.pdf")
plot(h)
dev.off()

cat("Done!\n")


