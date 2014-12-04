args = commandArgs(trailingOnly=T)
if(!length(args)){
  stop("Usage: Rscript makePeakMatrix.R <data-root-folder> <image-file-pattern> <output-file-prefix>\n
Example:\n\tRscript DataRelease_1/EachReplicate .+(Macro|Mega|CD4_Naive|Erythro|Mono)\\S+RDa$\n")
  
}

setwd(args[1])
files = list.files(args[1], args[2], recursive=T)
prefix = args[3]

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

cat("Replacing NAs and negative scores with zeros...\n")
for (i in 3:ncol(z)){
   z[is.na(z[,i]) | z[,i]<0, i] = 0
}

rdaname = paste0(prefix, ".Rda")
cat(paste0("Saving the result image as ", rdaname, "...\n")
save(z, file=rdaname)

txtname = paste0(prefix, ".txt")
cat(paste0("Writing out the result as ", txtname, "...\n")
write.table(z, txtname, quote = F, sep = "\t", col.names = T, row.names=F)

cat("Clustering samples based on 10000 random interactions...\n")

zsamp = z[sample(1:nrow(z), 10000),]
d = dist(t(zsamp[,3:ncol(zsamp)]))
h = hclust(d)

pdfname = paste0(prefix, "_tree.pdf")
cat(paste0("Saving the sample dendrogram as ", pdfname, "...\n")

pdf(pdfname)
plot(h)
dev.off()

cat("Done!\n")


