### Test changing the file

#Get code:
setwd("~/CHiCAGOtest/")
source("chicago.R")

#Overwrite file locations
fileDir = "/bi/group/sysgen/CHIC"
rmapfile= file.path(fileDir, "Digest_Human_HindIII.bed")
baitmapfile= file.path(fileDir, "Digest_Human_HindIII_baits_ID.bed")
nperbinfile = file.path(fileDir, "Digest_Human_HindIII_NperBin.txt")
nbaitsperbinfile = file.path(fileDir, "Digest_Human_HindIII_NbaitsPerBin.txt")
proxOEfile = file.path(fileDir, "proxOE_out.txt")


#Choose data to read in:
dataFileDir <- "/bi/group/sysgen/hESC"
datafile <- file.path(dataFileDir, "sample_hESC811_new/hESC811_new_bait_otherEnd_N_len_distSign.txt")
peakfile <- file.path(dataFileDir, "res_hESC_811/data/hESC_811.RData")


#Run CHiCAGO on data:

x <- readSample(datafile)
gc()
x <- normaliseBaits(x)
gc()
x <- normaliseOtherEnds(x)
gc()
f <- estimateDistFun(x)
gc()
x <- estimateTechnicalNoise(x)
gc()
#x <- estimateBrownianNoise(x, f) ##all
x <- estimateBrownianNoise(x, f, subset=1000)
gc()
x <- getPvals(x)
gc()
x <- getFDRs(x)
gc()
##FIXME call interactions

#save output somewhere...
save(x, f, file="results.Rda")

##All of the below was put in TestReport.Rmd

# ##a quick test (to make sure the table writes correctly)----------------
# #testfribble <- x[1:10,]
# #write.csv(, "")
# 
# ##compare against v1 output---------------------------------------------
# load(peakfile)
# 
# #baitinfo <- fread(baitmapfile)
# #fraginfo <- fread(rmapfile)
# 
# ##combine previous results into a nice(r) format
# cols <- c("baitID", "otherEndID", "isBait2bait")
# n1cat <- rbind(
#   as.data.table(n1pr[,c(cols,"logwp")]),
#   as.data.table(n1farpr[,c(cols,"logwp_bw")]),
#   use.names=FALSE ##because the colnames differ
#   )
# setnames(n1cat, c(cols,"logwp"))
# 
# x.sel <- as.data.table(x[,c("baitID", "otherEndID", "log.q")])
# 
# setkey(n1cat, baitID, otherEndID)
# setkey(x.sel, baitID, otherEndID)
# 
# venn(list(unique(x.sel$baitID), unique(n1cat$baitID)))
