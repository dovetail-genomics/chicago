### Example workflow

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

pi.rel = 1E5

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
x <- getScores(x, relAbundance = pi.rel)
gc()
##FIXME call interactions

#save output somewhere...
save(x, f, file="~/CHiCAGOtest/results.Rda")