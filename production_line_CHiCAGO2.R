### EDIT ME

### If running from RStudio interactively
# if(Sys.getenv("RSTUDIO_USER_IDENTITY")==""){
### or, simply, if clArgs is not defined:
if (!exists("clArgs")){
  clArgs = commandArgs(trailingOnly=T)
}
# else - expect to have clArgs defined explicitly in the local environment,
# such as clArgs = c(<input-file-name>, <output-folder>, <output-prefix>, <feature-folder>, <feature-list>)
# or clArgs = c(%<Nfiles>, <output-folder>, <output-prefix>, <feature-folder>, <feature-list>, <input-file-name-1>, ..., <input-file-name-N> <fileDirDigest> <headerDigest>)

print (clArgs)
print (length(clArgs))

if(length(clArgs)<5){
  stop ("For one sample, supply <input-file-name> <output-folder> <output-prefix> <feature-folder> <feature-list> as arguments\n
        for multiple samples, supply %<Nfiles> <output-folder> <output-prefix> <feature-folder> <feature-list> <input-file-1> ... <input-file-N> <fileDirDigest> <headerDigest>\n")
}

infname = clArgs[1]
outfolder = clArgs[2]
outprefix = paste0(outfolder, "/", clArgs[3])
featureFolder = paste0(clArgs[4],"/") # important as currently CompareSeqTotal requires featureFolder to end with a slash
featureList = clArgs[5]

features = read.table(featureList, stringsAsFactors=F)
files = features[,2]
names(files) = features[,1]




# CHiCAGOv2: calling interactions
scriptDir <- "/bi/apps/chicago/0.1.0.dev"
source(file.path(scriptDir, "chicago.R"))
# Testing enrichment of CHiCAGO peaks for genomic features of interest
source(file.path(scriptDir, "Functions_new_datatable.R"))
source(file.path(scriptDir, "run_peakEnrichment4Features.R"))
fileDir =  clArgs[length(clArgs)-1]
header= clArgs[length(clArgs)]
#fileDir = "/bi/group/sysgen/CHIC"
###

### Resource file locations
rmapfile= file.path(fileDir, paste0(header,".bed"))
baitmapfile= file.path(fileDir,  paste0(header,"_baits_ID.bed"))

nperbinfile = file.path(fileDir,  paste0(header,"_NperBin.txt"))
nbaitsperbinfile = file.path(fileDir, paste0(header,"_NbaitsPerBin.txt"))
proxOEfile = file.path(fileDir, paste0(header,"_proxOEout.txt"))
###############

### Fragment filtering and other settings
maxLBrownEst = 1.5e6 # maximum distance for Brownian noise estimation 
minFragLen = 150 # minimum allowed length of restriction fragments (otherwise ignored)
maxFragLen = 40000 # max allowed length of restriction fragments (otherwise ignored)
minNPerBait = 250 # minimum number of reads per bait (otherwise ignored)
binsize=20000 # bin size for parameter estimation 
removeAdjacent = TRUE # should restriction fragments immediately adjacent to baits be removed? 

### Free parameter: our expectation for the relative enrichment of true interactions close to the bait vs in trans
pi.rel = 1E5

### Score cutoff for writing out interaction files
outputCutoff=12



# 1. Read the input file(s)
# Note that for all downstream functions,
# a single replicate is supplied as a data frame
# and multiple replicates as a list of data frames.

cat("Reading data...\n")

if (length(grep("^%\\d+$", infname))){
  nfiles = as.numeric(gsub("%","", infname))
  cat ("Number of input files:" , nfiles, ".\n")
  if (length(clArgs)<3+nfiles){
    stop("Fewer files supplied than %<Nfiles>\n")
  }
  #infnames = clArgs[4:(3+nfiles)]
  infnames = clArgs[6:(5+nfiles)]
  cat ("Input files:\n")
  print (infnames)
  cat("\n*** Reading the input files... ***\n")
  x1 = vector("list")
  for (f in infnames){ 
    x1[[f]] = readSample(f)
  }
  cat("\n*** Running mergeSamples...\n")
  x = mergeSamples(x1)
}else{
  cat("Input file:", infname, "\n")
  cat("\n*** Reading the input file... ***\n")
  x = readSample(infname)
}

system(paste("mkdir -p", outfolder))

# 2. Running interaction calling

x <- chicagoPipeline(x, outprefix, pi.rel)

# 3. Saving image

cat("\n*** Saving image...\n")

save(x, file=paste0(outprefix,".RDa"))

# 4. Plotting examples

cat("\n*** Plotting examples...\n")

baits=plotBaits(x, outfile=paste0(outprefix, "_proxExamples.pdf"), xlim=c(-1e6,1e6))
plotBaits(x, baits=baits, outfile=paste0(outprefix, "_examples.pdf"))

# 5. Exporting results

cat("\n*** Exporting results...\n")
exportResults(x, outfile=paste0(outprefix,".txt"), cutoff=outputCutoff)
exportResults(x, outfile=paste0(outprefix,".ibed"), format="interBed", cutoff=outputCutoff)


# 6. Creating subfolders

cat("\n*** Sorting output files into subfolders...\n")

system(paste("mkdir -p", paste0(outfolder, "/diag_plots")))     
system(paste("mkdir -p", paste0(outfolder, "/examples")))
system(paste("mkdir -p", paste0(outfolder, "/data")))

system(paste0("mv ", outfolder, "/*.txt ", outfolder, "/data"))
system(paste0("mv ", outfolder, "/*.ibed ", outfolder, "/data"))
system(paste0("mv ", outfolder, "/*.RDa ", outfolder, "/data"))
system(paste0("mv ", outfolder, "/*xamples.pdf ", outfolder, "/examples"))
system(paste0("mv ", outfolder, "/*.pdf ", outfolder, "/diag_plots"))       

# 7. Producing overlaps with genomic features of interest - this can be for example ENCODE or BLUEPRINT chipseq files in bed format

# If a CompareSeqTotal run fails, try some combination of the following:
# - reducing the ncores parameter for multicore operations (currently set to 4, but the default is 8)
# - using only a sample of the negative set for analysis (negFraction<1, currently set to 0.3 for distal only)
# - using multicoreSampling=F (longer, but saves memory; currently set to TRUE for distal only; 
#                    note that in this case multicore will still be used in some other, less memory-hungry operations)
#
# We are working on further optimisations to this algorithm and will release them as soon as we are ready. 

# This warning message is expected and not a problem:
# Warning message:
# In `[.data.table`(bin_reads2, , `:=`(distbin3, udbin3)) :
#  Supplied 500 items to be assigned to xxx items of column 'distbin3' (xxx unused)

cat("\n*** Running compareSeq...\n")

resTable<-peakEnrichment4Features(x1=x[!is.na(x$distSign),], score=12, sample_number=100, no_bins=100, 
colname_score="score",folder=featureFolder, position_otherEnd_folder=fileDir, 
list_frag=files, filterB2B=TRUE,
colname_dist="distSign", beyond_dist=0, before_dist=1000000,
plot_name=paste0(outprefix, "_feature_overlaps_upto_1M.pdf"))

system(paste("mkdir -p", paste0(outfolder, "/overlap_data")))
write.table(resTable, quote=F, row.names= , col.names=T, file=   paste0(outprefix, "_feature_overlaps_upto1M.txt")) 
system(paste0("mv ", outfolder, "/*feature_overlaps*.* ", outfolder, "/overlap_data"))

cat("\n*** All done  ***\n")

