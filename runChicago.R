message("\n***runChicago.R\n\n")

message("Loading libraries...\n")

library(argparser)
library(Chicago)

dir.create.ifNotThere = function(path, ...){
  if(file.exists(path)){
    if(!file.info(path)$isdir){
      message("Cannot create directory ", path, " as it coincides with an existing file\n")
      return(0)
    }else{
      message("\nWarning: directory ", path, " exists and will be reused.")
      return(1)
    }
  }
  else{
    dir.create(path, ...)
  }
}

moveToFolder = function(pattern, where){
  files = list.files(".", pattern)
  if (!length(files)){
    message("No files found matching the pattern ", pattern)
    return(0)
  }
  res = vector("numeric", length(files))
  i=1
  for(f in files){
    res[i] = file.rename(f, file.path(where, f))
    if (!res[i]){
      message("Warning: did not succeed in moving ", f, " into ", where)
    }
    i = i+1
  }
  min(res)
}

args = commandArgs(trailingOnly=T)
spec = matrix(c("<input-files>", "Full path to the input file (or comma-separated list of files)", 
                "<output-prefix>", "Output file names prefix (cannot contain folders)"),  byrow=T, ncol=2)

p = arg.parser("Run Chicago from input files", name="Rscript runChicago.R")
p = add.argument(p, arg=spec[,1], help=spec[,2])

p = add.argument(p, arg="--settings-file", help = "Full path to Chicago settings file", default = NA)
p = add.argument(p, arg="--design-dir", 
                 help = "Folder with capture design files (note the settings file has priority over these)", default = "")

p = add.argument(p, arg="--print-memory", help = "Should chicagoPipeline print out memory use?", flag=T)

p = add.argument(p, arg="--cutoff", help = "Score cutoff for writing out peaks and testing feature enrichment", default = 5)
p = add.argument(p, arg="--export-format", 
                 help = "File format for writing out peaks: one or more of the following: seqMonk,interBed,washU_text,washU_track (comma-separated)", 
                 default = "washU_text")
p = add.argument(p, arg="--export-order", help = "Should the results be ordered by \"score\" or genomic \"position\"?", 
                 default = "position")

p = add.argument(p, arg="--rda", help = "Save the Chicago object as an RDa image (instead of the default RDS)", flag = T)
p = add.argument(p, arg="--save-df-only", help = "Save only the data part of the Chicago object, as a data frame (for compatibility)", flag = T)

p = add.argument(p, arg="--examples-prox-dist", help = "The distance limit for plotting \"proximal\" examples", default=1e6)
p = add.argument(p, arg="--examples-full-range", help = "Also plot examples for the full distance range", flag = T)

p = add.argument(p, arg="--output-dir", help = "The name of the output directory (can be a full path)", default="<output-prefix>")

p = add.argument(p, arg="--en-feat-files", 
                 help = "A comma-separated list of files with genomic feature coordinates for computing peaks' enrichment", 
                 default = NA)
p = add.argument(p, arg="--en-feat-list", 
                 help = "Same as above but the supplied file contains the feature names and 
                 locations of feature files (in the format <feature-name> <feature-file-location>", default = NA)
p = add.argument(p, arg="--en-feat-folder", 
                 help = "The folder, in which all feature files are located (if provided, --en-feature-file(s) don't need to list the full path)", 
                 default=NA)

p = add.argument(p, arg="--en-min-dist", help = "The lower distance limit for computing enrichment for features", default=0)
p = add.argument(p, arg="--en-max-dist", help = "The upper distance limit for computing enrichment for features", default=1e6)
p = add.argument(p, arg="--en-full-range", help = "Assess the enrichment for features for the full distance range (can be very slow!)", flag = T)
p = add.argument(p, arg="--en-sample-no", help = "The number of negative samples for computing enrichment for features", default=100)

opts = parse.args(p, args)

inputFiles = strsplit(opts[["<input-files>"]], "\\,")[[1]]
outPrefix_rel = opts[["<output-prefix>"]]

settingsFile = opts[["settings-file"]]
designDir = opts[["design-dir"]]

printMemory = opts[["print-memory"]]

cutoff = opts[["cutoff"]]
exportFormat = ifelse(is.na(opts[["export-format"]])[1], NA, strsplit(opts[["export-format"]], "\\,")[[1]])
exportOrder = opts[["export-order"]]

isRda = opts[["rda"]]
isDF = opts[["save-df-only"]]

proxLim = opts[["examples-prox-dist"]]
plotFull = opts[["examples-full-range"]]

outDir = ifelse(opts[["output-dir"]]=="<output-prefix>", outPrefix_rel, opts[["output-dir"]])

outPrefix = file.path(outDir, outPrefix_rel)

enSampleNumber = opts[["en-sample-no"]]
enMaxDist = opts[["en-max-dist"]]
if (opts[["en-full-range"]]){
  enMaxDist = NULL
  message("Warning: --en-full-range selected. Feature enrichment computation will be very slow\n")
}
enMinDist = opts[["en-min-dist"]]
enFeatFolder = opts[["en-feat-folder"]]
enFeatFiles = opts[["en-feat-files"]]
enFeatList = opts[["en-feat-list"]]

computeEnrichment = 1
if(is.na(enFeatFiles) & is.na(enFeatList)){
  message("Warning: neither --en-feat-files nor --en-feat-list provided. Feature enrichments will not be computed\n")
  computeEnrichment = 0
}

if(!is.na(enFeatFiles) & !is.na(enFeatList)){
  stop("Only one of --en-feat-files or --en-feat-list should be provided, not both\n")
}

if(!is.na(enFeatList)){
  featList = read.table(enFeatList, header=F, stringsAsFactors = F)
  if (ncol(featList)!=2){
    stop("--en-feat-list file should have two columns: <feature-name> <feature-file>")
  }
  enFeatFiles = unlist(featList[,2])
  names(enFeatFiles) = unlist(featList[,1])
}

if(is.na(settingsFile)){
  settingsFile = NULL
}

if (!all(exportFormat %in% c("seqMonk", "interBed", "washU_text", "washU_track"))) {
  stop("--exportFormat must be either seqMonk, interBed, washU_text or washU_track (or a comma-separated combination of these)\n")
}
if (!exportOrder %in% c("position", "score")) {
  stop("--exportOrder must be either position (default) or score\n")
}

message("Setting the experiment...\n")
cd = setExperiment(designDir = designDir, settingsFile = settingsFile)


if(length(inputFiles)>1){
  message("\nReading input files...\n")
  cd = readAndMerge(inputFiles, cd)
}else{
  message("\n")
  cd = readSample(inputFiles, cd)
}

if(!dir.create.ifNotThere(outDir, recursive = T)){
  stop(paste("Couldn't create folder", outDir, "\n"))
}


message("\n\nStarting chicagoPipeline...\n")
cd = chicagoPipeline(cd, outprefix = outPrefix, printMemory = printMemory)

if (isDF){
  message("\n\nSaving the image of Chicago output as a data frame...\n")
  y = cd@x
  setDF(y)
  if (isRda){ 
    save(y, file=paste0(outPrefix, "_df.RDa"))
  }else{
    saveRDS(y, paste0(outPrefix, "_df.Rds"))
  }
}else{
  message("\n\nSaving the Chicago object...\n")
  if (isRda){ 
    save(cd, file=paste0(outPrefix, ".RDa"))
  }else{
    saveRDS(cd, paste0(outPrefix, "_df.Rds"))
  }
}

message("\n\nPlotting examples...\n")
baits=plotBaits(cd, outfile=paste0(outPrefix, "_proxExamples.pdf"), xlim=c(-proxLim,proxLim))
if (plotFull){
  plotBaits(cd, baits=baits, outfile=paste0(outPrefix, "_examples.pdf"))
}

message("\n\nExporting peak lists...\n")
exportResults(cd, outfileprefix=outPrefix, cutoff=cutoff, format = exportFormat, order = exportOrder)

message("\n\nSorting output files into folders...\n")

curDir = getwd()
setwd(outDir)

if(!dir.create.ifNotThere("data")){
  stop("Couldn't create folder data\n")  
}
if(!dir.create.ifNotThere("diag_plots")){
  stop("Couldn't create folder diag_plots\n")
}
if(!dir.create.ifNotThere("examples")){
  stop("Couldn't create folder examples\n")      
}
if(!dir.create.ifNotThere("enrichment_data")){
  stop("Couldn't create folder enrichment_data\n")      
}

### TODO: if bgzip and tabix index files are going to be produced in washU-track mode, move them too. 
### Take into account that unlike other files, there may not be .gz / .gz.tbi created 
if (!moveToFolder("\\.txt", "data")){
  stop("Couldn't move txt files to data folder\n")
}

if("interBed" %in% exportFormat){
  if(!moveToFolder("\\.ibed", "data")){
    stop("Couldn't move ibed files to data folder\n")
  }
}

if(isRda){
  if(!moveToFolder("\\.RDa", "data")){
    stop("Couldn't move the RDa file to data folder\n")
  }
}else{
  if(!moveToFolder("\\.Rds", "data")){
    stop("Couldn't move the Rds file to data folder\n")
  }
}

if(!moveToFolder("xamples\\.pdf", "examples")){
  stop("Couldn't move pdf file(s) to example folder\n")
}

if(!moveToFolder("\\.pdf", "diag_plots")){
  stop("Couldn't move pdf files to diag_plots folder\n")
}

setwd(curDir)

if (computeEnrichment){
  message("\n\nComputing enrichment for features...\n")
  
  # enSampleNumber (100 by default)
  # enMaxDist (can be NULL if whole range; 1e6 by default)
  # enMinDist (0 by default)
  # enFeatFolder (can be NULL if full paths are provided)
  # enFeatFiles (can be a named vector but doesn't have to be) 
  # enOutputFile
  
  
  noBins = max(5, ceiling(ifelse(is.null(enMaxDist), max(abs(cd@x$distSign), na.rm = T), enMaxDist) - enMinDist)/1e4)

  enrichments = peakEnrichment4Features(cd, score=cutoff, sample_number=enSampleNumber, no_bins=noBins, 
                           colname_score="score", folder=enFeatFolder, list_frag=enFeatFiles, 
                           filterB2B=TRUE, min_dist=enMinDist, max_dist=enMaxDist,
                           plot_name=paste0(outPrefix, "_feature_overlaps.pdf"))
  
  enOutputFile = paste0(outPrefix, "_feature_overlaps.txt")
  cat(paste0("#\tmin_dist=", enMinDist, "\tmax_dist=", ifelse(is.null(enMaxDist), "whole_range", enMaxDist), "\n"), file = enOutputFile) 
  suppressWarnings(write.table(enrichments, quote=F, row.names= T, col.names=T, file= enOutputFile, append = T))
  setwd(outDir)
  if(!moveToFolder("feature_overlaps", "enrichment_data")){
     stop("Couldn't move feature_overlaps files to enrichment folder")
 }
}
setwd(curDir)
message("All done!\n")
