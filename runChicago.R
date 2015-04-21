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

p = add.argument(p, arg="--feature-files", 
                 help = "A comma-separated list of files with genomic feature coordinates for computing peaks' enrichment over these feature.", 
                 default = NA)
p = add.argument(p, arg="--feature-list", 
                 help = "Same as above but the supplied file contains the feature names and 
                 locations of feature files (in the format <feature-name> <feature-file-location>", default = NA)

p = add.argument(p, arg="--rda", help = "Save the Chicago object as an RDa image (instead of the default RDS)", flag = T)
p = add.argument(p, arg="--save-df-only", help = "Save only the data part of the Chicago object, as a data frame (for compatibility)", flag = T)

p = add.argument(p, arg="--examples-prox-dist", help = "The distance limit for plotting \"proximal\" examples", default=1e6)
p = add.argument(p, arg="--examples-full-range", help = "Also plot examples for the full distance range", flag = T)
p = add.argument(p, arg="--feat-max-dist", help = "The distance limit for computing enrichment for features", default=NA)

p = add.argument(p, arg="--output-dir", help = "The name of the output directory (can be a full path)", default="<output-prefix>")

opts = parse.args(p, args)

inputFiles = strsplit(opts[["<input-files>"]], "\\,")[[1]]
outPrefix_rel = opts[["<output-prefix>"]]

settingsFile = opts[["settings-file"]]
designDir = opts[["design-dir"]]

printMemory = opts[["print-memory"]]

cutoff = opts[["cutoff"]]
exportFormat = ifelse(is.na(opts[["export-format"]])[1], NA, strsplit(opts[["export-format"]], "\\,")[[1]])
exportOrder = opts[["export-order"]]

featureFiles = ifelse(is.na(opts[["feature-files"]])[1], NA, strsplit(opts[["feature-files"]], "\\,")[[1]])
featureList = opts[["feature-list"]]

isRda = opts[["rda"]]
isDF = opts[["save-df-only"]]

proxLim = opts[["examples-prox-dist"]]
plotFull = opts[["examples-full-range"]]

featDistLim = opts[["feat-max-dist"]]

outDir = ifelse(opts[["output-dir"]]=="<output-prefix>", outPrefix_rel, opts[["output-dir"]])

outPrefix = file.path(outDir, outPrefix_rel)

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
if(!dir.create.ifNotThere("enrichment-data")){
  stop("Couldn't create folder enrichment-data\n")      
}

### TODO: if bgzip and tabix index files are going to be produced in washU-track mode, move them too. 
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

if (!is.na(featureFiles)[1] | !is.na(featureList)){
  message("\n\nComputing enrichment for features...\n")
#   **TODO - UPDATE THIS**
#   peakEnrichment4Features(x1=x[!is.na(x$distSign),], score=12, sample_number=100, no_bins=100, 
#                           colname_score="score",folder=featureFolder, position_otherEnd_folder=fileDir, 
#                           list_frag=files, filterB2B=TRUE,
#                           colname_dist="distSign", beyond_dist=0, before_dist=1000000,
#                           plot_name=paste0(outprefix, "_feature_overlaps_upto_1M.pdf"))
#   
#   write.table(resTable, quote=F, row.names= , col.names=T, file=   paste0(outprefix, "_feature_overlaps_upto1M.txt")) 
#   setwd(outDir)
#   if(!moveToFolder("feature_overlaps", "enrichment")){
#     stop("Couldn't move feature_overlaps files to enrichment folder")
#   }
}
setwd(curDir)
message("All done!\n")