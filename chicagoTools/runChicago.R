message("\n***runChicago.R\n\n")

### Parsing the command line ###

library(argparser)
message("\n")

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

p = add.argument(p, arg="--en-min-dist", help = "The lower distance limit for computing enrichment for features", default="0")
p = add.argument(p, arg="--en-max-dist", help = "The upper distance limit for computing enrichment for features", default=1e6)
p = add.argument(p, arg="--en-full-cis-range", help = "Assess the enrichment for features for the full distance range [same chromosome only; use --en-trans in addition to include trans-interactions] (can be very slow!)", flag = T)
p = add.argument(p, arg="--en-sample-no", help = "The number of negative samples for computing enrichment for features", default=100)
p = add.argument(p, arg="--en-trans", help = "Include trans-interactions into enrichment analysis", flag=T)

p = add.argument(p, arg="--features-only", help = "Re-run feature enrichment analysis with Chicago output files. With this option, <input-files> must be either a single Rds file (must contain full Chicago objects) or '-', in which case the file location will be inferred automatically from <output-prefix> and files added to the corresponding folder.",  flag = T)

opts = parse.args(p, args)


### Auxilliary functions ###
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

### Main code ###

message("Loading the Chicago package and dependencies...\n")
library(Chicago)

featuresOnly = opts[["features-only"]]

inputFiles = strsplit(opts[["<input-files>"]], "\\,")[[1]]
outPrefix_rel = opts[["<output-prefix>"]]

settingsFile = opts[["settings-file"]]
designDir = opts[["design-dir"]]

printMemory = opts[["print-memory"]]

cutoff = opts[["cutoff"]]

exportFormat = NA
if(!is.na(opts[["export-format"]])[1]){
   exportFormat = strsplit(opts[["export-format"]], "\\,")[[1]]
}

exportOrder = opts[["export-order"]]

isRda = opts[["rda"]]
isDF = opts[["save-df-only"]]

proxLim = opts[["examples-prox-dist"]]
plotFull = opts[["examples-full-range"]]

outDir = ifelse(opts[["output-dir"]]=="<output-prefix>", outPrefix_rel, opts[["output-dir"]])

outPrefix = file.path(outDir, outPrefix_rel)

enSampleNumber = opts[["en-sample-no"]]
enMaxDist = opts[["en-max-dist"]]
if (opts[["en-full-cis-range"]]){
  enMaxDist = NULL
  message("Warning: --en-full-cis-range selected. Feature enrichment computation will be very slow\n")
}
enMinDist = opts[["en-min-dist"]]
enTrans = opts[["en-trans"]]

if (enMinDist == "NULL" & enTrans){
  message("Running enrichment analysis for trans interactions only.")
  enMinDist = NULL
  enMaxDist = NULL
}else{
  enMinDist = as.numeric(enMinDist)
}

enFeatFolder = opts[["en-feat-folder"]]
enFeatFiles = opts[["en-feat-files"]]
enFeatList = opts[["en-feat-list"]]


if(enTrans & !is.null(enMaxDist) & !is.null(enMinDist)){
  message("\nWarning: --en-trans specificed together with --en-max-dist and --en-min-dist (possibly keeping the default values), which was likely not intended. Running the enrichment analysis for the full range. To test trans only, rerun with --en-min-dist NULL.\n")
  enMaxDist = NULL
}
  
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

if(!featuresOnly){

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
}
  
if(!dir.create.ifNotThere(outDir, recursive = T)){
  stop(paste("Couldn't create folder", outDir, "\n"))
}

if(!featuresOnly){
  logfile = file(paste0(outPrefix, "_params.txt"), "w")
  cat("#  runChicago parameters:\n", file=logfile)
  for (arg in names(opts)){
          cat(paste(arg, opts[[arg]], sep="\t"), "\n", file=logfile)
  }
  close(logfile)
}else{
  logfile = file(paste0(outPrefix, "_params.txt"), "a")
  cat("#  runChicago featuresOnly parameters:\n", file=logfile)
  for (arg in names(opts)){
    cat(paste(arg, opts[[arg]], sep="\t"), "\n", file=logfile)
  }
  close(logfile)
}

if (!featuresOnly){
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
      saveRDS(cd, paste0(outPrefix, ".Rds"))
    }
  }
  
  logfile = file(paste0(outPrefix, "_params.txt"), "a")
  sink(logfile, append=T)
  cat("\n#  chicago pipeline settings (chicagoData@settings):\n")
  for (s in names(cd@settings)){
          cat(paste(s, cd@settings[[s]], sep="\t"), "\n")
  }
  cat("\n#  sessionInfo()\n")
  print(sessionInfo())
  sink(NULL)
  close(logfile)
  
  message("\n\nPlotting examples...\n")
  baits=plotBaits(cd, outfile=paste0(outPrefix, "_proxExamples.pdf"), xlim=c(-proxLim,proxLim))
  if (plotFull){
    plotBaits(cd, baits=baits, outfile=paste0(outPrefix, "_examples.pdf"))
  }
  
  message("\n\nExporting peak lists...\n")
  exportResults(cd, outfileprefix=outPrefix, cutoff=cutoff, format = exportFormat, order = exportOrder)
  
  message("\n\nSorting output files into folders...\n")
}

curDir = getwd()
setwd(outDir)


if(!dir.create.ifNotThere("data")){
  stop("Couldn't create folder data\n")  
}

if(!featuresOnly){
  if(!dir.create.ifNotThere("diag_plots")){
    stop("Couldn't create folder diag_plots\n")
  }
  if(!dir.create.ifNotThere("examples")){
    stop("Couldn't create folder examples\n")      
  }
}
  
if (!moveToFolder("\\.txt", "data")){
  stop("Couldn't move txt files to data folder\n")
}

if(!featuresOnly){
  if("interBed" %in% exportFormat){
    if(!moveToFolder("\\.ibed", "data")){
      stop("Couldn't move ibed files to data folder\n")
    }
  }

  ### If bgzip and tabix index files are going to be produced in washU-track mode, move them too. 
  ### Unlike other files, there may not be .gz / .gz.tbi created.
  if("washU_track" %in% exportFormat){
    if(!moveToFolder("\\.txt.gz", "data")){
      message("Couldn't move .txt.gz files (washU track) to data folder\n")
    }
    if(!moveToFolder("\\.txt.gz.tbi", "data")){
      message("Couldn't move .txt.gz.tbi files (washU track) to data folder\n")
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
}

if(!dir.create.ifNotThere("enrichment_data")){
  stop("Couldn't create folder enrichment_data\n")      
}

setwd(curDir)

if(featuresOnly){
  if (inputFiles[1]=="-"){
    dataPath = file.path(outDir, "data")
    cd = readRDS(file.path(dataPath, paste0(outPrefix_rel, ".Rds")))
  }else{
    cd = readRDS(inputFiles[1])
  }
}

if (computeEnrichment){
  message("\n\nComputing enrichment for features...\n")
  
  # enSampleNumber (100 by default)
  # enMaxDist (can be NULL if whole range; 1e6 by default)
  # enMinDist (0 by default)
  # enFeatFolder (can be NULL if full paths are provided)
  # enFeatFiles (can be a named vector but doesn't have to be) 
  # enOutputFile
  
  
  noBins = max(5, ceiling(ifelse(is.null(enMaxDist), max(abs(cd@x$distSign), na.rm = T), enMaxDist) - enMinDist)/1e4)
  
  plot_name = paste0(outPrefix, "_feature_overlaps.pdf")
  i=0
  raw.name = file.path(file.path(outDir, "enrichment_data"), paste0(outPrefix_rel, "_feature_overlaps.pdf"))
  # NB: since the file is getting moved post-hoc, the location we are checking is different from plot_name 
  while(file.exists(raw.name)){
    i = i+1
    plot_name = paste0(outPrefix, "_feature_overlaps", i, ".pdf")
    raw.name = file.path(file.path(outDir, "enrichment_data"), paste0(outPrefix_rel, "_feature_overlaps", i, ".pdf"))
  }
  
  if(i){
    message("Existing enrichment data found, so the data for this run will be saved as ", outPrefix_rel, "_feature_overlaps", i)  
  }
  
  enrichments = peakEnrichment4Features(cd, score=cutoff, sample_number=enSampleNumber, no_bins=noBins, 
                           colname_score="score", folder=enFeatFolder, list_frag=enFeatFiles, trans=enTrans,
                           filterB2B=TRUE, min_dist=enMinDist, max_dist=enMaxDist,
                           plot_name=plot_name)
  
  
  enOutputFile = paste0(outPrefix, "_feature_overlaps.txt")
  if(i){
    enOutputFile = paste0(outPrefix, "_feature_overlaps", i, ".txt")    
  }

  cat(paste0("#\tmin_dist=", enMinDist, "\tmax_dist=", ifelse(is.null(enMaxDist), "whole_range", enMaxDist), "\ttrans=", enTrans, "\n"), file = enOutputFile) 
  suppressWarnings(write.table(enrichments, quote=F, row.names= T, col.names=T, file= enOutputFile, append = T))
  setwd(outDir)
  if(!moveToFolder("feature_overlaps", "enrichment_data")){
     stop("Couldn't move feature_overlaps files to enrichment folder")
 }
}


if (!featuresOnly) { 
  setwd(curDir)
}
message("All done!\n")
