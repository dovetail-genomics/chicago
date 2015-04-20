library(argparser)
library(Chicago)

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
                 help = "File format for writing out peaks: one or more of the following: seqMonk,interBed,washU (comma-separated)", 
                 default = "washU")
p = add.argument(p, arg="--export-order", help = "Should the results be ordered by \"score\" or genomic \"position\"?", 
                 default = "position")

p = add.argument(p, arg="--feature-file", 
                 help = "File with genomic feature coordinates for computing peaks' enrichment over these feature.", 
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
outPrefix = opts[["<output-prefix>"]]

settingsFile = opts[["settings-file"]]
designDir = opts[["design-dir"]]

printMemory = opts[["print-memory"]]

cutoff = opts[["cutoff"]]
exportFormat = ifelse(is.na(opts[["export-format"]]), NA, strsplit(opts[["export-format"]], "\\,")[[1]])
exportOrder = opts[["export-order"]]

featureFile = opts[["feature-file"]]
featureList = opts[["feature-list"]]

isRda = opts[["rda"]]
isDF = opts[["save-df-only"]]

proxLim = opts[["examples-prox-dist"]]
plotFull = opts[["examples-full-range"]]

featDistLim = opts[["feat-max-dist"]]

outDir = ifelse(opts[["output-dir"]]=="<output-prefix>", outPrefix, opts[["output-prefix"]])

if(is.na(settingsFile)){
  settingsFile = NULL
}

if (!all(exportFormat %in% c("seqMonk", "interBed", "washU"))) {
  stop("--exportFormat must be either seqMonk, interBed or washU (or a comma-separated combination of these)\n")
}
if (!exportOrder %in% c("position", "score")) {
  stop("--exportOrder must be either position (default) or score\n")
}

message("Setting the experiment...\n")
cd = setExperiment(designDir = designDir, settingsFile = settingsFile)


if(length(inputFiles)>1){
  message("Reading input files...\n")
  cd = readAndMerge(inputFiles, cd)
}else{
  cd = readSample(inputFiles, cd)
}

if(system(paste("mkdir -p", outDir))){
  stop(paste("Couldn't create folder", outDir, "\n"))
}

curDir = getwd()
setwd(outDir)

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
if(system("mkdir -p data")){
  stop("Couldn't create folder data\n")  
}
if(system("mkdir -p diag_plots")){
  stop("Couldn't create folder diag_plots\n")
}
if(system("mkdir -p examples")){
  stop("Couldn't create folder examples\n")      
}
if(system("mkdir -p enrichment-data")){
  stop("Couldn't create folder enrichment-data\n")      
}

system("mv *.txt data/")

if("interBed" %in% exportFormat){
  system("mv *.ibed data/")
}

if(isRda){
  system("mv *.RDa data/")
}else{
  system("mv *.Rds data/")  
}

system("mv *xamples.pdf examples/")
system("mv *.pdf diag_plots/")       

if (!is.na(featureFile) | !is.na(featureList)){
  message("\n\nComputing enrichment for features...\n")
#   **TODO - UPDATE THIS**
#   peakEnrichment4Features(x1=x[!is.na(x$distSign),], score=12, sample_number=100, no_bins=100, 
#                           colname_score="score",folder=featureFolder, position_otherEnd_folder=fileDir, 
#                           list_frag=files, filterB2B=TRUE,
#                           colname_dist="distSign", beyond_dist=0, before_dist=1000000,
#                           plot_name=paste0(outprefix, "_feature_overlaps_upto_1M.pdf"))
#   
#   write.table(resTable, quote=F, row.names= , col.names=T, file=   paste0(outprefix, "_feature_overlaps_upto1M.txt")) 
#   system("mv *feature_overlaps*.* enrichment-data/")
}

setwd(curDir)
message("All done!\n")

