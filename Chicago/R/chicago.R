##This line prevents data.table syntax being flagged up as a NOTE in R CMD check
##("no visible binding for global variable")
utils::globalVariables(c("distSign", "isBait2bait", "otherEndID", "transLength",
                         "distbin2", "V1", "reads", "distbin3", "udbin2", "iTempVar",
                         "bin_reads", "otherEndID", "Bmean", "s_j", "s_i", "chr", "baitID", "ntotpb",
                         "distbin", "binCol", "ntot", "geomean", "bbm","s_iv","shape","rate","shapefit",
                         "ratefit","s_ivfilt", "J", "N", "refBinMean","tlb","tblb","baitChr",
                         "otherEndChr","log.p","Tmean","log.w","log.q","score","otherEndLen","nperbait",
                         "isAdjacent","isAllB2BProx", "samplefilename"))

## Settings functions -----------------

defaultSettings <- function()
{
  list(
    rmapfile= NA,
    baitmapfile= NA,
    nperbinfile = NA,
    nbaitsperbinfile = NA,
    proxOEfile = NA,
    Ncol = "N",
    baitmapFragIDcol=4,
    baitmapGeneIDcol=5,
    maxLBrownEst = 1.5e6,
    minFragLen = 150,
    maxFragLen = 40000,
    minNPerBait = 250,
    binsize=20000,
    removeAdjacent = TRUE,
    adjBait2bait=TRUE,
    tlb.filterTopPercent=0.01, 
    tlb.minProxOEPerBin=50000, 
    tlb.minProxB2BPerBin=2500,
    techNoise.minBaitsPerBin=1000, 
    brownianNoise.samples=5,
    brownianNoise.subset=1000,
    brownianNoise.seed=NA,
    baitIDcol = "baitID",
    otherEndIDcol = "otherEndID",
    otherEndLencol = "otherEndLen", 
    distcol = "distSign",
    weightAlpha = 34.1157346557331, 
    weightBeta = -2.58688050486759,
    weightGamma = -17.1347845819659,
    weightDelta = -7.07609245521541
  )
}

### This is where we now set the defaults
### Order of priority:
### settings override settings from settingsFile
### both override def.settings
setExperiment = function(designDir="", settings=list(), settingsFile=NULL,  
 def.settings=defaultSettings()){
  
  if(designDir == "" & identical(settings, list()) & is.null(settingsFile) & identical(def.settings, defaultSettings()))
  {
    stop("Design not specified. Please specify a design directory or design files.")
  }
  
  def.settings = .updateSettings(designDir, settings, settingsFile, def.settings) 
  cd = chicagoData(x=data.table(), params=list(), settings=def.settings)
}

modifySettings = function(cd, designDir=NULL, settings=list(), settingsFile=NULL){
  
  message("Warning: settings are not checked for consistency with the original ones.")
  cd@settings = .updateSettings(designDir, settings, settingsFile, def.settings=cd@settings, updateDesign=TRUE)   

  ##test validity of new object
  #validObject(cd)
  
  cd
}

.updateSettings = function(designDir, settings, settingsFile, def.settings, updateDesign=FALSE){
  modSettings = vector("list")
  
  if(!is.null(settingsFile)){
    
    message(paste0("Reading custom experimental settings from ", settingsFile, "..."))
    
    # http://stackoverflow.com/questions/6602881/text-file-to-list-in-r
    sf <- scan(settingsFile, what="", sep="\n")
    modSettings <- strsplit(sf, "[[:space:]]+")
    names(modSettings) <- sapply(modSettings, `[[`, 1)
    modSettings <- lapply(modSettings, `[`, -1)
    
    # convert numerically defined settings to numbers
    # do the same for logical settings
    # suppressing "NAs introduced by coercion" warnings
    suppressWarnings({
      modSettings <- lapply(modSettings, function(s){
        num_s = as.numeric(s);
        if(!is.na(num_s)) { return (num_s) };
        bool_s = as.logical(s);
        if(!is.na(bool_s)) { return(bool_s) };
        s 
      })
    })
  }

  for (s in names(settings)){
    modSettings[[s]] = settings[[s]]
  }
  
  for (s in names(modSettings)){
    def.settings[[s]] = modSettings[[s]]
  }
  
  if(is.na(def.settings[["baitmapfile"]]) | (updateDesign & !is.null(designDir))){
    def.settings[["baitmapfile"]] = locateFile("<baitmapfile>.baitmap", designDir, "\\.baitmap$")
  }else{
    if (!file.exists(def.settings[["baitmapfile"]])){
      stop(paste("No baitmap file found at the specified location", def.settings[["baitmapfile"]]))
    }
  }
  
  if(is.na(def.settings[["rmapfile"]]) | (updateDesign & !is.null(designDir))){
    def.settings[["rmapfile"]] = locateFile("<rmapfile>.rmap", designDir, "\\.rmap$")
  }else{
    if (!file.exists(def.settings[["rmapfile"]])){
      stop(paste("No rmap file found at the specified location", def.settings[["rmapfile"]]))
    }
  }
  
  if(is.na(def.settings[["nperbinfile"]]) | (updateDesign & !is.null(designDir))){
    def.settings[["nperbinfile"]] = locateFile("<nperbinfile>.npb", designDir, "\\.npb$")
  }else{
    if (!file.exists(def.settings[["nperbinfile"]])){
      stop(paste("No nperbin file found at the specified location", def.settings[["nperbinfile"]]))
    }
  }
  
  if(is.na(def.settings[["nbaitsperbinfile"]]) | (updateDesign & !is.null(designDir))){
    def.settings[["nbaitsperbinfile"]] = locateFile("<nbaitsperbinfile>.nbpb", designDir, "\\.nbpb$")
  }else{
    if (!file.exists(def.settings[["nbaitsperbinfile"]])){
      stop(paste("No nbaitsperbin file found at the specified location", def.settings[["nbaitsperbinfile"]]))
    }
  }
  
  if(is.na(def.settings[["proxOEfile"]]) | (updateDesign & !is.null(designDir))){
    def.settings[["proxOEfile"]] = locateFile("<proxOEfile>.poe", designDir, "\\.poe$")
  }else{
    if (!file.exists(def.settings[["proxOEfile"]])){
      stop(paste("No proxOE file found at the specified location", def.settings[["proxOEfile"]]))
    }
  }
 
  if (def.settings[["maxLBrownEst"]] %% def.settings[["binsize"]]){
      message("Warning: the supplied maxLBrownEst=", def.settings[["maxLBrownEst"]], " is not a multiple of binsize=", def.settings[["binsize"]], ". Will be truncated to the nearest bin boundary.\n")
      def.settings[["maxLBrownEst"]] = floor(def.settings[["maxLBrownEst"]]/def.settings[["binsize"]])*def.settings[["binsize"]]     
  }
 
  message("Checking the design files...")

  bm = .readBaitmap(def.settings)
  if(ncol(bm)<max(c(def.settings[["baitmapFragIDcol"]], def.settings[["baitmapGeneIDcol"]]))){
    stop("There are fewer columns in the baitmapfile than expected. Check that this file lists the genomic coordinates as well as both the IDs and names for each baited fragment,
and that the corresponding columns are specified in baitmapFragIDcol and baitmapGeneIDcol, respectively.")
  }
  rmap = .readRmap(def.settings)
  if (ncol(rmap)<4){
    stop("There are fewer columns in the rmap file than expected. This file should have 4 columns, listing the genomic coordinates and IDs for each restriction fragment.")
  }
  if (ncol(rmap)>4){
    stop("There are more columns in the rmap file than expected. This file should have 4 columns, listing the genomic coordinates and IDs for each restriction fragment. Check that rmapfile and baitmap files aren't swapped round\n")
  }
  
  if (any(duplicated(bm[[def.settings[["baitmapFragIDcol"]]]]))){
    stop(paste("Duplicated fragment IDs found in baitmapfile (listed below). 
               Check that the baitmapFragIDcol (usually 4) and baitmapGeneIDcol 
               (usually 5) aren't swapped round:\n", 
               paste(bm[[def.settings[["baitmapFragIDcol"]]]][duplicated(bm[[def.settings[["baitmapFragIDcol"]]]])], 
                     collapse=",")))
  }

  if (any(duplicated(rmap[[4]]))){
    stop(paste("Duplicated fragment IDs found in rmapfile (listed below):\n", 
               paste(rmap[[4]][duplicated(rmap[[4]])], collapse=",")))
  }

  setkeyv(rmap, names(rmap))
  setkeyv(bm, names(bm[c(1:3, def.settings[["baitmapFragIDcol"]])]))
  if (nrow(merge(rmap, bm))!=nrow(bm)){ 
     message("Error: Some entries of baitmap are not present in rmap (listed below).\n") 
     print(bm[!merge(bm, rmap)])
     stop()
  }

  message("Reading the settings from NPB file header...")
  header = readLines(def.settings[["nperbinfile"]], n=1)
  params = sapply(sapply(strsplit(header, "\t")[[1]],function(x)strsplit(x,"=")[[1]]), function(x)x[2])
  params = params[2:length(params)]
  names(params) = gsub("(\\S+)=.+", "\\1", names(params))
  minsize = as.numeric(params[["minFragLen"]])
  if (!is.null(settings[["minFragLen"]])){
    if (settings[["minFragLen"]] != minsize){
      stop("The minFragLen in the .npb file header is not equal to the custom-defined setting. Amend either setting (and if needed, generate a new .npb file) before running the analysis\n")
    }
  }
  if (!is.null(modSettings[["minFragLen"]])){
    if (modSettings[["minFragLen"]] != minsize){
      stop("The minFragLen in the .npb file header is not equal to the same setting defined in the settings file. Amend either setting (and if needed, generate a new .npb file) before running the analysis\n")
    }
  }
  if(def.settings[["minFragLen"]] != minsize){
    message("Amending the default maxFragLen setting from ", def.settings[["minFragLen"]], " to ", minsize, " specified in .npb file header.")
    def.settings[["minFragLen"]] = minsize
  }
  
  maxsize = as.numeric(params[["maxFragLen"]])
  if (!is.null(settings[["maxFragLen"]])){
    if (settings[["maxFragLen"]] != maxsize){
      stop("The maxFragLen in the .npb file header is not equal to the custom-defined setting. Amend either setting (and if needed, generate a new .npb file) before running the analysis\n")
    }
  }
  if (!is.null(modSettings[["maxFragLen"]])){
    if (modSettings[["maxFragLen"]] != maxsize){
      stop("The maxFragLen in the .npb file header is not equal to the same setting defined in the settings file. Amend either setting (and if needed, generate a new .npb file) before running the analysis\n")
    }
  }
  if(def.settings[["maxFragLen"]] != maxsize){
    message("Amending the default maxFragLen setting from ", def.settings[["maxFragLen"]], " to ", maxsize, " specified in .npb file header.")
    def.settings[["maxFragLen"]] = maxsize
  }
  
  maxl = as.numeric(params[["maxLBrownEst"]])
  if (!is.null(settings[["maxLBrownEst"]])){
    if (settings[["maxLBrownEst"]] != maxl){
      stop("The maxLBrownEst in the .npb file header is not equal to the custom-defined setting. Amend either setting (and if needed, generate a new .npb file) before running the analysis\n")
    }
  }
  if (!is.null(modSettings[["maxLBrownEst"]])){
    if (modSettings[["maxLBrownEst"]] != maxl){
      stop("The maxLBrownEst in the .npb file header is not equal to the same setting defined in the settings file. Amend either setting (and if needed, generate a new .npb file) before running the analysis\n")
    }
  }
  if(def.settings[["maxLBrownEst"]] != maxl){
    message("Amending the default maxLBrownEst setting from ", def.settings[["maxLBrownEst"]], " to ", maxl, " specified in .npb file header.")
    def.settings[["maxLBrownEst"]] = maxl
  }
  
  binsz = as.numeric(params[["binsize"]]) 
  if (!is.null(settings[["binsize"]])){
    if (settings[["binsize"]] != binsz){
      stop("The binsize in the .npb file header is not equal to the custom-defined setting. Amend either setting (and if needed, generate a new .npb file) before running the analysis\n")
    }
  }
  if (!is.null(modSettings[["binsize"]])){
    if (modSettings[["binsize"]] != binsz){
      stop("The binsize in the .npb file header is not equal to the same setting defined in the settings file. Amend either setting (and if needed, generate a new .npb file) before running the analysis\n")
    }
  }
  if(def.settings[["binsize"]] != binsz){
    message("Amending the default binsize setting from ", def.settings[["binsize"]], " to ", binsz, " specified in .npb file header.")
    def.settings[["binsize"]] = binsz
  }

  ra = params[["removeAdjacent"]]
  if (!is.null(settings[["removeAdjacent"]])){
    if ((settings[["removeAdjacent"]] == TRUE & toupper(ra)!="TRUE") | (settings[["removeAdjacent"]] == FALSE & toupper(ra)!="FALSE")){
      stop("The removeAdjacent setting in the .npb file header is not equal to the custom-defined setting. Amend either setting (and if needed, generate a new .npb file) before running the analysis\n")
    }
  }
  if (!is.null(modSettings[["removeAdjacent"]])){
    if ((modSettings[["removeAdjacent"]] == TRUE & toupper(ra)!="TRUE") | (modSettings[["removeAdjacent"]] == FALSE & toupper(ra)!="FALSE")){
      stop("The removeAdjacent in the .npb file header is not equal to the same setting defined in the settings file. Amend either setting (and if needed, generate a new .npb file) before running the analysis\n")
    }
  }
  if(def.settings[["removeAdjacent"]] == TRUE & toupper(ra)=="FALSE"){
    message("Amending the default removeAdjacent setting from ", def.settings[["removeAdjacent"]], " to FALSE specified in .npb file header.")
    def.settings[["removeAdjacent"]] = FALSE
  }
  if(def.settings[["removeAdjacent"]] == FALSE & toupper(ra)=="TRUE"){
    message("Amending the default removeAdjacent setting from ", def.settings[["removeAdjacent"]], " to TRUE specified in .npb file header.")
    def.settings[["removeAdjacent"]] = TRUE
  }


  if(basename(params[["rmapfile"]]) != basename(def.settings[["rmapfile"]])){
    stop("The basename of .rmap file used for generating the .npb file (according to the .npb header) and the one defined in experiment settings do not match. Amend either setting (and if needed, generate a new .npb file or rename your rmap file if sure it's for the same design) before running the analysis\n")
  }
  
  if(basename(params[["baitmapfile"]]) != basename(def.settings[["baitmapfile"]])){
    warning("The basename of .baitmap file used for generating the .npb file (according to the .npb header) and the one defined in experiment settings do not match. Please check this is intended.\n")
  }
  
  if (params[["removeb2b"]]!="True"){
    stop("The .npb file must be generated with removeb2b==True. Please generate a new file.\n")
  }

  message("Checking the integrity of the NPB file...")
  npb = .readNPBfile(def.settings)
  if(!all(bm[[4]] %in% npb[[1]]) | !all(npb[[1]] %in% bm[[4]]) ){
	stop("The lists of baits in baitmap and NPB files don't overlap fully") 
  }  
    
  message("Checking the integrity of the NBPB file...")
  nbpb = .readNbaitsPBfile(def.settings)
  if(!all(rmap[[4]] %in% nbpb[[1]]) | !all(nbpb[[1]] %in% rmap[[4]]) ){
        stop("The lists of restriction fragments in rmap and NBPB files don't overlap fully")
  }

  def.settings
}

.checkForIncompatibilities <- function(cd)
{
  ##Check for 1.1.4 or before
  if(cd@settings$tlb.minProxOEPerBin == 1000 | cd@settings$tlb.minProxB2BPerBin == 100)
  {
    message("WARNING: It looks like cd is an old chicagoData (version 1.1.4 or earlier).")
    message("If so, please update some of the settings - see the Chicago 1.1.5 entry in news(package='Chicago').")
    warning("tlb settings match Chicago 1.1.4 and before.")
  }
}

## Read-in functions ----------------

readSample = function(file, cd){
  
  message(paste("Reading", file))
  
  ##check for comment line at beginning of file
  testLine <- readLines(file, 1)
  if(substr(testLine, 1, 1) == "#")
  {
    x = fread(file, skip = 1L)
  }else{
    x = fread(file)
  }
    
  message("Processing input...")
  
  s = cd@settings
  
  if ( (! s$baitIDcol %in% names(x)) |   (! s$otherEndIDcol %in% names(x))  
       |  (! s$Ncol %in% names(x)) |  (! s$otherEndLencol %in% names(x))  | 
         (! s$distcol %in% names(x))){
    stop("Named columns baitIDcol = ", s$baitIDcol, ", otherEndIDcol = ", s$otherEndIDcol,
         ", Ncol = ", s$Ncol, ", otherEndLencol = ", s$otherEndLencol, " and distcol = ", s$distcol, 
         " must be present in the input file. Change these global parameters if names do not match\n")
  }
  
  setnames(x, s$baitIDcol, "baitID")
  setnames(x, s$otherEndIDcol, "otherEndID")
  setnames(x, s$Ncol, "N")
  setnames(x, s$distcol, "distSign")
  setnames(x, s$otherEndLencol, "otherEndLen")

  ## check that the file looks like it is produced with the specified array design
  bm = .readBaitmap(s)
  if (!all(x$baitID %in% bm[[s$baitmapFragIDcol]])){
    stop("Some entries of the input file have baitIDs not in the baitmapfile. Check that the specified design files are correct.\n")
  }
  rmap = .readRmap(s)
  if (!all(x$otherEndIDcol %in% rmap[[4]])){
    stop("Some entries of the input file have otherEndIDs not in the rmapfile. Check that the specified design files are correct.\n")
  }
  
  xlen = nrow(x)
  x = x[otherEndLen %between% c(s$minFragLen,s$maxFragLen)]
  message("minFragLen = ", s$minFragLen, " maxFragLen = ", s$maxFragLen)
  message("Filtered out ", xlen-nrow(x), " interactions involving other ends < minFragLen or > maxFragLen.")
  if(nrow(x) == 0) stop("All interactions have been filtered out.")
  
  setkey(x, baitID)

  ## remove self-ligation events
  oldlen = nrow(x)
  x = x[distSign!=0 | is.na(distSign)]
  if(oldlen>nrow(x)){
     message("Filtered out ", oldlen-nrow(x), " self-ligation events.\n")
  }
  
  ## remove baits that have no observations within the proximal range
  baitlen = length(unique(x$baitID)) 
  x = x[, nperbait:=sum(N), by=baitID]
  x = x[nperbait>=s$minNPerBait]
  message("minNPerBait = ", s$minNPerBait)
  message("Filtered out ", baitlen-length(unique(x$baitID)), " baits with < minNPerBait reads.\n")  
  set(x, NULL , "nperbait", NULL) # fast remove data.table column  
  if(nrow(x) == 0) stop("All interactions have been filtered out.")
  
  ## remove adjacent pairs
  if(s$removeAdjacent){
    x[, isAdjacent:=abs(baitID-otherEndID)==1, by=baitID]
    x = x[isAdjacent==FALSE]
    set(x, NULL, "isAdjacent", NULL)
    message("Removed interactions with fragments adjacent to baits.")
    if(nrow(x) == 0) stop("All interactions have been filtered out.")
  }

  ##remove baits without proximal non-bait2bait interactions
  baitlen = length(unique(x$baitID)) 
  x[, isBait2bait := FALSE]
  x[wb2b(otherEndID, s), isBait2bait:= TRUE] 
  x[, isAllB2BProx:={
    prox = abs(distSign)<s$maxLBrownEst & !is.na(distSign)
    if(!length(prox)){  TRUE  }
    else{ all(isBait2bait[prox]) }
  }, by=baitID]
  x = x[isAllB2BProx==FALSE]
  set(x, NULL, "isAllB2BProx", NULL)
  if(nrow(x) == 0) stop("All interactions have been filtered out.")
  
  message("Filtered out ", baitlen-length(unique(x$baitID)), " baits without proximal non-Bait2bait interactions\n")  

  chicagoData(x=x, settings = cd@settings, params=list()) # creating this object from scracth, so a single cd object could be passed to readSample as input for multiple replicates
}

readAndMerge = function(files, cd, ...){
  if (length(files)==1) readSample(files, cd)
  else mergeSamples(lapply(files, readSample, cd), ...)
}

getSkOnly <- function(files, cd)
{
  N <- length(files)
  if(any(!file.exists(files))) {stop("Could not find files: ", paste(files[!file.exists(files)], collapse=", "))}
  
  cdList <- vector("list", N)
  for(i in 1:N)
  {
    cdList[[i]] <- cd
    cdList[[i]] <- readSample(file = files[i],cd = cdList[[i]])
    sel <- !is.na(cdList[[i]]@x$distSign) & (abs(cdList[[i]]@x$distSign) < defaultSettings()$maxLBrownEst)
    cdList[[i]]@x <- cdList[[i]]@x[sel]
  }
  
  cdMerge <- mergeSamples(cdList)
  sk <- cdMerge@params$s_k
  names(sk) <- files
  sk
}

mergeSamples = function(cdl, normalise = TRUE, NcolOut="N", NcolNormPrefix="NNorm", 
                        mergeMethod=c("weightedMean", "mean")[1], repNormCounts = (mergeMethod=="mean") ){
  
  # Now takes a list of chicagoData classes as input and returns a single chicagoData class
  
  # If mergeMethod == "weightedMean", NcolOut is the weighted mean of the sample-wise counts
  # adjusted by the samples' respective scaling factors s_k
  # If mergeMethod == "normMean", sample-specific counts are first normalised by dividing by s_k
  # and NcolOut is computed as the mean of these normalised counts.
  
  if (! mergeMethod %in% c("weightedMean", "mean")){
    stop ("Unknown mergeMethod.\n")
  }
    
  attr = vector("list")
  
  # It's very inefficient to merge samples using all columns as keys. 
  # But if we just merge by the primary keys baitIDcol and otherEndIDcol,
  # the columns containing all other info (distance, length 
  # and anything else that may be present) will be duplicated, 
  # with NAs for all samples in which they are not detected. Parsing this is tedious as well.
  # So - still merge only by the primary keys 
  # but at the same time create a "backbone" of all the remaining fields for each combination 
  # of the primary keys, and then populate it with the counts for each sample.
  
  # Rename read count columns in each sample so they are uniquely identifiable 
    
  x = as.dataTableList(cdl)
  
  for (i in 1:length(x)){
    
    if(i>1){
      if(!identical(cdl[[i]]@settings, cdl[[i-1]]@settings)){
        stop("All samples to merge should have identical experiment settings")    
      }
    }
    
    Ncol = grep("^N$", names(x[[i]]), value=TRUE)
    if (!length(Ncol)){ 
      #In case names have already been changed to N.<k> - this being data.table
      Ncol = grep("^N\\.", names(x[[i]]), value=TRUE) 
      if (length(Ncol)==1){
        warning("Could not find column name \"N\" in element ", i, " of the input list. Using \"", Ncol, "\" instead\n")
      }
      else {
        stop(paste("Could not find column name \"N\" in element ", i, " of the input list.\n"))
      }
    }
    setnames(x[[i]], Ncol, paste0("N.", i))
        
    attr[[i]] = copy(x[[i]]) # So far I see no better way of doing this
    set(attr[[i]], NULL, paste0("N.", i), NULL) # remove the scores column
    
    remCols = names(x[[i]])[!names(x[[i]]) %in% c("baitID", "otherEndID", paste0("N.", i))]
    
    for (rc in remCols){
      set(x[[i]], NULL, rc, NULL)      
    }
    setkey(x[[i]], baitID, otherEndID)
    
  }
  
  attr = rbindlist(attr) 
  
  setkey(attr, baitID, otherEndID)
  # remove duplicates by key
  attr = unique(attr, by=key(attr))
  
  message("Merging samples...")
  
  xmerge = Reduce(function(...) merge(..., by=c("baitID", "otherEndID"), all=TRUE), x)
  
  setkey(xmerge, baitID, otherEndID)
  xmerge = attr[xmerge]
    
  for (i in 1:length(x)){
    iNcol = paste0("N.", i)    
    set(xmerge, which(is.na(xmerge[[iNcol]])), iNcol, 0) # an ugly but the most efficient way to replace NA's with zeros...
  }

  s = cdl[[1]]@settings
  
  message("Computing merged scores...\n")
  
  if (normalise){

    #   whichN = which(names(xmerge) %in% paste0("N.", 1:length(x))
    Nnames = names(xmerge) [ names(xmerge) %in% paste0("N.", 1:length(x))]
    
    s_ks = .getSampleScalingFactors(xmerge, s) 
            
    if(mergeMethod=="weightedMean"){
      # Sorry about this - but that's the only way to make it efficient...
      # Essentially, it's generating this (for the case of two replicates),
      # assuming NcolOut = "N"
      # N:=round((N.1*s_ks["N.1"]+N.2*s_ks["N.2"]+N.3*s_ks["N.3"])/sum(s_ks))
      xmerge[, eval(parse(text=paste0(NcolOut, ":=round((", paste0(Nnames, "*s_ks[\"", Nnames, "\"]", collapse="+"), ")/sum(s_ks))")))]
    }
    else{ # arithmetic mean
      # N:=round((N.1/s_ks["N.1"]+N.2/s_ks["N.2"]+N.3/s_ks["N.3"])/length(Nnames)
      xmerge[, eval(parse(text=paste0(NcolOut, ":=round((", paste0(Nnames, "/s_ks[\"", Nnames, "\"]", collapse="+"), ")/length(Nnames))")))]    
    }
    
    if(repNormCounts){
      # compute scaled per-sample quantities
      for (i in 1:length(x)){
        # e.g., NNorm.1:=round(N.1/s_ks["N.1"])
        xmerge[, eval(parse(text=paste0(NcolNormPrefix, ".", i,":=round(", Nnames[i], "/s_ks[\"", Nnames[i], "\"])")))]
      }
    }    
  }
  else{
    # N:=round((N.1+N.2+N.3)/length(Nnames))
    xmerge[, eval(parse(text=paste0(NcolOut, ":=round((", paste0(Nnames, collapse="+"), ")/length(Nnames))")))]
  }
  
  # Don't completely "cancel" interactions for which we have observed at least one read somewhere
  set(xmerge, which(!xmerge[[NcolOut]]), NcolOut, 1)
  
  chicagoData(x=xmerge, params=list(s_k=s_ks), settings=s)  # Sic! Now returns a chicagoData object, not a single table with attributes!
}

.getSampleScalingFactors = function(xs, s){
  
  # UPDATE: in computeNNorm=F (returnSKonly=T), will just return the s_k vector !!
  
  # Compute normalisation factors s_k based on the median number of reads per bait
  # within the distance range (0; maxLBrownEst)  
  # The normalisation factors themselves will be written as *attributes* of the output data frame. 
  # If compute NNorm==T, normalise sample-wise counts by dividing by 
  # the respective s_k, writing the results into NcolNormPrefix.<#sample> column. 
  
  # s is the current chicagoData object's settings list
  
  message("Computing sample scaling factors...\n")
  
  if (!is.data.frame(xs)){ 
    stop("xs must be a data table. If starting from a list of separate samples, use mergeSamples instead\n")
  }
  
  Ncols = grep("^N\\.(\\d+)", names(xs), value=TRUE) ##matches "N.[integer]"
  ns = as.numeric(gsub("N\\.(\\d+)", "\\1", Ncols)) ##collects [integer]s from the above
  n = max(ns)

  # message("n = ", n)
  
  # Only using the distance range  (0; maxLBrownEst) for normalisation
  # Saving the full table for output

  xs = xs[abs(distSign)<s$maxLBrownEst]
  
  setkey(xs, baitID)  
  Ncols = paste0("N.", 1:n)

  # Get the total number of other ends within the distance range (0; maxLBrownEst) for each bait 
  npb = .readNPBfile(s)   
  
  message("Computing normalisation scores...")
  whichN = names(npb)[2:ncol(npb)]
  npb = npb[, ntotpb := eval(parse(text = paste0(whichN, collapse="+"))) ] 
  for (nm in names(npb)[!names(npb)%in%c("baitID", "ntotpb")]){
    set(npb, NULL, nm, NULL)
  }
  setkey(npb, baitID)
  setkey(xs, baitID)
  xs = npb[xs]
  
  baitMeans = xs[, ntotpb[1],by=baitID]
  s_kjcols = paste0("s_", 1:n, "j")
  for (k in 1:n){
    ##for each sample, for each bait: take number of reads in proximal region, divide by number of fragments in proximal region (ntotpb)
    baitMeans[[ s_kjcols[k] ]] = xs[,sum(get(Ncols[k]))/ntotpb[1],by=baitID]$V1
  }
  
  baitMeans$geo_mean = apply(baitMeans[, s_kjcols, with=FALSE], 1,
                             function(x)if(any(x==0)){NA}else{geo_mean(x)})
  s_k = vector("numeric")
  for (k in 1:n){
    s_k [ k ] = median(baitMeans[[ s_kjcols[k] ]] / baitMeans [[  "geo_mean" ]], na.rm=TRUE)
  }

  names(s_k) = Ncols
  s_k
  
}

##Pipeline ------------

chicagoPipeline <- function(cd, outprefix=NULL, printMemory=FALSE)
{
  .checkForIncompatibilities(cd)
  
  message("\n*** Running normaliseBaits...\n")
  cd = normaliseBaits(cd)
  
  if(printMemory){
    print(gc(reset=TRUE))
  }
  
  message("\n*** Running normaliseOtherEnds...\n")
  cd = normaliseOtherEnds(cd,
                          outfile=ifnotnull(outprefix, paste0(outprefix, "_oeNorm.pdf"))
  )
  
  if(printMemory){
    print(gc(reset=TRUE))
  }
  
  message("\n*** Running estimateTechnicalNoise...\n")
  cd = estimateTechnicalNoise(cd,
                              outfile=ifnotnull(outprefix, paste0(outprefix, "_techNoise.pdf"))
  )
  
  if(printMemory){
    print(gc(reset=TRUE))
  }
  
  message("\n*** Running estimateDistFun...\n")
  
  ### Note that f is saved in cd@params
  cd = estimateDistFun(cd,
                       outfile=ifnotnull(outprefix, paste0(outprefix, "_distFun.pdf"))
  )
  
  if(printMemory){
    print(gc(reset=TRUE))
  }
  
  ### Note that f is saved as cd@params$f and  
  ### subset is saved as cd@settings$brownianNoise.subset
  message("\n*** Running estimateBrownianComponent...\n")
  cd = estimateBrownianComponent(cd)
  
  if(printMemory){
    print(gc(reset=TRUE))
  }  
  
  message("\n*** Running getPvals...\n")
  cd = getPvals(cd)
  
  if(printMemory){
    print(gc(reset=TRUE))
  }  
  
  message("\n*** Running getScores...\n")
  cd = getScores(cd)
  
  if(printMemory){
    print(gc(reset=TRUE))
  }  
  
  cd
}

##Pipeline functions ------------

normaliseBaits = function(cd, normNcol="NNb", shrink=FALSE, plot=TRUE, outfile=NULL, debug=FALSE){
  message("Normalising baits...")
  
  ##test to see if s_j, refBinMean columns exist, warn & delete if so
  replacedCols <- c("s_j", "refBinMean")
  sel <- replacedCols %in% colnames(cd@x)
  if(any(sel))
  {
    warning("Columns will be overwritten: ", paste(replacedCols[sel], collapse=", ")) 
    set(cd@x, j=replacedCols[sel], value=NULL) ##delete these columns
  }

  adjBait2bait = cd@settings$adjBait2bait
  
  # NON-URGENT TODO: An even more memory-efficient way of doing this would be to have .normaliseFragmentSets assign
  # the s_j, distbin and refBinMean columns by reference!
  
  if(debug){
    ##returns sbbm and not Chicago object! 
    .normaliseFragmentSets(x=cd@x, s=cd@settings, npb=.readNPBfile(s=cd@settings), viewpoint="bait", idcol="baitID", Ncol="N", adjBait2bait=adjBait2bait, shrink=shrink, refExcludeSuffix=NULL, plot=plot, outfile=outfile, debug=TRUE)
  }
  else{
    cd@x = .normaliseFragmentSets(x=cd@x, s=cd@settings, npb=.readNPBfile(s=cd@settings), viewpoint="bait", idcol="baitID", Ncol="N", adjBait2bait=adjBait2bait, shrink=shrink, refExcludeSuffix=NULL, plot=plot, outfile=outfile, debug=FALSE)
    
  }
  
  # sort by baitID, otherEndID and move distbin column to the end of the table 
  cd@x[, (normNcol):= pmax(1, round(N/s_j)) ] # do not completely "cancel" interactions that have one read 

  setkey(cd@x, baitID, otherEndID)

  othercols = names(cd@x)[!names(cd@x)%in% c("distbin", normNcol)]

  setcolorder(cd@x, c(othercols, "distbin", normNcol))
  
  cd
  
}

normaliseOtherEnds = function(cd, Ncol="NNb", normNcol="NNboe", plot=TRUE, outfile=NULL){
  
  # NON-URGENT TODO: instead of looking at bins defined by trans-counts & isB2B, look at bins defined by 
  # fragment length, mappability, GC content and isB2B; 
  # In this case, the call to .addTLB() will be substituted with something like .addLMGB() [to be written] 

  cd@x = .addTLB(cd)
  x = cd@x[abs(distSign)<=cd@settings$maxLBrownEst & is.na(distSign)==FALSE]
  
  message("Computing total bait counts...")
  nbpb = .readNbaitsPBfile(s=cd@settings)
  
  setkey(nbpb, otherEndID)
  setkey(x, otherEndID)
  
  nbpb = nbpb[x]
  
  # Compute the sums of observed baits per bin for pools of other ends - 
  # NB: we need to some only once for each other end, but each other end is present
  # more than once for each tlb
  setkeyv(nbpb, c("otherEndID", "distbin"))
  nbpb = unique(nbpb, by=key(nbpb))
  
  nbpbSum = nbpb[, sum(get(paste0("bin",as.integer(distbin)))), by=c("tlb","distbin")]
  setkeyv(nbpbSum, c("tlb", "distbin"))
  setkeyv(x, c("tlb", "distbin"))
  x = nbpbSum[x]
  setnames(x, "V1", "ntot")
  
  message("Computing scaling factors...")
  
  x = .normaliseFragmentSets(x, s = cd@settings, viewpoint="otherEnd", idcol="tlb", 
                             Ncol=Ncol, npb = NULL, shrink=FALSE, adjBait2bait=FALSE, refExcludeSuffix="B2B")
  
  setkey(x, tlb)
  x = unique(x, by=key(x)) # we don't need any other info than s_i for each tlb and it's one per tlb
  set(x, NULL, "NNb", NULL)
  set(x, NULL, "distbin", NULL)
  set(x, NULL, "ntot", NULL)
  
  if(plot){
    if (!is.null(outfile)){ pdf(outfile)}
    with(x,
         barplot(s_i, names.arg=tlb, col=sapply(tlb, function(x)ifelse(length(grep("B2B",x)), "darkblue", "red")),
                    xlab="tlb", ylab="s_i", main = "Brownian OE factors (s_i) estimated per OE pool")
         )
    legend("topleft", legend=c("non-B2B", "B2B"),fill=c("red", "darkblue"))
    if (!is.null(outfile)){ dev.off()}
  }

  message("Computing normalised counts...")
    
  setkey(cd@x, tlb)
  
  if("s_i" %in% colnames(cd@x))
  {
    set(cd@x, NULL, "s_i", NULL)
  }
  
  cd@x = merge(cd@x, x, all.x=TRUE, by="tlb")
  #setnames(xAll, "V1", "s_i")
  
  # if we can't estimate s_i robustly, assume it to be one
  set(cd@x, which(is.na(cd@x$s_i)), "s_i", 1)
  cd@x[, (normNcol):= pmax(1, round(get(Ncol)/s_i))]
        
  message("Post-processing...")
  setkey(cd@x, baitID, otherEndID) # to reorder
  
  othercols = names(cd@x)[!names(cd@x)%in%c("baitID", "otherEndID", "distbin", Ncol, normNcol)]
  
  setcolorder(cd@x, c("baitID", "otherEndID", "distbin", othercols, Ncol, normNcol))
  
  cd

}

##updates cd@params$distFunParams, can then be fed into .distFun()
estimateDistFun <- function (cd, method="cubic", plot=TRUE, outfile=NULL) {
  
  # Take the "refBinMean" column of the data x as f(d_b)
  # then interpolate & extrapolate to get f(d).
  # TODO output extra diagnostic information?
  
  if (!method == "cubic"){
    stop ("Unknown method.\n")
  }
  
  # Get f(d_b)
  setkey(cd@x, distbin)
  f.d <- unique(cd@x, by=key(cd@x))[is.na(refBinMean)==FALSE][, c("distbin", "refBinMean"), with=FALSE]

  setDF(f.d) # f.d is tiny, so no need to bother with it being a data.table
  f.d$midpoint <- seq(from=round(cd@settings$binsize/2), by=cd@settings$binsize, length.out=nrow(f.d))
  
  obs.min <- log(min(f.d$midpoint))
  obs.max <- log(max(f.d$midpoint))
  
  #   if(method == "lm") {
  #
  #     ##previously had arguments: n.obs.head=10, n.obs.tail=25
  #
  #     ##On log-scale, do a linear interpolation.
  #     ##Linear models applied to first (n.obs.head) observations, and to last (n.obs.tail) observations.
  #     
  #     ##Interpolation: Estimate f(d) (NB "rule" parameter = 1 forces NAs outside of range)
  #     log.f.obs <- approxfun(log(f.d$midpoint), log(f.d$refBinMean), rule=c(1,1))
  #     
  #     ##Extrapolation: Fit the "head" and "tail" of f using a linear model
  #     head.coef <- coefficients(lm(log(refBinMean)~log(midpoint), data = head(f.d, n.obs.head))) ##Fit for small d
  #     tail.coef <- coefficients(lm(log(refBinMean)~log(midpoint), data = tail(f.d, n.obs.tail))) ##Fit for large d
  #     
  #   }
  
  if(method == "cubic") {
    ##Spline - Cubic fit over observed interval, linear fit elsewhere, assume continuity of f(d) & f'(d).
    distFunParams <- list(method="cubic")
    
    ##cubic fit (quadratic not immensely different TBH)
    f.d.cubic <- lm(log(refBinMean) ~ log(midpoint) + I(log(midpoint)^2) + I(log(midpoint)^3), data = f.d)
    fit <- f.d.cubic$coefficients
    distFunParams[["cubicFit"]] <- fit
    
    ##Extrapolation: Fit the "head" and "tail" of f using continuity
    distFunParams[["obs.min"]] <- log(min(f.d$midpoint))
    distFunParams[["obs.max"]] <- log(max(f.d$midpoint))
    
    beta <- fit[2] + 2*fit[3]*c(obs.min, obs.max) + 3*fit[4]*(c(obs.min, obs.max)^2)
    alpha <- fit[1] + (fit[2] - beta)*c(obs.min, obs.max) + fit[3]*c(obs.min, obs.max)^2 + fit[4]*c(obs.min, obs.max)^3
    
    distFunParams[["head.coef"]] <- c(alpha[1], beta[1])
    distFunParams[["tail.coef"]] <- c(alpha[2], beta[2])
    
  }
  
  cd@params$distFunParams <- distFunParams
  
  if(plot)
  {
    if (!is.null(outfile)){ 
      pdf(outfile)
    }
    
    plotDistFun(cd)
#     my.log.d <- seq(from=obs.min, to=obs.max, length.out = 101)
#     my.d <- exp(my.log.d)
#     plot(my.log.d, log(.distFun(my.d, distFunParams)),
#          type="l",
#          main = "Distance function estimate",
#          xlab = "log distance",
#          ylab = "log f",
#          col = "Red")
     with(f.d, points(log(midpoint), log(refBinMean)))
     legend("topright", legend = c("Data", "Fit"), col = c("Black", "Red"), pch = c(1, NA), lty=c(0,1))
    if (!is.null(outfile)){
      dev.off()
    }
  }
  
  cd
}

.distFun <- function(d, distFunParams)
{
  ##d: distSign column of cd@x
  
  ##distFunParams: list of form
  ##list(method, ...)
  ##so list("cubic", head.coef, tail.coef, obs.min, obs.max, fit)
  ##in future, could support e.g. list("lm", ...)
  
  p <- distFunParams
  
  stopifnot(p[["method"]] == c("cubic")) ##compatibility with future "method"s
  
  ##Transfer parameters
  obs.max <- p[["obs.max"]]
  obs.min <- p[["obs.min"]]  
  head.coef <- p[["head.coef"]]
  tail.coef <- p[["tail.coef"]]
  fit <- p[["cubicFit"]]
  
  ##Put everything together to get the final function
  d <- log(d)
  
  out <- ifelse(d > obs.max,
                tail.coef[1] + d*tail.coef[2], ##Common case evaluated first (large d)
                ifelse(d < obs.min,
                       head.coef[1] + d*head.coef[2],
                       fit[1] + fit[2]*d + fit[3]*(d^2) + fit[4]*(d^3) ##2nd most common case evaluated second (small d)
                )
  )
  
  exp(out)
}

estimateBrownianComponent <- function(cd) {
  ##1) Reinstate zeros
  ##2) Add a "Bmean" column to x, giving expected Brownian component.
  ##3) Calculate dispersion by regressing against "Bmean", added to x as "dispersion" attribute
  ##subset: Since we don't need the entire data set, can just calculate based on a random subset of baits.
  ##!!NB!! Use set.seed to force subset analysis to be reproducible
  
  s = cd@settings
  adjBait2bait=s$adjBait2bait
  samples <- s$brownianNoise.samples
  subset <- s$brownianNoise.subset
  seed <- s$brownianNoise.seed
  maxLBrownEst = s$maxLBrownEst
  
  if (!is.na(seed)){
    set.seed(seed)
  }
  
  if(is.null(samples))
  {
    warning("brownianNoise.samples setting missing - cd was made using an old version of Chicago. Setting brownianNoise.samples to 5, which you can change using modifySettings()")
    cd@settings$brownianNoise.samples <- 5
  }
  
  siPresent <- "s_i" %in% colnames(cd@x)
  if(siPresent)
  {
    message("s_i factors found - estimating Brownian component...")
  } else {
    message("s_i factors NOT found - variance will increase, estimating Brownian component anyway...")
  }
  
  ##check if we are going to use the whole data set
  if(!is.na(subset))
  {
    if(!class(subset) %in% c("numeric","integer")) {stop("'subset' must be an integer.")}
    
    ##much faster than in situ data.frame calculation
    setkey(cd@x, baitID)
    if(nrow(cd@x[, .I[1], by=baitID]) <= subset){ 
      warning("subset > number of baits in data, so used the full dataset.\n")
      subset=NA
    }
  }
  ##if we're using the whole data set, there's no point repeatedly subsampling
  if(is.na(subset) & samples != 1){
    warning("We're using the whole data set to calculate dispersion. There's no reason to sample repeatedly in this case, so overriding brownianNoise.samples to 1.")
    samples <- 1
  }
  proxOE <- .readProxOEfile(s)
  cd@params$dispersion.samples <- replicate(samples, .estimateDispersion(cd, proxOE))
  message("Getting consensus dispersion estimate...")
  cd@params$dispersion <- mean(cd@params$dispersion.samples)
  
  cd@x <- .estimateBMean(cd@x, distFunParams=cd@params$distFunParams) ##NB: Different results from invocation of .estimateBMean() in .estimateDispersion.
  cd
}
estimateBrownianNoise <- estimateBrownianComponent

.estimateDispersion <- function(cd, proxOE)
{
  siPresent <- "s_i" %in% colnames(cd@x)

  s = cd@settings
  adjBait2bait=s$adjBait2bait
  samples <- s$brownianNoise.samples
  subset <- s$brownianNoise.subset
  maxLBrownEst = s$maxLBrownEst
  
  ##Pre-filtering: get subset of data, store as x
  ##---------------------------------------------
  
  if(!is.na(subset))
  {
    ##much faster than in situ data.frame calculation
    setkey(cd@x, baitID)
    if( nrow(cd@x[, .I[1], by=baitID])>subset){ 
      sel.sub <- sort(sample(unique(cd@x$baitID), subset))
      x <- cd@x[J(sel.sub)]
    }
    else{
      ##Use whole data set (a warning was triggered earlier)
      x <- cd@x
      subset=NA
    }
  }
  else{
    x <- cd@x
  }

  ##consider proximal region only...
  setkey(x, distSign)
  x = x[abs(distSign)<maxLBrownEst & is.na(distSign)==FALSE,] # will have NA for distal and trans interactions
  
  ##remove bait2bait...
  if (adjBait2bait){
    if (!"isBait2bait" %in% names(x)){
      x[, isBait2bait := FALSE]
      x[wb2b(otherEndID), isBait2bait:= TRUE] 
    }
    x = x[isBait2bait==FALSE]
  }
  
  ##1) Reinstate zeros:
  ##----------------
  ##1A) Choose some (uncensored) baits. Pick relevant proxOE rows. Note: censored fragments,
  ##   censored bait2bait pairs (etc...) already taken care of in pre-computation of ProxOE.  
  
  if(!is.na(subset)) {
    ## if we chose a subset of baits, restrict to that (none of these should be censored)
    sel.baits <- sel.sub
  } else {
    sel.baits <- unique(x$baitID)
  }
  
  setkey(proxOE, baitID)
  proxOE <- proxOE[J(sel.baits),]
  
  ##1B) Merge with our data, thus reinstating zero pairs.
  setkey(x, baitID, otherEndID)
  
  ##(make some lookup tables so we can get s_is, s_js later)
  # I don't think lookup tables are needed when we can subset and assign by reference
  # Probably more efficient to just use the original data table instead...
  
  # the lookup table can be generated on the fly x[, s_j[1], by=baitID]
  # but on the other hand, sjLookup is very small anyway
  sjLookup <- unique(x[,c("baitID","s_j"),with=FALSE])
  setkey(sjLookup, baitID)
  if(siPresent)
  {
    siLookup <- unique(x[,c("otherEndID","s_i"),with=FALSE])
    setkey(siLookup, otherEndID)
  }
  
  setkey(proxOE, baitID, otherEndID)
  
  if(siPresent)
  {
    x <- merge(x, proxOE, all.y=TRUE)[,c("baitID","otherEndID","s_i","s_j","N","distSign","dist"), with=FALSE]
  } else {
    x <- merge(x, proxOE, all.y=TRUE)[,c("baitID","s_j","N","distSign","dist"), with=FALSE]
  }
  ##Merging like this means that we are missing N, s_i, s_j information for most of the rows. So:
  ##1C) Repopulate table with information...
  
  # TODO: Recast following, avoid lookup tables, instead subset and assign by reference
  ## - 0s in Ncol
  x[is.na(N), N:=0]
  ## - s_js
  x[, s_j:=sjLookup[J(x$baitID)]$s_j ]
  ## - s_is (if present)
  if(siPresent)
  {
    x[, s_i:=siLookup[J(x$otherEndID)]$s_i ]
    if(any(is.na(x$s_i)))
    {
      ##If we don't have any information on a particular other end's s_i then...
      #       warning("Some other ends did not have s_i factors. Assuming s_i = 1 for these.")
      x[,s_i := ifelse(is.na(s_i), 1, s_i)]
    }
  }
  ## - distances
  ##Sanity check - the distances should agree (modulo rounding)
  if(any(removeNAs(abs(x$dist - abs(x$distSign))) > 1))
  {
    warning("estimateBrownianComponent: Distances in precomputed ProxOE file did not match distances supplied.")
  }
  
  x[, distSign:=dist]
  ##FIXME delete dist column
  
  ##2) Calculate Bmeans
  ##----------------
  x <- .estimateBMean(x, distFunParams=cd@params$distFunParams)
  
  ##3)Fit model
  ##---------
  message("Sampling the dispersion...")
  model <- glm.nb(formula= x$N ~ offset(log(x$Bmean)) + 0) 
  
  ##Construct Output
  ##----------------
  
  ##NB Parametrization: var = mu + (mu^2)/dispersion
  model$theta
}


.estimateBMean = function(x, distFunParams) {
  
  ##Adds a "Bmean" vector to a data.table x, giving expected value of Brownian component.
  ##NB updates by reference
  
  if("s_i" %in% colnames(x))
  {
    x[, Bmean:=s_j*s_i*.distFun(abs(distSign), distFunParams)]
  } else {
    warning("s_i factors NOT found in .estimateBMean - variance will increase, estimating means anyway...")
    x[, Bmean:=s_j*.distFun(abs(distSign), distFunParams)]
  }
  ##distcol == NA for trans-pairs. Thus set Bmean = 0.
  x[is.na(distSign), Bmean:=0]
  
  #out
  x
}

.addTLB = function(cd, adjBait2bait=TRUE){
  ##Assigns each fragment a "tlb" - a range containing the number of trans
  ##1) The bins are constructed based on reduced data - (outliers trimmed, )
  ##2) The bin endpoints are readjusted such that no fragments fall outside.
  ##3) These bins are then applied to the entire dataset.
  
  # cd is the current chicagoData object
  x = cd@x
  s = cd@settings
  
  filterTopPercent = s$tlb.filterTopPercent
  minProxOEPerBin = s$tlb.minProxOEPerBin
  minProxB2BPerBin = s$tlb.minProxB2BPerBin
  
  message("Preprocessing input...")
  
  # Checking whether in the input, we had distances at trans-interactions labeled as NA 
  # (as opposed to a dummy maximum distance)
  transNA = FALSE
  if(any(is.na(x$distSign))){
    transNA = TRUE
    transD = max(x[is.na(distSign)==FALSE]$distSign)+s$binsize
    x[is.na(distSign), distSign := transD]
  }
  else{
    transD = max(x$distSign)
    message("Warning: No NAs found in input. Assuming the max distance of ", transD, " is a dummy for trans-counts.")
  }
  
  message("Computing trans-counts...")
  
  if (adjBait2bait){
    if (!"isBait2bait" %in% names(x)){
      x[, isBait2bait := FALSE]
      x[wb2b(otherEndID), isBait2bait:= TRUE] 
    }
    setkey(x, otherEndID, distSign)
    # we are interested in the minimum distance for interactions that involve each other end
    # to then be able to check how many other ends we are pooling together for each tlb bin 
    # note that distSign[1] here means min(distSign) as it's keyed by this column
    # sum(distSign==transD) is a bit faster than length(.I[distSign==transD])
    
    # MEMORY-HUNGRY CODE
    # Watch this thread for possible solutions: http://stackoverflow.com/questions/29022185/how-to-make-this-r-data-table-code-more-memory-efficient?noredirect=1#comment46288765_29022185
    # Note that B2B interactions appear twice, once each way - this means that by=otherEndID is fine (no need for a by=baitID fudge).
    transLen = x[, list(sum(distSign==transD), isBait2bait[1], min(abs(distSign))), by=otherEndID]
    
    setnames(transLen, "V2", "isBait2bait")
    setnames(transLen, "V3", "distSign")
  }
  else{
    setkey(x, otherEndID, distSign)
    transLen = x[, list(length(.I[distSign==transD]), min(abs(distSign))), by=otherEndID]    
    setnames(transLen, "V2", "distSign")
  }
  
  setnames(transLen, "V1", "transLength")
  
  transLen0 = transLen
  transLen = transLen[transLength<=quantile(transLen$transLength,1-filterTopPercent/100)] 
  filteredLen = nrow(transLen0[!otherEndID %in% transLen$otherEndID])
  message("Filtering out ", filteredLen, " other ends with top ", filterTopPercent, "% number of trans-interactions")
  
  # first use cut2 to compute bin boundaries on the proximal range based on the desired minProxOEPerBin
  # (note cut2 doesn't guarantee that all bins will contain this min number of observations)
  # then use cut to split the whole dataset based on these bins
  # note we'll need the full dataset (and not only proximal interactions) assigned to tlb bins
  # when estimating technical noise
  # If adjBait2bait == TRUE, do this separately for bait2bait and non-bait2bait other ends.  
  
  if (adjBait2bait){
    transLen0 = transLen
    transLen = transLen0[isBait2bait==FALSE]
    transLenB2B = transLen0[isBait2bait==TRUE]
  }
  
  message("Binning...")
  
  cuts = cut2(transLen[abs(distSign)<= s$maxLBrownEst]$transLength, 
              m=minProxOEPerBin, onlycuts=TRUE)
  # for really depleted data sets, cuts is a single number, usually 0.
  if(length(cuts) == 1)
  {
    tlbClasses <- factor(rep(cuts, nrow(transLen)))
  } else {
    # If some other ends that do not feature in any proximal interactions have transLen's outside of the range
    # determined based on the proximal interactions, just move the boundaries of the first or last tlb bin accordingly...
    if (min(cuts)>min(transLen$transLength)){
      cuts[1] = min(transLen$transLength)
    }
    if (max(cuts)<max(transLen$transLength)){
      cuts[length(cuts)] = max(transLen$transLength)
    }
    tlbClasses <- cut(transLen$transLength, breaks=cuts, include.lowest=TRUE)
  }
  set(transLen, NULL, "tlb" , tlbClasses)
  
  if (adjBait2bait){
    
    cutsB2B = cut2(transLenB2B[abs(distSign)<= s$maxLBrownEst]$transLength, 
                   m=minProxB2BPerBin, onlycuts=TRUE)
    # for really depleted data sets, cutsB2B is a single number, usually 0.
    if(length(cutsB2B) == 1)
    {
      tlbClassesB2B <- factor(rep(cutsB2B, nrow(transLenB2B)))
    } else {
      # If some other ends that do not feature in any proximal interactions have transLen's outside of the range
      # determined based on the proximal interactions, just move the boundaries of the first or last tlb bin accordingly...
      if (min(cutsB2B)>min(transLenB2B$transLength)){
        cutsB2B[1] = min(transLenB2B$transLength)
      }
      if (max(cutsB2B)<max(transLenB2B$transLength)){
        cutsB2B[length(cutsB2B)] = max(transLenB2B$transLength)
      }
      tlbClassesB2B <- cut(transLenB2B$transLength, breaks=cutsB2B, include.lowest=TRUE)
    }
    set(transLenB2B, NULL, "tlb", tlbClassesB2B)
    levels(transLenB2B$tlb) = paste0(levels(transLenB2B$tlb), "B2B")
    
    transLen = rbind(transLen, transLenB2B)
  }      
  
  set(transLen, NULL, "transLength", NULL)
  set(transLen, NULL, "isBait2bait", NULL)
  set(transLen, NULL, "distSign", NULL)
  
  setkey(x, otherEndID)
  setkey(transLen, otherEndID)
  
  ##discard TLB if already present
  if("tlb" %in% colnames(x))
  {
    set(x, NULL, "tlb", NULL)
  }
  
  x = x[transLen] # note that if mode="even_filtered", we're not just merging, but also trimming x, 
  # removing the interactions with too "sticky" other ends and those mapping to very sparse bins
  
  if(transNA){
    x[distSign==max(x$distSign), distSign := NA]
  }
  
  x
}

estimateTechnicalNoise = function(cd, plot=TRUE, outfile=NULL){ 

# Estimate technical noise based on mean counts per bin, with bins defined based on trans-counts for baits _and_ other ends 
# Note we need raw read counts for this, as normalisation is done wrt Brownian component. 
  
# NB: filterTopPercent, minProxOEPerBin, minProxB2BPerBin
# are input parameters for .addTLB (the function for binning other ends) that is only called  
# if tlb's aren't already present in the input - and they will be present if other end normalisation
# has been applied 

  message("Estimating technical noise based on trans-counts...")
  .checkForIncompatibilities(cd)
  
  
  ##test to see if tblb, Tmean columns exist, warn & delete if so
  replacedCols <- c("tblb", "Tmean")
  sel <- replacedCols %in% colnames(cd@x)
  if(any(sel))
  {
    warning("Columns will be overwritten: ", paste(replacedCols[sel], collapse=", ")) 
    set(cd@x, j=replacedCols[sel], value=NULL) ##delete these columns
  }
  
  
  minBaitsPerBin = cd@settings$techNoise.minBaitsPerBin
  adjBait2bait = cd@settings$adjBait2bait
  
  ##TODO: considering turning trans counts into "Inf"

  if (!"tlb" %in% names(cd@x)){
    message("Binning other ends based on trans-counts...")
    cd@x = .addTLB(cd)
  }
  
  setkey(cd@x, baitID)
  
  message("Binning baits based on observed trans-counts...")
  
  transBaitLen = cd@x[, sum(is.na(distSign)), by=baitID] ##Number of trans counts per bait
  setnames(transBaitLen, "V1", "transBaitLen")
  
  transBaitLen$tblb = cut2(transBaitLen$transBaitLen, m=minBaitsPerBin, levels.mean=FALSE)
  
  setkey(transBaitLen, baitID)
  setkey(cd@x, baitID)
  cd@x = cd@x[transBaitLen]
  
  message("Defining interaction pools and gathering the observed numbers of trans-counts per pool...")
  
  # Getting the observed numbers of trans-counts 
  setkeyv(cd@x, c("tlb", "tblb"))
  Ntrans = cd@x[, { res=table(N[is.na(distSign)]) 
                    list(as.numeric(names(res)), as.numeric(res)) 
                  }, by=c("tlb", "tblb")]
  setnames(Ntrans, "V1", "N")
  setnames(Ntrans, "V2", "nobs")
  
  message("Computing the total number of possible interactions per pool...")
  message("Preparing the data...", appendLF = FALSE)
  
  # Now adding the zeros based on how many trans-interactions are possible in each (tlb, tblb) bin
  baitmap = .readBaitmap(cd@settings)
  rmap = .readRmap(cd@settings)
    
  setnames(rmap, "V1", "chr")
  setnames(baitmap, "V1", "chr")
  setnames(rmap, "V4", "otherEndID")
  setnames(baitmap, "V4", "baitID")
  
  message(".", appendLF=FALSE)
  
  baitmap = baitmap[,c("baitID", "chr"), with=FALSE]
  rmap = rmap[,c("otherEndID", "chr"), with=FALSE]
  
  setkey(baitmap, baitID)
  setkey(rmap, otherEndID)

  message(".", appendLF=FALSE)
  
  setkey(cd@x, baitID)
  cd@x = baitmap[cd@x]
  setnames(cd@x, "chr", "baitChr")
  setkey(cd@x, otherEndID)
  cd@x = rmap[cd@x]
  setnames(cd@x, "chr", "otherEndChr")
  setkey(cd@x, tlb, tblb)
  
  message("\nProcessing fragment pools", appendLF=FALSE)

  res = cd@x[, {
    message(".", appendLF=FALSE)
        
    baits = unique(baitID)
    oes = unique(otherEndID)
    
    bChr = baitChr[!duplicated(baitID)]
    oeChr = otherEndChr[!duplicated(otherEndID)]

    # computing the total number of pairs;
    # for bait2bait interactions, it's possible that some oe's will also be among the baits, in which case each such interaction between pairs of such baits will be counted twice, so take care of this 
    numPairs = sum(sapply(unique(bChr), function(this)length(bChr[bChr==this])*length(oeChr[oeChr!=this])))-length(baits[baits%in%oes]) 
    
    nTrans = sum(N[is.na(distSign)])
    Tmean = nTrans/numPairs
    
    list(numPairs=numPairs, nTrans=nTrans, Tmean=Tmean)
    
  }  , by=c("tlb", "tblb")]
  
  if(plot){ 
    message("\nPlotting...")
    if(!is.null(outfile)){ pdf(outfile)}
    par(mfrow=c(2,1))
    boxplot(Tmean~tblb, as.data.frame(res), main="Technical noise estimates per bait pool")
    boxplot(Tmean~tlb, as.data.frame(res), main="Technical noise estimates per other end pool")
    par(mfrow=c(1,1))
    if(!is.null(outfile)){dev.off()}
  }
  
  message("Post-processing the results...")
  
  setkeyv(cd@x, c("tlb", "tblb"))
  setkeyv(res, c("tlb", "tblb"))
  cd@x = res[cd@x]

  set(cd@x, NULL, "transBaitLen", NULL)
  set(cd@x, NULL, "numPairs", NULL)
  set(cd@x, NULL, "nTrans", NULL)
  set(cd@x, NULL, "baitChr", NULL)
  set(cd@x, NULL, "otherEndChr", NULL)
  
  setkey(cd@x, baitID, otherEndID) # re-sort this way
  othernames = names(cd@x)[!names(cd@x) %in% c("baitID", "otherEndID", "tlb", "tblb", "Tmean")]
  setcolorder(cd@x, c("baitID", "otherEndID", othernames, "tlb", "tblb", "Tmean")) ##tlb, tblb are the classes of the other ends, based on trans counts
  cd  
}

getPvals <- function(cd){
  ## - Calls p-values
  
  # No need for this anymore
  alpha = cd@params$dispersion
  x = cd@x
  if(is.null(alpha)) {stop("getPvals: 'dispersion' parameter of x not found.")}
  
  message("Calculating p-values...") 
  
  ##p-values:
  ##(gives P(X > x-1) = P(X >= x))
  ##Note that the cases Bmean = 0 and Bmean > 0 are considered separately.
  x[, log.p:=NA_real_]
  x[Bmean < .Machine$double.eps, log.p:=ppois(N - 1L, lambda=Tmean, lower.tail=FALSE, log.p=TRUE)]
  x[Bmean >= .Machine$double.eps, log.p:=pdelap(N - 1L, alpha, beta=Bmean/alpha, lambda=Tmean, lower.tail=FALSE, log.p=TRUE)]
  
  # Large N approximation
  # ---------------------
  
  ##In rare cases where pdelap returns Infs, estimate the p-value magnitude
  ##using an NB approximation, through method of moments argument
  ##NaNs occur when pdelap() thinks the p-value is negative (since can have 1 - 1 != 0),
  ##thus these are also approximated.
  sel <- which(is.infinite(x$log.p) | is.nan(x$log.p))
  if(length(sel) > 0)
  {
    message("Approximating ", length(sel), " very small p-values.")
    
    gamma <- x[sel,alpha*(1+Tmean/Bmean)^2] ##gamma is the "effective" dispersion
    
    ##in the case where Bmean << Tmean, gamma becomes Inf, so cap gamma above
    ##(should make very little difference since we are basically Poisson in this case)
    ##Case where Bmean >> Tmean causes no problem.
    gamma <- pmin(gamma, 1e10)
    
    x[sel,log.p := pnbinom(N - 1L, size=gamma, mu=Bmean+Tmean, lower.tail=FALSE, log.p=TRUE)]
    
    if(any(is.infinite(x[sel,log.p]))) {warning("Some log-p-values were infinite.")}
    if(any(is.nan(x[sel,log.p]))) {warning("Some log-p-values were NaNs.")}
  }
  if(any(is.na(x$log.p))) {warning("Some log-p-values were NA.")}
  cd
}

getScores <- function(cd, method="weightedRelative", includeTrans=TRUE, plot=TRUE, outfile=NULL)
{
  ## - If method="weightedRelative", we divide by weights (Genovese et al 2006)
  ##Note to self: Algebra is on P96 of my lab notebook.
  
  x <- cd@x
  set <- cd@settings
  avgFragLen <- .getAvgFragLength(cd) ##average fragment length
  
  if(!method %in% c("unweighted","weightedRelative")) {stop("method=",method," not recognized.")}
  
  if(!includeTrans)
  {
    x <- x[is.na(distSign)] ##Cannot delete row by reference yet?
  }
  
  if(method == "weightedRelative")
  {
    eta.bar <- .getEtaBar(cd)
    
    ##Get weights, weight p-values
    message("Calculating p-value weights...")
    x[, log.w:= .getWeights(abs(x$distSign), cd, eta.bar=eta.bar, includeTrans=includeTrans)]
    x[, log.q:= log.p - log.w] ##weighted p-val
    message("Calculating scores...")
    
    ##get score (more interpretable than log.q)
    minval <- .getWeights(0, cd, eta.bar=eta.bar, includeTrans=includeTrans) ##FIXME could be optimized a *lot*.
    x[,score := pmax(- minval - log.q, 0)]
    
  } else {
    x[,score := - log.p]
  }
  cd
}

.getAvgFragLength <- function(cd, rmapfile=NULL, excludeMT=TRUE)
{
  ##Normally, takes rmapfile from cd object.
  ##However, if rmapfile is specified, this overrides cd.
  
  if(!is.null(rmapfile))
  {
    rmap <- .readRmap(list(rmapfile=rmapfile))
  } else {
    rmap = .readRmap(cd@settings)
  }
  setnames(rmap, "V1", "chr")
  setnames(rmap, "V3", "end")
  if(excludeMT) {  
    chrMax <- rmap[chr != "MT",max(end),by="chr"] ##length of each chr
  }else{
    chrMax <- rmap[,max(end),by="chr"] ##length of each chr
  }
  sum(as.numeric(chrMax$V1))/nrow(rmap)
}

.getNoOfHypotheses <- function(cd, includeTrans=TRUE, excludeMT=TRUE)
{
  set <- cd@settings
  
  ##How many hypotheses are we testing? (algebra on p246 of JMC's lab notebook)
  rmap = .readRmap(set)
  setnames(rmap, "V1", "chr")
  setnames(rmap, "V3", "end")
  chrMax <- rmap[,max(end),by="chr"] ##length of each chr
  
  baitmap = .readBaitmap(set)
  nBaits <- table(baitmap$V1) ##number of baits on each chr
  
  avgFragLen <- .getAvgFragLength(cd) ##average fragment length
  
  chr <- chrMax$chr
  if(excludeMT && any(chr == "MT")) chr <- chr[chr != "MT"] ##remove mitochondria
  
  ##count # of hypotheses
  if(includeTrans)
  {
    Nhyp <- sum(nBaits)*(2*nrow(rmap) - sum(nBaits) - 1)/2L ##number of hypotheses being tested
  } else {
    temp <- chrMax[chrMax$chr != "MT",]
    temp$nBaits <- nBaits[as.character(temp$chr)]
    temp$nFrag <- as.integer(temp$V1/avgFragLen)
    temp$Nhyp <- with(temp, nBaits*(2*nFrag - nBaits - 1))
    Nhyp <- sum(temp$Nhyp) ##see p257 of JMC lab book
  }
  Nhyp
}

.getEtaBar <- function(cd, includeTrans=TRUE)
{
  set <- cd@settings
  
  ##1. Collect parameters
  alpha = set$weightAlpha
  beta = set$weightBeta
  gamma = set$weightGamma
  delta = set$weightDelta
  
  ##2. Get genomic/fragment map information
  rmap = .readRmap(set)
  setnames(rmap, "V1", "chr")
  setnames(rmap, "V3", "end")
  chrMax <- rmap[,max(end),by="chr"] ##length of each chr
  
  baitmap = .readBaitmap(set)
  nBaits <- table(baitmap$V1) ##number of baits on each chr
  
  chr <- as.character(names(nBaits))
  if(any(chr %in% c("MT", "chrMT")))
  {
    chr <- chr[chr != "MT"] ##no mitochondria
  }
  
  avgFragLen <- .getAvgFragLength(cd)
  
  ##count # of hypotheses
  Nhyp <- .getNoOfHypotheses(cd, includeTrans)
  
  ##3. Calculate eta.bar
  ##Loop, summing contributions of eta
  
  eta.sigma <- 0 
  for(c in chr)
  {
    ##length of chromosome
    d.c <- as.numeric(chrMax$V1[chrMax$chr == c])
    ##no of baits on chromosome
    n.c <- nBaits[c]
    
    for(i in 1:n.c) ##TODO nested for loop can be replaced with lapply & function
    {
      d = d.c*i/n.c
      d.near = min(d, d.c-d) ##dist to nearest chromosome end
      d.other <- seq(from=avgFragLen, to=max(avgFragLen,d.near), by = avgFragLen) ##locations of fragments
      d.other2 <- seq(from=d.near, to=d.c-d.near, by = avgFragLen)
      eta.sigma <- eta.sigma + 2*sum(expit(alpha + beta*log(d.other))) + sum(expit(alpha + beta*log(d.other2)))
      #eta[i] <- 2*sum(expit(alpha + beta*log(d.other))) + sum(expit(alpha + beta*log(d.other2)))
    }
  }
  
  eta.bar <- eta.sigma/Nhyp
  eta.bar
}

.getWeights <- function(dist, cd, eta.bar, includeTrans=TRUE)
{
  set <- cd@settings
  alpha = set$weightAlpha
  beta = set$weightBeta
  gamma = set$weightGamma
  delta = set$weightDelta
  
  ##4. Calculate weights
  eta <- expit(alpha + beta*log(naToInf(dist)))
  log.w <- log((expit(delta) - expit(gamma))*eta + expit(gamma)) -
    log((expit(delta) - expit(gamma))*eta.bar + expit(gamma))
  
  log.w
}

.normaliseFragmentSets = function(x, s, npb, viewpoint, idcol, Ncol, adjBait2bait=TRUE, shrink=TRUE, 
                          refExcludeSuffix=NULL, plot=TRUE, outfile=NULL, debug=FALSE){   #minPosBins = 5, 
  
  # The normalisation engine used for normaliseBaits and normaliseOtherEnds
  # "Viewpoint" will be used in the comments for either baits or sets of other ends,
  # depending on the direction of normalisation.
  
  # npb is a data table containing the number of reads per viewpoint per bin,
  # for viewpoint=="bait" it's just the table from the nperbinfile (read via .readNPBfile),
  # where npb's idcol should match x's idcol.
  # for viewpoit=="otherEnd" the table from the nbaitsperbin file needs preprocessing   
  # to sum over pools of other ends, and this is done in normaliseOtherEnds,
  # with the nbp column added to x before submitting to .normaliseFragmentSets.

  # s is the current chicagoData object's settings list
  
  bin=s$binsize
  
  if (!viewpoint %in% c("bait", "otherEnd")){
    stop("viewpoint must be either \"bait\" or \"otherEnd\"")
  }
  
  if(viewpoint=="bait"){
    scol = "s_j"
  
  	set(x, NULL, "distbin", cut(abs(x$distSign), seq(0, s$maxLBrownEst, bin))) # will have NA for distal and trans interactions

  	xAll = x # full data table with bait2bait and distal interactions
  
  	if (adjBait2bait){
      if (!"isBait2bait" %in% names(x)){
        x[, isBait2bait := FALSE]
        x[wb2b(otherEndID), isBait2bait:= TRUE] 
      }
    	x = x[isBait2bait==FALSE]
  	}
  
  	# x is the data table used to compute the scaling factors
	  x = x[is.na(distbin)==FALSE]    
  
  	setkeyv(x, idcol)
  	setkeyv(npb, idcol)
  	# compute the total number of valid other ends for each viewpoint and each bin
  	x = npb[x]
  
    # MEMORY-HUNGRY CODE
    # watch this thread for possible better solutions: http://stackoverflow.com/questions/29022185/how-to-make-this-r-data-table-code-more-memory-efficient?noredirect=1
    setkeyv(x, c(idcol, "distbin"))
    x[, binCol:=do.call(paste0, list("bin", as.integer(distbin))), by=distbin][, ntot:=get(binCol)[1],by=c(idcol, "distbin")]
    set(x, NULL, "binCol", NULL)
    
  }
  else{
	# for other ends, bait2bait adjustment and adding ntot is performed before calling this function 
  
    scol = "s_i"
  }

  message("Computing binwise means...")
  
  rmCol = names(x)[!names(x) %in% c(idcol, Ncol, "distbin", "ntot")]
  for (col in rmCol){
    set(x, NULL, col, NULL)    
  }
    
  # sbbm is the input matrix for normalisation that contains the mean number of reads per bin for each bait
  setkeyv(x, c(idcol, "distbin"))
  sbbm = x[, sum(get(Ncol))/ntot[1] , by=c(idcol, "distbin")]
  setkeyv(sbbm, "distbin")  
  setnames(sbbm, "V1", "bbm")
  setkey(sbbm, "distbin")
  
  # Compute the "virtual reference" profile to normalise against 
  if(viewpoint=="bait" | is.null(refExcludeSuffix)){ 
    sbbm[, geomean:=geo_mean(bbm), by="distbin"]
  } else{
    # for other end normalisation, not including pools containing b2b interactions into the virtual reference
    # computation, as they are a minority and we expect them to have much higher geo_means
    sbbm[, geomean:=geo_mean(bbm[-grep(paste0(".",refExcludeSuffix,"$"), get(idcol))]), by="distbin"]    
  }
  
  # DEseq-style normalisation
  if (!shrink | viewpoint=="otherEnd"){
    sbbm [, s_iv:=bbm/geomean]
    setkeyv(sbbm, idcol)    
    s_v = sbbm[, median(s_iv[!is.na(s_iv)]), by=idcol]
  }
  else{   # Same but with a Gamma shrinkage of binwise means - currently deprecated
    message("Computing shrunken means...")
        
    setkeyv(sbbm, "distbin")
    sbbm[, c("shape", "rate"):={r=fitdistr(bbm, "gamma"); 
                                list(r[[1]]["shape"], r[[1]]["rate"])}, by="distbin"]
    
    shr = sbbm[, list(shape[1], rate[1]), by="distbin"]
    setnames(shr, "V1", "shape")
    setnames(shr, "V2", "rate")
    setkeyv(shr, "distbin")
    shr[, shapefit:=predict(loess(shape~as.integer(distbin)))]
    shr[, ratefit:=predict(loess(rate~as.integer(distbin)))]

    if (plot){
      if (!is.null(outfile)){ pdf(outfile) }
      par(mfrow=c(1,2))
      plot(as.integer(shr$distbin), shr$shape)
      lines(as.integer(shr$distbin), shr$shapefit)
      plot(as.integer(shr$distbin), shr$rate)
      lines(as.integer(shr$distbin), shr$ratefit)
      if (!is.null(outfile)){ dev.off() }
    }
    
    setkeyv(x, "distbin")
    x = shr[, list(shapefit[1], ratefit[1]),by="distbin"][x]
    setnames(x, "V1", "shapefit")
    setnames(x, "V2", "ratefit")
    
    sbbm2 = x[, (sum(get(Ncol))+shapefit[1])/(ntot[1]+ratefit[1]), by=c(idcol, "distbin")]
    setkeyv(sbbm2, "distbin")  
    setnames(sbbm2, "V1", "bbm")
    setkey(sbbm2, "distbin")
    
    sbbm2[, geomean:=geo_mean(bbm), by="distbin"]
    sbbm2 [, s_iv:=bbm/geomean]
    
    # Be even more stringent: exclude bins where shrunken binmeans  
    # fall to a 5% percentile of the Gamma distribution  
    # This helps "drug" the bottom 5% up without affecting much else
    
    setkeyv(sbbm2, "distbin")
    sbbm2 = sbbm2[shr]
    sbbm2$p = pgamma(sbbm2$bbm, shape=sbbm2$shapefit, rate=sbbm2$ratefit)
    sbbm2$s_ivfilt = sbbm2$s_iv
    sbbm2$s_ivfilt[sbbm2$p<0.05] = NA
    s_v = sbbm2[, median(s_ivfilt[!is.na(s_ivfilt)]), by=idcol]
    sbbm = sbbm2
  }

  if(debug) {return(sbbm)}
  
  setnames(s_v, "V1", scol)  
  
  if(any(is.na(s_v[[scol]]))){
    message("The following viewpoints couldn't be robustly normalised (too sparse?) and will be removed:")
    print(s_v[is.na(get(scol))==TRUE][[idcol]])
    s_v = s_v[is.na(get(scol))==FALSE]
  }
    
  if(viewpoint=="otherEnd"){
  	xAll = x
  }
    
  # Assign s_v to each observation in viewpoint, but bait/OE refmeans only to the proximal ones 
  setkeyv(s_v, idcol)
  setkeyv(xAll, idcol) 
  xAll = s_v[xAll]# here we'll remove baits that couldn't be robustly normalised
  
  # Now, for viewpoint=="bait" bind geomeans to those bins for which they were computed
  # in getInteractionScores, they will be interpolated for exact distances, 
  # and extrapolated for theremaining more distal interactions  
  
  if (viewpoint=="bait"){
    gm = sbbm[, geomean[1], by="distbin"]
    setnames(gm, "V1", "refBinMean")
    setkeyv(gm, "distbin")
    setkeyv(xAll, "distbin")
    xAll = merge(xAll, gm, all.x=TRUE) 
  }
  
  xAll
  
}  

## Hidden "read" functions ----------------------

.readRmap = function(s){
  fread(s$rmapfile, colClasses = list(character=1))
}

.readBaitmap = function(s){
  fread(s$baitmapfile, colClasses = list(character=1))
}

.readNPBfile = function(s){
  
  # Reads a pre-made text file containing the numbers of fragments per bait per distance bin 
  # within the interval maxl, given binsize.
  # The file can be generated by countNperBin.py and its first line should start with # and 
  # contain the parameter settings used. In addition to maxl and binsize, it defines the 
  # filtering parameters for restriction fragment minsize, maxsize as well as a boolean
  # variable removeb2b specifying whether bait2bait interactions should also not be counted.
  # (In fact, they probably should never be).
  
  # s is the current chicagoData object's settings list
  
  message("Reading NPerBin file...")
  header = readLines(s$nperbinfile, n=1)
  params = sapply(sapply(strsplit(header, "\t")[[1]],function(x)strsplit(x,"=")[[1]]), function(x)x[2])
  params = params[2:length(params)]
  names(params) = gsub("(\\S+)=.+", "\\1", names(params))
  minsize = as.numeric(params[["minFragLen"]])
  if (minsize != s$minFragLen){
    stop("The minFragLen in the NPerBin file header is not equal to minFragLen defined in experiment settings. Amend either setting (and if needed, generate a new NPerBin file) before running the analysis\n")
  }
  maxsize = as.numeric(params[["maxFragLen"]])
  if (maxsize != s$maxFragLen){
    stop("The maxFragLen in the NPerBin file header is not equal to maxFragLen defined in experiment settings. Amend either setting (and if needed, generate a new NPerBin file) before running the analysis\n")
  }
  maxl = as.numeric(params[["maxLBrownEst"]])
  if (maxl != s$maxLBrownEst){
    stop("The maxLBrownEst in the NPerBin file header is not equal to maxLBrownEst defined in experiment settings. Amend either setting (and if needed, generate a new NPerBin file) before running the analysis\n")
  }
  binsz = as.numeric(params[["binsize"]]) 
  if (binsz != s$binsize){
    stop("The binsize in the NPerBin file header is not equal to binsize defined in experiment settings. Amend either setting (and if needed, generate a new NPerBin file) before running the analysis\n")
  }  
  if (params[["removeb2b"]]!="True"){
    stop("The NPerBin file must be generated with removeb2b==True. Please generate a new file.\n")
  }
  if ( (params[["removeAdjacent"]]=="True" & !s$removeAdjacent) | (params[["removeAdjacent"]]!="True" & s$removeAdjacent)  ){
    stop("The removeAdjacent parameter settings used for generating NPerBin file (according to its header) and defined in experiment settings do not match. Amend either setting (and if needed, generate a new NPerBin file) before running the analysis\n")
  }  
  if(basename(params[["rmapfile"]]) != basename(s$rmapfile)){
    stop("The .rmap file used for generating the NPerBin file (according to the NPerBin header) and the one defined in experiment settings do not match. Amend either setting (and if needed, generate a new NPerBin file) before running the analysis\n")
  }
  
## Not checking this for now as we have a mixup of _baits and _baits_ID files used at different times...
#   if(basename(params[["baitmapfile"]]) != basename(baitmapfile)){
#     stop("Bait files used for generating the NfragPerBin file and defined here do not match. 
#          Amend either setting before running the analysis\n")
#   }

  npb = fread(s$nperbinfile, skip=1L)
  setnames(npb, names(npb)[1], "baitID")
  for(i in 2:ncol(npb)){
    setnames(npb, names(npb)[i], paste0("bin", i-1))    
  }
  npb
}

.readNbaitsPBfile = function(s){
  
  # Reads a pre-made text file containing the numbers of baits per other end per distance bin 
  # within the interval maxl, given binsize.
  # The file can be generated by countNBaitsPerBin.py and its first line should start with # and 
  # contain the parameter settings used. 
  
  # s is the current chicagoData object's settings list
  
  message("Reading NBaitsPerBin file...")
  header = readLines(s$nbaitsperbinfile, n=1)
  params = sapply(sapply(strsplit(header, "\t")[[1]],function(x)strsplit(x,"=")[[1]]), function(x)x[2])
  params = params[2:length(params)]
  names(params) = gsub("(\\S+)=.+", "\\1", names(params))

  maxl = as.numeric(params[["maxLBrownEst"]])
  if (maxl != s$maxLBrownEst){
    stop("The maxLBrownEst in the NBaitsPerBin file header is not equal to maxLBrownEst defined in experiment settings. Amend either setting (and if needed, generate a new NBaitsPerBin file) before running the analysis\n")
  }

  # Currently binsize is called bin, but should correct this
  binsz = as.numeric(params[["binsize"]]) 
  if (binsz != s$binsize){
    stop("The binsize in the NBaitsPerBin file header is not equal to binsize defined in experiment settings. Amend either setting (and if needed, generate a new NBaitsPerBin file) before running the analysis\n")
  }  
  
  # Currently not in the file
#   if (params[["removeb2b"]]!="True"){
#     stop("The NfragPerBin file must be generated with removeb2b==True\n")
#   }
#   if ( (params[["removeAdjacent"]]=="True" & !removeAdjacent) | (params[["removeAdjacent"]]!="True" & removeAdjacent)  ){
#     stop("The removeAdjacent parameter settings used for generating NfragPerBin file and defined here do not match. 
#          Amend either setting before running the analysis\n")
#   }  
  
  if(basename(params[["rmapfile"]]) != basename(s$rmapfile)){
    stop("The .rmap file used for generating the NBaitsPerBin file (according to the NBaitsPerBin header) and the one defined in experiment settings do not match. Amend either setting (and if needed, generate a new NBaitsPerBin file) before running the analysis\n")
  }
  
  ## Not checking this for now as we have a mixup of _baits and _baits_ID files used at different times...
  #   if(basename(params[["baitmapfile"]]) != basename(baitmapfile)){
  #     stop("Bait files used for generating the NfragPerBin file and defined here do not match. 
  #          Amend either setting before running the analysis\n")
  #   }
  
  nbpb = fread(s$nbaitsperbinfile, skip=1L)
  setnames(nbpb, names(nbpb)[1], "otherEndID")
  for(i in 2:ncol(nbpb)){
    setnames(nbpb, names(nbpb)[i], paste0("bin", i-1))    
  }
  nbpb
}


.readProxOEfile <- function(s){
  
  # Reads a pre-computed text file that denotes which other ends are in the proximal 
  # range relative to each bait, and gives that distance. 
  # Note that fragments that are too small/too large have already been removed.
  # s is the current chicagoData object's settings list
  
  message("Reading ProxOE file...")
  header = readLines(s$proxOEfile, n=1)
  params = sapply(sapply(strsplit(header, "\t")[[1]],function(x)strsplit(x,"=")[[1]]), function(x)x[2])
  params = params[2:length(params)]
  names(params) = gsub("(\\S+)=.+", "\\1", names(params))
  minsize = as.numeric(params[["minFragLen"]])
  if (minsize != s$minFragLen){
    stop("The minFragLen specified in the ProxOE file header is not equal to minFragLen defined in experiment settings. Amend either parameter setting (and if needed, generate a new ProxOE file) before running the analysis\n")
  }
  maxsize = as.numeric(params[["maxFragLen"]])
  if (maxsize != s$maxFragLen){
    stop("The maxFragLen specified in the ProxOE file header is not equal to maxFragLen defined in experiment settings. Amend either parameter setting (and if needed, generate a new ProxOE file) before running the analysis\n")
  }
  maxl = as.numeric(params[["maxLBrownEst"]])
  if (maxl != s$maxLBrownEst){
    stop("The maxLBrownEst specified in the ProxOE file header is not equal to maxLBrownEst defined in experiment settings. Amend either parameter setting (and if needed, generate a new ProxOE file) before running the analysis\n")
  }
  binsz = as.numeric(params[["binsize"]]) 
  if (binsz != s$binsize){
    stop("The binsize specified in the ProxOE file header is not equal to binsize defined in experiment settigs. Amend either parameter setting (and if needed, generate a new ProxOE file) before running the analysis\n")
  }  
  if (params[["removeb2b"]]!="True"){
    stop("The ProxOE file must be generated with removeb2b==True. Please generate a new file.\n")
  }
  if ( (params[["removeAdjacent"]]=="True" & !s$removeAdjacent) | (params[["removeAdjacent"]]!="True" & s$removeAdjacent)  ){
    stop("The removeAdjacent parameter settings used for generating ProxOE file (according to its header) and defined in experiment settings do not match. Amend either setting (and if needed, generate a new ProxOE file) before running the analysis\n")
  }  
  if(basename(params[["rmapfile"]]) != basename(s$rmapfile)){
    stop("The .rmap files used for generating the ProxOE file (according to the ProxOE header) and the one defined in experiment settings do not match. Amend either setting (and if needed, generate a new ProxOE file) before running the analysis\n")
  }
  ## Not checking this for now as we have a mixup of _baits and _baits_ID files used at different times...
  #   if(basename(params[["baitmapfile"]]) != basename(baitmapfile)){
  #     stop("Bait files used for generating the ProxOE file and defined here do not match. 
  #          Amend either setting before running the analysis\n")
  #   }
  proxOE = fread(s$proxOEfile, skip=1L)
  setnames(proxOE, 1:3, c("baitID", "otherEndID", "dist"))
  proxOE
  }

## Export/plotting functions ----------------------

plotBaits=function(cd, pcol="score", Ncol="N", n=16, baits=NULL, plotBaitNames=TRUE, plotBprof=FALSE,plevel1 = 5, plevel2 = 3, outfile=NULL, removeBait2bait=TRUE, width=20, height=20, maxD=1e6, bgCol="black", lev2Col="blue", lev1Col="red", bgPch=1, lev1Pch=20, lev2Pch=20, ...)
{
  if(plotBaitNames){
    baitmap = .readBaitmap(cd@settings)
  }
  if (is.null(baits)){
    baits = sample(unique(cd@x$baitID),n)
  }
  else{
    n = length(baits)
  }
 
  if(plotBprof){
    disp = cd@params$dispersion 
  }
 
  if (!is.null(outfile)){ 
    pdf(outfile, width=width, height=height)
  }
  
  ncols = ceiling(n/4)
  nrows = ceiling(n/ncols)
  par(mfrow=c(nrows, ncols))
  
  setkey(cd@x, baitID)

  myArgs = list(...)
  xlimInArgs = ("xlim" %in% names(myArgs))

  minD = NULL
  if(xlimInArgs){
        minD = myArgs$xlim[1]
	maxD = myArgs$xlim[2]
  }else{
        if(!is.null(maxD)){
		minD = -maxD
	}
  }

  for(i in 1:n){

    this = cd@x[J(baits[i])]
    
    
    
    this = this[is.na(distSign)==FALSE]

    if (!is.null(maxD)){
       this = this[distSign<=maxD]
    }
    if(!is.null(minD)){
       this = this[distSign>=minD]
    }
     
    if (removeBait2bait){
       this = this[isBait2bait==FALSE]
    }

    setDF(this)
    this = this[order(this$distSign),]

    cols <- rep(bgCol, nrow(this))
    pchs <- rep(bgPch, nrow(this))
    sel1 <- this[,pcol] >=plevel1
    sel2 <- this[,pcol] >=plevel2
    cols[sel2] <- lev2Col ##less stringent first
    cols[sel1] <- lev1Col
    pchs[sel2] <- lev2Pch
    pchs[sel1] <- lev2Pch
    
    title = paste(baits[i], sep="")
    if(plotBaitNames){
         baitName = baitmap[baitmap$V4==baits[i]][, cd@settings$baitmapGeneIDcol, with=FALSE]
         if (length(grep(",",baitName))){
             baitName = gsub("(\\S+,).+","\\1", baitName)
             baitName = paste0(baitName, "...")
         }
         title = paste0(baitName, " (", title, ")")
    }    

    plot(this$distSign, this[,Ncol], xlab="Distance from viewpoint", ylab=Ncol, main=title, col=cols, pch=pchs, ...)
    abline(v=0, col="grey", lwd=1)

    if(plotBprof){
	      lines(this$distSign, this$Bmean, lwd=1, col="darkgrey")
        lines(this$distSign, this$Bmean+1.96*sqrt(this$Bmean+this$Bmean^2/disp), 
              lwd=1, lty=2, col="darkgrey")
    }
  }
  if (!is.null(outfile)){ 
    dev.off()
  }
  baits
}

plotDistFun <- function(cd, ...){
  ##TODO: alternative method where we get the observed values too
  params <- cd@params$distFunParams
  
  my.log.d <- seq(from = params$obs.min, to = params$obs.max, length.out = 101)
  my.d <- exp(my.log.d)
  plot(my.log.d, log(.distFun(my.d, params)), type = "l", 
       main = "Distance function estimate", xlab = "log distance", 
       ylab = "log f", col = "Red", ...)
}

exportResults <- function(cd, outfileprefix, scoreCol="score", cutoff=5, b2bcutoff=NULL,
                          format=c("seqMonk","interBed","washU_text"), order=c("position", "score")[1], removeMT=TRUE){
  
  if (any(c("rChr", "rStart", "rEnd", "rID", "bChr", "bStart", "bEnd", "bID") %in% colnames(cd@x))){
    stop ("Colnames x shouldn't contain rChr, rStart, rEnd, rID, bChr, bStart, bEnd, bSign, bID\n") 
  }
  if (!all(format %in% c("seqMonk","interBed", "washU_track", "washU_text"))){
    stop ("Format must be either seqMonk, interBed, washU_track or washU_text (or a vector containing several of these)\n")
  }
  if("washU_track" %in% format)
  {
    if (!requireNamespace("Rsamtools", quietly = TRUE)) {
      stop("Package Rsamtools required to export washU_track format.")
    }
  }

  if (! order %in% c("position","score")){
    stop ("Order must be either position (default) or score\n")
  }
  if (! removeMT %in% c(TRUE,FALSE)){
    stop("removeMT must be TRUE or FALSE")
  }
  
  message("Reading the restriction map file...")
  rmap = .readRmap(cd@settings)
  setnames(rmap, "V1", "rChr")
  setnames(rmap, "V2", "rStart")
  setnames(rmap, "V3", "rEnd")
  setnames(rmap, "V4", "otherEndID")
  
  message("Reading the bait map file...")
  baitmap = .readBaitmap(cd@settings)
  
  setnames(baitmap, "V1", "baitChr")
  setnames(baitmap, "V2", "baitStart")
  setnames(baitmap, "V3", "baitEnd")
  setnames(baitmap, cd@settings$baitmapFragIDcol, "baitID")
  setnames(baitmap, cd@settings$baitmapGeneIDcol, "promID")
  
  message("Preparing the output table...")
  
  if (is.null(b2bcutoff)){
    x = cd@x[ get(scoreCol)>=cutoff ]
  }
  else{
    x = cd@x[ (isBait2bait==TRUE & get(scoreCol)>=b2bcutoff ) | 
                ( isBait2bait==FALSE & get(scoreCol)>=cutoff )]
  }
  
  x = x[, c("baitID", "otherEndID", "N", scoreCol), with=FALSE]
  
  setkey(x, otherEndID)
  setkey(rmap, otherEndID)
  
  x = merge(x, rmap, by="otherEndID", allow.cartesian = TRUE)
  setkey(x, baitID)
  
  setkey(baitmap, baitID)  
  x = merge(x, baitmap, by="baitID", allow.cartesian = TRUE)
  
  # note that baitmapGeneIDcol has been renamed into "promID" above 
  bm2 = baitmap[,c ("baitID", "promID"), with=FALSE]
  
  setDF(x)
  setDF(bm2)
  
  # this way we can be sure that the new column will be called promID.y  
  out = merge(x, bm2, by.x="otherEndID", by.y="baitID", all.x=TRUE, all.y=FALSE, sort=FALSE)
  out[is.na(out$promID.y), "promID.y"] = "."
  
  out = out[,c("baitChr", "baitStart", "baitEnd", "promID.x", "rChr", "rStart", "rEnd", "otherEndID", scoreCol, "N", "promID.y")]
  
  names(out) = c("bait_chr", "bait_start", "bait_end", "bait_name", "otherEnd_chr", "otherEnd_start", "otherEnd_end", "otherEnd_ID", "score", "N_reads", "otherEnd_name")
  
  out$N_reads [ is.na(out$N_reads) ] = 0
  out$score = round(out$score,2)
  
  if (order=="position"){
    out = out[order(out$bait_chr, out$bait_start, out$bait_end, out$otherEnd_chr, out$otherEnd_start, out$otherEnd_end), ]
  }
  if (order=="score"){
    out = out[order(out$score, decreasing=TRUE), ]
  }
  
  if(removeMT)
  {
    ##Remove mitochondrial DNA
    selMT <- tolower(out$bait_chr) == c("chrmt")
    if(any(selMT))
    {
      out <- out[!selMT,]
    }
  }
  
  out0=out
  
  if ("seqMonk" %in% format){
    message("Writing out for seqMonk...")
    out[,"bait_name"] = gsub(",", "|", out[,"bait_name"], fixed=TRUE)
    
    out$newLineOEChr = paste("\n",out[,"otherEnd_chr"], sep="")    
    out = out[,c("bait_chr", "bait_start", "bait_end", "bait_name", "N_reads", "score", "newLineOEChr", "otherEnd_start", "otherEnd_end", "otherEnd_name", "N_reads", "score")]
    
    write.table(out, paste0(outfileprefix,"_seqmonk.txt"), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
  }  
  if ("interBed" %in% format){
    message("Writing out interBed...")
    out = out0[,c("bait_chr", "bait_start", "bait_end", "bait_name", 
                  "otherEnd_chr", "otherEnd_start", "otherEnd_end", "otherEnd_name", 
                  "N_reads", "score")]
    write.table(out, paste0(outfileprefix,".ibed"), sep="\t", quote=FALSE, row.names=FALSE)
  }
  
  if(any(c("washU_text", "washU_track") %in% format)){
    message("Preprocessing for WashU outputs...")
    ##WashU formats
    
    out = out0[,c("bait_chr", "bait_start", "bait_end", 
                  "otherEnd_chr", "otherEnd_start", "otherEnd_end", "otherEnd_name",
                  "score")]
    
    if(!(tolower(substr(out0$bait_chr[1], 1, 3)) == "chr"))
    {
      out$bait_chr <- paste0("chr", out$bait_chr)
      out$otherEnd_chr <- paste0("chr", out$otherEnd_chr)
    }
    
    ##Bait to bait interactions can be asymmetric in terms of score. Here, we find asymmetric interactions and delete the minimum score
    setDT(x)
    setkey(x, baitID, otherEndID)
    x$ReversedInteractionScore <- x[J(x$otherEndID, x$baitID), get(scoreCol)]
    sel <- x[,get(scoreCol)] > x$ReversedInteractionScore
    sel <- ifelse(is.na(sel), TRUE, sel) ##"FALSE" entries in sel should correspond to minima
    setDF(x)
    x$ReversedInteractionScore <- NULL

    if(removeMT)
    {
      sel <- sel[!selMT] ##re-remove mitochondrial interactions (still present in x)
    }
    out <- out[sel,]
    
    if("washU_text" %in% format)
    {
      message("Writing out text file for WashU browser upload...")
      res = paste0(out$bait_chr, ",", out$bait_start, ",", out$bait_end, "\t", out$otherEnd_chr,",", out$otherEnd_start, ",", out$otherEnd_end, "\t", out$score) 
      writeLines(res, con=paste0(outfileprefix,"_washU_text.txt"))
    }
    
    if("washU_track" %in% format)
    {
      message("Writing out track for WashU browser...")
    
      ##this format requires a duplicate of each row, with bait/otherEnd reversed
      out$i = seq(1,nrow(out)*2,2)
      appendOut <- out[,c("otherEnd_chr", "otherEnd_start", "otherEnd_end", "bait_chr", "bait_start", "bait_end", "otherEnd_name", "score", "i")]
      names(appendOut) <- names(out)
      out <- rbind(out, appendOut)
      sel <- order(out$bait_chr, out$bait_start)
      out <- out[sel,]
      
      res = apply(out, 1, function(x){
        lines = paste0(x["bait_chr"], "\t", x["bait_start"], "\t", x["bait_end"], "\t", x["otherEnd_chr"],":", x["otherEnd_start"], "-", x["otherEnd_end"], ",", x["score"],"\t", x["i"], "\t", ".") 
        lines
      })
      res = gsub(" ", "", res)
      writeLines(res, con=paste0(outfileprefix,"_washU_track.txt"))

      Rsamtools::bgzip(paste0(outfileprefix,"_washU_track.txt"),
          dest=paste0(outfileprefix,"_washU_track.txt.gz"))
      Rsamtools::indexTabix(paste0(outfileprefix,"_washU_track.txt.gz"), format="bed")
    }
  }
}

exportToGI <- function(cd, scoreCol="score", cutoff=5, b2bcutoff=NULL,
                       order=c("position", "score")[1], removeMT=TRUE)
{
  if (any(c("rChr", "rStart", "rEnd", "rID", "bChr", "bStart", "bEnd", "bID") %in% colnames(cd@x))){
    stop ("Colnames x shouldn't contain rChr, rStart, rEnd, rID, bChr, bStart, bEnd, bSign, bID\n") 
  }

  if (! order %in% c("position","score")){
    stop ("Order must be either position (default) or score\n")
  }
  if (! removeMT %in% c(TRUE,FALSE)){
    stop("removeMT must be TRUE or FALSE")
  }
  
  message("Reading the restriction map file...")
  rmap = .readRmap(cd@settings)
  setnames(rmap, "V1", "rChr")
  setnames(rmap, "V2", "rStart")
  setnames(rmap, "V3", "rEnd")
  setnames(rmap, "V4", "otherEndID")
  
  message("Reading the bait map file...")
  baitmap = .readBaitmap(cd@settings)
  
  setnames(baitmap, "V1", "baitChr")
  setnames(baitmap, "V2", "baitStart")
  setnames(baitmap, "V3", "baitEnd")
  setnames(baitmap, cd@settings$baitmapFragIDcol, "baitID")
  setnames(baitmap, cd@settings$baitmapGeneIDcol, "promID")
  
  message("Preparing the output table...")
  
  if (is.null(b2bcutoff)){
    x = cd@x[ get(scoreCol)>=cutoff ]
  }
  else{
    x = cd@x[ (isBait2bait==TRUE & get(scoreCol)>=b2bcutoff ) | 
                ( isBait2bait==FALSE & get(scoreCol)>=cutoff )]
  }
  
  x = x[, c("baitID", "otherEndID", "N", scoreCol), with=FALSE]
  
  setkey(x, otherEndID)
  setkey(rmap, otherEndID)
  
  x = merge(x, rmap, by="otherEndID", allow.cartesian = TRUE)
  setkey(x, baitID)
  
  setkey(baitmap, baitID)  
  x = merge(x, baitmap, by="baitID", allow.cartesian = TRUE)
  
  # note that baitmapGeneIDcol has been renamed into "promID" above 
  bm2 = baitmap[,c ("baitID", "promID"), with=FALSE]
  
  setDF(x)
  setDF(bm2)
  
  # this way we can be sure that the new column will be called promID.y  
  out = merge(x, bm2, by.x="otherEndID", by.y="baitID", all.x=TRUE, all.y=FALSE, sort=FALSE)
  out[is.na(out$promID.y), "promID.y"] = "."
  
  out = out[,c("baitChr", "baitStart", "baitEnd", "baitID", "promID.x", "rChr", "rStart", "rEnd", "otherEndID", scoreCol, "N", "promID.y")]
  
  names(out) = c("bait_chr", "bait_start", "bait_end", "bait_ID", "bait_name", "otherEnd_chr", "otherEnd_start", "otherEnd_end", "otherEnd_ID", "score", "N_reads", "otherEnd_name")
  
  out$N_reads [ is.na(out$N_reads) ] = 0
  out$score = round(out$score,2)
  
  if (order=="position"){
    out = out[order(out$bait_chr, out$bait_start, out$bait_end, out$otherEnd_chr, out$otherEnd_start, out$otherEnd_end), ]
  }
  if (order=="score"){
    out = out[order(out$score, decreasing=TRUE), ]
  }
  
  if(removeMT)
  {
    ##Remove mitochondrial DNA
    selMT <- tolower(out$bait_chr) == c("chrmt")
    if(any(selMT))
    {
      out <- out[!selMT,]
    }
  }
  
  #out
  
  ##convert out to a GI
  anchor.one = with(out, GenomicRanges::GRanges(as.character(bait_chr), IRanges::IRanges(start=bait_start, end=bait_end)))
  anchor.two = with(out, GenomicRanges::GRanges(as.character(otherEnd_chr), IRanges::IRanges(start=otherEnd_start, end=otherEnd_end)))
  GenomicInteractions::GenomicInteractions(anchor.one, anchor.two, experiment_name="CHiCAGO calls",
                      counts=out$N_reads, baitName=out$bait_name, otherEndName=out$otherEnd_name,
                      score= out$score)
}

copyCD <- function(cd)
{
  if(class(cd) != "chicagoData") {stop ("cd is not a chicagoData object")}
  newCD <- cd
  newCD@x <- copy(cd@x)
  newCD
}

## Misc functions ---------------------------

wb2b = function(oeID, s, baitmap=NULL){
  # s is the current chicagoData object's settings list
  if (is.null(baitmap)){
    baitmap = .readBaitmap(s)
  }
  which(oeID %in% baitmap[[s$baitmapFragIDcol]])
}

whichbait2bait = function(x, baitmap=NULL){
  stop("whichbait2bait is deprecated. Use wb2b instead")
}

as.dataTableList <- function(cd){ 
  # takes a list of chicagoData objects cd and returns a list of data tables 
  # from the respective cd@x slots
  lapply(cd, function(cdi)cdi@x)
}

geo_mean <- function(data){    
  log_data <- log(data);    
  gm <- exp(mean(log_data[is.finite(log_data)]));    
  return(gm) 
} # http://stackoverflow.com/questions/2602583/geometric-mean-is-there-a-built-in

logit <- function(p){log(p/(1-p))}

expit <- function(x){1/(1+exp(-x))}

removeNAs <- function(x) {x[!is.na(x)]}

naToInf <- function(x)
{
  ifelse(is.na(x), Inf, x) ##Convert NAs to infs.
}

ifnotnull = function(var, res){ if(!is.null(var)){res}}

locateFile = function(what, where, pattern){
  message("Locating ", what, " in ", where, "...")
  filename = list.files(where, pattern)
  
  if (length(filename)!=1){
    stop(paste0("Could not unambigously locate a ", what, " in ", where, ". Please specify explicitly in settings\n"))
  }
  
  message("Found ", filename)
  file.path(where, filename)
}
