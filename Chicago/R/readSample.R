readSample <-
function(file, cd){
  
  message(paste("Reading", file))
  
  x = fread(file)
  
  message("Processing input...")
  
  s = cd@settings
  
  if ( (! s$baitIDcol %in% names(x)) |   (! s$otherEndIDcol %in% names(x))  
       |  (! s$Ncol %in% names(x)) |  (! s$otherEndLencol %in% names(x))  | 
         (! s$distcol %in% names(x))){
    stop("Named columns baitIDcol = ", s$baitIDcol, ", otherEndIDcol = ", s$otherEndIDcol,
         ", Ncol = ", s$Ncol, ", otherEndLencol = ", s$otherEndLencol, " and distcol = ", s$distcol, 
         " must be present in the input file. Change these global parameters if names do not macth\n")
  }
  
  setnames(x, s$baitIDcol, "baitID")
  setnames(x, s$otherEndIDcol, "otherEndID")
  setnames(x, s$Ncol, "N")
  setnames(x, s$distcol, "distSign")
  setnames(x, s$otherEndLencol, "otherEndLen")
    
  xlen = nrow(x)
  x = x[otherEndLen %between% c(s$minFragLen,s$maxFragLen)]
  message("minFragLen = ", s$minFragLen, " maxFragLen = ", s$maxFragLen)
  message("Filtered out ", xlen-nrow(x), " interactions involving other ends < minFragLen or > maxFragLen.")

  setkey(x, baitID)
    
  ## remove baits that have no observations within the proximal range
  baitlen = length(unique(x$baitID)) 
  x = x[, nperbait:=sum(N), by=baitID]
  x = x[nperbait>=s$minNPerBait]
  message("minNPerBait = ", s$minNPerBait)
  message("Filtered out ", baitlen-length(unique(x$baitID)), " baits with < minNPerBait reads.\n")  
  set(x, NULL , "nperbait", NULL) # fast remove data.table column  

  ## remove adjacent pairs
  if(s$removeAdjacent){
    x[, isAdjacent:=abs(baitID-otherEndID)==1, by=baitID]
    x = x[isAdjacent==FALSE]
    set(x, NULL, "isAdjacent", NULL)
    message("Removed interactions with fragments adjacent to baits.")
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
  
  message("Filtered out ", baitlen-length(unique(x$baitID)), " baits without proximal non-Bait2bait interactions\n")  

  chicagoData(x=x, settings = cd@settings, params=list()) # creating this object from scracth, so a single cd object could be passed to readSample as input for multiple replicates
}
