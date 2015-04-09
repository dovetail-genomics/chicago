mergeSamples <-
function(cdl, normalise = TRUE, NcolOut="N", NcolNormPrefix="NNorm", 
                        mergeMethod=c("weightedMean", "mean")[1]){
  
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
    
    Ncol = grep("^N$", names(x[[i]]), value=T)
    if (!length(Ncol)){ 
      #In case names have already been changed to N.<k> - this being data.table
      Ncol = grep("^N\\.", names(x[[i]]), value=T) 
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
  attr = unique(attr)
  
  message("Merging samples...")
  
  xmerge = Reduce(function(...) merge(..., by=c("baitID", "otherEndID"), all=T), x)
  
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
