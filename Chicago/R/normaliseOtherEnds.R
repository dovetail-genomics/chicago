normaliseOtherEnds <-
function(cd, Ncol="NNb", normNcol="NNboe", plot=TRUE, outfile=NULL){

  # NON-URGENT TODO: instead of looking at bins defined by trans-counts & isB2B, look at bins defined by 
  # fragment length, mappability, GC content and isB2B; 
  # In this case, the call to .addTLB() will be substituted with something like .addLMGB() [to be written] 

  cd@x = .addTLB(cd)
  x = cd@x[abs(distSign)<=cd@settings$maxLBrownEst & is.na(distSign)==F]
  
  message("Computing total bait counts...")
  nbpb = .readNbaitsPBfile(s=cd@settings)
  
  setkey(nbpb, otherEndID)
  setkey(x, otherEndID)
  
  nbpb = nbpb[x]
  
  # Compute the sums of observed baits per bin for pools of other ends - 
  # NB: we need to some only once for each other end, but each other end is present
  # more than once for each tlb
  setkeyv(nbpb, c("otherEndID", "distbin"))
  nbpb = unique(nbpb)
  
  nbpbSum = nbpb[, sum(get(paste0("bin",as.integer(distbin)))), by=c("tlb","distbin")]
  setkeyv(nbpbSum, c("tlb", "distbin"))
  setkeyv(x, c("tlb", "distbin"))
  x = nbpbSum[x]
  setnames(x, "V1", "ntot")
  
  message("Computing scaling factors...")
  
  x = .normaliseFragmentSets(x, s = cd@settings, viewpoint="otherEnd", idcol="tlb", Ncol=Ncol, npb = NULL, shrink=FALSE, adjBait2bait=FALSE, refExcludeSuffix="B2B")
  
  setkey(x, tlb)
  x = unique(x) # we don't need any other info than s_i for each tlb and it's one per tlb
  set(x, NULL, "NNb", NULL)
  set(x, NULL, "distbin", NULL)
  set(x, NULL, "ntot", NULL)
  
  if(plot){
    if (!is.null(outfile)){ pdf(outfile)}
    with(x, barplot(s_i, names.arg=tlb, col=sapply(tlb, function(x)ifelse(length(grep("B2B",x)), "darkblue", "red")), xlab="tlb", ylab="s_i"))
    legend("topleft", legend=c("non-B2B", "B2B"),fill=c("red", "darkblue"))
    if (!is.null(outfile)){ dev.off()}
  }

  message("Computing normalised counts...")
    
  setkey(cd@x, tlb)
  
  cd@x = merge(cd@x, x, all.x=T, by="tlb")
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
