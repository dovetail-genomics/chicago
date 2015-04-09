estimateTechnicalNoise <-
function(cd, plot=TRUE, outfile=NULL){ 

# Estimate technical noise based on mean counts per bin, with bins defined based on trans-counts for baits _and_ other ends 
# Note we need raw read counts for this, as normalisation is done wrt Brownian noise. 
  
# NB: filterTopPercent, minProxOEPerBin, minProxB2BPerBin
# are input parameters for .addTLB (the function for binning other ends) that is only called  
# if tlb's aren't already present in the input - and they will be present if other end normalisation
# has been applied 
    
  message("Estimating technical noise based on trans-counts...")

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
  message("Preparing the data...", appendLF = F)
  
  # Now adding the zeros based on how many trans-interactions are possible in each (tlb, tblb) bin
  baitmap = fread(cd@settings$baitmapfile)
  rmap = fread(cd@settings$rmapfile)
    
  setnames(rmap, "V1", "chr")
  setnames(baitmap, "V1", "chr")
  setnames(rmap, "V4", "otherEndID")
  setnames(baitmap, "V4", "baitID")
  
  message(".", appendLF=F)
  
  baitmap = baitmap[,c("baitID", "chr"), with=FALSE]
  rmap = rmap[,c("otherEndID", "chr"), with=FALSE]
  
  setkey(baitmap, baitID)
  setkey(rmap, otherEndID)

  message(".", appendLF=F)
  
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
    if(!is.null(outfile)){dev.off()}
    par(mfrow=c(1,1))
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
