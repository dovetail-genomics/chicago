estimateBrownianNoise <-
function(cd, distFun, Ncol="N", reEstimateMean=FALSE) {
  ##1) Reinstate zeros
  ##2) Add a "Bmean" column to x, giving expected Brownian noise.
  ##3) Calculate dispersion by regressing against "Bmean", added to x as "dispersion" attribute
  ##subset: Since we don't need the entire data set, can just calculate based on a random subset of baits.
  ##!!NB!! Use set.seed to force subset analysis to be reproducible
  
  s = cd@settings
  adjBait2bait=s$adjBait2bait
  subset=s$brownianNoise.subset
  seed = s$brownianNoise.seed
  maxLBrownEst = s$maxLBrownEst
  
  if (!is.null(seed)){
    set.seed(seed)
  }
  
  if(reEstimateMean) {stop("reEstimateMean=TRUE not implemented yet.")}
  
  siPresent <- "s_i" %in% colnames(cd@x)
  if(siPresent)
  {
    message("s_i factors found - estimating Brownian noise...")
  } else {
    message("s_i factors NOT found - variance will increase, estimating Brownian noise anyway...")
  }
  
  
  ##Pre-filtering: get subset of data, store as x
  ##---------------------------------------------
  
  if(!is.null(subset))
  {
    if(!class(subset) %in% c("numeric","integer")) {stop("'subset' must be an integer.")}
    
    ##much faster than in situ data.frame calculation
    setkey(cd@x, baitID)
    if( nrow(cd@x[, .I[1], by=baitID])>subset){ 
      sel.sub <- sort(sample(unique(cd@x$baitID), subset))
      x <- cd@x[J(sel.sub)]
    }
    else{
      x <- cd@x
      subset=NULL
      warning("subset > number of baits in data, so use the full dataset.\n")
    }
  }
  else{
    x <- cd@x
  }
  
  ##consider proximal region only...
  setkey(x, distSign)
  x = x[abs(distSign)<maxLBrownEst & is.na(distSign)==F,] # will have NA for distal and trans interactions
  
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
  ##1A) Grab precomputed data.table: baitID, otherends in range, distance
  proxOE <- .readProxOEfile(s)
  
  ##1B) Choose some (uncensored) baits. Pick relevant proxOE rows. Note: censored fragments,
  ##   censored bait2bait pairs (etc...) already taken care of in pre-computation of ProxOE.  
  
  if(!is.null(subset)) {
    ## if we chose a subset of baits, restrict to that (none of these should be censored)
    sel.baits <- sel.sub
  } else {
    sel.baits <- unique(x$baitID)
  }
  
  setkey(proxOE, baitID)
  proxOE <- proxOE[J(sel.baits),]
  
  ##1C) Merge with our data, thus reinstating zero pairs.
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
  ##1D) Repopulate table with information...
  
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
      warning("Some other ends did not have s_i factors. Assuming s_i = 1 for these.")
      x[,s_i := ifelse(is.na(s_i), 1, s_i)]
    }
  }
  ## - distances
  ##Sanity check - the distances should agree (modulo rounding)
  if(any(removeNAs(abs(x$dist - abs(x$distSign))) > 1))
  {
    warning("estimateBrownianNoise: Distances in precomputed ProxOE file did not match distances supplied.")
  }
  
  x[, distSign:=dist]
  ##FIXME delete dist column
  
  ##2) Calculate Bmeans
  ##----------------
  x <- .estimateBMean(x, distFun=cd@params$f)
  
  ##3)Fit model
  ##---------
  message("Calculating dispersion...")
  model <- glm.nb(formula= x$N ~ offset(log(x$Bmean)) + 0) 
  
  ##Construct Output
  ##----------------
  
  ##NB Parametrization: var = mu + (mu^2)/dispersion
  cd@params$dispersion <- model$theta
  
  if(reEstimateMean)
  {
    stop("Not implemented yet")
    ##Basically you need to grab the means from the x object - reestimating means on xAll doesn't work due to 0 truncation.
  } else {
    cd@x <- .estimateBMean(cd@x, cd@params$f)
  }
  cd
}
