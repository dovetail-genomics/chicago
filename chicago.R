# (c) Mikhail Spivakov, Paula Freire-Pritchett, Jonathan Cairns

library(data.table)
library(matrixStats)
library(MASS)
library(Hmisc)
library(Delaporte)

# fileDir = "../Resources"
fileDir = "/bi/group/sysgen/CHIC"
rmapfile= file.path(fileDir, "Digest_Human_HindIII.bed")
baitmapfile= file.path(fileDir, "Digest_Human_HindIII_baits_ID.bed")
nperbinfile = file.path(fileDir, "Digest_Human_HindIII_NperBin.txt")
nbaitsperbinfile = file.path(fileDir, "Digest_Human_HindIII_NbaitsPerBin.txt")
proxOEfile = file.path(fileDir, "proxOE_out.txt")

baitIDcol = "baitID"
otherEndIDcol = "otherEndID"
otherEndLencol = "otherEndLen"
distcol = "distSign"
bincol = "distbin"
Ncol = "N"

baitmapFragIDcol=4 
baitmapGeneIDcol=5

maxLBrownEst = 1.5e6
minFragLen = 150
maxFragLen = 40000
minNPerBait = 250 
binsize=20000
removeAdjacent = TRUE


chicagoPipeline <- function(x, outprefix, pi.rel)
{
  message("\n*** Running normaliseBaits...\n")
  x = normaliseBaits(x)
  
  message("\n*** Running normaliseOtherEnds...\n")
  x = normaliseOtherEnds(x, outfile=paste0(outprefix, "_oeNorm.pdf"))
  
  message("\n*** Running estimateTechnicalNoise...\n")
  x = estimateTechnicalNoise(x, outfile=paste0(outprefix, "_techNoise.pdf"))
  
  message("\n*** Running estimateDistFun...\n")
  f = estimateDistFun(x, outfile=paste0(outprefix, "_distFun.pdf"))
  
  message("\n*** Running estimateBrownianNoise...\n")
  x = estimateBrownianNoise(x, f, subset=1000)
  
  message("\n*** Running getPvals...\n")
  x = getPvals(x)
  
  message("\n*** Running getScores...\n")
  x = getScores(x, relAbundance = pi.rel)
  
  invisible(x)
}

readSample = function(file){
  
  x = fread(file)
  
  if ( (! baitIDcol %in% names(x)) |   (! otherEndIDcol %in% names(x))  
       |  (! Ncol %in% names(x)) |  (! otherEndLencol %in% names(x))  | 
         (! distcol %in% names(x))){
    stop("Named columns baitIDcol = ", baitIDcol, ", otherEndIDcol = ", otherEndIDcol,
         ", Ncol = ", Ncol, ", otherEndLencol = ", otherEndLencol, " and distcol = ", distcol, 
         " must be present in the input file. Change these global parameters if names do not macth\n")
  }
  
  xlen = nrow(x)
  x = x[get(otherEndLencol)>=minFragLen & get(otherEndLencol)<=maxFragLen]
  message("minFragLen = ", minFragLen, " maxFragLen = ", maxFragLen)
  message("Filtered out ", xlen-nrow(x), " interactions involving other ends < minFragLen or > maxFragLen.")

  setkeyv(x, baitIDcol)
    
  ## remove baits that have no observations within the proximal range
  baitlen = length(unique(x[[baitIDcol]])) 
  x = x[, nperbait:=sum(get(Ncol)), by=baitIDcol]
  x = x[nperbait>=minNPerBait]
  message("minNPerBait = ", minNPerBait)
  message("Filtered out ", baitlen-length(unique(x[[baitIDcol]])), " baits with < minNPerBait reads.\n")  
  x$nperbait = NULL

  ## remove adjacent pairs
  if(removeAdjacent){
    x[, isAdjacent:=abs(get(baitIDcol)-get(otherEndIDcol))==1, by=baitIDcol]
    x = x[isAdjacent==FALSE]
    x$isAdjacent = NULL
    message("Removed interactions with fragments adjacent to baits.")
  }
  
  ##remove baits without proximal non-bait2bait interactions
  baitlen = length(unique(x[[baitIDcol]])) 
  x$isBait2bait = FALSE
  x$isBait2bait[whichbait2bait(as.data.frame(x))] = TRUE
  x[, isAllB2BProx:={
    if(!length(.I[abs(distSign)<maxLBrownEst & !is.na(distSign)])){  TRUE  }
    else{ all(isBait2bait[abs(distSign)<maxLBrownEst & !is.na(distSign)]) }
  }, by=baitIDcol]
  x = x[isAllB2BProx==FALSE]
  x$isAllB2BProx = NULL
  
  message("Filtered out ", baitlen-length(unique(x[[baitIDcol]])), " baits without proximal non-Bait2bait interactions\n")  

  invisible(as.data.frame(x))
}

mergeSamples = function(x, normalise = TRUE, NcolOut="N", NcolNormPrefix="NNorm",
                        mergeMethod=c("weightedMean", "normMean")[1]){
  
  # If mergeMethod == "weightedMean", NcolOut is the weighted mean of the sample-wise counts
  # adjusted by the samples' respective scaling factors s_k
  # If mergeMethod == "normMean", sample-specific counts are first normalised by dividing by s_k
  # and NcolOut is computed as the mean of these normalised counts.
  
  if (! mergeMethod %in% c("weightedMean", "normMean")){
    stop ("Unknown mergeMethod.\n")
  }
  
  message("Preprocessing samples...")
  
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
  for (i in 1:length(x)){  
    names(x[[i]])[which(names(x[[i]])==Ncol)] = paste(Ncol, i, sep=".")
    
    x[[i]] = data.table(x[[i]])
    
    attr[[i]] = x[[i]]
    # remove the read counts from the "backbone" - the only thing we expect to be different across the samples.  
    attr[[i]][[paste(Ncol, i, sep=".")]] = NULL
    
    x[[i]] = x[[i]][, c(baitIDcol, otherEndIDcol, paste(Ncol, i, sep=".")), with=F]
    setkeyv(x[[i]], c(baitIDcol, otherEndIDcol))
    
  }
  
  attr = rbindlist(attr)
  setkeyv(attr, c(baitIDcol, otherEndIDcol))
  # remove duplicates by key
  attr = unique(attr)
  
  message("Merging samples...")
  xmerge = Reduce(function(...) merge(..., by=c(baitIDcol, otherEndIDcol), all=T), x)
  
  setkeyv(xmerge, c(baitIDcol, otherEndIDcol))
  xmerge = attr[xmerge]
  
  xmerge = as.data.frame(xmerge)    
  
  for (i in 1:length(x)){
    iNcol = paste(Ncol,i, sep=".")    
    xmerge[is.na(xmerge[,iNcol]), iNcol] = 0
  }

  whichN = which(names(xmerge) %in% paste(Ncol,1:length(x), sep=".") )
  
  if (normalise){
    xmerge = normaliseSamples(xmerge, computeNNorm=(mergeMethod=="normMean"), 
                              NcolNormPrefix=NcolNormPrefix) 
  
    if(mergeMethod=="weightedMean"){
      names_sk = paste0("s_k",1:length(x))
      s_ks = sapply(attributes(xmerge)[names_sk], function(x)x[[1]])
      for (k in 1:length(x)){
        xmerge[, paste0("Nsk.",k)] = xmerge[, whichN[k]] * s_ks[k]
      }
      whichNsk = which(names(xmerge) %in% paste("Nsk",1:length(x), sep=".") )
      xmerge[, NcolOut] = round(rowSums(xmerge[,whichNsk])/sum(s_ks))   
      for (k in 1:length(x)){
        xmerge[, paste0("Nsk.",k)] = NULL
      }
    }
    else{
      for (i in 1:length(x)){
        normNcol = paste(NcolNormPrefix,i, sep=".")    
        xmerge[is.na(xmerge[,normNcol]), normNcol] = 0
      }
      whichNNorm =  which( names(xmerge) %in% paste(NcolNormPrefix,1:length(x), sep=".") )  
      xmerge[,NcolOut] = round(rowMeans(xmerge[,whichNNorm]))      
    }
  }
  else{
    xmerge[,NcolOut] = round(rowMeans(xmerge[,whichN]))
  }
  
  # Don't completely "cancel" interactions for which we have observed at least one read somewhere
  xmerge[!xmerge[,NcolOut],NcolOut] = 1
  
  invisible(xmerge)
}

normaliseSamples = function(xs, computeNNorm = T, NcolNormPrefix="NNorm"){
  
  # Compute normalisation factors s_k based on the median number of reads per bait
  # within the distance range (0; maxLBrownEst)  
  # The normalisation factors themselves will be written as *attributes* of the output data frame. 
  # If compute NNorm==T, normalise sample-wise counts by dividing by 
  # the respective s_k, writing the results into NcolNormPrefix.<#sample> column. 
  
  message("Normalising samples...")
  
  if (!is.data.frame(xs)){ 
    stop("xs must be a data frame. If starting from a list of separate samples, use mergeSamples instead\n")
  }
  
  Ncols = grep(paste0("^", Ncol, "\\.(\\d+)"), names(xs), value=T)
  ns = as.numeric(gsub(paste0(Ncol, "\\.(\\d+)"), "\\1", Ncols))
  n = max(ns)

  # cat("n =", n, "\n")
  
  xs = data.table(xs)
  
  # Only using the distance range the distance range (0; maxLBrownEst) for normalisation
  # Saving the full table for output
  xs0 = xs
  xs = xs[abs(get(distcol))<maxLBrownEst]
  
  setkeyv(xs, baitIDcol)  
  Ncols = paste(Ncol, 1:n, sep=".")

  # Get the total number of other ends within the distance range (0; maxLBrownEst) for each bait 
  npb = .readNPBfile()   
  
  message("Computing normalisation scores...")
  npb$ntotpb = rowSums(npb[2:ncol(npb)])
  npb = npb[,c(baitIDcol, "ntotpb")]
  npb = data.table(npb)
  setkeyv(npb, baitIDcol)
  setkeyv(xs, baitIDcol)
  xs = npb[xs]
  
  baitMeans = xs[, ntotpb[1],by=baitID]
  s_kjcols = paste0("s_", 1:n, "j")
  for (k in 1:n){
    baitMeans[[ s_kjcols[k] ]] = xs[,sum(get(Ncols[k]))/ntotpb[1],by=baitID]$V1
  }
  
  baitMeans$geo_mean = apply(baitMeans[, s_kjcols, with=F], 1,
                             function(x)if(any(x==0)){NA}else{geo_mean(x)})
  s_k = vector("numeric")
  for (k in 1:n){
    s_k [ k ] = median(baitMeans[[ s_kjcols[k] ]] / baitMeans [[  "geo_mean" ]], na.rm=T)
  }

  s_kcols = paste0("s_k", 1:n)
  for (k in 1:n){
    attributes(xs0)[s_kcols[k]] = s_k[k]
  }
    
  if (computeNNorm){
    for (k in 1:n){
      incol = paste(Ncol, k, sep=".")
      outcol = paste(NcolNormPrefix, k, sep=".")
      # note these normalised values won't be used if merge method = "weightedMean"
      xs0[, outcol] = round(xs0[[incol]]/s_k[k])
      # do not "cancel" an interaction even if its normalised value rounds down to zero 
      set(xs0, j=outcol, i=which(!xs0[[outcol]] & xs0[[incol]]), value=1)      
    } 
  }
  
  xs0 = as.data.frame(xs0)
  invisible(xs0)
}

normaliseBaits = function(x, normNcol="NNb", adjBait2bait=TRUE, shrink=FALSE, plot=TRUE, outfile=NULL){
  message("Normalising baits...")
  alpha <- attributes(x)$dispersion ##store dispersion (if applicable)
  x = .normaliseFragmentSets(x=x, npb=.readNPBfile(), 
                            viewpoint="bait", idcol=baitIDcol, Ncol=Ncol, adjBait2bait=adjBait2bait, 
                            shrink=shrink, refExcludeSuffix=NULL, plot=plot, outfile=outfile)
  # sort by baitID, otherEndID and move distbin column to the end of the table 
  x[, normNcol] = round(x[,Ncol]/x$s_j)
  x[x[,normNcol]==0, normNcol] = 1 # do not completely "cancel" interactions that have one read 
  x = x[order(x[,baitIDcol], x[, otherEndIDcol]),]
  distbincol = which(names(x)=="distbin")
  othercols = which(names(x)!="distbin")
  x = x[, c(othercols, distbincol)]
  attributes(x)$dispersion <- x ##return dispersion (if applicable)
  invisible(x) 
}

normaliseOtherEnds = function(x, Ncol="NNb", normNcol="NNboe", adjBait2bait=TRUE, 
                              filterTopPercent=0.01, minProxOEPerBin=1000, minProxB2BPerBin=minProxOEPerBin/10,
                              plot=TRUE, outfile=NULL){

  alpha <- attributes(x)$dispersion ##store dispersion (if applicable)
  
  # NON-URGENT TODO: instead of looking at bins defined by trans-counts & isB2B, look at bins defined by 
  # fragment length, mappability, GC content and isB2B; 
  # In this case, the call to .addTLB() will be substituted with something like .addLMGB() [to be written] 
  x = .addTLB(x, adjBait2bait=adjBait2bait, filterTopPercent=filterTopPercent, 
              minProxOEPerBin=minProxOEPerBin, minProxB2BPerBin=minProxB2BPerBin)
  
  xAll = x
  
  x = x[abs(x[,distcol])<=maxLBrownEst & !is.na(x[,distcol]),]
  x = data.table(x)
  
  message("Computing total bait counts...")
  nbpb = .readNbaitsPBfile()
  nbpb = data.table(nbpb)
  
  setkeyv(nbpb, otherEndIDcol)
  setkeyv(x, otherEndIDcol)
  
  nbpb = nbpb[x]
  
  # Compute the sums of observed baits per bin for pools of other ends - 
  # NB: we need to some only once for each other end, but each other end is present
  # more than once for each tlb
  setkeyv(nbpb, c(otherEndIDcol, "distbin"))
  nbpb = unique(nbpb)
  
  nbpbSum = nbpb[, sum(.SD[,as.integer(distbin)+1,with=F]), by=c("tlb","distbin")]
  setkeyv(nbpbSum, c("tlb", "distbin"))
  setkeyv(x, c("tlb", "distbin"))
  x = nbpbSum[x]
  setnames(x, "V1", "ntot")
  
  x = as.data.frame(x)
  xt = x[,c("tlb", baitIDcol, Ncol, "distbin", "ntot")]  
  
  message("Computing scaling factors...")
  
  xtnorm=.normaliseFragmentSets(xt, viewpoint="otherEnd", idcol="tlb", Ncol=Ncol, npb = NULL, shrink=FALSE, 
                                adjBait2bait=FALSE, refExcludeSuffix="B2B")
  
  xtnorm = xtnorm[, c("tlb", "s_i")]  
  xtnorm = data.table(xtnorm)
  setkeyv(xtnorm, "tlb")
  xtnorm = unique(xtnorm) # we don't need any other info than s_i for each tlb and it's one per tlb

  if(plot){
    if (!is.null(outfile)){ pdf(outfile)}
    with(xtnorm, barplot(s_i, names.arg=tlb, 
                         col=sapply(tlb, function(x)ifelse(length(grep("B2B",x)), "darkblue", "red")), xlab="tlb", ylab="s_i"))
    legend("topleft", legend=c("non-B2B", "B2B"),fill=c("red", "darkblue"))
    if (!is.null(outfile)){ dev.off()}
  }

  message("Computing normalised counts...")
    
  xAll = data.table(xAll)
  setkeyv(xAll, "tlb")
  
  xAll = merge(xAll, xtnorm, all.x=T)
  #setnames(xAll, "V1", "s_i")
  
  # if we can't estimate s_i robustly, assume it to be one
  xAll$s_i[is.na(xAll$s_i)] = 1
  
  xAll[[normNcol]] = round(xAll[[Ncol]]/xAll$s_i)
  xAll[[normNcol]][xAll[[normNcol]]==0] = 1
  
  message("Post-processing...")
  setkeyv(xAll, c(baitIDcol, otherEndIDcol)) # to reorder
  
  x = as.data.frame(xAll) # note - redefining x
  othercols = names(x)[(!names(x)%in%c(baitIDcol, otherEndIDcol))]
  x = x[, c(baitIDcol, otherEndIDcol, othercols)]  
  
  attributes(x)$dispersion <- x ##return dispersion (if applicable)
  
  invisible(x)
}

estimateDistFun <- function (x, method="cubic", n.obs.head=10, n.obs.tail=25, logScale=FALSE, outfile=NULL) {
  
  # Take the "refBinMean" column of the data x as f(d_b)
  # then interpolate & extrapolate to get f(d).
  # TODO output extra diagnostic information?
  # TODO optimize with data.table
  
  if (!method %in% c("lm", "cubic")){
    stop ("Unknown method.\n")
  }
  
  # Get f(d_b)
  f.d <- unique(x[!is.na(x$refBinMean),c("distbin", "refBinMean")]) ##delete rows with NAs from baits that are too far away
  f.d <- f.d[order(f.d$refBinMean, decreasing=TRUE),]
  f.d$midpoint <- seq(from=round(binsize/2), by=binsize, length.out=nrow(f.d))
  
  obs.min <- log(min(f.d$midpoint))
  obs.max <- log(max(f.d$midpoint))
  
  if(method == "lm") {
    ##On log-scale, do a linear interpolation.
    ##Linear models applied to first (n.obs.head) observations, and to last (n.obs.tail) observations.
    
    ##Interpolation: Estimate f(d) (NB "rule" parameter = 1 forces NAs outside of range)
    log.f.obs <- approxfun(log(f.d$midpoint), log(f.d$refBinMean), rule=c(1,1))
    
    ##Extrapolation: Fit the "head" and "tail" of f using a linear model
    head.coef <- coefficients(lm(log(refBinMean)~log(midpoint), data = head(f.d, n.obs.head))) ##Fit for small d
    tail.coef <- coefficients(lm(log(refBinMean)~log(midpoint), data = tail(f.d, n.obs.tail))) ##Fit for large d
    
    log.f.head <- function(x, head.coef) {head.coef[1] + x*head.coef[2]}
    log.f.tail <- function(x, tail.coef) {tail.coef[1] + x*tail.coef[2]} ##explicitly stated in case of later change
    
  }
  
  if(method == "cubic") {
    ##Spline - Cubic fit over observed interval, linear fit elsewhere, assume continuity of f(d) & f'(d).
    
    ##cubic fit (quadratic not immensely different TBH)
    f.d.cubic <- lm(log(refBinMean) ~ log(midpoint) + I(log(midpoint)^2) + I(log(midpoint)^3), data = f.d)
    fit <- f.d.cubic$coefficients
    
    ##Interpolation: Estimate f(d) from cubic (NB see "rule" parameter for what to do outside range)
    log.f.obs <- function(x, fit. = fit) {fit.[1] + fit.[2]*x + fit.[3]*(x^2) + fit.[4]*(x^3)}
    
    ##Extrapolation: Fit the "head" and "tail" of f using continuity
    obs.min <- log(min(f.d$midpoint))
    obs.max <- log(max(f.d$midpoint))
    
    beta <- fit[2] + 2*fit[3]*c(obs.min, obs.max) + 3*fit[4]*(c(obs.min, obs.max)^2)
    alpha <- fit[1] + (fit[2] - beta)*c(obs.min, obs.max) + fit[3]*c(obs.min, obs.max)^2 + fit[4]*c(obs.min, obs.max)^3
    
    head.coef <- c(alpha[1], beta[1])
    tail.coef <- c(alpha[2], beta[2])
    
    log.f.head <- function(x, head.coef.=head.coef) {head.coef.[1] + x*head.coef.[2]}
    log.f.tail <- function(x, tail.coef.=tail.coef) {tail.coef.[1] + x*tail.coef.[2]} ##explicitly stated in case of later change
    
  }
  
  ##Put everything together to get the final function
  ##All these appended dots (e.g. "head.conf.") are to avoid errors of form "promise already under evaluation"
  log.f <- function(x, head.coef.=head.coef, tail.coef.=tail.coef, obs.min.=obs.min, obs.max.=obs.max)
  {
    ifelse(x > obs.max.,
           log.f.tail(x, tail.coef.), ##Common case evaluated first
           ifelse(x < obs.min.,
                  log.f.head(x, head.coef.),
                  log.f.obs(x))
    )
  }
  if(logScale)
  {
    f <- log.f
  } else {
    f <- function(x) exp(log.f(log(x)))
  }

  if (!is.null(outfile)){ 
    pdf(outfile)
  }
    curve(log.f.obs, obs.min, obs.max,
          main = paste0("Distance function (points = obs, line = ", method, " fit)"),
          xlab = "log(distance)",
          ylab = "log(f(d))")
    with(f.d, points(log(midpoint), log(refBinMean)))
  if (!is.null(outfile)){ 
    dev.off()
  }
  
  f
}

estimateBrownianNoise <- function(x, distFun, Ncol="N", adjBait2bait=TRUE, subset=NULL, reEstimateMean=FALSE) {
  ##1) Reinstate zeros
  ##2) Add a "Bmean" column to x, giving expected Brownian noise.
  ##3) Calculate dispersion by regressing against "Bmean", added to x as "dispersion" attribute
  ##subset: Since we don't need the entire data set, can just calculate based on a random subset of pairs.
  ##!!NB!! Use set.seed to force subset analysis to be reproducible
  
  if(reEstimateMean) {stop("reEstimateMean=TRUE not implemented yet.")}
  
  siPresent <- "s_i" %in% colnames(x)
  if(siPresent)
  {
    message("s_i factors found - estimating Brownian noise...")
  } else {
    message("s_i factors NOT found - variance will increase, estimating Brownian noise anyway...")
  }
  xAll <- x ##duplicate before we start reinstating zeros
  
  ##Pre-filtering
  ##-------------
  
  if(!is.null(subset))
  {
    if(!class(subset) %in% c("numeric","integer")) {stop("'subset' must be an integer.")}
    
    ##much faster than in situ data.frame calculation
    x <- data.table(x)
    setkeyv(x, baitIDcol)
    if( nrow(x[, .I[1], by=baitIDcol])>subset){ 
       sel.sub <- sort(sample(unique(x[[baitIDcol]]), subset))
       x <- as.data.frame(x[J(sel.sub),])
    }
    else{
       x <- as.data.frame(x)
       subset=NULL
       warning("subset > number of baits in data, so use the full dataset.\n")
    }
  }
  
  ##consider proximal region only...
  sel <- abs(x[,distcol]) < maxLBrownEst ##be careful to omit NAs
  sel <- ifelse(is.na(sel), FALSE, sel)
  x = x[sel,] # will have NA for distal and trans interactions
  
  ##remove bait2bait...
  if (adjBait2bait){
    x$isBait2bait = FALSE
    x$isBait2bait[whichbait2bait(x)] = TRUE
    x = x[!x$isBait2bait,]
  }
  
  ##Reinstate zeros:
  ##----------------
  ##1) Grab precomputed data.table: baitID, otherends in range, distance
  proxOE <- .readProxOEfile()
  
  ##2) Delete censored baits (i.e. those with too few observations). Note: censored fragments,
  ##   censored bait2bait pairs (etc...) already taken care of in pre-computation.
  
  if(!is.null(subset)) {
    ## if we chose a subset of baits, restrict to that (none of these should be censored)
    sel.baits <- sel.sub
  } else {
    sel.baits <- unique(x[,baitIDcol])
  }
    
  setkeyv(proxOE, baitIDcol)
  proxOE <- proxOE[J(sel.baits),]
      
  ##3) Merge with our data, thus reinstating zero pairs.
  x <- as.data.table(x)
  setkeyv(x, c(baitIDcol, otherEndIDcol))
  
  ##(make some lookup tables so we can get s_is, s_js later)
  sjLookup <- unique(x[,c(baitIDcol,"s_j"),with=FALSE])
  setkeyv(sjLookup, baitIDcol)
  if(siPresent)
  {
    siLookup <- unique(x[,c(otherEndIDcol,"s_i"),with=FALSE])
    setkeyv(siLookup, otherEndIDcol)
  }
  
  setkeyv(proxOE, c(baitIDcol, otherEndIDcol))
  
  if(siPresent)
  {
    x <- merge(x, proxOE, all.y=TRUE)[,c(baitIDcol,otherEndIDcol,"s_i","s_j",Ncol,distcol,"dist"), with=FALSE]
  } else {
    x <- merge(x, proxOE, all.y=TRUE)[,c(baitIDcol,"s_j",Ncol,distcol,"dist"), with=FALSE]
  }
  
  ##4) Repopulate table with information...
  ## - 0s in Ncol
  x[[Ncol]] <- ifelse(is.na(x[[Ncol]]), 0, x[[Ncol]])
  ## - s_js
  x$s_j <- sjLookup[J(x[[baitIDcol]])]$s_j
  ## - s_is (if present)
  if(siPresent)
  {
    x$s_i <- siLookup[J(x[[otherEndIDcol]])]$s_i
    if(any(is.na(x$s_i)))
    {
      ##If we don't have any information on a particular other end's s_i then...
      warning("Some other ends did not have s_i factors. Assuming s_i = 1 for these.")
      x$s_i <- ifelse(is.na(x$s_i), 1, x$s_i)
    }
  }
  ## - distances
  ##Sanity check - the distances should agree (modulo rounding)
  if(any(removeNAs(abs(x$dist - abs(x$distSign))) > 1))
  {
    warning("Distances in precomputed ProxOE file did not match distances supplied. Using precomputed distances.")
  }
  x$distSign <- x$dist
  
  ##Calculate Bmeans
  ##----------------
  x$Bmean <- .estimateBMean(x, distFun, reEstimateMean)
  
  ##Fit model
  ##---------
  message("Calculating dispersion...")
  model <- glm.nb(formula= x[[Ncol]] ~ offset(log(x$Bmean)) + 0) 
  
  ##Construct Output
  ##----------------
  
  ##NB Parametrization: var = mu + (mu^2)/dispersion
  attributes(xAll)$dispersion <- model$theta
  
  if(reEstimateMean)
  {
    stop("Not implemented yet")
    ##Basically you need to grab the means from the x object - reestimating means on xAll doesn't work due to 0 truncation.
  } else {
    xAll$Bmean <- .estimateBMean(xAll, distFun, reEstimateMean=FALSE)
  }
  invisible(xAll)
}


.estimateBMean = function(x, distFun, reEstimateMean=FALSE) {
  ##1) Gives a "Bmean" vector of length nrow(x), giving expected Brownian noise.
  
  if(reEstimateMean) {stop("reEstimateMean=TRUE not implemented yet.")}
  ##this is a little tricky because we need to reinstate zeros to do it
  ##Two strategies - grab restriction fragment information & calculate explicitly, or cleverly use Mikhail's precomputed tables somehow
  
  
  ##Calculate means
  ##---------------
  
  if("s_i" %in% colnames(x))
  {
    #message("s_i factors found - estimating means...")
    out <- with(x, s_j*s_i*distFun(abs(distSign)))
  } else {
    #message("s_i factors NOT found - variance will increase, estimating means anyway...")
    out <- with(x, s_j*distFun(abs(distSign)))
  }
  ##distcol == NA for trans-pairs. Thus Bmean = 0.
  out <- ifelse(is.na(x$distSign), 0, out)
  
  if(reEstimateMean)
  {
    ##TODO Reserved for future use
  }
  
  out
}


estimateTechnicalNoise = function(x, Ncol="N", filterTopPercent=0.01, minBaitsPerBin=1000, 
                                  minProxOEPerBin=1000, minProxB2BPerBin=minProxOEPerBin/40, adjBait2bait=T,
                                  plot=TRUE, outfile=NULL){ 

# Estimate technical noise based on mean counts per bin, with bins defined based on trans-counts for baits _and_ other ends 
# Note we need raw read counts for this, as normalisation is done wrt Brownian noise. 
  
# NB: filterTopPercent, minProxOEPerBin, minProxB2BPerBin
# are input parameters for .addTLB (the function for binning other ends) that is only called  
# if tlb's aren't already present in the input - and they will be present if other end normalisation
# has been applied 
    
  message("Estimating technical noise based on trans-counts...")

  alpha <- attributes(x)$dispersion ##store dispersion (if applicable)
  ##TODO: considering turning trans counts into "Inf"

  if (!"tlb" %in% names(x)){
    message("Binning other ends based on trans-counts...")
    x = .addTLB(x, adjBait2bait=adjBait2bait, filterTopPercent=filterTopPercent, 
                minProxOEPerBin=minProxOEPerBin, minProxB2BPerBin=minProxB2BPerBin)
  }
  
  x = data.table(x)
  setkeyv(x, baitIDcol)
  
  message("Binning baits based on observed trans-counts...")
  
  transBaitLen = x[, sum(is.na(get(distcol))), by=baitIDcol] ##Number of trans counts per bait
  setnames(transBaitLen, "V1", "transBaitLen")
  
  transBaitLen$tblb = cut2(transBaitLen$transBaitLen, m=minBaitsPerBin, levels.mean=FALSE)
  
  setkeyv(transBaitLen, baitIDcol)
  setkeyv(x, baitIDcol)
  x = x[transBaitLen]
  
  message("Defining interaction pools and gathering the observed numbers of trans-counts per pool...")
  
  # Getting the observed numbers of trans-counts 
  setkeyv(x, c("tlb", "tblb"))
  Ntrans = x[, { res=table(N[is.na(get(distcol))]) 
                    list(as.numeric(names(res)), as.numeric(res)) 
  }, by=c("tlb", "tblb")]
  setnames(Ntrans, "V1", "N")
  setnames(Ntrans, "V2", "nobs")
  
  message("Computing the total number of possible interactions per pool...")
  message("Preparing the data...", appendLF = F)
  
  # Now adding the zeros based on how many trans-interactions are possible in each (tlb, tblb) bin
  baitmap = fread(baitmapfile)
  rmap = fread(rmapfile)
    
  setnames(rmap, "V1", "chr")
  setnames(baitmap, "V1", "chr")
  setnames(rmap, "V4", otherEndIDcol)
  setnames(baitmap, "V4", baitIDcol)
  
  message(".", appendLF=F)
  
  baitmap = baitmap[,c(baitIDcol, "chr"), with=FALSE]
  rmap = rmap[,c(otherEndIDcol, "chr"), with=FALSE]
  
  setkeyv(baitmap, baitIDcol)
  setkeyv(rmap, otherEndIDcol)

  message(".", appendLF=F)
  
  setkeyv(x, baitIDcol)
  x = baitmap[x]
  setnames(x, "chr", "baitChr")
  setkeyv(x, otherEndIDcol)
  x = rmap[x]
  setnames(x, "chr", "otherEndChr")
  setkey(x, tlb, tblb)
  
  message(".", appendLF=F)
      
#   expand.grid.chr = function(seq1,seq2, chr1, chr2) { 
#     # like expand.grid, expands the combinations of seq1 and seq2, 
#     # but instead outputs the corresponding values of chr1 and chr2.
#     data.table(as.data.frame(
#       cbind(rep.int(chr1, length(seq2)),
#           c(t(matrix(rep.int(chr2, length(seq1)), nrow=length(seq2))))),
#       stringsAsFactors=F), key=c("V1", "V2"))    
#   }

  message("\nProcessing fragment pools", appendLF=FALSE)

  res = x[, {
#     message("(", tlb, ",", tblb, ")\t", appendLF=FALSE)
    message(".", appendLF=FALSE)
        
    baits = unique(get(baitIDcol))
    oes = unique(get(otherEndIDcol))
    
    whichBaitFirstOccur = (1:length(get(baitIDcol)))[!duplicated(get(baitIDcol))] 
    bChr = baitChr[whichBaitFirstOccur]
    
    whichOEFirstOccur = (1:length(get(otherEndIDcol)))[!duplicated(get(otherEndIDcol))] 
    oeChr = otherEndChr[whichOEFirstOccur]

    numPairs=sum(sapply(unique(bChr), function(x)length(bChr[bChr==x])*length(oeChr[oeChr!=x])))-
      length(baits[baits%in%oes]) # for bait2bait interactions, it's possible that some oe's will also be among the baits,  
                                  # in which case each such interaction between pairs of such baits will be counted twice 
    nTrans=sum(get(Ncol)[is.na(get(distcol))])
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
  
  x = data.table(x)
  setkeyv(x, c("tlb", "tblb"))
  setkeyv(res, c("tlb", "tblb"))
  x = res[x]
  
  setkeyv(x, c(baitIDcol, otherEndIDcol)) # re-sort this way
  x = as.data.frame(x)
  
  othernames = names(x)[!names(x) %in% c("tlb", "tblb", "Tmean")]
  x = x[, c(othernames, "tlb", "tblb", "Tmean")] ##tlb, tblb are the classes of the other ends, based on trans counts
    
  attributes(x)$dispersion <- alpha ##return dispersion (if applicable)

  invisible(x)
}

getPvals <- function(x, Ncol="N", outcol="log.p", plot=TRUE){
  ## - Calls p-values

  alpha = attributes(x)$dispersion
  if(is.null(alpha)) {stop("getPvals: 'dispersion' attribute of x not found.")}
  
  ##p-values:
  ##(gives P(X > x-1) = P(X >= x))
  message("Calculating p-values...") 
  
  ##The "ifelse" is because pdelap cannot deal with beta=0.
  ##TODO can probably optimize this:
  x[,outcol] <- ifelse(
    x$Bmean < .Machine$double.eps,
    ppois(x[,Ncol] - 1L, lambda=x$Tmean, lower.tail=FALSE, log.p=TRUE),
    pdelap(x[,Ncol] - 1L, alpha, beta=x$Bmean/alpha, lambda=x$Tmean, lower.tail=FALSE, log.p=TRUE)
  )

  # Large N approximation ---------------------------------------------------
  
  ##In rare cases where pdelap returns Infs, estimate the p-value magnitude
  ##using an NB approximation, through method of moments argument
  ##NaNs occur when pdelap() thinks the p-value is negative (since can have 1 - 1 != 0),
  ##thus these are also approximated.
  sel <- which(is.infinite(x[,outcol]) | is.nan(x[,outcol]))
  if(length(sel) > 0)
  {
    message("Approximating ", length(sel), " very small p-values.")
    x.inf <- x[sel,]
  
    gamma <- with(x.inf, alpha*(1+Tmean/Bmean)^2)
    ##in the case where Bmean << Tmean, gamma becomes Inf, so cap gamma above
    ##(should make very little difference since we are basically Poisson in this case)
    ##Case where Bmean >> Tmean causes no problem.
    gamma <- pmin(gamma, 1e10)
    x.inf[,outcol] <- with(x.inf,
                           pnbinom(x.inf[,Ncol] - 1L, size=gamma, mu=Bmean+Tmean, lower.tail=FALSE, log.p=TRUE)
    )
    x[sel,outcol] <- x.inf[,outcol]
    if(any(is.infinite(x.inf[,outcol]))) {warning("Some log-p-values were infinite.")}
    if(any(is.nan(x.inf[,outcol]))) {warning("Some log-p-values were NaNs.")}
  }
  if(any(is.na(x[,outcol]))) {warning("Some log-p-values were NA.")}
  invisible(x)
}

getScores <- function(x, method="weightedRelative",
                    relAbundance=1E5, includeBait2Bait=TRUE, plot=T, outfile=NULL)
{
  ## - If method="weightedRelative", we divide by weights (Genovese et al 2006)
  ## - Then, use Benjamini-Hochberg to calculate FDRs
  
  ##relAbundance is the "relative abundance" i.e. (pi_1^B)/(pi_1^T)
  ##where pi_1^B = probability of interaction in Brownian-dominated bin
  ##and   pi_1^T = probability of interaction in Technical-dominated bin
  
  ##Note to self: Algebra is on P96 of my lab notebook.
  
  if(!method %in% c("unweighted","weightedRelative")) {stop("method=",method," not recognized.")}
  col <- switch(method, unweighted="log.p", weightedRelative="log.q")
  
  if(!includeBait2Bait)
  {
    message("Removing bait2bait interactions...")
    sel <- whichbait2bait(x)
    if(length(sel) > 0)
    {
      x <- x[-sel,] ##remove all bait2bait interactions here
    }
  }
  
  if(method == "weightedRelative")
  {
    message("Calculating q-values...")
    if(is.null(relAbundance))
    {
      stop("relAbundance estimation not supported yet.")
      ##Tech: look at p-values with distance in range Tdom ---------------------------------------------------------------
      sel.T <- abs(x[,distcol]) > Trange[1] & abs(x[,distcol]) < Trange[2]
      sel.T[is.na(sel.T)] <- FALSE ##could be changed to allow trans counts
      log.p.T <- x$log.p[sel.T]
      
      ##FIXME reinstate zeros?
      #extraZeros.T <- 0; warning("Technical estimate needs fixing")
      
      ##Brownian: look at p-values with distance in range Bdom, but ONLY for baits with "large enough" s_j. --------------
      Tbar = mean(x$Tmean, na.rm=TRUE)
      x.sel <- x[!is.na(x[,distcol]),] ##kill NA distances (i.e. trans pairs)
      x.sel <- x.sel[abs(x.sel[,distcol]) > Brange[1] & abs(x.sel[,distcol]) < Brange[2],] ##restrict to range of interest
      x.sel <- as.data.table(x.sel)
      setkeyv(x.sel, baitIDcol)
      
      ##collect list of baits for which Brownian noise still dominates at large distance
      baits <- x.sel[,c(baitIDcol, "s_j"), with=FALSE]
      baits <- baits[!duplicated(baits$baitID),]
      sel.baits <- baits$baitID[baits$s_j > Tbar/distfun(Brange[2])]
      ProxOE <- .readProxOEfile() ##we'll need this later to reinstate zeros
      
      if(length(sel.baits) < nrow(baits)) ##if any baits fail to pass our criterion...
      {
        ##cut them out
        x.sel <- x.sel[J(sel.baits),] ##observed
        setkeyv(ProxOE, baitIDcol) ##possible
        ProxOE <- ProxOE[J(sel.baits),]
      }
      if(length(sel.baits) == 0) {stop("relAbundance estimation failed. Either f(d) is too low, or all s_j are too low.")}
      
      log.p.B <- x.sel$log.p
      
      ##Calculate the number of zeros to reinstate
      ##(total no of possible observed pairs minus no of obs pairs)
      #extraZeros.B <- sum(ProxOE$dist > Brange[1] & ProxOE$dist < Brange[2]) - nrow(x.sel)
      
      ##Plot histograms -------------------------------------------------------------------------------------------------
      if(plot)
      {
        if (!is.null(outfile)){ pdf(outfile) }
        hist(log.p.B, 1000, main="log(p-values), Brownian zone (except zeros)")
        hist(log.p.T, 1000, main="log(p-values), Technical zone (except zeros)")
        hist(exp(log.p.B), 1000, main="p-values, Brownian zone (except zeros)")
        hist(exp(log.p.T), 1000, main="p-values, Technical zone (except zeros)")
        if (!is.null(outfile)){ dev.off() }        
      }
      
      ##Get abundances --------------------------------------------------------------------------------------------------
      #return(list(B=log.p.B, T=log.p.T)) #debug
      
      relAbundance <- estimateRelAbundance(log.p.T, log.p.B, extraZeros.T, extraZeros.B, cutoff=10000)
      message("Calculated relAbundance=", relAbundance)
    }
    
    ##Calculate eta based on logistic regression approach
    B <- logit(0.99) ##eta(d=0)
    d.B <- maxLBrownEst ##distance where eta=0.5. default: d.b=1.5E6
    A <- B/d.B
    eta <- expit(B - A*abs(x[,distcol]))
    
    eta[is.na(x[,distcol])] <- 0
    
    ##calculate eta.bar:
    ##get genome size
    rmap = fread(rmapfile)
    setkey(rmap, V1)
    temp <- rmap[, max(V3), by="V1"]
    G <- sum(as.numeric(temp[[2]]))
    
    eta.sigma <- 2*sum(1/(1+exp(4000*A*(1:10000) - B))) ##accurate to 4dp ##FIXME hardcoded 6-cutter here!
    eta.bar <- eta.sigma*4000/G
    
    ##Get weights, weight p-values
    x$log.w <- log(1 + eta*(relAbundance - 1)) - log(1 + eta.bar*(relAbundance - 1)) ##weight
    x$log.q <- x$log.p - x$log.w ##weighted p-val
    x$score <- -x$log.q ##final score (may omit log.q in final release)
  }
  
#   ##calculate FDR
#   ##How many hypotheses are we testing? Depends how many fragments we are considering. (algebra on p129 of JMC's lab notebook)
#   N.frag <- nrow(fread(rmapfile))
#   
#   if(includeBait2Bait)
#   {
#     N.hyp <- (N.frag^2 - N.frag)/2
#   }else{
#     N.bait <- nrow(fread(baitmapfile))
#     N.hyp <- (N.frag^2 - N.frag + N.bait - N.bait^2)/2
#   }
# 
#   ##sort pvals and derive threshold
#   pvals <- x[,col]
#   if(any(is.na(pvals)))
#   {
#     warning(sum(is.na(pvals))," ", col," values were NA - assuming that means they are -Inf.") ##i.e. underflow in pdelap()
#     pvals[is.na(pvals)] <- (-Inf)
#   }
#   #message("Correcting ", col," values for FDR...")
#  
#   temp.FDR <- pvals + log(N.hyp) - log(rank(pvals, ties.method="max"))
#   sel <- order(pvals, na.last=FALSE) ##"NA" means -Inf, thus these should be first
#   x$log.FDR[sel] <- rev(cummin(rev(temp.FDR[sel]))) ##for final version, pmin(..., 0) to prevent FDR > 1

  
  x
}

.normaliseFragmentSets = function(x, npb, viewpoint, idcol, Ncol, adjBait2bait=TRUE, shrink=TRUE, 
                          refExcludeSuffix=NULL, plot=TRUE, outfile=NULL){   #minPosBins = 5, 
  
  # The normalisation engine used for normaliseBaits and normaliseOtherEnds
  # "Viewpoint" will be used in the comments for either baits or sets of other ends,
  # depending on the direction of normalisation.
  
  # npb is a data frame containing the number of reads per viewpoint per bin,
  # for viewpoint=="bait" it's just the table from the nperbinfile (read via .readNPBfile),
  # where npb's idcol should match x's idcol.
  # for viewpoit=="otherEnd" the table from the nbaitsperbin file needs preprocessing   
  # to sum over pools of other ends, and this is done in normaliseOtherEnds,
  # with the nbp column added to x before submitting to .normaliseFragmentSets.

  bin=binsize
  
  if (!viewpoint %in% c("bait", "otherEnd")){
    stop("viewpoint must be either \"bait\" or \"otherEnd\"")
  }
  
  if(viewpoint=="bait"){
    scol = "s_j"
  
  	x$distbin = cut(abs(x[,distcol]), seq(0, maxLBrownEst, bin)) # will have NA for distal and trans interactions

  	xAll = data.table(x) # full data table with bait2bait and distal interactions
  
  	if (adjBait2bait){
    	x$isBait2bait = FALSE
    	x$isBait2bait[whichbait2bait(x)] = TRUE
    	x = x[!x$isBait2bait,]
  	}
  
  	x = data.table(x)
  
  	# x is the data table used to compute the scaling factors
	  x = x[is.na(distbin)==FALSE]    
  
  	npb = data.table(npb)
  
  	setkeyv(x, idcol)
  	setkeyv(npb, idcol)
  	# compute the total number of valid other ends for each viewpoint and each bin
  	x = npb[x]
  	setkeyv(x, c(idcol, "distbin"))
	  x[, ntot:=get(paste0("bin", as.integer(distbin)))[1], by=c(idcol, "distbin")]

  }
  else{
	# for other ends, bait2bait adjustment and adding ntot is performed before calling this function 
  	x = data.table(x)
    scol = "s_i"
  }

  message("Computing binwise means...")
  
  x = x[, c(idcol, Ncol, "distbin", "ntot"), with=FALSE]
    
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
    # computation, as they are a minority of them and we expect them to have much higher geo_means
    sbbm[, geomean:=geo_mean(bbm[-grep(paste0(".",refExcludeSuffix,"$"), get(idcol))]), by="distbin"]    
  }
  
  # DEseq-style normalisation
  if (!shrink | viewpoint=="otherEnd"){
    sbbm [, s_iv:=bbm/geomean]
    setkeyv(sbbm, idcol)    
    s_v = sbbm[, median(s_iv[!is.na(s_iv)]), by=idcol]
  }
  else{   # Same but with a Gamma shrinkage of binwise means  
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
    xAll = merge(xAll, gm, all.x=T) 
  }
  
  xAll = as.data.frame(xAll)
  invisible(xAll)
  
}  

.readNPBfile = function(){
  
  # Reads a pre-made text file containing the numbers of fragments per bait per distance bin 
  # within the interval maxl, given binsize.
  # The file can be generated by countNperBin.py and its first line should start with # and 
  # contain the parameter settings used. In addition to maxl and binsize, it defines the 
  # filtering parameters for restriction fragment minsize, maxsize as well as a boolean
  # variable removeb2b specifying whether bait2bait interactions should also not be counted.
  # (In fact, they probably should never be).
  
  message("Reading NperBin file...")
  header = readLines(nperbinfile, n=1)
  params = sapply(sapply(strsplit(header, "\t")[[1]],function(x)strsplit(x,"=")[[1]]), function(x)x[2])
  params = params[2:length(params)]
  names(params) = gsub("(\\S+)=.+", "\\1", names(params))
  minsize = as.numeric(params[["minFragLen"]])
  if (minsize != minFragLen){
    stop("The minFragLen in the NfragPerBin file is not equal to minFragLen defined here. 
         Amend either parameter setting before running the analysis\n")
  }
  maxsize = as.numeric(params[["maxFragLen"]])
  if (maxsize != maxFragLen){
    stop("The maxFragLen in the NfragPerBin file is not equal to maxFragLen defined here. 
         Amend either parameter setting before running the analysis\n")
  }
  maxl = as.numeric(params[["maxLBrownEst"]])
  if (maxl != maxLBrownEst){
    stop("The maxLBrownEst in the NfragPerBin file is not equal to maxLBrownEst defined here. 
         Amend either parameter setting before running the analysis\n")
  }
  binsz = as.numeric(params[["binsize"]]) 
  if (binsz != binsize){
    stop("The binsize in the NfragPerBin file is not equal to binsize defined here. 
         Amend either parameter setting before running the analysis\n")
  }  
  if (params[["removeb2b"]]!="True"){
    stop("The NfragPerBin file must be generated with removeb2b==True\n")
  }
  if ( (params[["removeAdjacent"]]=="True" & !removeAdjacent) | (params[["removeAdjacent"]]!="True" & removeAdjacent)  ){
    stop("The removeAdjacent parameter settings used for generating NfragPerBin file and defined here do not match. 
         Amend either setting before running the analysis\n")
  }  
  if(basename(params[["rmapfile"]]) != basename(rmapfile)){
    stop("Rmap files used for generating the NfragPerBin file and defined here do not match. 
         Amend either setting before running the analysis\n")
  }
## Not checking this for now as we have a mixup of _baits and _baits_ID files used at different times...
#   if(basename(params[["baitmapfile"]]) != basename(baitmapfile)){
#     stop("Bait files used for generating the NfragPerBin file and defined here do not match. 
#          Amend either setting before running the analysis\n")
#   }
  npb = as.data.frame(fread(nperbinfile))
  names(npb)[1] = baitIDcol
  names(npb)[2:ncol(npb)] = paste0("bin", 1:(ncol(npb)-1))
  npb
}

.readNbaitsPBfile = function(){
  
  # Reads a pre-made text file containing the numbers of baits per other end per distance bin 
  # within the interval maxl, given binsize.
  # The file can be generated by countNBaitsPerBin.py and its first line should start with # and 
  # contain the parameter settings used. 
  
  message("Reading NbaitsPerBin file...")
  header = readLines(nbaitsperbinfile, n=1)
  params = sapply(sapply(strsplit(header, "\t")[[1]],function(x)strsplit(x,"=")[[1]]), function(x)x[2])
  params = params[2:length(params)]
  names(params) = gsub("(\\S+)=.+", "\\1", names(params))

  maxl = as.numeric(params[["maxLBrownEst"]])
  if (maxl != maxLBrownEst){
    stop("The maxLBrownEst in the NfragPerBin file is not equal to maxLBrownEst defined here. 
         Amend either parameter setting before running the analysis\n")
  }

  # Currently binsize is called bin, but should correct this
  binsz = as.numeric(params[["binsize"]]) 
  if (binsz != binsize){
    stop("The binsize in the NfragPerBin file is not equal to binsize defined here. 
         Amend either parameter setting before running the analysis\n")
  }  
  
  # Currently not in the file
#   if (params[["removeb2b"]]!="True"){
#     stop("The NfragPerBin file must be generated with removeb2b==True\n")
#   }
#   if ( (params[["removeAdjacent"]]=="True" & !removeAdjacent) | (params[["removeAdjacent"]]!="True" & removeAdjacent)  ){
#     stop("The removeAdjacent parameter settings used for generating NfragPerBin file and defined here do not match. 
#          Amend either setting before running the analysis\n")
#   }  
  
  if(basename(params[["rmapfile"]]) != basename(rmapfile)){
    stop("Rmap files used for generating the NfragPerBin file and defined here do not match. 
         Amend either setting before running the analysis\n")
  }
  
  ## Not checking this for now as we have a mixup of _baits and _baits_ID files used at different times...
  #   if(basename(params[["baitmapfile"]]) != basename(baitmapfile)){
  #     stop("Bait files used for generating the NfragPerBin file and defined here do not match. 
  #          Amend either setting before running the analysis\n")
  #   }
  
  nbpb = as.data.frame(fread(nbaitsperbinfile))
  names(nbpb)[1] = otherEndIDcol
  names(nbpb)[2:ncol(nbpb)] = paste0("bin", 1:(ncol(nbpb)-1))
  nbpb
}


.readProxOEfile <- function(){
  
  # Reads a pre-computed text file that denotes which other ends are in the proximal range relative to each
  # bait, and gives that distance. Note that fragments that are too small/too large have already been removed.
  
  message("Reading ProxOE file...")
  header = readLines(proxOEfile, n=1)
  params = sapply(sapply(strsplit(header, "\t")[[1]],function(x)strsplit(x,"=")[[1]]), function(x)x[2])
  params = params[2:length(params)]
  names(params) = gsub("(\\S+)=.+", "\\1", names(params))
  minsize = as.numeric(params[["minFragLen"]])
  if (minsize != minFragLen){
    stop("The minFragLen in the ProxOE file is not equal to minFragLen defined here. 
         Amend either parameter setting before running the analysis\n")
  }
  maxsize = as.numeric(params[["maxFragLen"]])
  if (maxsize != maxFragLen){
    stop("The maxFragLen in the ProxOE file is not equal to maxFragLen defined here. 
         Amend either parameter setting before running the analysis\n")
  }
  maxl = as.numeric(params[["maxLBrownEst"]])
  if (maxl != maxLBrownEst){
    stop("The maxLBrownEst in the ProxOE file is not equal to maxLBrownEst defined here. 
         Amend either parameter setting before running the analysis\n")
  }
  binsz = as.numeric(params[["binsize"]]) 
  if (binsz != binsize){
    stop("The binsize in the ProxOE file is not equal to binsize defined here. 
         Amend either parameter setting before running the analysis\n")
  }  
  if (params[["removeb2b"]]!="True"){
    stop("The ProxOE file must be generated with removeb2b==True\n")
  }
  if ( (params[["removeAdjacent"]]=="True" & !removeAdjacent) | (params[["removeAdjacent"]]!="True" & removeAdjacent)  ){
    stop("The removeAdjacent parameter settings used for generating ProxOE file and defined here do not match. 
         Amend either setting before running the analysis\n")
  }  
  if(basename(params[["rmapfile"]]) != basename(rmapfile)){
    stop("Rmap files used for generating the ProxOE file and defined here do not match. 
         Amend either setting before running the analysis\n")
  }
  ## Not checking this for now as we have a mixup of _baits and _baits_ID files used at different times...
  #   if(basename(params[["baitmapfile"]]) != basename(baitmapfile)){
  #     stop("Bait files used for generating the ProxOE file and defined here do not match. 
  #          Amend either setting before running the analysis\n")
  #   }
  proxOE = fread(proxOEfile)
  setnames(proxOE, 1:3, c(baitIDcol, otherEndIDcol, "dist"))
  #setkeyv(proxOE, c(baitIDcol, otherEndIDcol))
  proxOE
  }

.addTLB = function(x, adjBait2bait=TRUE, filterTopPercent=0.01, minProxOEPerBin=1000, minProxB2BPerBin=100){
    
  message("Preprocessing input...")
  
  # Checking whether in the input, we had distances at trans-interactions labeled as NA 
  # (as opposed to a dummy maximum distance)
  transNA = FALSE
  if(any(is.na(x[,distcol]))){
    transNA = TRUE
    transD = max(x[!is.na(x[,distcol]),distcol])+binsize
    x[is.na(x[,distcol]), distcol] = transD
  }
  else{
    transD = max(x[,distcol])
    message("Warning: No NAs found in input. Assuming the max distance of ", transD, " is a dummy for trans-counts.")
  }
  
  message("Computing trans-counts...")
  
  if (adjBait2bait){ # if we want special bins for bait2bait interactions 
    if(!"isBait2bait" %in% names(x)){
      x$isBait2bait = FALSE
      x$isBait2bait[whichbait2bait(x)] = TRUE
    }
    x = data.table(x)
    setkeyv(x, otherEndIDcol)
    # we are interested in the minimum distance for interactions that involve each other end
    # to then be able to check how many other ends we are pooling together for each tlb bin   
    transLen = x[, list(length(.I[get(distcol)==transD]), isBait2bait[1], min(get(distcol))), by=otherEndIDcol]
    setnames(transLen, "V2", "isBait2bait")
    setnames(transLen, "V3", distcol)
  }
  else{
    x = data.table(x)
    setkeyv(x, otherEndIDcol)
    transLen = x[, list(length(.I[get(distcol)==transD]), min(get(distcol))), by=otherEndIDcol]    
    setnames(transLen, "V2", distcol)
  }
  
  setnames(transLen, "V1", "transLen")
  
  transLen = as.data.frame(transLen)
  transLen0 = transLen
  transLen = transLen[transLen$transLen<quantile(transLen$transLen,1-filterTopPercent/100), ] 
  # apparently data.tables don't like columns with the same name as themselves, so convert back to a data frame
  
  filteredLen = nrow(transLen0[!transLen0[,otherEndIDcol] %in% transLen[,otherEndIDcol],])
  message("Filtering out ", filteredLen, " other ends with top ", filterTopPercent, "% number of trans-interactions")
  
  # first use cut2 to compute bin boundaries on the proximal range based on the desired minProxOEPerBin
  # (note cut2 doesn't guarantee that all bins will contain this min number of observations)
  # then use cut to split the whole dataset based on these bins
  # note we'll need the full dataset (and not only proximal interactions) assigned to tlb bins
  # when estimating technical noise
  # If adjBait2bait == TRUE, do this separately for bait2bait and non-bait2bait other ends.  
  
  if (adjBait2bait){
    transLen0 = transLen
    transLen = transLen0[!transLen0$isBait2bait,]
    transLenB2B = transLen0[transLen0$isBait2bait,]
    
  }
  
  message("Binning...")
  
  cuts = cut2(transLen$transLen[abs(transLen[,distcol])<= maxLBrownEst], 
              m=minProxOEPerBin, onlycuts=T)
  # If some other ends that do not feature in any proximal interactions have transLen's outside of the range
  # determined based on the proximal interactions, just move the boundaries of the first or last tlb bin accordingly...
  if (min(cuts)>min(transLen$transLen)){
    cuts[1] = min(transLen$transLen)
  }
  if (max(cuts)<max(transLen$transLen)){
    cuts[length(cuts)] = max(transLen$transLen)
  }
  transLen$tlb = cut(transLen$transLen, breaks=cuts, include.lowest=T)
  
  if (adjBait2bait){
    
    cutsB2B = cut2(transLenB2B$transLen[abs(transLenB2B[,distcol])<= maxLBrownEst], 
                   m=minProxB2BPerBin, onlycuts=T)
    # If some other ends that do not feature in any proximal interactions have transLen's outside of the range
    # determined based on the proximal interactions, just move the boundaries of the first or last tlb bin accordingly...
    if (min(cutsB2B)>min(transLenB2B$transLen)){
      cutsB2B[1] = min(transLenB2B$transLen)
    }
    if (max(cutsB2B)<max(transLenB2B$transLen)){
      cutsB2B[length(cutsB2B)] = max(transLenB2B$transLen)
    }
    transLenB2B$tlb = cut(transLenB2B$transLen, breaks=cutsB2B, include.lowest=T)
    levels(transLenB2B$tlb) = paste0(levels(transLenB2B$tlb), "B2B")
    
    transLen = rbind(transLen, transLenB2B)
  }      
  
  transLen$transLen = NULL
  transLen$isBait2bait = NULL
  transLen[,distcol] = NULL
  
  transLen = data.table(transLen)
  
  setkeyv(x, otherEndIDcol)
  setkeyv(transLen, otherEndIDcol)
  
  x = x[transLen] # note that if mode="even_filtered", we're not just merging, but also trimming x, 
  # removing the interactions with too "sticky" other ends and those mapping to very sparse bins
  
  x = as.data.frame(x)
  
  if(transNA){
    x[x[,distcol]==max(x[,distcol]), distcol] = NA
  }
  
  invisible(x)
}

plotBaits=function(x, pcol="score", Ncol="N", n=16, baits=NULL, plotBaitNames=TRUE, plotBprof=FALSE,      
                              plevel1 = 12, plevel2 =10, outfile=NULL, removeBait2bait=TRUE, 
                              width=20, height=20, maxD=NULL, ...){
  if(plotBaitNames){
    baitmap = fread(baitmapfile)
  }
  if (is.null(baits)){
    baits = sample(unique(x[,baitIDcol]),n)
  }
  else{
    n = length(baits)
  }
 
  if(plotBprof){
    disp = attributes(x)$dispersion 
  }
 
  if (!is.null(outfile)){ 
    pdf(outfile, width=width, height=height)
  }
  if(n>=4){
    par(mfrow=c(4, ceiling(n/4)))
  }
  else{
    par(mfrow=c(1, n))
  }

  x = data.table(x)
  setkeyv(x, baitIDcol)

  for(i in 1:n){

    this = x[get(baitIDcol)==baits[i]]
    this = this[is.na(distSign)==FALSE]

    if (!is.null(maxD)){
       this = this[abs(distSign)<=maxD]
    }
     
    if (removeBait2bait){
       this = this[isBait2bait==FALSE]
    }
    this = as.data.frame(this)
    this = this[order(this[,distcol]),]

    cols <- rep("Black", nrow(this))
    pchs <- rep(1, nrow(this))
    sel1 <- this[,pcol] >=plevel1
    sel2 <- this[,pcol] >=plevel2
    cols[sel2] <- "Blue" ##less stringent first
    cols[sel1] <- "Red"
    pchs[sel1 | sel2] <- 20
    
    title = paste("baitID=", baits[i], sep="")
    if(plotBaitNames){
         title = paste(title, baitmap$V5[baitmap$V4==baits[i]])
    }    

    plot(this[,distcol], this[,Ncol], xlab=distcol, ylab=Ncol, main=title, col=cols, pch=pchs, ...)
    abline(v=0, col="grey", lwd=1)

    if(plotBprof){
	lines(this[,distcol], this$Bmean, lwd=1, col="darkgrey")
        lines(this[,distcol], this$Bmean+1.96*sqrt(this$Bmean+this$Bmean^2/disp), lwd=1, lty=2, col="darkgrey")
    }
  }
  if (!is.null(outfile)){ 
    dev.off()
  }
  baits
}

exportResults = function(x, outfileprefix, pcol="score", cutoff, format=c("seqMonk","interBed","washU"), order=c("position", "score")[1],
                                 b2bcutoff=NULL, abscutoff=T, rmap=NULL, baitmap=NULL, #logp=T, 
                                 b2bcol = "isBait2bait"){
  
  #stop("Not yet supported")
  
  if (any(c("rChr", "rStart", "rEnd", "rID", "bChr", "bStart", "bEnd", "bID") %in% colnames(x))){
    stop ("Colnames x shouldn't contain rChr, rStart, rEnd, rID, bChr, bStart, bEnd, bSign, bID\n") 
  }
  if (! format %in% c("seqMonk","interBed", "washU")){
    stop ("Format must be either seqMonk, interBed or washU (or a vector containing several of these)\n")
  }
  if (! order %in% c("position","score")){
    stop ("Order must be either position (default) or score\n")
  }
  
  if (is.null(rmap)){
    cat("Reading the restriction map file...\n")
    rmap = as.data.frame(read.table(rmapfile))
  }
  names(rmap) = c("rChr", "rStart", "rEnd", "otherEndID")
  
  if (is.null(baitmap)){
    cat("Reading the bait map file...\n")
    baitmap = as.data.frame(fread(baitmapfile))
  }
  names(baitmap)[1:3] = c("bChr", "bStart", "bEnd") 
  names(baitmap)[baitmapGeneIDcol] = "bID"
  
  #if (baitIDformat=="chr_st_end"){
  #  if (!"index" %in% names(baitmap)){
  #    cat("Generating baitmap index...\n")
  #    baitmap[,1] = as.character(baitmap[,1])
  #    baitmap$index = paste(baitmap[,1], baitmap[,2], baitmap[,3], sep="_")  
  #  }
  #}
  cat("Preparing the output table...\n")
  # just in case
  baitmap = baitmap[!duplicated(baitmap$V4),]
  rmap = rmap[!duplicated(rmap$otherEndID),]
  if (abscutoff) { ps = abs(x[,pcol]) }
  else { ps = x[,pcol] }
  
  if (is.null(b2bcutoff)){
    x = x[ ps>=cutoff, ]
  }
  else{
    x = x[ ( (!x[,b2bcol]) & ps>=cutoff ) | ( (x[,b2bcol]) & ps>=b2bcutoff ) , ]
  }
  x = data.table(x)
  rmap = data.table(rmap)
  setnames(x, otherEndIDcol, "otherEndID")
  setkey(x, otherEndID)
  setkey(rmap, otherEndID)
  
  xa = merge(x,rmap, all.x=T, all.y=F, allow.cartesian=T)
  baitmap = data.table(baitmap)
  setkey(xa, baitID)
  
  #if (baitIDformat =="chr_st_end") { 
  #  setnames(baitmap, "index", names(baitmap)[baitmapFragIDcol], "baitID")
  #  setkey(baitmap, baitID)  
  #  xb = merge(xa, baitmap, all.x=T, all.y=F, allow.cartesian=T)
  #}
  #else{
  setnames(baitmap, names(baitmap)[baitmapFragIDcol], "baitID")
  setkey(baitmap, baitID)  
  xb = merge(xa, baitmap, all.x=T, all.y=F, allow.cartesian=T)
  #}
  
  out = as.data.frame(xb)
  
  out = merge(out, as.data.frame(baitmap)[,c(baitmapFragIDcol, baitmapGeneIDcol)], # note that baitmapGeneIDcol has been renamed into "bID" above 
              by.x=otherEndIDcol, by.y=1, all.x=T, all.y=F, sort=F)  # this way we can be sure that the new column will be called bID.y
  
  if (any (is.na(out$bID.y))) { 
    out$bID.y[is.na(out$bID.y)] = "."
  }
  #out = out[order(out$bID.x, out[,otherEndIDcol]),]
  
  out = out[,c("bChr", "bStart", "bEnd", "bID.x", "rChr", "rStart", "rEnd", otherEndIDcol, pcol, Ncol, "bID.y"),]
  names(out) = c("bait_chr", "bait_start", "bait_end", "bait_name", "otherEnd_chr", 
                 "otherEnd_start", "otherEnd_end", "otherEnd_ID", "score", "N_reads", "otherEnd_name")
  
  out$N_reads [ is.na(out$N_reads) ] = 0
  
  out$score = round(out$score,2)
  if (order=="position"){
    out = out[order(out$bait_chr, out$bait_start, out$bait_end, out$otherEnd_chr, out$otherEnd_start, out$otherEnd_end), ]
  }
  if (order=="score"){
    out = out[order(out$score, decreasing=T), ]
  }
  out0=out
   
  if ("seqMonk" %in% format){
    cat("Writing out for seqMonk...\n")
    out[,"bait_name"] = gsub(",", "|", out[,"bait_name"], fixed=T)
    
    #out$star = "*"
    #out$starNewLineOEChr = paste("*\n",out[,"otherEnd_chr"], sep="")
    
    out$newLineOEChr = paste("\n",out[,"otherEnd_chr"], sep="")
    
    out = out[,c("bait_chr", "bait_start", "bait_end", "bait_name", "N_reads", "score", "newLineOEChr", 
                 "otherEnd_start", "otherEnd_end", "otherEnd_name", "N_reads", "score")]
    cat("Writing out for seqMonk...\n")		
    write.table(out, paste0(outfileprefix,"_seqmonk.txt"), sep="\t", quote=F, row.names=F, col.names=F)
    
  }	
  if ("interBed" %in% format){
    cat("Writing out interBed...\n")
    out = out0[,c("bait_chr", "bait_start", "bait_end", "bait_name", 
                 "otherEnd_chr", "otherEnd_start", "otherEnd_end", "otherEnd_name", 
                 "N_reads", "score")]
    write.table(out, paste0(outfileprefix,".ibed"), sep="\t", quote=F, row.names=F)	
  }
  if("washU" %in% format){
   cat("Writing out for washU browser...\n")
   out = out0[,c("bait_chr", "bait_start", "bait_end", "bait_name", 
                 "otherEnd_chr", "otherEnd_start", "otherEnd_end", "otherEnd_name", 
                 "N_reads", "score")]
   out$i = seq(1,nrow(out)*2,2)
   res = apply(out,1,function(x){
         lines = paste0(x["bait_chr"], "\t", x["bait_start"], "\t", x["bait_end"], "\t", x["otherEnd_chr"],":", x["otherEnd_start"], "-", x["otherEnd_end"], ",", x["score"],"\t", x["i"], "\t", ".") 
         if(x["otherEnd_name"]=="."){ # for bait2bait interactions the second line will be created anyway as the results are duplicated in the original dataframe
            lines = paste0(lines,  "\n",
                paste0(x["otherEnd_chr"], "\t", x["otherEnd_start"], "\t", x["otherEnd_end"], "\t", x["bait_chr"],":", x["bait_start"], "-", x["bait_end"], ",", x["score"],"\t", as.numeric(x["i"])+1, "\t", "."))
         } 
         lines
    })
    res = gsub(" ", "", res)
    writeLines(res, con=paste0(outfileprefix,"_washU.txt"))
   }
  
}


whichbait2bait = function(x, baitmap=NULL){
  if (is.null(baitmap)){
    baitmap = as.data.frame(fread(baitmapfile))
  }
  which(x[,otherEndIDcol] %in% baitmap[,baitmapFragIDcol])
}

distbinToDist = function(db, abs=TRUE){
  distSt = gsub("^\\(", "", db)
  distSt = gsub("\\,.+$", "", distSt)
  ifelse(abs, abs(as.numeric(distSt)), as.numeric(distSt))
}
distbinToDistVec = Vectorize(distbinToDist, "db")

geo_mean <-function(data){    
  log_data <- log(data);    
  gm <- exp(mean(log_data[is.finite(log_data)]));    
  return(gm) 
} # http://stackoverflow.com/questions/2602583/geometric-mean-is-there-a-built-in

logit <- function(p){log(p/(1-p))}

expit <- function(x){1/(1+exp(-x))}

removeNAs <- function(x) {x[!is.na(x)]}
