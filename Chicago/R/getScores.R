getScores <-
function(cd, method="weightedRelative", includeTrans=TRUE, plot=T, outfile=NULL)
{
  ## - If method="weightedRelative", we divide by weights (Genovese et al 2006)
  ##Note to self: Algebra is on P96 of my lab notebook.
  
  x <- cd@x
  set <- cd@settings
  avgFragLen <- .getAvgFragLength(cd) ##average fragment length
  
  if(!method %in% c("unweighted","weightedRelative")) {stop("method=",method," not recognized.")}
  
  if(!includeTrans)
  {
    x <- x[!is.na(x$distSign),] ##Cannot delete row by reference yet?
  }
  
  if(method == "weightedRelative")
  {
    ##Get weights, weight p-values
    message("Calculating p-value weights...")
    x[, log.w:= .getWeights(abs(x$distSign), cd, includeTrans=includeTrans)]
    x[, log.q:= log.p - log.w] ##weighted p-val
    message("Calculating scores...")
    
    ##get score (more interpretable than log.q)
    minval <- .getWeights(0, cd, includeTrans=includeTrans) ##FIXME could be optimized a *lot*.
    x[,score := pmax(- minval - log.q, 0)]
    
  } else {
    stop("Only method='weightedRelative' available currently.")
  }
  cd
}
