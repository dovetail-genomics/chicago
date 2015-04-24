peakEnrichment4Features <- function(x1=NULL, score="score", colname_score, min_dist=NULL, max_dist=NULL, 
                                    folder=NULL, list_frag=NULL, sep="\t", filterB2B=TRUE, b2bcol="isBait2bait",  
                                    unique=TRUE, no_bins, sample_number,
                                    plot_name=NULL, position_otherEnd= NULL,colname_dist=NULL) {
  # Extract significant interactions
  # Be aware that you can trim for a specific window
  if (is.null(colname_score)){
    stop("colname_score needs to be specified...\n")
  }
  
  # Check if list_frag is a named vector and gives it names if not:
  if(is.null(names(list_frag))){ names(list_frag) = paste0("Feature", 1:length(list_frag)) }
  
  if (is.data.frame(x1)) {
    cat("Input is not a CHiCAGO object but a dataframe.\n
        Warning: This will require the explicit specification of parameters that otherwise would be contained\n
        in the CHiCAGO object settings...\n")
    if (!any(c(is.null(position_otherEnd),is.null(colname_dist)))){
      cat(paste0("All parameters were specified:\n
                 position_otherEnd= ",position_otherEnd,"\n
                 colname_dist = ",colname_dist,"\n
                 colname_score = ",colname_score,"\n)"))
    }
    else {
      if (is.null(position_otherEnd)){
        stop("position_otherEnd needs to be specified...\n")
      }
      if (is.null(colname_dist)){
        stop("colname_dist needs to be specified...\n")
      }
    } 
    setDT(x1)
    } else {
      cat("Input is a CHiCAGO object...\n")
      position_otherEnd <- x1@settings$rmapfile
      colname_dist <- "distSign"
      x1 = x1@x   
    }
  
  if(filterB2B){
    cat("Filtering out bait2bait interactions...\n")
    x1 <- x1[get(b2bcol)==FALSE]
  }

  # negFraction a numeric value indicating that a fraction of non-significant interactions which will be used to draw samples.
  # The default value is 1 in order for the function to use all the non-significant interactions to draw samples. 
  # This parameter is particularly useful when the user wants to query a large interval of distances and
  # the number of non-significant interactions is very large.
  # By setting this parameter to a value smaller than one, the user may reduce the computing time of this function.
#   if (negFraction<1){
#     cat("Taking a fraction of the negative set...\n")
#     x1pos <- x1[colname_score>=score]
#     x1neg <- x1[colname_score<score]
#     if(!distal){
#       negLen = nrow(x1neg)
#       x1neg <- x1neg[sample(1:negLen, ceiling(negLen*negFraction))]
#     }
#     else{
#       setkeyv(x1neg, c("baitID", "distSign"))
#       pairs = x1neg[, .I[1], by=c("baitID", "distSign")]
#       pairs$V1=NULL
#       samp = pairs[sample(1:nrow(pairs), ceiling(nrow(pairs)*negFraction))]
#       setkeyv(samp, c("baitID", "distSign"))  
#       x1neg <- x1neg[samp]   
#     }
#     x1 = rbind(x1pos, x1neg)
#   }
  
  cat("Overlap our reads with Features...\n")
  x1<- overlapFragWithFeatures(x = x1, position_otherEnd=position_otherEnd, folder = folder, list_frag = list_frag)
  cat("Separate significant interactions from non-significant interactions...\n")
  x1 <- .splitCHiC(x1=x1, threshold=score, colname_score=colname_score, colname_dist=colname_dist, min_dist=min_dist, max_dist=max_dist)
  if (unique){
    cat("Removing duplicated other-ends from significant interactions (same will happen with samples)...\n")
    x1[[1]] <- x1[[1]][!duplicated(otherEndID)]
  }

  cat("Bin non-significant interactions according to distance from bait before drawing random samples...\n")
  x1[[2]] <- .binning(sign=x1[[1]], no_bins=no_bins, x1_nonsign=x1[[2]])

  cat("Draw random samples...\n")
  result_3 <- .drawSamples(x1_nonsign=x1[[2]], sample_number=sample_number,unique = unique)
  
  cat("Sum number of overlaps with feature in our significant interactions and in our samples...\n")
  result_5<-.plotNumberOL(x_sign = x1[[1]], s=result_3, files = list_frag,  plot_name=plot_name)
  return(result_5) 
}
