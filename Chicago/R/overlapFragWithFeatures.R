overlapFragWithFeatures <- function(x=NULL,folder=NULL, list_frag, sep="\t", position_otherEnd=NULL) {
  if (is.null(x)) {
    stop("Missing sample")
  }
  if (is.data.frame(x)) {
    if (is.null(position_otherEnd)){
      stop("position_otherEnd needs to be specified...\n")
    } 
    if (!is.data.table(x)){setDT(x)}
  } else {
    position_otherEnd <- x@settings$rmapfile
    x <- x@x   
  }
  
  Digest <- .readBedList(folder=NULL, list_frag = c(Digest=position_otherEnd), sep=sep, rm.MT = T)[[1]]
  setnames(Digest, names(Digest)[4], "otherEndID")
  
  # Get Features to overlap
  features <- .readBedList(folder=folder, list_frag=list_frag, rm.MT = T)
  
  featuresMapped2Digest<-lapply(features, function(feat) {
    setkeyv(feat, names(feat)[1:3])
    foverlaps(Digest, feat, by.x=names(Digest)[1:3], by.y=names(feat)[1:3], nomatch = 0, mult = "first") # nomatch = 0 equivalent to merge(..., all=F) 
                                                                                      # nomatch = NA is equivalent to merge(..., all=T)
  })
  
  for ( i in 1:length(featuresMapped2Digest)) {
    x[, names(featuresMapped2Digest)[i]:= otherEndID %in% featuresMapped2Digest[[i]]$otherEndID]
  }
  return(x)  
}
