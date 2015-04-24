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
  
  Digest <- .convertBedFormat2GR(folder=NULL, list_frag = c(Digest=position_otherEnd), sep=sep, rm.MT = T)[[1]]
  
  # Get Features to overlap
  features <- .convertBedFormat2GR(folder=folder, list_frag=list_frag, rm.MT = T)
  names(features)<-names(list_frag)
  
  featuresMapped2Digest<-lapply(features, function(i) {
    i2=subsetByOverlaps(Digest,i,ignore.strand=TRUE)
    i2<-as.data.frame(i2)
    setDT(i2)
    setnames(x = i2,old = "Feature...4.ncol.Feature..",new="otherEndID")
    return(i2)
  })
  
  for ( i in 1:length(featuresMapped2Digest)) {
    x[,names(featuresMapped2Digest)[i]:= otherEndID %in%  featuresMapped2Digest[[i]]$otherEndID]
  }
  return(x)  
}
