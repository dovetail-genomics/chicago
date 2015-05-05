peakEnrichment4Features <- function(x1=NULL, score=5, colname_score="score",
                                    min_dist=NULL, max_dist=NULL, folder=NULL,  list_frag, sep="\t", filterB2B=TRUE, 
                                    b2bcol="isBait2bait", unique=TRUE, no_bins, sample_number,
                                    plot_name=NULL, position_otherEnd= NULL,colname_dist=NULL,trans=FALSE) {
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

  if(!trans){
    cat("Removing trans interactions from analysis...\n")
    x1 <- x1[!is.na(get(colname_dist))]
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
# This function bins results assigns probabilities to bins depending on their distance from bait
.binning <- function(sign, no_bins, x1_nonsign) {
  if(!is.data.table(sign)){setDT(sign)}
  # Bin distances from bait in sign - 100 bins
  if (is.null(sign$dist)) {
    sign[,dist := abs(distSign)]
  }
  
  sign[,distbin2 := cut(dist, breaks=(no_bins))]
  # Re-order sign so ensure that the bins between the pools of sign and non-sign interactions will match!
  sign <- sign[order(dist)]
  
  # Calculate how many other-ends in this bin
  bin_reads2 <- sign[,length(dist), by="distbin2"]
  setnames(bin_reads2,"distbin2","udbin2")
  setnames(bin_reads2,"V1","bin_reads")
  
  # Bin distances from bait in x1_nonsign
  if(!is.data.table(x1_nonsign)){setDT(x1_nonsign)}
  
  if (is.null(x1_nonsign$dist)) {
    x1_nonsign[,dist:=abs(distSign)]
  }
  
  x1_nonsign <- x1_nonsign[order(dist)]
  
  x1_nonsign[,distbin3:= cut(dist, breaks=(no_bins))]
  udbin3<-unique(x1_nonsign$distbin3)
  udbin3<-udbin3[order(udbin3)]
  
  bin_reads2[,distbin3:=udbin3]
  
  # Assign to each bin, how many other-ends should be sampled
  setkey(x1_nonsign, distbin3)
  setkey(bin_reads2, distbin3)
  
  x1_nonsign<-x1_nonsign[bin_reads2[,udbin2:=NULL],allow.cartesian=T]
 
  x1_nonsign[is.na(bin_reads),bin_reads:=0]
  
  # Provide correct indexing for non-sign paired-end reads
  x1_nonsign[,i:=seq(1,nrow(x1_nonsign))]
  return(x1_nonsign)
}
.drawSamples <- function(x1_nonsign, sample_number, unique=TRUE) {
  sample_NP <- list()
  if(!is.data.table(x1_nonsign)){setDT(x1_nonsign)}  
  setkey(x1_nonsign,distbin3)
  sample_NP <-  lapply(1:sample_number, function(j) {
    b <- x1_nonsign[,.I[sample(1:length(.I),bin_reads[1],replace=TRUE)],by="distbin3"]
    s1 <- x1_nonsign[b$V1]
    
    if(unique){
      s1<-s1[!duplicated(otherEndID)]
    }
    return(s1)
  })
  if (length(sample_NP)<sample_number){
    cat("Warning: The Number of samples generated is smaller than the number requested.
        This may cause troubles in the downstream processing.")
  }
  return(sample_NP)
}
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
.plotNumberOL <- function(x_sign,s, files, plot_name=NULL) {
  x_sign[,distbin2:=NULL]
  x_sign[,dist:=NULL]
  x_sign<-colSums(x_sign[,(ncol(x_sign)-length(files)+1):ncol(x_sign),with=FALSE],na.rm = T)
  
  sample_number<- length(s)
  featureSumsMatrix <- matrix(rep(0),length(files)*sample_number,nrow=sample_number,ncol=length(files))
  for (i in 1:sample_number){
    x<-s[[i]]
    x[,dist:=NULL]
    x[,distbin3:=NULL]
    x[,bin_reads:=NULL]
    x[,i:=NULL]
    featureSums <- colSums(x[,(ncol(x)-length(files)+1):ncol(x),with=FALSE],na.rm = T)
    featureSumsMatrix[i,]<-featureSums
  }
  colnames(featureSumsMatrix)<-names(files)
  
  # Calculate Mean, SD, EB_low and EB_high for each row of the dataframe.
  # Store results for all features in a matrix.
  
  Mean <- colMeans(featureSumsMatrix)
  SD <- apply(featureSumsMatrix,2,sd)
  
  EB_high <- Mean + 1.96 *SD
  EB_low <- Mean - 1.96 *SD
  result3 <- data.frame(Mean, SD, EB_high, EB_low)
  
  # Plot results
  cat("Plot barplot number of overlaps for features and samples...\n")
  if(!is.null(plot_name)) {pdf(paste0(plot_name), width=15, height=15)}
  
  data <- cbind(x_sign, result3)
  d <- as.matrix(data[,c(1,2)])
  
  toplot <- barplot(t(d), beside=TRUE, main = "Number of interactions in our samples that map to a GF", col=c("lightyellow","lightblue"),
                    legend = c("Significant Reads", "Random Samples"), names.arg=rownames(data), ylab = c("Number of Overlaps with Feature") )
  arrows(toplot[2,], data$Mean, toplot[2,], data$EB_high, length=0.1, angle=90)
  arrows(toplot[2,], data$Mean, toplot[2,], data$EB_low, length= 0.1, angle=90)
  
  if(!is.null(plot_name)) {
    dev.off()
  }
  
  # Return Matrix with Number of overlaps for ou significant interactions dataset and our samples
  colnames(data)<-c("OLwithSI","MeanOLwithSamples", "SDOLwithSample","HigherCI", "LowerCI")
  cat("Return Table with results...\n\n")
  return(data[,c("OLwithSI","MeanOLwithSamples", "SDOLwithSample", "LowerCI","HigherCI")])
}
.readBedList <- function(folder=NULL, list_frag=NULL, sep="\t", header=FALSE, rm.MT = FALSE) {
  
  if (is.null(list_frag) ) {
    stop("Please provide list with files of Genomic features to overlap")
  }
  list_names <- names(list_frag) 
  if (!is.null(folder) ) { list_frag = file.path(folder, list_frag)} # note list_frag names are lost
  
  result <- list()
  mt <- c("M", "MT")
  j <- 1
  for (bed in list_frag) {
    Feature <- fread(input=bed, sep=sep, header=header,stringsAsFactors=FALSE)
    firstCol <- names(Feature)[1]
    Feature[, c(firstCol):=gsub("chr", "", get(firstCol))]
    if (rm.MT){
      Feature <- Feature[ ! get(firstCol) %in% mt ]
    }
    result[[list_names[j]]] <- Feature
    j <- j+1
  }
  return(result)
}

.splitCHiC <-  function(x1=NULL, filename=NULL, threshold, colname_score, colname_dist=NULL, min_dist=NULL, max_dist=NULL) {
  if (is.null(x1) & is.null(filename)) {
    stop("Please provide file with paired-end reads")
  }
  else if (!is.null(filename)) {
    x1 <- read.table(samplefilename, header=TRUE)
  }  
  if(!is.data.table(x1)){setDT(x1)}
  if (!is.null(colname_dist)) {
    if (is.null(max_dist) & is.null(min_dist)) {
      cat("No distance from bait to trim sample was provided...\n")
    }
    else {
      x1 <- x1[,dist := abs(get(colname_dist))]
      if (!is.null(max_dist)) {x1<-x1[dist<=max_dist]}
      if (!is.null(min_dist)) {x1<-x1[dist>=min_dist]}
    }
  }
  result <- list(x1[get(colname_score)>=threshold],
                 x1[get(colname_score)<threshold])
  return(result)
}
