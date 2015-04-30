library(GenomicRanges)
library(parallel)
library(data.table)

Extract <- function(x1=NULL, filename=NULL, score, colname_score, colname_dist=NULL, beyond_dist=NULL, before_dist=NULL, significant=TRUE) {
  if (is.null(x1) & is.null(filename)) {
    stop("Please provide file with paired-end reads")
  }
  else if (!is.null(filename)) {
    x1 <- read.table(samplefilename, header=TRUE)
  }  
  if (significant) {
    # Significant according to test
    result <- x1[x1[,colname_score]>score,]
  }
  else {
    # NOT Significant according to test
    result <- x1[x1[,colname_score]<score,]
  }
  result$dist <- abs(result$distSign)
  if (!is.null(colname_dist)) {
    if (is.null(before_dist) & is.null(beyond_dist)) {
      stop("Please provide distance from bait to trim sample")
    }
    else {
      if (is.null(result$dist)) {
        result$dist<-abs(result$distSign)
      }
      if (!is.null(before_dist)) {
        result<-result[result$dist<=before_dist,]
      }
      else {
        result<-result[result$dist>=beyond_dist,]
      }
    }
  }  
  return(result)
}

.splitCHiC <-  function(x1=NULL, filename=NULL, threshold, colname_score, colname_dist=NULL, beyond_dist=NULL, before_dist=NULL) {
  if (is.null(x1) & is.null(filename)) {
    stop("Please provide file with paired-end reads")
  }
  else if (!is.null(filename)) {
    x1 <- read.table(samplefilename, header=TRUE)
  }  
  if(!is.data.table(x1)){setDT(x1)}
  if (!is.null(colname_dist)) {
    if (is.null(before_dist) & is.null(beyond_dist)) {
      cat("No distance from bait to trim sample was provided...\n")
    }
    else {
      x1 <- x1[,dist := abs(get(colname_dist))]
      if (!is.null(before_dist)) {x1<-x1[dist<=before_dist]}
      if (!is.null(beyond_dist)) {x1<-x1[dist>=beyond_dist]}
    }
  }
  result <- list(x1[get(colname_score)>=threshold],
                 x1 <- x1[get(colname_score)<threshold])
  return(result)
}


.convertBedFormat2GR <- function(folder=NULL, list_frag=NULL, sep="\t", header=TRUE, rm.MT = FALSE) {
  if (is.null(folder) ) {
    stop("Please provide location for samples")
  }
  if (is.null(list_frag) ) {
    stop("Please provide list with files of Genomic features to overlap")
  }
  
  list_names <- names(list_frag)
  j <- 1
  result <- list()
  for (i in list_frag) {
    i <- paste0(folder,i)
    Feature <- read.table(file=i,sep=sep, header=header,stringsAsFactors=FALSE)
    chrM <- grep(Feature[,1],pattern="M")
    if (rm.MT & length(chrM)>0){
      Feature[-chrM,]->Feature      
    }
    if(length(grep("chr",Feature[1,1]))==0){
      Feature[,1]=paste0("chr",Feature[,1])
    }
    if (ncol(Feature)>3){
      gr <- GRanges( seqnames = Rle(Feature[,1]),
                     ranges=IRanges( start = Feature[,2], end = Feature[,3]),
                     strand= Rle(strand("*"),nrow(Feature)),
                     Feature[,4:ncol(Feature)])
    }
    else {
      gr <- GRanges( seqnames = Rle(Feature[,1]),
                     IRanges( start = Feature[,2], end = Feature[,3]))
      
      
    }
    result[[j]] <- (assign(list_names[j],gr))
    j <- j+1
  }
  return(result)
}

# This function bins results assigns probabilities to bins depending on their distance from bait
.binning <- function(sign, no_bins, x1_nonsign) {
#.binning <- function(sign, no_bins, x1_nonsign, distal) {
  if(!is.data.table(sign)){setDT(sign)}
  # Bin distances from bait in sign - 100 bins
  if (is.null(sign$dist)) {
    # sign$dist<-abs(sign$distSign)
    sign[,dist := abs(distSign)]
  }
  
  # sign$distbin2 <- cut(sign$dist, breaks=(no_bins))
  sign[,distbin2 := cut(dist, breaks=(no_bins))]
  
  # Calculate how many other-ends in this bin
  # sign <- data.table(sign)
  bin_reads2 <- sign[,length(dist), by="distbin2"]
  setnames(bin_reads2,"distbin2","udbin2")
  setnames(bin_reads2,"V1","bin_reads")
  # bin_reads <- tapply(sign$dist, sign$distbin2, length)
  # bin_reads2 <- data.frame(udbin2=names(bin_reads),bin_reads=as.vector(bin_reads))
  
  # Bin distances from bait in x1_nonsign
  if(!is.data.table(x1_nonsign)){setDT(x1_nonsign)}

  if (is.null(x1_nonsign$dist)) {
    x1_nonsign[,dist:=abs(distSign)]
  }
#   if (distal) {
#     x1_nonsign <- x1_nonsign[dist>=min(dist) & dist<=max(dist),]
#   }

  x1_nonsign[,distbin3:= cut(dist, breaks=(no_bins))]
  udbin3<-unique(x1_nonsign$distbin3)
  udbin3<-udbin3[order(udbin3)]
  
  # bin_reads2$udbin3<-udbin3
  bin_reads2[,distbin3:=udbin3]
  
  # Assign to each bin, how many other-ends should be sampled
  # x1_nonsign$bin_reads <- mclapply(x1_nonsign$distbin3, function(x) {bin_reads2$bin_reads[bin_reads2$udbin3==x]}, mc.cores=8)
  # x1_nonsign <- data.table(x1_nonsign)
  setkey(x1_nonsign, distbin3)
  setkey(bin_reads2, distbin3)
  
  x1_nonsign<-x1_nonsign[bin_reads2[,udbin2:=NULL],allow.cartesian=T]
  # x1_nonsign <- as.data.frame(x1_nonsign)
  
                            
  # x1_nonsign$bin_reads[is.na(x1_nonsign$bin_reads)]=0
  x1_nonsign[is.na(bin_reads),bin_reads:=0]
  
  # Provide correct indexing for non-sign paired-end reads
  # x1_nonsign$i <- seq(1,nrow(x1_nonsign))
  x1_nonsign[,i:=seq(1,nrow(x1_nonsign))]
  sign[,distbin2:=NULL]
  return(x1_nonsign)
}


overlapFragWithFeatures <- function(x=NULL,folder=NULL, position_otherEnd=NULL, list_frag, sep="\t", header=TRUE) {
  if (is.null(x)) {
    stop("Missing sample")
  }
  if (is.null(position_otherEnd)) {
    stop("Missing file with other end start and end coordinates")
  }
  if (!is.data.table(x)){setDT(x)}
    
  HindIII <- .convertBedFormat2GR(folder="", list_frag = c(HindIII=position_otherEnd), sep=sep, header=FALSE,rm.MT = T)[[1]]
  
  # Get Features to overlap
  features <- .convertBedFormat2GR(folder=folder, list_frag=list_frag, sep=sep, header=header,rm.MT = T)
  names(features)<-names(list_frag)
  
  featuresMapped2HindIII<-lapply(features, function(i) {
    i2=subsetByOverlaps(HindIII,i,ignore.strand=TRUE)
    #return(as.data.frame(i2))
    i2<-as.data.frame(i2)
    setDT(i2)
    setnames(x = i2,old = "Feature...4.ncol.Feature..",new="otherEndID")
    return(i2)
  })
  n <-names(x)
  for ( i in 1:length(featuresMapped2HindIII)) {
    x[,names(featuresMapped2HindIII)[i]:= otherEndID %in%  featuresMapped2HindIII[[i]]$otherEndID]
  }
  return(x)  
}


.drawSamples <- function(x1_nonsign, sample_number, unique=TRUE) {
  sample_NP <- list()
  if(!is.data.table(x1_nonsign)){setDT(x1_nonsign)}  
  setkey(x1_nonsign,distbin3)
  sample_NP <-  lapply(1:sample_number, function(j) {
    b <- x1_nonsign[,.I[sample(1:length(.I),bin_reads[1],replace=TRUE)],by="distbin3"]
    
    #s1<-as.data.frame(x1_nonsign)[b$V1,]
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

.plotNumberOL <- function(x_sign,s, files, plot_name=NULL) {
#   x_sign<-as.data.frame(x_sign)
#   x_sign$dist<-NULL
#   x_sign<-colSums(x_sign[,(ncol(x_sign)-length(files)+1):ncol(x_sign)],na.rm = T)
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
#     x$dist <- NULL
#     x$distbin3<-NULL
#     x$bin_reads <- NULL
#     x$i <- NULL
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
    cat(paste0("Plot saved under the name ",plot_name," in your working directory...\n"))
  }
  
  # Return Matrix with Number of overlaps for ou significant interactions dataset and our samples
  colnames(data)<-c("OLwithSI","MeanOLwithSamples", "SDOLwithSample","HigherCI", "LowerCI")
  cat("Return Table with results...\n\n")
  return(data[,c("OLwithSI","MeanOLwithSamples", "SDOLwithSample", "LowerCI","HigherCI")])
}


peakEnrichment4Features <- function(x1=NULL, score, colname_score, colname_dist=NULL, 
                                    beyond_dist=NULL, before_dist=NULL, no_bins, sample_number, folder=NULL, 
                                    position_otherEnd= NULL,filename=NULL, list_frag=NULL, 
                                    sep="\t", header=TRUE, plot_name=NULL, distal=FALSE, coldist=NULL, 
                                    unique=TRUE, filterB2B=FALSE, b2bcol="isBait2bait", negFraction = 1) {
  # Extract significant interactions
  # Be aware that you can trim for a specific window
  if (is.null(colname_score)){
    stop("colname_score needs to be specified...\n")
  } 
  
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
    colname_dist <- x1@settings$distcol
    x1 = x1@x   
  }
    
  if(filterB2B){
    cat("Filtering out bait2bait interactions...\n")
    x1 <- x1[b2bcol==FALSE]
  }
  
  #### This part of the code is deprecated ###
  
  if (negFraction<1){
    cat("Taking a fraction of the negative set...\n")
    x1pos <- x1[colname_score>=score]
    x1neg <- x1[colname_score<score]
    if(!distal){
      negLen = nrow(x1neg)
      x1neg <- x1neg[sample(1:negLen, ceiling(negLen*negFraction))]
    }
    else{
      setkeyv(x1neg, c("baitID", coldist))
      pairs = x1neg[, .I[1], by=c("baitID", coldist)]
      pairs$V1=NULL
      samp = pairs[sample(1:nrow(pairs), ceiling(nrow(pairs)*negFraction))]
      setkeyv(samp, c("baitID", coldist))  
      x1neg <- x1neg[samp]	 
      # setDF(x1neg)  ### fixed from as.data.table
    }
    x1 = rbind(x1pos, x1neg)
  }
  
  ######
  
  cat("Overlap our reads with Features...\n")
  x1<- overlapFragWithFeatures(x = x1, folder = featureFolder, position_otherEnd = position_otherEnd,list_frag = list_frag)
  cat("Separate significant interactions from non-significant interactions...\n")
  x1 <- .splitCHiC(x1=x1, filename=filename, threshold=score, colname_score=colname_score, colname_dist=colname_dist, beyond_dist=beyond_dist, before_dist=before_dist)
  if (unique){
    cat("Removing duplicated other-ends from significant interactions (same will happen with samples)...\n")
    x1[[1]] <- x1[[1]][!duplicated(otherEndID)]
  }
  # Bin non-significant interactions according to distance from bait before drawing random samples
  cat("Bin non-significant interactions according to distance from bait before drawing random samples...\n")
  x1[[2]] <- .binning(sign=x1[[1]], no_bins=no_bins, x1_nonsign=x1[[2]], distal=distal)
  # Draw random samples
  cat("Draw random samples...\n")
  result_3 <- .drawSamples(x1_nonsign=x1[[2]], sample_number=sample_number,unique = unique)
  cat("Sum number of overlaps with feature in our significant interactions and in our samples...\n")
  result_5<-plotNumberOL(x_sign = x1[[1]], s=result_3, files = list_frag,  plot_name=plot_name)
  return(result_5)
  
}
