library(GenomicRanges)
library(multicore)
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


convertBedFormat2GR <- function(folder=NULL, list_frag=NULL, sep="\t", header=TRUE) {
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
    Feature <- read.table(file=i,sep=sep, header=header,stringsAsFactors=F)
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


overlapFragWithFeatures <- function(x=NULL,folder=NULL, position_otherEnd_folder=NULL, position_otherEnd_file =NULL, list_frag, sep="\t", header=TRUE) {
  if (is.null(x)) {
    stop("Missing sample")
  }
  if (is.null(position_otherEnd_file)) {
    stop("Missing file with other end start and end coordinates")
  }
  
  # Get HindIII fragments coordinates and IDs
  #position_otherEnd <- read.table(position_otherEnd, header=FALSE)
  #colnames(position_otherEnd)<- c("chr","start","end","ID")
  
  HindIII <- convertBedFormat2GR(folder=position_otherEnd_folder, list_frag = c(HindIII=position_otherEnd_file), sep=sep, header=F)[[1]]
  
  # Get Features to overlap
  features <- convertBedFormat2GR(folder=folder, list_frag=list_frag, sep=sep, header=header)
  names(features)<-names(list_frag)
  
  featuresMapped2HindIII<-lapply(features, function(i) {
    i2=subsetByOverlaps(HindIII,i,ignore.strand=TRUE)
    return(as.data.frame(i2))
  })
  
  n <-names(x)
  for ( i in 1:length(featuresMapped2HindIII)) {
    x[,ncol(x)+1] <- x$otherEndID %in%  featuresMapped2HindIII[[i]][,ncol(featuresMapped2HindIII[[i]])]
  }
  names(x)<-c(n,names(featuresMapped2HindIII))
  return(x)
  
}


drawSamples <- function(x1_nonsign, sample_number, unique=T) {
  sample_NP <- list()
  x1_nonsign<-data.table(x1_nonsign)
  setkey(x1_nonsign,distbin3)
  sample_NP <-  lapply(1:sample_number, function(j) {
    b <- x1_nonsign[,.I[sample(1:length(.I),bin_reads[1],replace=T)],by="distbin3"]
    s1<-as.data.frame(x1_nonsign)[b$V1,]
    if(unique){
      s1<-s1[!duplicated(s1$otherEndID),]
    }
    return(s1)
  })
  if (length(sample_NP)<sample_number){
    browser()
  }
  return(sample_NP)
}

plotNumberOL <- function(x_sign,s, files, plot_name=NULL) {
  x_sign$dist<-NULL
  x_sign<-colSums(x_sign[,(ncol(x_sign)-length(files)+1):ncol(x_sign)])
  sample_number<- length(s)
  featureSumsMatrix <- matrix(rep(0),length(files)*sample_number,nrow=sample_number,ncol=length(files))
  for (i in 1:sample_number){
    x<-s[[i]]
    x$dist <- NULL
    x$bin_reads <- NULL
    x$i <- NULL
    featureSums <- colSums(x[,(ncol(x)-length(files)+1):ncol(x)])
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


peakEnrichment4Features <- function(x1=NULL, filename=NULL, score, colname_score, colname_dist=NULL, beyond_dist=NULL, before_dist=NULL,no_bins, sample_number, 
                                    restriction_enz_file=NULL, folder_samples=NULL, generic_name, folder=NULL, position_otherEnd_folder = "/bi/group/sysgen/CHIC/",
                                    position_otherEnd_file = "Digest_Human_HindIII.bed",list_frag=NULL, sep="\t", header=TRUE, 
                                    plot_name=NULL, distal=FALSE, coldist=NULL, unique=TRUE, filterB2B=FALSE,
                                    b2bcol="isBait2bait", negFraction = 1) {
  # Extract significant interactions
  # Be aware that you can trim for a specific window
  if (any(c("V1", "V2", "V3", "V4") %in% names(x1))){
    stop ("x1 column names cannot be called V1, V2, V3, V4 - sorry...\n")
  }
  
  position_otherEnd <- paste0(position_otherEnd_folder,position_otherEnd_file)
  
  if(filterB2B){
    cat("Filtering out bait2bait interactions...\n")
    x1 <- x1[! x1[,b2bcol], ]
  }
  if (negFraction<1){
    cat("Taking a fraction of the negative set...\n")
    x1pos <- x1[x1[,colname_score]>=score,]
    x1neg <- x1[x1[,colname_score]<score,]
    if(!distal){
      negLen = nrow(x1neg)
      x1neg <- x1neg[sample(1:negLen, ceiling(negLen*negFraction)),]
    }
    else{
      x1neg <- data.table(x1neg)
      setkeyv(x1neg, c("baitID", coldist))
      pairs = x1neg[, .I[1], by=c("baitID", coldist)]
      pairs$V1=NULL
      samp = pairs[sample(1:nrow(pairs), ceiling(nrow(pairs)*negFraction))]
      setkeyv(samp, c("baitID", coldist))  
      x1neg <- x1neg[samp]	 
      x1neg = as.data.frame(x1neg)  ### fixed from as.data.table
    }
    x1 = rbind(x1pos, x1neg)
  }
  
  cat("Overlap our reads with Features")
  x1<- overlapFragWithFeatures(x = x1, folder = featureFolder, position_otherEnd_folder = "/bi/group/sysgen/CHIC/",
                               position_otherEnd_file = "Digest_Human_HindIII.bed",list_frag = list_frag)
  cat("Extract significant interactions...\n")
  result_1 <- Extract(x1=x1, filename=filename, score=score, colname_score=colname_score, colname_dist=colname_dist, beyond_dist=beyond_dist, before_dist=before_dist, significant=TRUE)
  if (unique){
    cat("Removing duplicated other-ends from significant interactions (same will happen with samples)...\n")
    result_1 <- result_1[!duplicated(result_1$otherEndID),]
  }
  # Extract non-significant interactions
  # You are trimming for the same window that was specified for significant interactions
  cat("Extract non-significant interactions...\n")
  result_2 <- Extract(x1=x1, filename=filename, score=score, colname_score=colname_score, colname_dist=colname_dist, beyond_dist=beyond_dist, before_dist=before_dist, significant=FALSE)
  # Bin non-significant interactions according to distance from bait before drawing random samples
  cat("Bin non-significant interactions according to distance from bait before drawing random samples...\n")
  result_2 <- Binning(sign=result_1, no_bins=no_bins, x1_nonsign=result_2, distal=distal)
  # Draw random samples
  cat("Draw random samples...\n")
  result_3 <- drawSamples(x1_nonsign=result_2, sample_number=sample_number)
  cat("Sum number of overlaps with feature in our significant interactions and in our samples...\n")
  result_5<-plotNumberOL(x_sign = result_1, s=result_3, files = list_frag,  plot_name=plot_name)
  return(result_5)
  
}
