peakEnrichment4Features <- function(x1, folder=NULL,  list_frag, no_bins=NULL, sample_number, position_otherEnd= NULL,colname_dist=NULL,
                                    score=5, colname_score="score",min_dist=0, max_dist=NULL,  sep="\t", filterB2B=TRUE, 
                                    b2bcol="isBait2bait", unique=TRUE,plot_name=NULL, trans=FALSE, plotPeakDensity=FALSE) {
  # Check that all features have different names
  if (any(duplicated(names(list_frag)))){
    stop(paste0("Feature(s) ", unique(names(list_frag)[duplicated(names(list_frag))])," appear more than once. Please make sure that each feature appears only once or is given a different name.\n", collapse=","))
  }
  # Extract significant interactions
  # Be aware that you can trim for a specific window
  if (is.null(colname_score)){
    stop("colname_score needs to be specified...\n")
  }
  
  # Check if list_frag is a named vector and gives it names if not:
  if(is.null(names(list_frag))){ names(list_frag) = paste0("Feature", 1:length(list_frag)) }
  
  if (is.data.frame(x1)) {
    warning("Input is not a CHiCAGO object but a dataframe. This requires the explicit specification of parameters that otherwise would be contained in the CHiCAGO object settings...")
    if (!any(c(is.null(position_otherEnd),is.null(colname_dist)))){
      message(paste0("All parameters were specified:\n
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
    message("Input is a CHiCAGO object...")
    position_otherEnd <- x1@settings$rmapfile
    colname_dist <- "distSign"
    x1 = x1@x   
  }
  
  if(filterB2B){
    message("Filtering out bait2bait interactions...")
    x1 <- x1[get(b2bcol)==FALSE]
  }
  
  if(!trans){
    message("Removing trans interactions from analysis...")
    x1 <- x1[!is.na(get(colname_dist))]
    
  } else {
    if(is.null(min_dist)){
      message("Enrichment for features computed for trans interactions only...")
      x1 <- x1[is.na(get(colname_dist))]
    } else {
      message("Enrichment for features computed for trans and cis interactions...")
    }
  }
  
  message("Overlap our reads with Features...")
  x1<- overlapFragWithFeatures(x = x1, position_otherEnd=position_otherEnd, folder = folder, list_frag = list_frag)
  message("Separate significant interactions from non-significant interactions...")
  x1 <- .splitCHiC(x1=x1, threshold=score, colname_score=colname_score, colname_dist=colname_dist, min_dist=min_dist, max_dist=max_dist, trans=trans)
  if (unique){
    message("Removing duplicated other-ends from significant interactions (same will happen with samples)...")
    x1[[1]] <- x1[[1]][!duplicated(otherEndID)]
  }
  
  if (!trans) {
    message("Bin non-significant interactions according to distance from bait before drawing random samples...")
    x1[[2]] <- .binning(sign=x1[[1]], no_bins=no_bins, x1_nonsign=x1[[2]], min_dist=min_dist, max_dist=max_dist)
    message("Draw random samples...")
    if (plotPeakDensity){     
      plot(density(abs(x1[[1]]$distSign)))
    }
    result_3 <- .drawSamples(x1_nonsign=x1[[2]], sample_number=sample_number,unique = unique)
    if (plotPeakDensity){
      j=1
      for (i in c(1,round(sample_number/2),sample_number)){
        lines(density(abs(result_3[[i]]$distSign)),col=c("red","magenta","blue")[j])
        j=j+1
      }
    }   
  } 
  else {
    if(is.null(min_dist)){
      x1[[2]][,distbin3:="trans"]
      x1[[2]][,bin_reads:=nrow(x1[[1]])]
      result_3 <- .drawSamples(x1_nonsign=x1[[2]], sample_number=sample_number,unique = unique)
    } else {
      # trans
      x1[[3]] <- x1[[2]][is.na(x1[[2]]$distSign)]
      x1[[3]][,distbin3:="trans"]
      x1[[3]][,bin_reads:=nrow(x1[[1]][is.na(x1[[1]]$distSign)])] # We want the same number as significant trans-interactions.   
      
      # cis
      x1[[4]] <- .binning(sign=x1[[1]][!is.na(x1[[1]]$distSign)], no_bins=no_bins, x1_nonsign=x1[[2]][!is.na(x1[[2]]$distSign)], min_dist=min_dist, max_dist=max_dist)
      
      # combine trans and cis
      x1[[2]] <- rbindlist(list(x1[[3]],x1[[4]][,names(x1[[3]]),with=FALSE]))
      x1[[3]] <- NULL
      x1[[4]] <- NULL
      
      if (plotPeakDensity){     
        plot(density(abs(x1[[1]][!is.na(x1[[1]]$distSign)]$distSign)))
      }
      result_3 <- .drawSamples(x1_nonsign=x1[[2]], sample_number=sample_number,unique = unique)
      if (plotPeakDensity){
        j=1
        for (i in c(1,round(sample_number/2),sample_number)){
          lines(density(abs(result_3[[i]][!is.na(result_3[[i]]$distSign)]$distSign)),col=c("red","magenta","blue")[j])
          j=j+1
        }
      }   
    }
    
  }
  
  message("Sum number of overlaps with feature in our significant interactions and in our samples...")
  result_5<-.plotNumberOL(x_sign = x1[[1]], s=result_3, files = list_frag,  plot_name=plot_name)
  return(result_5) 
}
# This function bins results assigns probabilities to bins depending on their distance from bait
.binning <- function(sign, no_bins=NULL, x1_nonsign,min_dist=NULL, max_dist=NULL) {
  if(!is.data.table(sign)){setDT(sign)}
  # Bin distances from bait in sign - 100 bins
  
  if (is.null(sign$dist)) {
    sign[,dist := abs(distSign)]
  }
  
  if(is.null(min_dist)){
    min_dist=min(sign$dist)
  }
  
  if(is.null(max_dist)){
    max_dist=max(sign$dist)
  }
  
  if(is.null(no_bins)){
    cat("no_bins not specified.\nParameter will be specified so that bin size is approximately 10kb.\nIn addition, if the computed number falls below 5 bins, no_bins will be set to 5.\n")
    no_bins <- max(5, ceiling((max_dist - min_dist)/1e4))
    cat(paste0("no_bins was set to ",no_bins,".\n"))
  }
  
  
  #sign[,distbin2 := cut(dist, breaks=(no_bins))]
  sign[,distbin2 := cut(dist, breaks=seq(min_dist,max_dist,length.out=no_bins+1),include.lowest = TRUE)]
  # Re-order sign so ensure that the bins between the pools of sign and non-sign interactions will match!
  sign <- sign[order(dist)]
  
  # Calculate how many other-ends in this bin
  bin_reads <- data.table(distbin2=factor(x=levels(sign$distbin2),levels=levels(sign$distbin2)),reads=rep(0,length(levels(sign$distbin2))),key="distbin2")
  
  bin_reads2 <- sign[,length(dist), by="distbin2"]
  setkey(bin_reads2,distbin2)
  bin_reads2[bin_reads]->bin_reads2
  bin_reads2[,bin_reads:=sum(V1,reads,na.rm = TRUE),by=distbin2]
  bin_reads2[,reads:=NULL]
  bin_reads2[,V1:=NULL]

  
  setnames(bin_reads2,"distbin2","udbin2")
  
  # Bin distances from bait in x1_nonsign
  if(!is.data.table(x1_nonsign)){setDT(x1_nonsign)}
  
  if (is.null(x1_nonsign$dist)) {
    x1_nonsign[,dist:=abs(distSign)]
  }
  
  x1_nonsign <- x1_nonsign[order(dist)]
  
  if(min(x1_nonsign$dist)<min_dist){
    x1_nonsign <- x1_nonsign[dist>=min_dist]
  }
  
  if(max(x1_nonsign$dist)>max_dist){
    x1_nonsign <- x1_nonsign[dist<=max_dist]
  }
  
  
  x1_nonsign[,distbin3:= cut(dist, breaks=seq(min_dist,max_dist,length.out=no_bins+1),include.lowest=TRUE)]
  udbin3<-unique(x1_nonsign$distbin3)
  
  # Make sure that udbin2 and udbin3 have the same arguments
  bin_reads2 <- bin_reads2[bin_reads2$udbin2 %in% udbin3]
  
  # Ready to merge           
  bin_reads2[,distbin3:=udbin3]
  
  # Assign to each bin, how many other-ends should be sampled
  setkey(x1_nonsign, distbin3)
  setkey(bin_reads2, distbin3)
  
  x1_nonsign<-x1_nonsign[bin_reads2[,udbin2:=NULL],allow.cartesian=TRUE]
 
  x1_nonsign[is.na(bin_reads),bin_reads:=0]
  
  # Provide correct indexing for non-sign paired-end reads
  x1_nonsign[,iTempVar:=seq(1,nrow(x1_nonsign))]
  
  return(x1_nonsign)
}
.drawSamples <- function(x1_nonsign, sample_number, unique=TRUE) {
  sample_NP <- list()
  if(!is.data.table(x1_nonsign)){setDT(x1_nonsign)}  
  setkey(x1_nonsign,distbin3)
  sample_NP <-  lapply(1:sample_number, function(j) {
    bTempVar <- x1_nonsign[,.I[sample(1:length(.I),bin_reads[1],replace=TRUE)],by="distbin3"]
    s1 <- x1_nonsign[bTempVar$V1]
    
    if(unique){
      s1<-s1[!duplicated(otherEndID)]
    }
    return(s1)
  })
  if (length(sample_NP)<sample_number){
    warning("The Number of samples generated is smaller than the number requested.
        This may cause troubles in the downstream processing.")
  }
  return(sample_NP)
}
overlapFragWithFeatures <- function(x=NULL,folder=NULL, list_frag, position_otherEnd=NULL, sep="\t") {
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
  
  Digest <- .readBedList(folder=NULL, list_frag = c(Digest=position_otherEnd), sep=sep, rm.MT = TRUE,is.Digest=TRUE)[[1]]
  setnames(Digest, names(Digest)[4], "otherEndID")
  
  # Get Features to overlap
  features <- .readBedList(folder=folder, list_frag=list_frag, rm.MT = TRUE)
  
  featuresMapped2Digest<-lapply(features, function(feat) {
    setkeyv(feat, names(feat)[1:3])
    foverlaps(Digest, feat, by.x=names(Digest)[1:3], by.y=names(feat)[1:3], nomatch = 0, mult = "first") # nomatch = 0 equivalent to merge(..., all=F) 
                                                                                      # nomatch = NA is equivalent to merge(..., all=T)
  })
  
  for ( iTempIndex in 1:length(featuresMapped2Digest)) {
    x[, names(featuresMapped2Digest)[iTempIndex]:= otherEndID %in% featuresMapped2Digest[[iTempIndex]]$otherEndID]
  }
  return(x)  
}
.plotNumberOL <- function(x_sign,s, files, plot_name=NULL) {
 
  if("distbin2" %in% names(x_sign)){x_sign[,distbin2:=NULL]}
  if("dist" %in% names(x_sign)){x_sign[,dist:=NULL]}

  x_sign<-colSums(x_sign[,(ncol(x_sign)-length(files)+1):ncol(x_sign),with=FALSE],na.rm = TRUE)
  
  sample_number<- length(s)
  featureSumsMatrix <- matrix(rep(0),length(files)*sample_number,nrow=sample_number,ncol=length(files))
  for (k in 1:sample_number){
    x<-s[[k]]
    if("dist" %in% names(x)){x[,dist:=NULL]}
    if("distbin3" %in% names(x)){x[,distbin3:=NULL]}
    if("bin_reads" %in% names(x)){x[,bin_reads:=NULL]}
    if("iTempVar" %in% names(x)){x[,iTempVar:=NULL]}
    featureSums <- colSums(x[,(ncol(x)-length(files)+1):ncol(x),with=FALSE],na.rm = TRUE)
    featureSumsMatrix[k,]<-featureSums
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
  message("Plotting barplot number of overlaps for features and samples...")
  if(!is.null(plot_name)) {pdf(paste0(plot_name), width=15, height=15)}
  
  data <- cbind(x_sign, result3)
  d <- as.matrix(data[,c(1,2)])
  
  if(any(!d)){
    message("Some (or all) overlaps are zero - warnings will be generated.")
  }
  
  toplot <- barplot(t(d), beside=TRUE, main = "Number of interactions in our samples that map to a GF", col=c("lightyellow","lightblue"),
                    legend = c("Significant Reads", "Random Samples"), names.arg=rownames(data), ylab = c("Number of Overlaps with Feature") )
  arrows(toplot[2,], data$Mean, toplot[2,], data$EB_high, length=0.1, angle=90)
  arrows(toplot[2,], data$Mean, toplot[2,], data$EB_low, length= 0.1, angle=90)
  
  if(!is.null(plot_name)) {
    dev.off()
  }
  
  # Return Matrix with Number of overlaps for ou significant interactions dataset and our samples
  colnames(data)<-c("OLwithSI","MeanOLwithSamples", "SDOLwithSample","HigherCI", "LowerCI")
  message("Returning table with results...\n")
  return(data[,c("OLwithSI","MeanOLwithSamples", "SDOLwithSample", "LowerCI","HigherCI")])
}
.readBedList <- function(folder=NULL, list_frag=NULL, sep="\t", header=FALSE, rm.MT = FALSE, is.Digest=FALSE) {
  
  if (is.null(list_frag) ) {
    stop("Please provide list with files of Genomic features to overlap")
  }
  list_names <- names(list_frag) 
  if (!is.null(folder) ) { list_frag = file.path(folder, list_frag)} # note list_frag names are lost
  
  result <- list()
  mt <- c("M", "MT")
  j <- 1
  for (bed in list_frag) {
    if(is.Digest) {Feature <- fread(input=bed, sep=sep, header=header,stringsAsFactors=FALSE)}
    else {Feature <- fread(input=bed, sep=sep, header=header,stringsAsFactors=FALSE, select=1:3)}
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

.splitCHiC <-  function(x1=NULL, filename=NULL, threshold, colname_score, colname_dist=NULL, min_dist=NULL, max_dist=NULL, trans=FALSE) {
  if (is.null(x1) & is.null(filename)) {
    stop("Please provide file with paired-end reads")
  }
  else if (!is.null(filename)) {
    x1 <- read.table(samplefilename, header=TRUE)
  }  
  if (is.null(min_dist)){min_dist=0}
  if(!is.data.table(x1)){setDT(x1)}
  if (!is.null(colname_dist)) {
    if (is.null(max_dist) & min_dist==0) {
      message("No distance from bait to trim sample was provided...")
    }
    else {
      x1 <- x1[,dist := abs(get(colname_dist))]
      if (!is.null(max_dist)) {x1<-x1[dist<=max_dist]}
      if (!(min_dist==0)) {x1<-x1[dist>=min_dist]}
    }
  }
  result <- list(x1[get(colname_score)>=threshold],
                 x1[get(colname_score)<threshold])
  return(result)
}
