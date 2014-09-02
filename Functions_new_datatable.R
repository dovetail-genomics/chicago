# Functions to be used in Mikhail's pipeline

library(GenomicRanges)
library(multicore)
library(data.table)
source('~/pipeline/write2.R')

# Functions to filter poor baits after 
filterPoorBaits <- function(x,Nreads) {
  # Calculate number of reads per bait
  rpb <- tapply(x$N, x$baitID, sum)
  rpb <- data.frame(baitID=names(rpb), Ntotal=rpb)
  # Select poor baits
  # Nreads is the threshold you establish for the minimum number of reads per bait
  poorBaits <- rpb$baitID[rpb$Ntotal<=Nreads]
  x <- x[!x$baitID %in% poorBaits,]
  return(x)
}

     
     
     
#This function will extract the reads which fit above or below a score (the score is in -log)
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

# This function creates a file with the paired-end fragments scored significant according to the test (WB or PL)
Save_sign_reads <- function(sign,restriction_enz_file=NULL, reads_filename=NULL) {
    if (is.null(restriction_enz_file) | is.null(reads_filename) ) {
        stop("Missing file")
        }

    a <- merge(sign,chr1_res, by.x ="otherEndID", by.y="V4")
    chr <- paste0("chr",a$V1)
   aux <- cbind(chr, a[,c("V2","V3","baitID")])
    names(aux)<- c("chr","start","end","baitID")
    write2(aux, reads_filename)
    }

# This function bins results assigns probabilities to bins depending on their distance from bait
Binning <- function(sign, no_bins, x1_nonsign, distal) {
  # Bin distances from bait in sign - 100 bins
  if (is.null(sign$dist)) {
    sign$dist<-abs(sign$distSign)
  }
  
  sign$distbin2 <- cut(sign$dist, breaks=(no_bins))
  
  # Calculate how many other-ends in this bin
  sign <- data.table(sign)
  bin_reads2 <- sign[,length(dist), by="distbin2"]
  setnames(bin_reads2,"distbin2","udbin2")
  setnames(bin_reads2,"V1","bin_reads")
  # bin_reads <- tapply(sign$dist, sign$distbin2, length)
  # bin_reads2 <- data.frame(udbin2=names(bin_reads),bin_reads=as.vector(bin_reads))
  
  # Bin distances from bait in x1_nonsign
  if (is.null(x1_nonsign$dist)) {
    x1_nonsign$dist<-abs(x1_nonsign$distSign)
  }
  if (distal) {
    x1_nonsign <- x1_nonsign[x1_nonsign$dist>=min(sign$dist) & x1_nonsign$dist<=max(sign$dist),]
  }
  
  x1_nonsign$distbin3 <- cut(x1_nonsign$dist, breaks=(no_bins))
  udbin3<-unique(x1_nonsign$distbin3)
  udbin3<-udbin3[order(udbin3)]
  
  # bin_reads2$udbin3<-udbin3
  bin_reads2[,distbin3:=udbin3]

  # Assign to each bin, how many other-ends should be sampled
  # x1_nonsign$bin_reads <- mclapply(x1_nonsign$distbin3, function(x) {bin_reads2$bin_reads[bin_reads2$udbin3==x]}, mc.cores=8)
  x1_nonsign <- data.table(x1_nonsign)
  setkey(x1_nonsign, distbin3)
  setkey(bin_reads2, distbin3)

  x1_nonsign<-x1_nonsign[bin_reads2[,udbin2:=NULL],allow.cartesian=T]
  x1_nonsign <- as.data.frame(x1_nonsign)

  # x1_nonsign$bin_reads <- unlist(x1_nonsign$bin_reads)                              
  x1_nonsign$bin_reads[is.na(x1_nonsign$bin_reads)]=0
  # Provide correct indexing for non-sign paired-end reads
  x1_nonsign$i <- seq(1,nrow(x1_nonsign))
  return(x1_nonsign)
}



# Draw samples
Draw_samples <- function(x1_nonsign, sample_number, restriction_enz_file=NULL, ncores=8, multicoreSampling=T) {
  sample_NP <- list()
  if (is.null(restriction_enz_file) ) {
    stop("Missing file")
  }
  chr1_res <- read.table(restriction_enz_file)
  x1_nonsign<-data.table(x1_nonsign)
  setkey(x1_nonsign,distbin3)
  if (multicoreSampling){
    sample_NP <-  mclapply(1:sample_number, function(j) {
      b <- x1_nonsign[,.I[sample(1:length(.I),bin_reads[1],replace=T)],by="distbin3"]
      s1<-as.data.frame(x1_nonsign)[b$V1,]
      a2 <- merge(s1, chr1_res,  by.x ="otherEndID", by.y="V4")
          chr <- paste0("chr",a2$V1)
          aux2<-data.frame(chr=chr, start=a2$V2, end=a2$V3, baitID=a2$baitID, dist=a2$dist)
          return(aux2)
        }, mc.cores = ncores)
  }
  else{
    sample_NP <-  lapply(1:sample_number, function(j) {
      b <- x1_nonsign[,.I[sample(1:length(.I),bin_reads[1],replace=T)],by="distbin3"]
      s1<-as.data.frame(x1_nonsign)[b$V1,]
      a2 <- merge(s1, chr1_res,  by.x ="otherEndID", by.y="V4")
      chr <- paste0("chr",a2$V1)
      aux2<-data.frame(chr=chr, start=a2$V2, end=a2$V3, baitID=a2$baitID, dist=a2$dist)
      return(aux2)
    })
  }
  if (length(sample_NP)<sample_number){
    browser()
  }
  return(sample_NP)
}


# Save samples which were drawn
# Variable "samples" has to be a list!
Save_samples <- function(folder=NULL, samples) {
    if (is.null(folder) ) {
        stop("Please provide name of folder to store samples")
        }
    for (j in 1:length(samples)) {
        aux2 <- samples[[j]]
        waux2 <- paste0(folder,"ReadsToOverlap_sNP",j,".txt")
        write2(aux2, waux2)
        }
    }

######

# This function takes a list of files with fragments and converts those files into GRanges objects
library("GenomicRanges")
ConvertGR <- function(folder=NULL, list_frag=NULL, sep="\t", header=TRUE) {
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
        Feature <- read.table(file=i,sep=sep, header=header)
        gr <- GRanges( seqnames = Rle(Feature[,1]),
            IRanges( start = Feature[,2], end = Feature[,3]))

        result[[j]] <- (assign(list_names[j],gr))
        j <- j+1
        }
    return(result)
    }

# This Function will give gaps between the Genomic Ranges Object initially created
# This Function is useful when you want to know the segments of the genome that do not contain a given feature
ConvertGR_FindGaps <- function(folder=NULL, list_frag=NULL, sep="\t", header=TRUE) {
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
    Feature <- read.table(file=i,sep=sep, header=header)
    gr <- GRanges( seqnames = Rle(Feature[,1]),
                   IRanges( start = Feature[,2], end = Feature[,3]))
    gr<-gaps(gr) 
    
    result[[j]] <- (assign(list_names[j],gr))
    j <- j+1
  }
  return(result)
}



# This function takes a list of files with fragments and returns how many fragments are in each file
LengthsGR <- function(folder=NULL, list_frag, sep="\t", header=TRUE) {
  list_names <- names(list_frag)
    aux <- ConvertGR(folder=folder, list_frag=list_frag, sep=sep, header=header)
  
    # Can you separate number of features per chromosome?
  
    result <- c()
    for (i in aux) {
        result <- cbind(result,length(i))
        }
    result <- as.vector(result)
    result <- data.frame(list_names,result)
    return(result)
    }


# This function will count overlaps and subsets between given genomic features and a file containing coord of reads
# In this function, the user has to chose which sort of comparison needs to be done

Compare_seq <- function(sample=NULL, folder=NULL, position_otherEnd=NULL, list_frag, sep="\t", header=TRUE, CountOL=FALSE, Subset_GtoF=FALSE, Subset_FtoG=FALSE, distal=FALSE, unique_reads=TRUE) {
    if (is.null(sample)) {
        stop("Missing sample")
        }
    if (is.null(position_otherEnd)) {
        stop("Missing file with other end start and end coordinates")
        }
        position_otherEnd <- read.table(position_otherEnd, header=FALSE)
        colnames(position_otherEnd)<- c("chr","start","end","ID")
    if (is.character(sample)) {
        S <- read.table(sample, sep=sep, header=header)
        }
    else {
        S <- sample
        }
    
    if (distal) {
      S <- merge(S,position_otherEnd, by.x="otherEndID", by.y="ID")
      S<-data.table(S)
      setkeyv(S, coldist)
      start <- S[,min(start),by=coldist]
      start <- start[,V1]
      
      end <- S[,max(end),by=coldist]
      end <- end[,V1]
        
      
      chr <- position_otherEnd$chr[position_otherEnd$start %in% start]
      
      S<-data.frame(chr=chr,start=start,end=end)
    }
    
    # Convert sample into the right format for using GRanges, if necessary
    if (ncol(S)>5) {
        S <- merge(S, position_otherEnd,by.x="otherEndID",by.y="ID")
        S <- S[,c("chr","start","end","baitID")]
        S$chr <- paste0("chr",S$chr)
    }
    # Sample or significant other-ends to overlap
    S <- GRanges( seqnames = Rle(S[,1]),
                  IRanges( start = S[,2], end = S[,3]))
    if (unique_reads){
      S <- unique(S)
    }
    # Genomic features to overlap
    list_names <- names(list_frag)
    x <- ConvertGR(folder=folder, list_frag=list_frag, sep=sep, header=header)
    Lengths <- LengthsGR(folder=folder,list_frag=list_frag, sep=sep, header=header)
    Overlaps<- NULL
    Subsets_GtoF <- NULL
    Subsets_FtoG <- NULL
    # Conditions to determine which sort of overlap is required
    if (CountOL==TRUE) {
        Nu_Overlaps <- sapply(x, function(i) {
            sum(countOverlaps(i,S,ignore.strand=TRUE))
            })
        Overlaps <- data.frame(Lengths, Nu_Overlaps)
        colnames(Overlaps) <- c("Genomic_Features","Total_Number", "Number_OL")
        }
    if (Subset_GtoF==TRUE) {
        Nu_subsets <- sapply(x, function(i) {
            length(subsetByOverlaps(i,S,ignore.strand=TRUE))
            })
        Subsets_GtoF <- data.frame(Lengths, Nu_subsets)
        colnames(Subsets_GtoF) <- c("Genomic_Features","Total_Number", "Number_Subsets_GFvsSR")
        }
    if (Subset_FtoG==TRUE) {
        Nu_subsets2 <- sapply(x, function(i) {
            length(subsetByOverlaps(S,i,ignore.strand=TRUE))
            })
        Subsets_FtoG <- data.frame(Lengths, Nu_subsets2)
        colnames(Subsets_FtoG) <- c("Genomic_Features","Total_Number", "Number_Subsets_SRvsGF")
        }
    result <- list(Overlaps, Subsets_GtoF, Subsets_FtoG)
    if (is.null(result)){
      cat("Trying to return NULL from CompareSeq\n")
      result <- NA
    } 
    return(result)
    }

# This function creates a list with all our samples of non-significant paired-end reads
Read_samples <- function(folder=NULL,no_samples, generic_name, sep="\t", header=TRUE) {
    if (is.null(folder)) {
        stop("Folder with samples is missing")
        }
    samples <- mclapply(1:no_samples, function (i) {
        filetoadd <- paste0( folder,generic_name,i,".txt")
        read.table(filetoadd,sep=sep,header=header)
        })
    return(samples)
    }

# This function is an optimized version of Compare_seq
# In this version, we can get which fragments have overlaps. We are no longer restricted to just numbers.
Compare_seq_vs2 <- function(sample=NULL, folder=NULL, position_otherEnd=NULL, list_frag, sep="\t", header=TRUE, Interested_in_GF=TRUE,CountOL=FALSE, Subset_GtoF=FALSE, Subset_FtoG=FALSE, Counts=TRUE, Which=FALSE) {
  if (is.null(sample)) {
    stop("Missing sample")
  }
  if (is.null(position_otherEnd)) {
    stop("Missing file with other end stat and end coordinates")
  }
  position_otherEnd <- read.table(position_otherEnd, header=FALSE)
  colnames(position_otherEnd)<- c("chr","start","end","ID")
  if (is.character(sample)) {
    S <- read.table(sample, sep=sep, header=header)
  }
  else {
    S <- sample
  }
  # Convert sample into the right format for using GRanges, if necessary
  if (ncol(S)>5) {
    S <- merge(S, position_otherEnd,by.x="otherEndID",by.y="ID")
    S <- S[,c("chr","start","end","baitID")]
    S$chr <- paste0("chr",S$chr)
  }
  # Sample or significant other-ends to overlap
  S <- GRanges( seqnames = Rle(S[,1]),
                IRanges( start = S[,2], end = S[,3]))
  # Genomic features to overlap
  if(Interested_in_GF){
    list_names <- names(list_frag)
    x <- ConvertGR(folder=folder, list_frag=list_frag, sep=sep, header=header)
  }
  else {
    list_names <- names(list_frag)
    list_names<-paste0("None_",list_names)
    x <- ConvertGR_FindGaps(folder=folder, list_frag=list_frag, sep=sep, header=header)
  }
  
  # Lengths <- LengthsGR(folder=folder,list_frag=list_frag, sep=sep, header=header)
  which_OL<-NULL
  which_subset<-NULL
  which_subset2<-NULL
  Overlaps<- NULL
  Subsets_GtoF <- NULL
  Subsets_FtoG <- NULL
  # Conditions to determine which sort of overlap is required
  if (CountOL==TRUE) {
    which_OL<- lapply(x, function(i) {
      countOverlaps(i,S,ignore.strand=TRUE)
    })
    names(which_OL)<- names(list_frag)
    Nu_Overlaps <- sapply(which_OL, function(i) {sum(i)})
    Overlaps <- data.frame(Number_OL=Nu_Overlaps, row.names=names(list_frag))
  }
  if (Subset_GtoF==TRUE) {
    which_subset <- lapply(x, function(i) {
      subsetByOverlaps(i,S,ignore.strand=TRUE)
    })
    names(which_subset)=names(list_frag)
    Nu_subsets <- sapply(which_subset, function(i) {length(i)})
    Subsets_GtoF <- data.frame(Number_Subsets_GFvsSR=Nu_subsets, row.names=names(list_frag))
  }
  if (Subset_FtoG==TRUE) {
    which_subset2 <- lapply(x, function(i) {
      subsetByOverlaps(S,i,ignore.strand=TRUE)
    })
    names(which_subset2)=names(list_frag)
    Nu_subsets2 <- sapply(which_subset2, function(i) {length(i)})
    Subsets_FtoG <- data.frame(Number_Subsets_SRvsGF=Nu_subsets2, row.names=names(list_frag))
  }
  if(Counts) {
    result <- list(Overlaps, Subsets_GtoF, Subsets_FtoG)}
  else if(Which){ result<- list(which_OL, which_subset, which_subset2)}
  return(result)
}






########

# This function takes a list of samples and a list of genomic feature, and overlaps each sample with each genomic feature.
# Then, it calculates the mean of the number overlaps and the standard deviation.
# This function can use countOverlaps and subsetByOverlaps
# This function returns a list with a data frame for each option(countOverlaps or subsetByOverlaps)
Compare_seq2 <- function(folder_samples=NULL,no_samples, generic_name, position_otherEnd=NULL, 
                         samples=NULL, folder, list_frag, sep="\t", header=TRUE, CountOL=FALSE, 
                         Subset_GtoF=FALSE, Subset_FtoG=FALSE,ncores=8, multicoreSampling=T) {
  if (is.null(samples)& is.null(folder_samples)) {
    stop("Missing file")
  }
  if (!is.null(folder_samples)) {
    samples <- Read_samples(folder=folder_samples,no_samples=no_samples, generic_name=generic_name, sep=sep, header=header)
  }
  
  list_names <- names(list_frag)
  if(multicoreSampling){
    cat("First mclappy... \n")
    aux<- mclapply(samples, function(i) {
      Compare_seq(sample=i, folder=folder, position_otherEnd=position_otherEnd, list_frag=list_frag, sep=sep, header=header, CountOL=CountOL, Subset_GtoF=Subset_GtoF, Subset_FtoG=Subset_FtoG)
    }, mc.cores=ncores)
  }  
  else{
    cat("Lapply... \n")
    aux<- lapply(samples, function(i) {
      Compare_seq(sample=i, folder=folder, position_otherEnd=position_otherEnd, list_frag=list_frag, sep=sep, header=header, CountOL=CountOL, Subset_GtoF=Subset_GtoF, Subset_FtoG=Subset_FtoG)
    })
  }
  flag=FALSE
  index<-c() 
  for (i in 1:length(samples)) {
    if (is.null(aux[[i]])) {
      cat("sample ",i, "did not work the first time...\n")
      index <- c(index,i)
      flag=TRUE
    }
  }
  
  if(flag) {
    # Assume this never happens to lapply, only to mclapply 
    cat("Second mclapply... \n")
    aux_n <- mclapply(samples[index], function(i) { 
      Compare_seq(sample=i, folder=folder, position_otherEnd=position_otherEnd, list_frag=list_frag, sep=sep, header=header, CountOL=CountOL, Subset_GtoF=Subset_GtoF, Subset_FtoG=Subset_FtoG)
    }, mc.cores=ncores)
    aux <- c(aux,aux_n)
  }
 
  aux1 <- sapply(aux,function(x) { x[[1]][,3] })
  aux2 <- sapply(aux,function(x) { x[[2]][,3] })
  aux3 <- sapply(aux,function(x) { x[[3]][,3] })
  
  result1<-NULL
  result2<-NULL
  result3<-NULL

  if (CountOL==TRUE) {
    # Check how many Genomic Features you are overlapping
    if (length(list_frag)==1) {
        Mean <- mean(aux1)
        SD <- sd(aux1)
    }
    else {
        Mean <- as.vector(tapply(aux1, row(aux1), mean))
        SD <- as.vector(tapply(aux1, row(aux1), sd))
    }
   
    EB_high <- Mean + 1.96 *SD
    EB_low <- Mean - 1.96 *SD
    result1 <- data.frame(list_names, Mean, SD, EB_high, EB_low)
  }
  if (Subset_GtoF==TRUE) {
    if (length(list_frag)==1) {
      Mean <- mean(aux2)
      SD <- sd(aux2)
    }
    else {
        Mean <- as.vector(tapply(aux2, row(aux2), mean))
        SD <- as.vector(tapply(aux2, row(aux2), sd))
    }
    EB_high <- Mean + 1.96 *SD
    EB_low <- Mean - 1.96 *SD
    result2 <- data.frame(list_names, Mean, SD, EB_high, EB_low)
  }
  if (Subset_FtoG==TRUE) {
    if (length(list_frag)==1) {
        Mean <- mean(aux2)
        SD <- sd(aux2)
    }
    else {
        Mean <- as.vector(tapply(aux3, row(aux3), mean))
        SD <- as.vector(tapply(aux3, row(aux3), sd))
    }
    EB_high <- Mean + 1.96 *SD
    EB_low <- Mean - 1.96 *SD
    result3 <- data.frame(list_names, Mean, SD, EB_high, EB_low)
  }
  result<-list(result1, result2, result3)
  return(result)
}



Show_OL <- function(signOL, sampleOL, CountOL=FALSE, Subset_GtoF=FALSE, Subset_FtoG=FALSE) {
    k <- CountOL + Subset_GtoF + Subset_FtoG
    par(mfrow=c(k,1))
    if (CountOL==TRUE) {
        data <- data.frame(signOL[[1]][,3],sampleOL[[1]])
        d <- as.matrix(data[,c(1,3)])
        toplot <- barplot(t(d), beside=TRUE, main = "Total number of Overlaps", col=c("lightyellow","lightblue"), 
            legend = c("Significant Reads", "Random Samples"), names.arg=data[,2], ylab = c("Number of Overlaps with Feature") )
        arrows(toplot[2,], data$Mean, toplot[2,], data$EB_high, length=0.1, angle=90)
        arrows(toplot[2,], data$Mean, toplot[2,], data$EB_low, length= 0.1, angle=90)
        }
    if (Subset_GtoF==TRUE) {
        data <- data.frame(signOL[[2]][,3],sampleOL[[2]])
        d <- as.matrix(data[,c(1,3)])
        toplot <- barplot(t(d), beside=TRUE, main = "Number of GF that map to an interaction in our samples", col=c("lightyellow","lightblue"), 
            legend = c("Significant Reads", "Random Samples"), names.arg=data[,2], ylab = c("Number of Overlaps with Feature") )
        arrows(toplot[2,], data$Mean, toplot[2,], data$EB_high, length=0.1, angle=90)
        arrows(toplot[2,], data$Mean, toplot[2,], data$EB_low, length= 0.1, angle=90)
        }
    if (Subset_FtoG==TRUE) {
        data <- data.frame(signOL[[3]][,3],sampleOL[[3]])
        d <- as.matrix(data[,c(1,3)])
        toplot <- barplot(t(d), beside=TRUE, main = "Number of interactions in our samples that map to a GF", col=c("lightyellow","lightblue"), 
            legend = c("Significant Reads", "Random Samples"), names.arg=data[,2], ylab = c("Number of Overlaps with Feature") )
        arrows(toplot[2,], data$Mean, toplot[2,], data$EB_high, length=0.1, angle=90)
        arrows(toplot[2,], data$Mean, toplot[2,], data$EB_low, length= 0.1, angle=90)
        }
    }




CompareSeqTotal <- function(x1=NULL, filename=NULL, score, colname_score, colname_dist=NULL, beyond_dist=NULL, before_dist=NULL,no_bins, sample_number, 
  restriction_enz_file=NULL, folder_samples=NULL, generic_name, folder=NULL, position_otherEnd=NULL, list_frag=NULL, sep="\t", header=TRUE, 
  CountOL=FALSE, Subset_GtoF=FALSE, Subset_FtoG=TRUE, plot_overlaps=FALSE, plot_name=NULL, distal=FALSE, coldist=NULL, unique_reads=TRUE, filterB2B=FALSE,
  b2bcol="isBait2bait", ncores = 8, multicoreSampling=T, negFraction = 1) {
    # Extract significant interactions
    # Be aware that you can trim for a specific window
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
		x1neg = as.data.table(x1neg)
      }
      x1 = rbind(x1pos, x1neg)
    }
    cat("Extract significant interactions...\n")
    result_1 <- Extract(x1=x1, filename=filename, score=score, colname_score=colname_score, colname_dist=colname_dist, beyond_dist=beyond_dist, before_dist=before_dist, significant=TRUE)
    # Extract non-significant interactions
    # You are trimming for the same window that was specified for significant interactions
    cat("Extract non-significant interactions...\n")
    result_2 <- Extract(x1=x1, filename=filename, score=score, colname_score=colname_score, colname_dist=colname_dist, beyond_dist=beyond_dist, before_dist=before_dist, significant=FALSE)
    # Bin non-significant interactions according to distance from bait before drawing random samples
    cat("Bin non-significant interactions according to distance from bait before drawing random samples...\n")
    result_2 <- Binning(sign=result_1, no_bins=no_bins, x1_nonsign=result_2, distal=distal)
    # Draw random samples
    cat("Draw random samples...\n")
    result_3 <- Draw_samples(x1_nonsign=result_2, sample_number=sample_number, restriction_enz_file=position_otherEnd, ncores=ncores, multicoreSampling=multicoreSampling)
    cat("Overlap significant interactions with Genomic Features...\n")
    list_names <- names(list_frag)
    result_4 <- Compare_seq(sample=result_1, folder=folder, position_otherEnd=position_otherEnd, list_frag=list_frag, sep=sep, 
                            header=header, CountOL=CountOL, Subset_GtoF=Subset_GtoF, Subset_FtoG=Subset_FtoG, unique_reads=unique_reads)
    cat("Overlap samples with Genomic Features...\n")
    result_5 <- Compare_seq2(samples=result_3, folder=folder, position_otherEnd=position_otherEnd, list_frag=list_frag, sep=sep, 
                             header=header, CountOL=CountOL, Subset_GtoF=Subset_GtoF, Subset_FtoG=Subset_FtoG, 
                             folder_samples=folder_samples,no_samples=sample_number, generic_name=generic_name, ncores=ncores, multicoreSampling=multicoreSampling)
    cat("Plot results...\n")
    if (plot_overlaps==TRUE) {
      if(is.null(plot_name)) {
        stop("Please provide name for plot.")
      } 
      pdf(paste0(plot_name), width=15, height=15)
    }
    Show_OL(signOL=result_4, sampleOL=result_5, CountOL=CountOL, Subset_GtoF=Subset_GtoF, Subset_FtoG=Subset_FtoG)
    if (plot_overlaps==TRUE) {
      dev.off()
        }  

}

#######################
### Example logwp=6 ###
# CompareSeqTotal(x1=z12_allG, score=6, sample_number=10, no_bins=100, colname_score="logwp",
#                 folder="/bi/home/paulafp/interactions/EncodeData/",
#                 position_otherEnd="/bi/group/sysgen/CHIC/Digest_Human_HindIII.bed", list_frag=files_GF_h1)
#######################


# Function for computing the overlaps (number and nature) present in significant interactions from a single bait 
OLSingleBait <- function(x, bait, GF, score=6, position_oe="/bi/group/sysgen/CHIC/Digest_Human_HindIII.bed") {
  
  position_oe <- read.table("/bi/group/sysgen//CHIC/Digest_Human_HindIII.bed")
  names(position_oe)<-c("chr","start","end","ID")
  
  Bait_int <- x[x$baitID==bait & x$logwp>=score,]

  res <- c()
    
  for (i in 1:nrow(Bait_int)) {
    res <- rbind(res,position_oe[position_oe$ID==Bait_int$otherEndID[i],])
  }
  
  Bait_int<- res
  

  Bait_int_GR <- GRanges( seqnames = Rle(paste0("chr",Bait_int[,1])),
                          IRanges( start = Bait_int[,2], end = Bait_int[,3]),
                          baitID=Bait_int[,4])
  
  overlap<-c()
  for (i in 1:length(Bait_int_GR)) {
    for (j in 1:length(GF)) {
      CountOL <-countOverlaps(Bait_int_GR[i],GF[[j]])
      newrow <- c(otherEndID = Bait_int_GR$baitID[[i]], CountOL=CountOL, GenomicFeature=names(GF)[[j]])
      overlap <- rbind(overlap,newrow)
    }
  }
  rownames(overlap)<- NULL
  overlap<- as.data.frame(overlap)
  
  return(overlap)
}


# Function copied from the internet to convert data frames into Granges objects
# https://stat.ethz.ch/pipermail/bioconductor/2011-November/042333.html

# keepColumns indicate whether additional data frame columns should be
# Sput into the GRanges.  ignoreStrand makes the strand of the
# constructed GRanges equal to *.  It assumes the input data.frame has
# columns chr/seqnames, start, end.

data.frame2GRanges <- function(df, keepColumns = FALSE, ignoreStrand = FALSE) {
  stopifnot(class(df) == "data.frame")
  stopifnot(all(c("start", "end") %in% names(df)))
  stopifnot(any(c("chr", "seqnames") %in% names(df)))
  if("seqnames" %in% names(df))
    names(df)[names(df) == "seqnames"] <- "chr"
  if(!ignoreStrand && "strand" %in% names(df)) {
    if(is.numeric(df$strand)) {
      strand <- ifelse(df$strand == 1, "+", "*")
      strand[df$strand == -1] <- "-"
      df$strand <- strand
    }
    gr <- GRanges(seqnames = df$chr,
                  ranges = IRanges(start = df$start, end = df$end),
                  strand = df$strand)
  } else {
    gr <- GRanges(seqnames = df$chr,
                  ranges = IRanges(start = df$start, end = df$end))
  }
  if(keepColumns) {
    dt <- as(df[, setdiff(names(df), c("chr", "start", "end", "strand"))],
             "DataFrame")
    elementMetadata(gr) <- dt
  }
  names(gr) <- rownames(df)
  gr
}

