.convertBedFormat2GR <- function(folder=NULL, list_frag=NULL, sep="\t", header=FALSE, rm.MT = FALSE) {
  if (is.null(list_frag) ) {
    stop("Please provide list with files of Genomic features to overlap")
  }
  list_names <- names(list_frag) 
  if (!is.null(folder) ) { list_frag = file.path(folder, list_frag)} # note list_frag names are lost

  j <- 1
  result <- list()
  for (i in list_frag) {
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
