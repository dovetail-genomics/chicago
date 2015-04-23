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
