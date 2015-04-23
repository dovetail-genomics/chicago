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
