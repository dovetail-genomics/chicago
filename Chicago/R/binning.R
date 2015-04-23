# This function bins results assigns probabilities to bins depending on their distance from bait
.binning <- function(sign, no_bins, x1_nonsign) {
  if(!is.data.table(sign)){setDT(sign)}
  # Bin distances from bait in sign - 100 bins
  if (is.null(sign$dist)) {
    sign[,dist := abs(distSign)]
  }
  
  sign[,distbin2 := cut(dist, breaks=(no_bins))]
  
  # Calculate how many other-ends in this bin
  bin_reads2 <- sign[,length(dist), by="distbin2"]
  setnames(bin_reads2,"distbin2","udbin2")
  setnames(bin_reads2,"V1","bin_reads")
  
  # Bin distances from bait in x1_nonsign
  if(!is.data.table(x1_nonsign)){setDT(x1_nonsign)}
  
  if (is.null(x1_nonsign$dist)) {
    x1_nonsign[,dist:=abs(distSign)]
  }
  
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
  sign[,distbin2:=NULL]
  return(x1_nonsign)
}
