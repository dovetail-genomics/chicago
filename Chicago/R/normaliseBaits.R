normaliseBaits <-
function(cd, normNcol="NNb", shrink=FALSE, plot=TRUE, outfile=NULL, debug=FALSE){
  message("Normalising baits...")

  adjBait2bait = cd@settings$adjBait2bait
  
  # NON-URGENT TODO: An even more memory-efficient way of doing this would be to have .normaliseFragmentSets assign
  # the s_j, distbin and refBinMean columns by reference!
  
  if(debug){
    ##returns sbbm and not Chicago object! 
    .normaliseFragmentSets(x=cd@x, s=cd@settings, npb=.readNPBfile(s=cd@settings), viewpoint="bait", idcol="baitID", Ncol="N", adjBait2bait=adjBait2bait, shrink=shrink, refExcludeSuffix=NULL, plot=plot, outfile=outfile, debug=T)
  }
  else{
    cd@x = .normaliseFragmentSets(x=cd@x, s=cd@settings, npb=.readNPBfile(s=cd@settings), viewpoint="bait", idcol="baitID", Ncol="N", adjBait2bait=adjBait2bait, shrink=shrink, refExcludeSuffix=NULL, plot=plot, outfile=outfile, debug=F)
    
  }
  
  # sort by baitID, otherEndID and move distbin column to the end of the table 
  cd@x[, (normNcol):= pmax(1, round(N/s_j)) ] # do not completely "cancel" interactions that have one read 

  setkey(cd@x, baitID, otherEndID)

  othercols = names(cd@x)[!names(cd@x)%in% c("distbin", normNcol)]

  setcolorder(cd@x, c(othercols, "distbin", normNcol))
  
  cd
  
}
