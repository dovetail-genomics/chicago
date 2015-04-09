chicagoPipeline <-
function(cd, outprefix, printMemory=FALSE)
{
  message("\n*** Running normaliseBaits...\n")
  cd = normaliseBaits(cd)

  if(printMemory){
    print(gc(reset=T))
  }
  
  message("\n*** Running normaliseOtherEnds...\n")
  cd = normaliseOtherEnds(cd, outfile=paste0(outprefix, "_oeNorm.pdf"))
  
  if(printMemory){
    print(gc(reset=T))
  }
  
  message("\n*** Running estimateTechnicalNoise...\n")
  cd = estimateTechnicalNoise(cd, outfile=paste0(outprefix, "_techNoise.pdf"))
  
  if(printMemory){
    print(gc(reset=T))
  }
  
  message("\n*** Running estimateDistFun...\n")
  
  ### Note that f is now saved in cd@params
  cd = estimateDistFun(cd, outfile=paste0(outprefix, "_distFun.pdf"))
  
  if(printMemory){
    print(gc(reset=T))
  }  
  
  ### Note that f is now saved as cd@params$f and  
  ### subset is saved as cd@settings$brownianNoise.subset
  message("\n*** Running estimateBrownianNoise...\n")
  cd = estimateBrownianNoise(cd)

  if(printMemory){
    print(gc(reset=T))
  }  
  
  message("\n*** Running getPvals...\n")
  cd = getPvals(cd)
  
  if(printMemory){
    print(gc(reset=T))
  }  
  
  message("\n*** Running getScores...\n")
  cd = getScores(cd)
  
  if(printMemory){
    print(gc(reset=T))
  }  
  
  cd
}
