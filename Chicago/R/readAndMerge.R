readAndMerge <-
function(files, cd, ...){
  mergeSamples(lapply(files, readSample, cd), ...)
}
