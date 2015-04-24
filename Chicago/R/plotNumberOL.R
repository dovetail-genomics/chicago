.plotNumberOL <- function(x_sign,s, files, plot_name=NULL) {
  x_sign[,dist:=NULL]
  x_sign<-colSums(x_sign[,(ncol(x_sign)-length(files)+1):ncol(x_sign),with=FALSE],na.rm = T)
  
  sample_number<- length(s)
  featureSumsMatrix <- matrix(rep(0),length(files)*sample_number,nrow=sample_number,ncol=length(files))
  for (i in 1:sample_number){
    x<-s[[i]]
    x[,dist:=NULL]
    x[,distbin3:=NULL]
    x[,bin_reads:=NULL]
    x[,i:=NULL]
    featureSums <- colSums(x[,(ncol(x)-length(files)+1):ncol(x),with=FALSE],na.rm = T)
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
#     cat(paste0("Plot saved under the name ",plot_name," in your working directory...\n"))
  }
  
  # Return Matrix with Number of overlaps for ou significant interactions dataset and our samples
  colnames(data)<-c("OLwithSI","MeanOLwithSamples", "SDOLwithSample","HigherCI", "LowerCI")
  cat("Return Table with results...\n\n")
  return(data[,c("OLwithSI","MeanOLwithSamples", "SDOLwithSample", "LowerCI","HigherCI")])
}
