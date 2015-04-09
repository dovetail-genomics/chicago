plotBaits <-
function(cd, pcol="score", Ncol="N", n=16, baits=NULL, plotBaitNames=TRUE, plotBprof=FALSE,plevel1 = 5, plevel2 = 3, outfile=NULL, removeBait2bait=TRUE, width=20, height=20, maxD=NULL, bgCol="black", lev2Col="blue", lev1Col="red", bgPch=1, lev1Pch=20, lev2Pch=20, ...)
{
  if(plotBaitNames){
    baitmap = fread(cd@settings$baitmapfile)
  }
  if (is.null(baits)){
    baits = sample(unique(cd@x$baitID),n)
  }
  else{
    n = length(baits)
  }
 
  if(plotBprof){
    disp = cd@params$dispersion 
  }
 
  if (!is.null(outfile)){ 
    pdf(outfile, width=width, height=height)
  }
  if(n>=4){
    par(mfrow=c(4, ceiling(n/4)))
  }
  else{
    par(mfrow=c(n, 1))
  }

  setkey(x@x, baitID)

  for(i in 1:n){

    this = cd@x[baitID==baits[i]]
    
    
    
    this = this[is.na(distSign)==FALSE]

    if (!is.null(maxD)){
       this = this[abs(distSign)<=maxD]
    }
     
    if (removeBait2bait){
       this = this[isBait2bait==FALSE]
    }

    setDF(this)
    this = this[order(this$distSign),]

    cols <- rep(bgCol, nrow(this))
    pchs <- rep(bgPch, nrow(this))
    sel1 <- this[,pcol] >=plevel1
    sel2 <- this[,pcol] >=plevel2
    cols[sel2] <- lev2Col ##less stringent first
    cols[sel1] <- lev1Col
    pchs[sel2] <- lev2Pch
    pchs[sel1] <- lev2Pch
    
    title = paste(baits[i], sep="")
    if(plotBaitNames){
         baitName = baitmap[baitmap$V4==baits[i]][, cd@settings$baitmapGeneIDcol, with=F]
         if (length(grep(",",baitName))){
             baitName = gsub("(\\S+,).+","\\1", baitName)
             baitName = paste0(baitName, "...")
         }
         title = paste0(baitName, " (", title, ")")
    }    

    plot(this$distSign, this[,Ncol], xlab=distcol, ylab=Ncol, main=title, col=cols, pch=pchs, ...)
    abline(v=0, col="grey", lwd=1)

    if(plotBprof){
	      lines(this[,distcol], this$Bmean, lwd=1, col="darkgrey")
        lines(this[,distcol], this$Bmean+1.96*sqrt(this$Bmean+this$Bmean^2/disp), 
              lwd=1, lty=2, col="darkgrey")
    }
  }
  if (!is.null(outfile)){ 
    dev.off()
  }
  baits
}
