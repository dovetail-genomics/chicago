exportResults <-
function(cd, outfileprefix, scoreCol="score", cutoff, b2bcutoff=NULL, format=c("seqMonk","interBed","washU"), order=c("position", "score")[1]){
    
  if (any(c("rChr", "rStart", "rEnd", "rID", "bChr", "bStart", "bEnd", "bID") %in% colnames(x))){
    stop ("Colnames x shouldn't contain rChr, rStart, rEnd, rID, bChr, bStart, bEnd, bSign, bID\n") 
  }
  if (!all(format %in% c("seqMonk","interBed", "washU"))){
    stop ("Format must be either seqMonk, interBed or washU (or a vector containing several of these)\n")
  }
  if (! order %in% c("position","score")){
    stop ("Order must be either position (default) or score\n")
  }
  
  message("Reading the restriction map file...")
  rmap = fread(cd@settings$rmapfile)
  setnames(rmap, "V1", "rChr")
  setnames(rmap, "V2", "rStart")
  setnames(rmap, "V3", "rEnd")
  setnames(rmap, "V4", "otherEndID")
  
  message("Reading the bait map file...")
  baitmap = fread(cd@settings$baitmapfile)

  setnames(baitmap, "V1", "baitChr")
  setnames(baitmap, "V2", "baitStart")
  setnames(baitmap, "V3", "baitEnd")
  setnames(baitmap, cd@settings$baitmapFragIDcol, "baitID")
  setnames(baitmap, cd@settings$baitmapGeneIDcol, "promID")
  
  message("Preparing the output table...")

  if (is.null(b2bcutoff)){
    x = cd@x[ get(scoreCol)>=cutoff ]
  }
  else{
    x = cd@x[ (isBait2bait==T & get(scoreCol)>=b2bcutoff ) | 
                ( isBait2bait==F & get(scoreCol)>=cutoff )]
  }

  x = x[, c("baitID", "otherEndID", "N", scoreCol), with=F]
  
  setkey(x, otherEndID)
  setkey(rmap, otherEndID)
  
  x = merge(x, rmap, by="otherEndID", allow.cartesian = T)
  setkey(x, baitID)
  
  setkey(baitmap, baitID)  
  x = merge(x, baitmap, by="baitID", allow.cartesian = T)
        
  # note that baitmapGeneIDcol has been renamed into "promID" above 
  bm2 = baitmap[,c ("baitID", "promID"), with=F]

  setDF(x)
  setDF(bm2)
  
  # this way we can be sure that the new column will be called promID.y  
  out = merge(x, bm2, by.x="otherEndID", by.y="baitID", all.x=T, all.y=F, sort=F)
  out[is.na(out$promID.y), "promID.y"] = "."
  
  out = out[,c("baitChr", "baitStart", "baitEnd", "promID.x", "rChr", "rStart", "rEnd", "otherEndID", scoreCol, "N", "promID.y")]
  
  names(out) = c("bait_chr", "bait_start", "bait_end", "bait_name", "otherEnd_chr", "otherEnd_start", "otherEnd_end", "otherEnd_ID", "score", "N_reads", "otherEnd_name")
  
  out$N_reads [ is.na(out$N_reads) ] = 0
  out$score = round(out$score,2)

  if (order=="position"){
    out = out[order(out$bait_chr, out$bait_start, out$bait_end, out$otherEnd_chr, out$otherEnd_start, out$otherEnd_end), ]
  }
  if (order=="score"){
    out = out[order(out$score, decreasing=T), ]
  }
  out0=out
   
  if ("seqMonk" %in% format){
    message("Writing out for seqMonk...")
    out[,"bait_name"] = gsub(",", "|", out[,"bait_name"], fixed=T)
    
    out$newLineOEChr = paste("\n",out[,"otherEnd_chr"], sep="")    
    out = out[,c("bait_chr", "bait_start", "bait_end", "bait_name", "N_reads", "score", "newLineOEChr", "otherEnd_start", "otherEnd_end", "otherEnd_name", "N_reads", "score")]
  
    write.table(out, paste0(outfileprefix,"_seqmonk.txt"), sep="\t", quote=F, row.names=F, col.names=F)
  }	
  if ("interBed" %in% format){
    message("Writing out interBed...")
    out = out0[,c("bait_chr", "bait_start", "bait_end", "bait_name", 
                 "otherEnd_chr", "otherEnd_start", "otherEnd_end", "otherEnd_name", 
                 "N_reads", "score")]
    write.table(out, paste0(outfileprefix,".ibed"), sep="\t", quote=F, row.names=F)	
  }
  if("washU" %in% format){
   message("Writing out for washU browser...")
   out = out0[,c("bait_chr", "bait_start", "bait_end", 
                 "otherEnd_chr", "otherEnd_start", "otherEnd_end", "otherEnd_name",
                 "score")]
   
   ##Bait to bait interactions can be asymmetric in terms of score. Here, we find asymmetric interactions and delete the minimum score
   setDT(x)
   setkey(x, baitID, otherEndID)
   x$ReversedInteractionScore <- x[J(x$otherEndID, x$baitID), get(scoreCol)]
   sel <- x[,get(scoreCol)] > x$ReversedInteractionScore
   sel <- ifelse(is.na(sel), TRUE, sel) ##"FALSE" entries in sel should correspond to minima
   setDF(x)
   x$ReversedInteractionScore <- NULL
   out <- out[sel,]
   
   out$i = seq(1,nrow(out)*2,2)
   res = apply(out,1,function(x){
      lines = paste0(x["bait_chr"], "\t", x["bait_start"], "\t", x["bait_end"], "\t", x["otherEnd_chr"],":", x["otherEnd_start"], "-", x["otherEnd_end"], ",", x["score"],"\t", x["i"], "\t", ".") 
      lines = paste0(lines,  "\n",
        paste0(x["otherEnd_chr"], "\t", x["otherEnd_start"], "\t", x["otherEnd_end"], "\t", x["bait_chr"],":", x["bait_start"], "-", x["bait_end"], ",", x["score"],"\t", as.numeric(x["i"])+1, "\t", "."))
         lines
    })
    res = gsub(" ", "", res)
    writeLines(res, con=paste0(outfileprefix,"_washU.txt"))
   }
  
}
