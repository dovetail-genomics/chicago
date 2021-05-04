library(data.table)
library(argparser)
options("scipen"=100, "digits"=4)

args = commandArgs(trailingOnly=T)

p = arg_parser("Bin restriction fragments from rmap", name="Rscript makeBins.R")
p = add_argument(p, arg="<rmap>", 
                 help="Full path to rmap file (bed file containing coordinates of the restriction fragments). 
                       Rmap should be coordinate sorted.")
p = add_argument(p, arg="--baitmap", 
                 help="Full path to baitmap file (bed file containing coordinates of baited restriction fragments). 
                       Baitmap should be coordinate sorted. If not provided, it will look for .baitmap file matching rmap file header.",
                 default=NA)
p = add_argument(p, arg="--output_prefix", 
                 help="Prefix for resulting new rmap with binning incorporated. 
                       If none provided, it will default to rmap header.", default=NA)
p = add_argument(p, arg="--binsize", 
                 help = "Size of bins in bp. Default: 5000.", default = 5000)
p = add_argument(p, arg="--start", 
                 help = "Column name of start coordinate from rmap. By default, it assumes that rmap does not contain a header, thus column will be V2.", default = "V2")
p = add_argument(p, arg="--end",
                 help = "Column name of end coordinate from rmap. By default, it assumes that rmap does not contain a header, thus column will be V3.", default = "V3")
p = add_argument(p, arg="--include_baits",
                 help = "Flag specifying whether to include baited fragments in the bins. By default, the baited fragments will be left 'solitary' as separate bins.", flag=TRUE)
p = add_argument(p, arg="--verbose",
                 help = "Flag specifying whether to print process steps.", flag=TRUE)


opts = parse_args(p, args)

rmap = opts[["<rmap>"]]
baitmap = opts[["baitmap"]]
output_prefix = opts[["output_prefix"]]
bin_size = opts[["binsize"]]
chr_start = opts[["start"]]
chr_end = opts[["end"]]
include_baits = opts[["include_baits"]]
verb = opts[["verbose"]]

if(is.na(baitmap)){
  baitmap <- gsub(".rmap",".baitmap",rmap)
  message("Using baitmap file: ", baitmap) 
}
baitmap_dt = fread(baitmap)
baitmap_dt[,V1:=as.character(V1)]

if (is.na(output_prefix)){
  output_prefix = gsub(".rmap","",rmap)
}


rmap_dt = fread(rmap)
setnames(rmap_dt, names(rmap_dt)[1], "chr")

makeBins = function(test, binsize=5000, st, end, verbose){

  bins = vector("numeric", nrow(test))
  previ=1
  prevchr = test$chr[1]
  i=1
  k=1
  while (i <= nrow(test)){
    while(test[[end]][i]-test[[st]][previ]<binsize & test$chr[i]==prevchr & i<nrow(test)){
      prevchr=test$chr[i]
      i=i+1
    }
    # We have to treat reaching chr end and exceeding bin size differently,
    # because for chr end we draw the bin boundary before we overshoot,
    # whereas for bin size we draw the bin boundary *after* we overshoot so we don't end up with very small bins preceding large frags
    if( test$chr[i]!=prevchr){
      if (verbose){
        message("previ=", previ, " i=",i, " k=", k, " test[[end]][i-1]-test[[st]][previ]=", test[[end]][i-1]-test[[st]][previ], " test$chr[i]=", test$chr[i], " prevchr=", prevchr)
      }else{
        message("*") # this happens at the end of the old chr
      }
      bins[previ:(i-1)]=k
    } else {
      if (verbose){
        message("previ=", previ, " i=",i, " k=", k, " test[[end]][i]-test[[st]][previ]=", test[[end]][i]-test[[st]][previ], " test$chr[i]=", test$chr[i], " prevchr=", prevchr)
      }else{
        if(i>1) if(prevchr!=test$chr[i-1]) message(".") # this happens at the start of the new chr
      }
      bins[previ:i] = k
      i=i+1
    }
    k=k+1
    previ=i
    prevchr=test$chr[i]
  }
  bins

}

makeBins2 = function(test, bm, binsize=5000, st="V2", end="V3", verbose=F){
  test[, isBaited:=V4 %in% bm$V4]
  setnames(test, names(test)[1], "chr")
  bins = vector("numeric", nrow(test))
  previ=1
  prevchr = test$chr[1]
  i=1
  k=1
  while (i <= nrow(test)){
    while(test[[end]][i]-test[[st]][previ]<binsize & test$chr[i]==prevchr & i<nrow(test) & !test$isBaited[i]){
      prevchr=test$chr[i]
      i=i+1
    }
    # We have to treat reaching chr end and exceeding bin size differently,
    # because for chr end we draw the bin boundary before we overshoot,
    # whereas for bin size we draw the bin boundary *after* we overshoot so we don't end up with very small bins preceding large frags
    if( test$chr[i]!=prevchr | test$isBaited[i]){
      if (verbose){
        message("previ=", previ, " i=",i, " k=", k, " test[[end]][i-1]-test[[st]][previ]=", test[[end]][i-1]-test[[st]][previ], " test$chr[i]=", test$chr[i], " prevchr=", prevchr)   
      }
      bins[previ:(i-1)]=k
      if(test$isBaited[i]){
        if (verbose){
          message("previ=", previ, " i=",i, " k=", k, " test[[end]][i]-test[[st]][previ]=", test[[end]][i]-test[[st]][previ], " test$chr[i]=", test$chr[i], " prevchr=", prevchr, " isBaited=", test$isBaited[i])   
        }
        bins[i]=k+1
        k=k+1
        i=i+1
      }
    }
    else{
      if (verbose){
        message("previ=", previ, " i=",i, " k=", k, " test[[end]][i]-test[[st]][previ]=", test[[end]][i]-test[[st]][previ], " test$chr[i]=", test$chr[i], " prevchr=", prevchr)   
      }else{
        if(i>1) if(prevchr!=test$chr[i-1]) message(".") # this happens at the start of the new chr
      }
      bins[previ:i] = k
      i=i+1
    }
    k=k+1
    previ=i
    prevchr=test$chr[i]
  }
  bins
}

binRmap <- function(test,bins,binsize,agg_baits=TRUE){
  if (agg_baits == TRUE){
    tag = paste0("_",binsize/1000,"kb.rmap")
  } else {
    tag = paste0("_",binsize/1000,"kb_sol_baits.rmap")
  }
  #if (is.na(output_prefix)){
  #  output_prefix = gsub(".rmap","",rmap)
  #}
  brmp_filename = paste0(output_prefix,tag)
  
  test[,bin:=allBins]
  
  binned_rmap = test[, .(first(chr), first(V2), last(V3)), by="bin"]
  setcolorder(binned_rmap, c(2:4,1))
  
  # Sanity Check
  if(any(binned_rmap$V2>binned_rmap$V3)){
    message("ERROR: Chromosome start should not be greater than chromosome end! Exiting...")
    quit(save = "no", status = 1, runLast = FALSE)
  } else {
    write.table(binned_rmap,brmp_filename,quote = FALSE,sep = "\t",row.names = FALSE,col.names = FALSE)
    msg = paste0("Rmap has been binned with ", bin_size/1000,
                 "kb bins. Results saved in ", brmp_filename, " .")
    message(msg)
  }
  return(binned_rmap)
}

editBaitmap <- function(bm,brmp,binsize,agg_baits=TRUE){
  tag = paste0("_",binsize/1000,"kb")
  setkey(bm, V1,V2,V3)
  setkey(brmp, V1, V2, V3)
  bm2 = foverlaps(binned_rmap, bm, by.x=c("V1", "V2", "V3"), by.y=c("V1", "V2", "V3"), nomatch = 0, mult="all")
  if (agg_baits==TRUE) {
    # aggregating over multiple baits such that each bin is present only once
    tag = paste0(tag,".baitmap")
  } else {
    # “leave baited fragments unbinned” approach
    tag = paste0(tag,"_sol_baits.baitmap")
  }
  bm2 <- bm2[, .(V1[1],i.V2[1],i.V3[1], paste(V5, collapse=",")), by="bin"]
  setcolorder(bm2, c(2:4,1,5))
  bm2_filename = paste0(output_prefix,tag)
  write.table(bm2, bm2_filename, col.names=F, sep="\t", quote=F, row.names=F)
  return(bm2)
}


#######################
message(paste0("Creating ",bin_size,"kb bins..."))
if (include_baits == TRUE){
  message("All restrictions fragments will be binned, including baits...")
  allBins <- makeBins(test=rmap_dt, binsize=bin_size, st=chr_start, end=chr_end, verbose=verb)
} else {
  message("All restrictions fragments except baits will be binned...")
  allBins <- makeBins2(test=rmap_dt, bm=baitmap_dt, binsize=bin_size, st=chr_start, end=chr_end, verbose=verb)
}
message("Binning rmap...")
binned_rmap  <- binRmap(test = rmap_dt, bins = allBins, binsize = bin_size, agg_baits=include_baits)
message("Editing baitmap...")
save.image("testing.RData")
new_baitmap <- editBaitmap(bm = baitmap_dt,brmp=binned_rmap,binsize=bin_size, agg_baits=include_baits)
save.image("testing.RData")
message("Done!")


