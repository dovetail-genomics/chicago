#!/usr/bin/env Rscript

suppressMessages(library(argparser))
args = commandArgs(trailingOnly=T)

p <- arg_parser("Plotting data sparsity within the CHiCAGO background model", name="Rscript plotBackgroundDropouts.R", hide.opts = TRUE)
p <- add_argument(p, arg="<Rds>",
                  help="Full path to Rds file(s), comma separated e.g. file1.rds,file2.rds...")
p = add_argument(p, arg="--output_prefix",
                 help="Prefix of resulting plot", default="CHiCAGOplot")
p = add_argument(p, arg="--outdir",
                 help="Full path to output directory", default=".")
p = add_argument(p, arg="--verbose",
                 help = "Flag specifying whether to print process steps.", flag=TRUE)

opts = parse_args(p, args)

my_rds = opts[["<Rds>"]]
boxplotname = opts[["output_prefix"]]
outdir = opts[["outdir"]]
verb = opts[["verbose"]]

suppressMessages(library(Chicago))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(data.table))


getSparsityDistBin <- function(CHiCAGO.Rds) {
  # Get files
  chicdata <- readRDS(CHiCAGO.Rds)
  samplename <- basename(CHiCAGO.Rds)
  int <- chicdata@x
  poe <- fread(chicdata@settings$proxOEfile, skip = 1)
  names(poe) = c("baitID", "otherEndID", "poe_dist")
  # Restrict Rds to maxLBrownEst
  LBD <- as.numeric(chicdata@settings$maxLBrownEst)
  binsize <- as.numeric(chicdata@settings$binsize)
  int[, absDist := abs(distSign)]
  int_short <- int[absDist <= LBD]
  # Remove bait2bait as these are not used in background estimation
  int_short_nob2b <- int_short[isBait2bait == FALSE]
  # Join with the poe file - which is already restricted to maxLBrownEst
  setkey(poe, baitID, otherEndID)
  int_poe <- int_short_nob2b[poe, on = c("baitID", "otherEndID")] 
  # int_poe has all interacting fragments regardless of whether they have reads
  # Remove baits where every other end is zero as they do not contribute to the model
  int_poe$N[is.na(int_poe$N)] <- 0
  temp <- as.data.table(int_poe %>% group_by(baitID) %>% summarize(total_count = sum(N)))
  baits_to_remove1 <- temp[total_count == 0]
  baits_to_remove <- as.data.table(baits_to_remove1$baitID)
  names(baits_to_remove)[1] = "baitID"
  new_int_poe <- int_poe[!baits_to_remove, on = "baitID"]
  # bin poe_dist with breaks corresponding to the binsize parameter
  my_breaks <- seq(from = 0, to = LBD, by = binsize)
  int_poe_bins <- as.data.table(new_int_poe %>% group_by (Distance = cut(poe_dist, breaks = my_breaks, include.lowest = TRUE)))
  # for each bait within each distance bin, find the proportion of dropouts (Other ends with count N of 0)
  total <- as.data.table(int_poe_bins %>% group_by(baitID, Distance) %>% tally())
  names(total)[3] = "n_totalOE"
  exist <- int_poe_bins[N != 0]
  total_exist <- as.data.table(exist %>% group_by(baitID, Distance) %>% tally())
  names(total_exist)[3] = "n_existOE"
  miss <- int_poe_bins[N == 0]
  total_miss <- as.data.table(miss %>% group_by(baitID, Distance) %>% tally())
  names(total_miss)[3] = "n_missOE"
  # Join total, exist and miss together
  all_bins1 <- total_exist[total, on = c("baitID", "Distance")]
  all_bins <- total_miss[all_bins1, on = c("baitID", "Distance")]
  all_bins[is.na(all_bins)] <- 0
  # Get proportion missing
  all_bins[, prop_miss := n_missOE/n_totalOE]
  all_bins$Sample <- samplename
  return(all_bins) # Table containing dropout rate per bait, per distance bin
}


########### Get sparsity per distance bin (binsize) within the Brownian component (maxLBrownEst)

message("Reading input files and calculating data sparsity in distance bins...")
files <- as.list(scan(text=my_rds, what="", sep = ","))
bintable <- data.table()
for(file in files) {
  result <- getSparsityDistBin(file)
  bintable <- rbind(bintable, result, fill = TRUE)
}
pdf(file = paste(outdir, "/", boxplotname, "_backgroundSparsity.pdf", sep = ""), width =10, height = 8)
message("Plotting graph...")
f <- ggplot(bintable, (aes(x = Distance, y = prop_miss)))
f + geom_boxplot(aes(fill = Sample), outlier.shape = NA) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~Sample, scales = "free") +
  scale_x_discrete(breaks = levels(bintable$Distance)[c(T, rep(F, 9))]) + # for simple viewing shows on x axis 1 bin in every 10
  ylab("Proportion of missing other ends") +
  xlab("Distance bin")
invisible(dev.off())
############## 

message("...done!")
