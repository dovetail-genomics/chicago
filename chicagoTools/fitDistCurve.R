#! /usr/bin/env Rscript
# Determine the weighting function based on call overlap across replicates

# Functions used to fit a bounded logistic regression model.

# Essentially, p, the probability of an interaction between two fragments, is a function of d_i, the distance between those fragments:

# $p = (expit \delta - expit \gamma)\times expit(\alpha + \beta\log{d_i}) + expit(\gamma)$

##Inputs

library(argparser)

if (packageVersion("argparser") < 0.3) {
  stop("argparser version (", packageVersion("argparser"), ") is out of date - 0.3 or later is required. Please open R and run install.packages('argparser') to update.")
}

args = commandArgs(trailingOnly=TRUE)

p = arg_parser("Get the parameters for the p-value weighting curve (alpha through delta). Specify '--inputs' OR '--summaryInput'.", name="Rscript fitDistCurve.R")

p = add_argument(p, arg="<output-prefix>", help = "Output files: settings file for use in CHiCAGO, summary object, and plot.", default = "fitDistCurve.settingsFile")
p = add_argument(p, arg="--inputs", help = "Comma-separated list: locations of saved CHiCAGO objects (either .Rda or .Rds).")
p = add_argument(p, arg="--summaryInput", help = "An .Rda file of summary information -- the max P-val for each putative interaction, and the location of the .rmap file. This file will be generated if it wasn't provided.")
p = add_argument(p, arg="--threshold", help = "Threshold applied to -log(p) values (NB: not the CHiCAGO score!).", type="numeric", default = -10L)
p = add_argument(p, arg="--subsets", help = "Number of subsets to partition the data into. Parameters estimated on subsets, median taken.", type="numeric", default = 5L)
p = add_argument(p, arg="--largeBinSize", help = "Largest bin size to consider", type="numeric", default = 1E6L)
p = add_argument(p, arg="--binNumber", help = "Number of large bins", type="numeric", default = 16L)
p = add_argument(p, arg="--halfNumber", help = "First bin is subdivided into halves - the number of times to do this", type = "numeric", default = 5L)


opts = parse_args(p, args)

print(opts)

if(packageVersion("argparser") < 0.4)
{
  names(opts) <- gsub("-", "_", names(opts))
}

outputPrefix = opts[["<output_prefix>"]]
inputs <- opts[["inputs"]]
summaryInput <- opts[["summaryInput"]]
threshold <- opts[["threshold"]]
subsets <- opts[["subsets"]]

largeBinSize = opts[["largeBinSize"]]
binNumber = opts[["binNumber"]]
halfNumber = opts[["halfNumber"]]


##number of subsets to estimate parameters from, individually
Nsub <- subsets

##bins
breaks <- c(0, (2^((-halfNumber):0))*largeBinSize)
if(binNumber > 1)
{
  breaks <- c(breaks, (2L:binNumber)*largeBinSize)
}

bins <- data.frame(
  start=breaks + 1,
  end=c(breaks[-1], Inf)
)
bins$mid <- with(bins, 0.5*(end + start - 1))


##--------------------------------------------------

library(Chicago)

expit <- function(x)
{
  exp(x)/(1+exp(x))
}

##parse input files
validInput <- xor(is.na(inputs), is.na(summaryInput))
if(!validInput) {stop("Please specify --inputs or --summaryInput (but not both)")}

##check for output files already existing?
##FIXME


##############-------------------FIXME FROM HERE ON IN---------------################

message("Bins used:")
print(bins)

##location of restriction fragment map (rmapfile)
#rmapfile <- "/bi/group/sysgen/CHIC/Digest_Human_HindIII.bed" FIXME can get from CHiCAGO objects


##---------------------------------------------------------------
##Mode where inputs are specified

if(!is.na(inputs))
{
  samples <- strsplit(inputs, split = ",")[[1]]
  
  invalidFileNames <- !seq_along(samples) %in% grep(".[Rr][Dd][AaSs]$", samples)
  if(any(invalidFileNames))
  {
    stop("Some of the submitted files do not end with .Rda or .Rds.")
  }
  print(samples)
  if(any(!file.exists(samples)))
  {
    stop("Files not found: ", paste(samples[!file.exists(samples)],collapse=", "))
  }
  
    ##loop through the results, grabbing data
    for(s in samples)
    {
      ##Load sample
      gc()
      
      if(identical(grep(".[Rr][Dd][Ss]$", s), 1L))
      {
        message("Found a .Rds file: ", s)
        ##it's an Rds/rds, load it straight in
        x <- readRDS(s)
      } else {
        message("Found a .Rda file:", s)
        ##it's an Rda/rda, find the chicagoData object
        attach(s)
        myObjects <- ls(paste0("file:", s))
        detach(paste0("file:", s), character.only=TRUE)
        
        load(s)
        ##which objects are chicagoData objects?
        sel <- sapply(myObjects, function(x)
          {
          eval(parse(text=paste0("'chicagoData' %in% class(",x,")")))
        })
        if(length(sel) == 0)
        {
          warning("File ", s, " contained no objects at all. Skipping this file.")
          rm(list=myObjects)
          next
        }
        if(sum(sel) > 1) {stop("File ", s, " contained multiple chicagoData objects; currently not supported by this script. Aborting.")}
        if(sum(sel) == 0) {
          warning("File ", s, " did not contain a chicagoData object. Skipping this file.")
          rm(list=myObjects)
          next
        }
        
        chicagoObjectName <- names(sel)[sel]
        message("Found unique chicagoData object: ", chicagoObjectName)
        
        ##make sure the chicagoData object is copied to "x"
        eval(parse(text=paste0("x <- ", chicagoObjectName)))
        
        ##remove all objects, but not "x"
        if("x" %in% myObjects){
          sel <- "x" == myObjects
          myObjects <- myObjects[!sel]
        }
        rm(list=myObjects)
      }
      
      ##log which map file is being used
      if(!exists("rmapfile"))
      {
        rmapfile <- x@settings$rmapfile
      } else {
        if (rmapfile != x@settings$rmapfile)
        {
          warning("rmapfile location did not match across input files. The rmapfile specified first takes precedence; please check that your data were analysed using the same rmapfile.")
        }
      }

      x <- x@x
      requiredCols <- c("baitID","otherEndID","distSign","log.p")
      x <- x[,requiredCols, with=FALSE]
      
      #put x in merge-conducive form
      setnames(x, "log.p", "new.log.p")
      setkey(x, baitID, otherEndID, distSign)
      
      gc()
      
      message("Updating maxPvals...")
      if(exists("maxPvals"))
      {
        ##clever merging - 
        
        ##do merge, retaining only the highest p-val
        maxPvals <- merge(x, maxPvals, all=TRUE)
        
        maxPvals$log.p[is.na(maxPvals$log.p)] <- 0 ##since NA means p-val=1
        maxPvals$new.log.p[is.na(maxPvals$new.log.p)] <- 0 ##since NA means p-val=1
        
        maxPvals$log.p <- pmax(maxPvals$log.p, maxPvals$new.log.p)
        maxPvals$new.log.p <- NULL
        
        ##we CANNOT discard 0s
        
        gc()
        
      } else {
        maxPvals <- data.table(x)
        setnames(maxPvals, "new.log.p", "log.p")
        setkey(maxPvals, baitID, otherEndID, distSign)
        gc()
      }
      rm(x)
    }
    ##now we can discard 1s from maxPvals
    maxPvals <- maxPvals[maxPvals$log.p < 0,]


    rdaFile <- paste0(outputPrefix, "_summaryInput.Rda")
    message("Saving file: ", rdaFile)
    save(maxPvals, rmapfile, file=rdaFile)
  
  gc()
}

##---------------------------------------------------------------
##Mode where maxPvals is specified

if(!is.na(summaryInput))
{
  if(!identical(grep(".[Rr][Dd][Aa]$", summaryInput), 1L))
  {
    stop("summaryInput file is not an .Rda")
  }
  if(!file.exists(summaryInput))
  {
    stop("Could not find summaryInput file: ", summaryInput)
  }
  message("Reading summaryInput file: ", summaryInput)
  load(summaryInput) ##should contain maxPvals and rmapfile (can change code to check for this?)
}

##---------------------------------------------------------------
##All modes

#Data input: summaryInput object, data.table with cols: baitID, otherendID, distSign, log.p

#(log.p can be omitted if threshold <- NA)

#Replace following code chunk to supply your own maxPvals object:

#To calculate the number of pairs expected at each distance, we need to be a bit clever. (derivation: p215-217, JMC lab book)

##Obtain some information about the fragment map

rmap = fread(rmapfile)
setnames(rmap, "V1", "chr")
setnames(rmap, "V3", "end")
rmap <- rmap[! chr %in% c("MT", "chrM"),]
chrMax <- rmap[,max(end),by="chr"] ##length of each chr

##Chromsome lengths
C = as.numeric(sort(chrMax$V1))
##Genome size
G = sum(C)
##Number of baits
Nbaits <- length(unique(maxPvals$baitID))
##average fragment length
Lfrag <- Chicago:::.getAvgFragLength(rmapfile=rmapfile)

##calculate number of hypotheses at given distances
bins$N <- as.numeric(NA)
for(i in 1:(nrow(bins)-1))
{
  bins$N[i] <- with(bins,
              (2*(end[i] - start[i])*Nbaits/(G*Lfrag))*sum(pmax(0, C - 0.5*end[i] - 0.5*start[i]))
               )
}
bins$N[nrow(bins)] <- with(bins,
              (Nbaits/(G*Lfrag))*sum((pmax(0, C - start[nrow(bins)]))^2)
               )
#number of trans pairs is:
Ntrans <- (G^2 - sum(C^2))*Nbaits/(G*Lfrag)


#bins$N.old <- with(bins, 2*Nbaits*(end - start)/4000) ##can be improved

#label each interactions by its distance bin.

maxPvals$bin <- cut(abs(maxPvals$distSign), c(breaks, Inf))

##Split data into subsets

#Idea is to split data into $N_{sub}$ pieces, and somehow combine these intelligently to avoid outliers.

#We will get a matrix "obsJack". Each column contains one subset.

maxPvals$BIDbin <- cut(maxPvals$baitID, Nsub)
table(maxPvals$BIDbin)

outJack <- vector("list", Nsub)
obsJack <- matrix(nrow = nrow(bins)+1,ncol=Nsub)

for(i in 1:Nsub)
{
  lev = levels(maxPvals$BIDbin)[i]

  ##get subset of matrix
  temp <- maxPvals[maxPvals$BIDbin == lev,]
  
  ##calculate "obsmax"
  if(is.na(threshold))
    {
      obsJack[,i] <- as.matrix(with(temp, table(bin, useNA="always")))
    } else {
      obsJack[,i] <- as.matrix(with(temp, table(bin[log.p < threshold], useNA="always")))
    }
}

print(obsJack)


##fit model

#1. Construct likelihood function to maximise:

f <- function(params)
{
  if(length(params) != 4) stop("f got passed an invalid parameter vector.")
  alpha <- params[1]
  beta <- params[2]
  gamma <- params[3]
  delta <- params[4]
  
  d <- log(data.lr$mid)
  x <- data.lr$Ints
  N <- data.lr$N
  
  if(is.nan(beta*Inf))
    {
      p <- (expit(delta)-expit(gamma))*expit(alpha) + expit(gamma)
    } else {
      p <- (expit(delta)-expit(gamma))*expit(alpha + beta*d) + expit(gamma)
    }
  
  summands <- x*log(p) + (N-x)*log(1-p)
  #print(p)
  -sum(summands)
}

# 2. Construct data, optim

for(i in 1:Nsub)
  {
    data.lr <- cbind(rbind(bins, c(Inf,Inf,Inf,Ntrans)), Ints=obsJack[,i])
    data.lr$N <- data.lr$N %/% Nsub ##CRUCIAL!
    data.lr$notInts <- with(data.lr, N - Ints)
    #data.lr$mid <- with(data.lr, 0.5*(end + start - 1))
    
    start <- c(4*14.5, -4, -14.6, -8)
    #start <- c(0,-0.0000000001,-10, -1)
    out <- optim(start, f)
    out2 <- optim(out$par, f)
    outJack[[i]] <- out2
  }

outJackPar <- do.call(cbind,lapply(outJack, function(x){x$par}))
rownames(outJackPar) <- c("alpha","beta","gamma","delta")
outJackPar <- rbind(outJackPar, d0=-outJackPar["alpha",]/outJackPar["beta",])
outJackPar <- rbind(outJackPar, expd0=exp(outJackPar["d0",]))
outJackPar <- rbind(outJackPar, logPiRel=outJackPar["delta",] - outJackPar["gamma",])
outJackPar <- rbind(outJackPar, piRel=exp(outJackPar["logPiRel",]))

medianJackPar <- apply(outJackPar, 1, median)

# Plotting

p.fit <- function(d, params)
{
  alpha <- params[1]
  beta <- params[2]
  gamma <- params[3]
  delta <- params[4]
  
  if(is.nan(beta*Inf))
    {
      p <- (expit(delta)-expit(gamma))*expit(alpha) + expit(gamma)
    } else {
      p <- (expit(delta)-expit(gamma))*expit(alpha + beta*d) + expit(gamma)
    }
  p
}

##lines from fit to each subset:
pdf(paste0(outputPrefix, "_mediancurveFit.pdf"))
plot(log(data.lr$mid), log(p.fit(log(data.lr$mid), outJack[[1]]$par)), type="l", ylim = c(-19,-4), col=rainbow(5)[1], main="Data subsetting", xlab="log_e(distance)", ylab="log_e(probability of interaction)")
for(i in 2:Nsub)
{
  lines(log(data.lr$mid), log(p.fit(log(data.lr$mid), outJack[[i]]$par)), col=rainbow(5)[i])
}

##median fit
lines(log(data.lr$mid), log(p.fit(log(data.lr$mid), medianJackPar[1:4])), col="black", lty=2, pch=7, type="o")
legend("topright", legend = "Median parameters", col = "Black", pch = 7, lty=2)
dev.off()


##full fit
##get table for all observations
data.lr <- cbind(rbind(bins, c(Inf,Inf,Inf,Ntrans)), Ints=rowSums(obsJack))
data.lr$N <- data.lr$N
data.lr$notInts <- with(data.lr, N - Ints)
xlims <- range(log(data.lr$mid))
ylims <- range(log(p.fit(log(data.lr$mid), medianJackPar[1:4])))

pdf(paste0(outputPrefix, "_curveFit.pdf"))
plot(log(data.lr$mid), log(p.fit(log(data.lr$mid), medianJackPar[1:4])), type="l", ylim = c(-19,-4), col="Red", main="P-value weighting curve: fit to data", xlab="log_e(distance)", ylab="log_e(probability of interaction)")
points(log(data.lr$mid), log(data.lr$Ints) - log(data.lr$N), col="black")
legend("topright", legend = c("Data", "Fit"), col = c("Black", "Red"), pch = c(1, NA), lty=c(0,1))
dev.off()


## FINAL OUTPUT
outputTable <- cbind(c("weightAlpha", "weightBeta", "weightGamma", "weightDelta"), medianJackPar[1:4])
write.table(outputTable, paste0(outputPrefix,".settings"), quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
message("Output saved.")
