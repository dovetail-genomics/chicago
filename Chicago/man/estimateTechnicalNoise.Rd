\name{estimateTechnicalNoise}
\alias{estimateTechnicalNoise}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimate Technical Noise
}
\description{
Calculates the expected technical noise based on trans read pairs.
}
\usage{
estimateTechnicalNoise(cd, plot = TRUE, outfile = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cd}{A \code{chicagoData} object.}
  \item{plot}{
Logical - if \code{TRUE}, then output a diagnostic plot.
}
  \item{outfile}{\code{NULL}, or a character string. If \code{NULL}, the diagnostic plot is outputted to the current plotting device. If a character string, e.g. \code{outfile="tech.pdf"}, then the plot will be outputted to that file.}
}
\warning{
The object \code{cd} is updated by reference. Thus, \code{cd} will be altered. See vignette for further information.
}
\details{

}
\value{
An object of class \code{chicagoData}, with additional columns "tlb", "tblb", "Tmean".
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Jonathan Cairns
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{chicagoPipeline}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (cd, plot = TRUE, outfile = NULL) 
{
    message("Estimating technical noise based on trans-counts...")
    minBaitsPerBin = cd@settings$techNoise.minBaitsPerBin
    adjBait2bait = cd@settings$adjBait2bait
    if (!"tlb" \%in\% names(cd@x)) {
        message("Binning other ends based on trans-counts...")
        cd@x = .addTLB(cd)
    }
    setkey(cd@x, baitID)
    message("Binning baits based on observed trans-counts...")
    transBaitLen = cd@x[, sum(is.na(distSign)), by = baitID]
    setnames(transBaitLen, "V1", "transBaitLen")
    transBaitLen$tblb = cut2(transBaitLen$transBaitLen, m = minBaitsPerBin, 
        levels.mean = FALSE)
    setkey(transBaitLen, baitID)
    setkey(cd@x, baitID)
    cd@x = cd@x[transBaitLen]
    message("Defining interaction pools and gathering the observed numbers of trans-counts per pool...")
    setkeyv(cd@x, c("tlb", "tblb"))
    Ntrans = cd@x[, {
        res = table(N[is.na(distSign)])
        list(as.numeric(names(res)), as.numeric(res))
    }, by = c("tlb", "tblb")]
    setnames(Ntrans, "V1", "N")
    setnames(Ntrans, "V2", "nobs")
    message("Computing the total number of possible interactions per pool...")
    message("Preparing the data...", appendLF = F)
    baitmap = fread(cd@settings$baitmapfile)
    rmap = fread(cd@settings$rmapfile)
    setnames(rmap, "V1", "chr")
    setnames(baitmap, "V1", "chr")
    setnames(rmap, "V4", "otherEndID")
    setnames(baitmap, "V4", "baitID")
    message(".", appendLF = F)
    baitmap = baitmap[, c("baitID", "chr"), with = FALSE]
    rmap = rmap[, c("otherEndID", "chr"), with = FALSE]
    setkey(baitmap, baitID)
    setkey(rmap, otherEndID)
    message(".", appendLF = F)
    setkey(cd@x, baitID)
    cd@x = baitmap[cd@x]
    setnames(cd@x, "chr", "baitChr")
    setkey(cd@x, otherEndID)
    cd@x = rmap[cd@x]
    setnames(cd@x, "chr", "otherEndChr")
    setkey(cd@x, tlb, tblb)
    message("\nProcessing fragment pools", appendLF = FALSE)
    res = cd@x[, {
        message(".", appendLF = FALSE)
        baits = unique(baitID)
        oes = unique(otherEndID)
        bChr = baitChr[!duplicated(baitID)]
        oeChr = otherEndChr[!duplicated(otherEndID)]
        numPairs = sum(sapply(unique(bChr), function(this) length(bChr[bChr == 
            this]) * length(oeChr[oeChr != this]))) - length(baits[baits \%in\% 
            oes])
        nTrans = sum(N[is.na(distSign)])
        Tmean = nTrans/numPairs
        list(numPairs = numPairs, nTrans = nTrans, Tmean = Tmean)
    }, by = c("tlb", "tblb")]
    if (plot) {
        message("\nPlotting...")
        if (!is.null(outfile)) {
            pdf(outfile)
        }
        par(mfrow = c(2, 1))
        boxplot(Tmean ~ tblb, as.data.frame(res), main = "Technical noise estimates per bait pool")
        boxplot(Tmean ~ tlb, as.data.frame(res), main = "Technical noise estimates per other end pool")
        if (!is.null(outfile)) {
            dev.off()
        }
        par(mfrow = c(1, 1))
    }
    message("Post-processing the results...")
    setkeyv(cd@x, c("tlb", "tblb"))
    setkeyv(res, c("tlb", "tblb"))
    cd@x = res[cd@x]
    set(cd@x, NULL, "transBaitLen", NULL)
    set(cd@x, NULL, "numPairs", NULL)
    set(cd@x, NULL, "nTrans", NULL)
    set(cd@x, NULL, "baitChr", NULL)
    set(cd@x, NULL, "otherEndChr", NULL)
    setkey(cd@x, baitID, otherEndID)
    othernames = names(cd@x)[!names(cd@x) \%in\% c("baitID", "otherEndID", 
        "tlb", "tblb", "Tmean")]
    setcolorder(cd@x, c("baitID", "otherEndID", othernames, "tlb", 
        "tblb", "Tmean"))
    cd
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
