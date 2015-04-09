\name{normaliseOtherEnds}
\alias{normaliseOtherEnds}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
normaliseOtherEnds(cd, Ncol = "NNb", normNcol = "NNboe", plot = TRUE, outfile = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cd}{
%%     ~~Describe \code{cd} here~~
}
  \item{Ncol}{
%%     ~~Describe \code{Ncol} here~~
}
  \item{normNcol}{
%%     ~~Describe \code{normNcol} here~~
}
  \item{plot}{
%%     ~~Describe \code{plot} here~~
}
  \item{outfile}{
%%     ~~Describe \code{outfile} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (cd, Ncol = "NNb", normNcol = "NNboe", plot = TRUE, 
    outfile = NULL) 
{
    cd@x = .addTLB(cd)
    x = cd@x[abs(distSign) <= cd@settings$maxLBrownEst & is.na(distSign) == 
        F]
    message("Computing total bait counts...")
    nbpb = .readNbaitsPBfile(s = cd@settings)
    setkey(nbpb, otherEndID)
    setkey(x, otherEndID)
    nbpb = nbpb[x]
    setkeyv(nbpb, c("otherEndID", "distbin"))
    nbpb = unique(nbpb)
    nbpbSum = nbpb[, sum(get(paste0("bin", as.integer(distbin)))), 
        by = c("tlb", "distbin")]
    setkeyv(nbpbSum, c("tlb", "distbin"))
    setkeyv(x, c("tlb", "distbin"))
    x = nbpbSum[x]
    setnames(x, "V1", "ntot")
    message("Computing scaling factors...")
    x = .normaliseFragmentSets(x, s = cd@settings, viewpoint = "otherEnd", 
        idcol = "tlb", Ncol = Ncol, npb = NULL, shrink = FALSE, 
        adjBait2bait = FALSE, refExcludeSuffix = "B2B")
    setkey(x, tlb)
    x = unique(x)
    set(x, NULL, "NNb", NULL)
    set(x, NULL, "distbin", NULL)
    set(x, NULL, "ntot", NULL)
    if (plot) {
        if (!is.null(outfile)) {
            pdf(outfile)
        }
        with(x, barplot(s_i, names.arg = tlb, col = sapply(tlb, 
            function(x) ifelse(length(grep("B2B", x)), "darkblue", 
                "red")), xlab = "tlb", ylab = "s_i"))
        legend("topleft", legend = c("non-B2B", "B2B"), fill = c("red", 
            "darkblue"))
        if (!is.null(outfile)) {
            dev.off()
        }
    }
    message("Computing normalised counts...")
    setkey(cd@x, tlb)
    cd@x = merge(cd@x, x, all.x = T, by = "tlb")
    set(cd@x, which(is.na(cd@x$s_i)), "s_i", 1)
    cd@x[, `:=`((normNcol), pmax(1, round(get(Ncol)/s_i)))]
    message("Post-processing...")
    setkey(cd@x, baitID, otherEndID)
    othercols = names(cd@x)[!names(cd@x) \%in\% c("baitID", "otherEndID", 
        "distbin", Ncol, normNcol)]
    setcolorder(cd@x, c("baitID", "otherEndID", "distbin", othercols, 
        Ncol, normNcol))
    cd
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
