\name{normaliseOtherEnds}
\alias{normaliseOtherEnds}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Normalise Other Ends
}
\description{
Compute \code{s_i} normalisation factors for other ends.
}
\usage{
normaliseOtherEnds(cd, Ncol = "NNb", normNcol = "NNboe", plot = TRUE, outfile = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cd}{A \code{chicagoData} object.}
  \item{Ncol}{
The name of the column in \code{cd} that contains normalised counts.
}
  \item{normNcol}{
FIXME%%     ~~Describe \code{normNcol} here~~
}
  \item{plot}{
If \code{TRUE}, output a diagnostic plot.
}
  \item{outfile}{\code{NULL}, or a character string. If \code{NULL}, the diagnostic plot is outputted to the current plotting device. If a character string, e.g. \code{outfile="tech.pdf"}, then the plot will be outputted to that file.}
}
\section{Warning}{
The object \code{cd} is updated by reference. Thus, \code{cd} will be altered. See vignette for further information.
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
An object of class \code{chicagoData}.
}
\author{
Mikhail Spivakov, Jonathan Cairns, Paula Freire Pritchett
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
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
