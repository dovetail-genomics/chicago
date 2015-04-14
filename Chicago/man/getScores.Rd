\name{getScores}
\alias{getScores}
\title{
Get CHiCAGO scores.
}
\description{
Converts p-values into a CHiCAGO score, using p-value weighting.
}
\usage{
getScores(cd, method = "weightedRelative", 
      includeTrans = TRUE, plot = TRUE, outfile = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cd}{A \code{chicagoData} object.}
  \item{method}{
  Either "weightedRelative" (recommended), or "unweighted".
}
  \item{includeTrans}{
  If \code{FALSE}, trans interactions are discounted.
}
  \item{plot}{
  Plot a diagnostic plot.
}
  \item{outfile}{
  A string containing a .pdf file location to write to.
}
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
\references{
Genovese, C. R., Roeder, K., and Wasserman, L. (2006). False discovery control with p-value weighting. Biometrika, 93, 509-524. doi:10.1093/biomet/93.3.509
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
function (cd, method = "weightedRelative", includeTrans = TRUE, 
    plot = T, outfile = NULL) 
{
    x <- cd@x
    set <- cd@settings
    avgFragLen <- .getAvgFragLength(cd)
    if (!method \%in\% c("unweighted", "weightedRelative")) {
        stop("method=", method, " not recognized.")
    }
    if (!includeTrans) {
        x <- x[!is.na(x$distSign), ]
    }
    if (method == "weightedRelative") {
        message("Calculating p-value weights...")
        x[, `:=`(log.w, .getWeights(abs(x$distSign), cd, includeTrans = includeTrans))]
        x[, `:=`(log.q, log.p - log.w)]
        message("Calculating scores...")
        minval <- .getWeights(0, cd, includeTrans = includeTrans)
        x[, `:=`(score, pmax(-minval - log.q, 0))]
    }
    else {
        stop("Only method='weightedRelative' available currently.")
    }
    cd
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
