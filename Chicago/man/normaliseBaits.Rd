\name{normaliseBaits}
\alias{normaliseBaits}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
normaliseBaits(cd, normNcol = "NNb", shrink = FALSE, plot = TRUE, outfile = NULL, debug = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cd}{A \code{chicagoData} object.}
  \item{normNcol}{
%%     ~~Describe \code{normNcol} here~~
}
  \item{shrink}{
%%     ~~Describe \code{shrink} here~~
}
  \item{plot}{
%%     ~~Describe \code{plot} here~~
}
  \item{outfile}{
%%     ~~Describe \code{outfile} here~~
}
  \item{debug}{
%%     ~~Describe \code{debug} here~~
}
}
\warning{
The object \code{cd} is updated by reference. Thus, \code{cd} will be altered. See vignette for further information.
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
An object of class \code{chicagoData}.
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
function (cd, normNcol = "NNb", shrink = FALSE, plot = TRUE, 
    outfile = NULL, debug = FALSE) 
{
    message("Normalising baits...")
    adjBait2bait = cd@settings$adjBait2bait
    if (debug) {
        .normaliseFragmentSets(x = cd@x, s = cd@settings, npb = .readNPBfile(s = cd@settings), 
            viewpoint = "bait", idcol = "baitID", Ncol = "N", 
            adjBait2bait = adjBait2bait, shrink = shrink, refExcludeSuffix = NULL, 
            plot = plot, outfile = outfile, debug = T)
    }
    else {
        cd@x = .normaliseFragmentSets(x = cd@x, s = cd@settings, 
            npb = .readNPBfile(s = cd@settings), viewpoint = "bait", 
            idcol = "baitID", Ncol = "N", adjBait2bait = adjBait2bait, 
            shrink = shrink, refExcludeSuffix = NULL, plot = plot, 
            outfile = outfile, debug = F)
    }
    cd@x[, `:=`((normNcol), pmax(1, round(N/s_j)))]
    setkey(cd@x, baitID, otherEndID)
    othercols = names(cd@x)[!names(cd@x) \%in\% c("distbin", normNcol)]
    setcolorder(cd@x, c(othercols, "distbin", normNcol))
    cd
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
