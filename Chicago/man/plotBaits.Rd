\name{plotBaits}
\alias{plotBaits}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot Baits
}
\description{
Plot the read counts around baits.
}
\usage{
plotBaits(cd, pcol = "score", Ncol = "N", n = 16, baits = NULL, plotBaitNames = TRUE, plotBprof = FALSE, plevel1 = 5, plevel2 = 3, outfile = NULL, removeBait2bait = TRUE, width = 20, height = 20, maxD = NULL, bgCol = "black", lev2Col = "blue", lev1Col = "red", bgPch = 1, lev1Pch = 20, lev2Pch = 20, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cd}{A \code{chicagoData} object.}
  \item{pcol}{
The name of the column that contains the score.
}
  \item{Ncol}{
The name of the column that contains counts.
}
  \item{n}{
The number of baits to plot (ignored if \code{baits} is specified).
}
  \item{baits}{
The IDs of the baits to plot.
}
  \item{plotBaitNames}{
If \code{TRUE}, the names of the baits, rather than their IDs, will appear in the plot.
}
  \item{plotBprof}{
If \code{TRUE}, display a line representing the expected Brownian noise at each distance.
}
  \item{plevel1, plevel2}{
%%     ~~Describe \code{plevel1} here~~
}
%  \item{plevel2}{
%%     ~~Describe \code{plevel2} here~~
%}
  \item{outfile}{
%%     ~~Describe \code{outfile} here~~
}
  \item{removeBait2bait}{
%%     ~~Describe \code{removeBait2bait} here~~
}
  \item{width}{
%%     ~~Describe \code{width} here~~
}
  \item{height}{
%%     ~~Describe \code{height} here~~
}
  \item{maxD}{
%%     ~~Describe \code{maxD} here~~
}
  \item{bgCol}{
%%     ~~Describe \code{bgCol} here~~
}
  \item{lev2Col}{
%%     ~~Describe \code{lev2Col} here~~
}
  \item{lev1Col}{
%%     ~~Describe \code{lev1Col} here~~
}
  \item{bgPch}{
%%     ~~Describe \code{bgPch} here~~
}
  \item{lev1Pch}{
%%     ~~Describe \code{lev1Pch} here~~
}
  \item{lev2Pch}{
%%     ~~Describe \code{lev2Pch} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
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
  function(FIXME) {FIXME}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
