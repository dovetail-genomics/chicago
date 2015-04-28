\name{mergeSamples}
\alias{mergeSamples}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Merge samples together.
}
\description{
Merge a number of \code{chicagoData} objects together, summarising their counts into a normalised value.
}
\usage{
mergeSamples(cdl, normalise = TRUE, NcolOut = "N",
      NcolNormPrefix = "NNorm", mergeMethod = c("weightedMean", "mean")[1])
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cdl}{A \code{list} of \code{chicagoData} objects.}
  \item{normalise}{
If \code{TRUE}, use a normalisation procedure, specified by \code{mergeMethod}, to arrive at a normalised count. If \code{FALSE}, take the mean number of reads.
}
  \item{NcolOut}{
The column to store the normalised counts in.
}
  \item{NcolNormPrefix}{
Each sample gains a normalised count column, that begins with this prefix.
}
  \item{mergeMethod}{
If \code{mergeMethod == "weightedMean"}, then NcolOut is the weighted mean of the sample-wise counts adjusted by the samples' respective scaling factors \code{s_k}. If \code{mergeMethod == "normMean"}, then sample-specific counts are first normalised by dividing by s_k, and NcolOut is computed as the mean of these normalised counts.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
An object of class \code{chicagoData}.
}
\author{
Mikhail Spivakov, Jonathan Cairns, Paula Freire Pritchett
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{readAndMerge}}
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
