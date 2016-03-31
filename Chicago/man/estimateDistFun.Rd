\name{estimateDistFun}
\alias{estimateDistFun}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimate the Distance Function
}
\description{
Estimates the function that models how the expected number of counts decreases with increasing distance.
}
\usage{
estimateDistFun(cd, method = "cubic", plot = TRUE, outfile = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cd}{A \code{chicagoData} object.}
  \item{method}{
Choice of method: "cubic" is currently the only allowed option, which fits a cubic function with linear extrapolation, on a log-log scale.
}
%  \item{n.obs.head}{
%If \code{method="lm"}, the number of observations in the head (see details).
%}
%  \item{n.obs.tail}{
%If \code{method="lm"}, the number of observations in the tail (see details).
%}
%  \item{logScale}{
%If \code{TRUE}, the function \code{f()} will output on the log.scale. Note that \code{f()}'s input is always assumed to be on the log scale.
%}
  \item{plot}{
  Output a diagnostic plot.
}
  \item{outfile}{
  If \code{NULL}, plot to current device. Otherwise, plot to the .pdf file \code{outfile}.
}
}
\details{
By default, we look in 75 distance bins, and a cubic fit is used. For distances that lie outside of the bin boundaries, it is assumed that the function is log-linear, with continuity of \code{f} and its first derivative on the log-scale. %If \code{method="lm"}, then it is assumed that the distance function consists of two segments. The first \code{n.obs.head} are used to estimate the first segment, and the last \code{n.obs.tail} are used to estimate the second segment.
}
\value{
An object of class \code{chicagoData}, with the parameters of the distance function present as params(cd)$distFunParams.
}
\author{
Mikhail Spivakov, Jonathan Cairns, Paula Freire Pritchett
}

\seealso{
\code{\link{chicagoPipeline}}, \code{\link{plotDistFun}}
}
\examples{
data(cdUnitTest)
estimateDistFun(cdUnitTest)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
