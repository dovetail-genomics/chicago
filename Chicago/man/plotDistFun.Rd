\name{plotDistFun}
\alias{plotDistFun}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot the Distance Function
}
\description{
Estimates the function that models how the expected number of counts decreases with increasing distance.
}
\usage{
plotDistFun(cd, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cd}{A \code{chicagoData} object.}
  \item{\dots}{Further arguments passed to \code{\link{plot}}.}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
A plot.
}
\author{
Jonathan Cairns
}

\seealso{
\code{\link{estimateDistFun}}
}
\examples{
data(cdUnitTest)
plotDistFun(cdUnitTest)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
