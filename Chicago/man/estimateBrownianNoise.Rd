\name{estimateBrownianNoise}
\alias{estimateBrownianNoise}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimate Brownian Noise
}
\description{
Estimates the dispersion, and adds a a \code{Bmean} column giving the mean Brownian noise.
}
\usage{
estimateBrownianNoise(cd)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cd}{A \code{chicagoData} object.}
}
\section{Warning}{
The object \code{cd@x} is updated by reference. Thus, \code{cd@x} will be altered. See vignette for further information.
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
\code{\link{chicagoPipeline}}
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
