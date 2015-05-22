\name{modifySettings}
\alias{modifySettings}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Modify Settings
}
\description{
Modify the settings in a \code{chicagoData} object.
}
\usage{
modifySettings(cd, settings)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cd}{A \code{chicagoData} object.}
  \item{settings}{A named list.}
}
\details{
\code{cd}'s settings are updated. Setting \code{names(settings)[1]} is set to (settings)[[1]], and so on. For a list of available settings, see \code{\link{defaultSettings}}.
}
\value{
An object of class \code{chicagoData}.
}
\author{
Mikhail Spivakov, Jonathan Cairns, Paula Freire Pritchett
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{defaultSettings}}
}
\examples{
designDir <- file.path(system.file("extdata", package="Chicago"), "unitTestDesign")

cd <- setExperiment(designDir)

##Suppose I want to increase the zone in which Brownian noise is assumed to dominate
##(by default, 1500000):
print(cd@settings$maxLBrownEst)

cd <- modifySettings(cd, list(maxLBrownEst = 3000000L))
print(cd@settings$maxLBrownEst)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
