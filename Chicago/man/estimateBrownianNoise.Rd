\name{estimateBrownianNoise}
\alias{estimateBrownianNoise}
\title{
Estimate Brownian Noise
}
\description{
Estimates the dispersion, and adds a a \code{Bmean} column giving the mean Brownian noise.
}
\usage{
estimateBrownianNoise(cd)
}
\arguments{
  \item{cd}{A \code{chicagoData} object.}
}
\section{Warning}{
The object \code{intData(x)} is updated by reference. Thus, \code{intData(cd)} will be altered. See vignette for further information.
}
\value{
An object of class \code{chicagoData}.
}
\author{
Mikhail Spivakov, Jonathan Cairns, Paula Freire Pritchett
}

\seealso{
\code{\link{chicagoPipeline}}
}
\examples{

data(cdUnitTest)

##modifications to cdUnitTest, ensuring it uses correct design directory
designDir <- file.path(system.file("extdata", package="Chicago"), "unitTestDesign")
cdUnitTest <- modifySettings(cd=cdUnitTest, designDir=designDir)

##make cdUnitTest use the full subset of baits
cdUnitTest <- modifySettings(cd=cdUnitTest, settings=list(brownianNoise.subset=NA))

cdUnitTest <- estimateBrownianNoise(cdUnitTest)
}
