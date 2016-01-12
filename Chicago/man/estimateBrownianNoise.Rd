\name{estimateBrownianNoise}
\alias{estimateBrownianNoise}
\title{
Estimate Brownian Noise
}
\description{
Estimates the dispersion, and adds a a \code{Bmean} column giving the mean Brownian noise.

Usually, the dispersion is not calculated on the full dataset - rather, a subsample of baits is taken, and the dispersion is calculated on that. The number of baits used is taken from \code{brownianNoise.subset} (with an \code{NA} value meaning that the entire dataset is used, and no subsampling is performed).
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
