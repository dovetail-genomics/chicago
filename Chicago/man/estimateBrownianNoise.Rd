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

data(cdUnitTest)

##modifications to cdUnitTest, ensuring it uses correct design directory
designDir <- file.path(system.file("extdata", package="Chicago"), "unitTestDesign")
new.settings <- paste0(file.path(designDir, "unitTest"), c(".rmap", ".baitmap", ".npb", ".nbpb", ".poe"))
new.settings <- as.list(new.settings)
names(new.settings) <- c("rmapfile", "baitmapfile", "nperbinfile", "nbaitsperbinfile", "proxOEfile")
cdUnitTest <- modifySettings(cd = cdUnitTest, new.settings)

##make cdUnitTest use the full subset of baits
cdUnitTest <- modifySettings(cd = cdUnitTest, list(brownianNoise.subset=NA))

cdUnitTest <- estimateBrownianNoise(cdUnitTest)
}

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
