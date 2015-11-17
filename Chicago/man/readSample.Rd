\name{readSample}
\alias{readSample}
\title{
Read Sample
}
\description{
This function reads input data from a file, into a \code{chicagoData} object.
}
\usage{
readSample(file, cd)
}
\arguments{
  \item{file}{
The location of an input file FIXME more details!
}
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
\examples{
filesDir <- file.path(system.file("extdata", package="Chicago"), "unitTestData")
file <- file.path(filesDir, dir(filesDir))[1]
print(file) ##we will read in this file

designDir <- file.path(system.file("extdata", package="Chicago"), "unitTestDesign")

cd <- setExperiment(designDir=designDir)
cd <- readAndMerge(file, cd)
}
