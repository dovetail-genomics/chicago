\name{chicagoPipeline}
\alias{chicagoPipeline}
\title{
CHiCAGO pipeline function
}
\description{
This function runs data through the CHiCAGO pipeline.
}
\usage{
chicagoPipeline(cd, outprefix = NULL, printMemory = FALSE)
}
\arguments{
  \item{cd}{A \code{chicagoData} object.}
  \item{outprefix}{\code{NULL}, or a character string. If \code{NULL}, diagnostic plots are outputted to the current plotting device. If a character string, then pdfs will be generated for a series of diagnostic plots, in files of form "[outprefix]_[plotname].pdf". For example, \code{outprefix="experiment1"} leads to files experiment1_oeNorm.pdf, etc...}
  \item{printMemory}{Set to \code{TRUE} for memory diagnostics.}
}
\section{Warning}{
The object \code{intData(cd)} is updated by reference. Thus, \code{intData(cd)} will be altered. See vignette for further information.
}
\details{
This pipeline runs the following functions in order:

- \code{\link{normaliseBaits}}

- \code{\link{normaliseOtherEnds}}

- \code{\link{estimateTechnicalNoise}}

- \code{\link{estimateDistFun}}

- \code{\link{estimateBrownianNoise}}

- \code{\link{getPvals}}

- \code{\link{getScores}}

It does not export the output. Use \code{\link{exportResults}} for this.
}
\value{
An object of class \code{chicagoData}.
}
\author{
Mikhail Spivakov, Jonathan Cairns, Paula Freire Pritchett
}

\seealso{
\code{\link{exportResults}}
}
\examples{
  ##Read in some raw data
  filesDir <- file.path(system.file("extdata", package="Chicago"), "unitTestData")
  file <- file.path(filesDir, dir(filesDir))[1]
  print(file) ##we will read in this file

  designDir <- file.path(system.file("extdata", package="Chicago"), "unitTestDesign")
  
  ##Add a setting specific to the unit test data! Do not use in practice!
  if(!interactive()) {
    settings <- list(brownianNoise.samples=1)
  } else {
    settings <- NULL
  }

  cd <- setExperiment(designDir=designDir, settings=settings)
  cd <- readAndMerge(file, cd)
}
