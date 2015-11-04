\name{estimateTechnicalNoise}
\alias{estimateTechnicalNoise}
\title{
Estimate Technical Noise
}
\description{
Calculates the expected technical noise based on trans read pairs.
}
\usage{
estimateTechnicalNoise(cd, plot = TRUE, outfile = NULL)
}
\arguments{
  \item{cd}{A \code{chicagoData} object.}
  \item{plot}{
Logical - if \code{TRUE}, then output a diagnostic plot.
}
  \item{outfile}{\code{NULL}, or a character string. If \code{NULL}, the diagnostic plot is outputted to the current plotting device. If a character string, e.g. \code{outfile="tech.pdf"}, then the plot will be outputted to that file.}
}
\section{Warning}{
The object \code{intData(cd)} is updated by reference. Thus, \code{intData(cd)} will be altered. See vignette for further information.
}
\value{
An object of class \code{chicagoData}, with additional columns "tlb", "tblb", "Tmean".
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

cdUnitTest <- estimateTechnicalNoise(cdUnitTest)
}
