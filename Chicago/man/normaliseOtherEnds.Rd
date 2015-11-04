\name{normaliseOtherEnds}
\alias{normaliseOtherEnds}
\title{
Normalise Other Ends
}
\description{
Compute \code{s_i} normalisation factors for other ends, and normalised counts.
}
\usage{
normaliseOtherEnds(cd, Ncol = "NNb", normNcol = "NNboe", plot = TRUE, outfile = NULL)
}
\arguments{
  \item{cd}{A \code{chicagoData} object.}
  \item{Ncol}{
The name of an input column in \code{intData(cd)} that contains counts normalised by bait (i.e. it is output from \code{\link{normaliseBaits}}.
}
  \item{normNcol}{
The name of an output column that will contain counts normalised by other ends (in addition to any normalisation already performed on the \code{Ncol} column). Useful for plotting.
}
  \item{plot}{
If \code{TRUE}, output a diagnostic plot.
}
  \item{outfile}{\code{NULL}, or a character string. If \code{NULL}, the diagnostic plot is outputted to the current plotting device. If a character string, e.g. \code{outfile="tech.pdf"}, then the plot will be outputted to that file.}
}
\section{Warning}{
The object \code{intData(cd)} is updated by reference. Thus, \code{intData(cd)} will be altered. See vignette for further information.
}
\details{
A \code{chicagoData} object: \code{intData(cd)} gains new columns \code{s_i}, and normalised output \code{NNboe} (unless the \code{normNcol} parameter is altered).
}
\value{
An object of class \code{chicagoData}.
}
\author{
Mikhail Spivakov, Jonathan Cairns, Paula Freire Pritchett
}

\examples{
##FIXME: example can be run by loading data package if it is installed, once it exists

if("PCHiCdata" \%in\% rownames(installed.packages()))
{
  library(PCHiCdata)
  data(smESC)
  
  ##modifiy smESC to use correct design directory
  designDir <- file.path(system.file("extdata", package="PCHiCdata"), "mm9TestDesign")
  smESC <- modifySettings(cd=smESC, designDir=designDir)
  
  ##normalise here...
  normaliseOtherEnds(smESC)
  
} else {
  warning("Please install the PCHiCdata package to run this example.")
}

}
