\name{normaliseBaits}
\alias{normaliseBaits}
\title{
Normalise Baits
}
\description{
Calculate normalisation factors \code{s_j} for each bait.
}
\usage{
normaliseBaits(cd, normNcol = "NNb", shrink = FALSE,
      plot = TRUE, outfile = NULL, debug = FALSE)
}
\arguments{
  \item{cd}{A \code{chicagoData} object.}
  \item{normNcol}{
The name of the column in \code{cd} that contains normalised counts.
}
  \item{shrink}{
Deprecated.
}
  \item{plot}{
If \code{TRUE}, output a diagnostic plot.
}
  \item{outfile}{\code{NULL}, or a character string. If \code{NULL}, the diagnostic plot is outputted to the current plotting device. If a character string, e.g. \code{outfile="tech.pdf"}, then the plot will be outputted to that file.}
  \item{debug}{
Deprecated.
}
}
\section{Warning}{
The object \code{intData(cd)} is updated by reference. Thus, \code{intData(cd)} will be altered. See vignette for further information.
}
\details{
A \code{chicagoData} object: \code{intData(cd)} gains a new column \code{s_j}, and normalised output \code{NNb} (unless the \code{normNcol} parameter is altered).
}
\value{
An object of class \code{chicagoData}.
}
\author{
Mikhail Spivakov, Jonathan Cairns, Paula Freire Pritchett
}

\examples{
data(cdUnitTest)

##modifications to cdUnitTest, ensuring it uses correct design directory
designDir <- file.path(system.file("extdata", package="Chicago"), "unitTestDesign")
cdUnitTest <- modifySettings(cd=cdUnitTest, designDir=designDir)

cdUnitTest <- normaliseBaits(cdUnitTest)

}
