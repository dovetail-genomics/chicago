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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (cd, outprefix, printMemory = FALSE) 
{
    message("\n*** Running normaliseBaits...\n")
    cd = normaliseBaits(cd)
    if (printMemory) {
        print(gc(reset = TRUE))
    }
    message("\n*** Running normaliseOtherEnds...\n")
    cd = normaliseOtherEnds(cd, outfile = paste0(outprefix, "_oeNorm.pdf"))
    if (printMemory) {
        print(gc(reset = TRUE))
    }
    message("\n*** Running estimateTechnicalNoise...\n")
    cd = estimateTechnicalNoise(cd, outfile = paste0(outprefix, 
        "_techNoise.pdf"))
    if (printMemory) {
        print(gc(reset = TRUE))
    }
    message("\n*** Running estimateDistFun...\n")
    cd = estimateDistFun(cd, outfile = paste0(outprefix, "_distFun.pdf"))
    if (printMemory) {
        print(gc(reset = TRUE))
    }
    message("\n*** Running estimateBrownianNoise...\n")
    cd = estimateBrownianNoise(cd)
    if (printMemory) {
        print(gc(reset = TRUE))
    }
    message("\n*** Running getPvals...\n")
    cd = getPvals(cd)
    if (printMemory) {
        print(gc(reset = TRUE))
    }
    message("\n*** Running getScores...\n")
    cd = getScores(cd)
    if (printMemory) {
        print(gc(reset = TRUE))
    }
    cd
  }
}
