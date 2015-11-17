\name{getScores}
\alias{getScores}
\title{
Get CHiCAGO scores.
}
\description{
Converts p-values into a CHiCAGO score, using p-value weighting.
}
\usage{
getScores(cd, method = "weightedRelative", 
      includeTrans = TRUE, plot = TRUE, outfile = NULL)
}
\arguments{
  \item{cd}{A \code{chicagoData} object.}
  \item{method}{
  Either "weightedRelative" (recommended), or "unweighted".
}
  \item{includeTrans}{
  If \code{FALSE}, trans interactions are discounted.
}
  \item{plot}{
  Plot a diagnostic plot.
}
  \item{outfile}{
  A string containing a .pdf file location to write to.
}
}
\section{Warning}{
The object \code{intData(cd)} is updated by reference. Thus, \code{intData(cd)} will be altered. See vignette for further information.
}
\details{
Weighting is performed using the parameters \code{weightAlpha}, \code{weightBeta}, \code{weightGamma}, \code{weightDelta}. Briefly, this function calculates weights \code{w} that decrease with increasing distance. Then, we construct weighted p-values \code{p/w}. As a result, the significance of long-range interactions is upweighted, and the significance of short-range interactions is downweighted.

Finally, the output score is calculated as \code{-log(p/w) - log(w_max)}, where \code{w_max} is the highest attainable weight, and provided the score is positive (otherwise it is set to 0).

Please see the CHiCAGO paper and its supplementary for full details.
}
\value{
An object of class \code{chicagoData}.
}
\references{
Genovese, C. R., Roeder, K., and Wasserman, L. (2006). False discovery control with p-value weighting. Biometrika, 93, 509-524. doi:10.1093/biomet/93.3.509
}
\author{
Jonathan Cairns
}

\seealso{
\code{\link{chicagoPipeline}}
}
\examples{
data(cdUnitTest)

##modifications to cdUnitTest, ensuring it uses correct design directory
designDir <- file.path(system.file("extdata", package="Chicago"), "unitTestDesign")
cdUnitTest <- modifySettings(cd=cdUnitTest, designDir=designDir)

cdUnitTest <- getScores(cdUnitTest)
}
