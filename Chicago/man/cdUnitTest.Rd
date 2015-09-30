\name{cdUnitTest}
\alias{cdUnitTest}
\docType{data}
\title{
ChicagoData object for unit testing
}
\description{
This data set is used for unit testing - it is too small to run all of the steps of CHiCAGO. For a toy data set that is large enough, please see the data package. (Note that cdUnitTest is a subset of those data.)
}
\usage{data("cdUnitTest")}

\details{
The data are derived from mouse ESCs. They are a subset of the object \code{\link[PCHiCdata]{smESC}} (from the \code{PCHiCdata} package)
}

\value{
A \code{\link{chicagoData}} object.
}

\source{
Schoenfelder, S. et al. "The pluripotent regulatory circuitry connecting promoters to their long-range interacting elements." Genome research 25.4 (2015): 582-597.
}
\seealso{
  \code{\link[PCHiCdata]{smESC}}, \code{\link{chicagoData}}
}
\examples{
data(cdUnitTest)

##modifications to cdUnitTest, ensuring it uses correct design directory
designDir <- file.path(system.file("extdata", package="Chicago"), "unitTestDesign")
cdUnitTest <- modifySettings(cd=cdUnitTest, designDir=designDir)

print(cdUnitTest)
}
\keyword{datasets}
