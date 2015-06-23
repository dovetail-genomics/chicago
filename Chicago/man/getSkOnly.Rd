\name{getSkOnly}
\alias{getSkOnly}

\title{
Get S_k factors from multiple replicates
}
\description{
Finds s_k scaling factors for a (potentially large) number of samples. Typically, these factors are used as library size factors in some sort of differential count algorithm (DESeq, EdgeR, baySeq, ...) to find differential binding events between samples.
}
\usage{
getSkOnly(files, cd)
}

\arguments{
  \item{files}{
  Character vector containing the locations of the .chinput files to read in.
}
  \item{cd}{
  A blank \code{chicagoData} object for reference, usually created with {\code{\link{setExperiment}}}.
}
}
\value{
Numeric vector of s_k factors.
}
\author{
Jonathan Cairns, Paula Freire Pritchett, Mikhail Spivakov
}

\seealso{
readAndMerge
}
\examples{
filesDir <- file.path(system.file("extdata", package="Chicago"), "unitTestData")
files <- file.path(filesDir, dir(filesDir))
print(files) ##we will read in and merge these files

designDir <- file.path(system.file("extdata", package="Chicago"), "unitTestDesign")

cd <- setExperiment(designDir=designDir)
s_k <- getSkOnly(files, cd)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
