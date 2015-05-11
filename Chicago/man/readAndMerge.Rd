\name{readAndMerge}
\alias{readAndMerge}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Read And Merge
}
\description{
A wrapper that calls readSample() on a number of files, then mergeSamples().
}
\usage{
readAndMerge(files, cd, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{files}{
  Character vector containing the locations of the files to read in.
}
  \item{cd}{
  A \code{chicagoData} object, usually created with {\code{\link{setExperiment}}}.
}
  \item{...}{
  Further arguments passed to \code{\link{mergeSamples}}.
}
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
\code{\link{readSample}}, \code{\link{mergeSamples}}
}
\examples{
filesDir <- file.path(system.file("extdata", package="Chicago"), "unitTestData")
files <- file.path(filesDir, dir(filesDir))
print(files) ##we will read in and merge these files

designDir <- file.path(system.file("extdata", package="Chicago"), "unitTestDesign")

cd <- setExperiment(designDir=designDir)
cd <- readAndMerge(files, cd)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
