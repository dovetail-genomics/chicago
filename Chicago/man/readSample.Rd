\name{readSample}
\alias{readSample}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Read Sample
}
\description{
This function reads input data from a file, into a \code{chicagoData} object.
}
\usage{
readSample(file, cd)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file}{
The location of an input file FIXME more details!
}
  \item{cd}{A \code{chicagoData} object.}
}
\section{Warning}{
The object \code{cd@x} is updated by reference. Thus, \code{cd@x} will be altered. See vignette for further information.
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

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
filesDir <- file.path(system.file("extdata", package="Chicago"), "unitTestData")
file <- file.path(filesDir, dir(filesDir))[1]
print(file) ##we will read in this file

designDir <- file.path(system.file("extdata", package="Chicago"), "unitTestDesign")

cd <- setExperiment(designDir=designDir)
cd <- readAndMerge(file, cd)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
