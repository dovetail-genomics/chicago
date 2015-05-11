\name{estimateTechnicalNoise}
\alias{estimateTechnicalNoise}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimate Technical Noise
}
\description{
Calculates the expected technical noise based on trans read pairs.
}
\usage{
estimateTechnicalNoise(cd, plot = TRUE, outfile = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cd}{A \code{chicagoData} object.}
  \item{plot}{
Logical - if \code{TRUE}, then output a diagnostic plot.
}
  \item{outfile}{\code{NULL}, or a character string. If \code{NULL}, the diagnostic plot is outputted to the current plotting device. If a character string, e.g. \code{outfile="tech.pdf"}, then the plot will be outputted to that file.}
}
\section{Warning}{
The object \code{cd@x} is updated by reference. Thus, \code{cd@x} will be altered. See vignette for further information.
}
%\details{

%}
\value{
An object of class \code{chicagoData}, with additional columns "tlb", "tblb", "Tmean".
}
\author{
Mikhail Spivakov, Jonathan Cairns, Paula Freire Pritchett
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{chicagoPipeline}}
}
\examples{
data(cdUnitTest)

##modifications to cdUnitTest, ensuring it uses correct design directory
designDir <- file.path(system.file("extdata", package="Chicago"), "unitTestDesign")
rmapfileLoc <- file.path(designDir, "unitTest.rmap")
baitmapfileLoc <- file.path(designDir, "unitTest.rmap")
cdUnitTest <- modifySettings(cd = cdUnitTest, list(rmapfile=rmapfileLoc, baitmapfile=baitmapfileLoc))

cdUnitTest <- estimateTechnicalNoise(cdUnitTest)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
