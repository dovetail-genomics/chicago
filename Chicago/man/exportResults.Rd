\name{exportResults}
\alias{exportResults}

\title{
Export Results
}
\description{
Export the results from a \code{chicagoData} object to disk.
}
\usage{
exportResults(cd, outfileprefix, scoreCol = "score", cutoff = 5, b2bcutoff = NULL,
        format = c("seqMonk", "interBed", "washU_track", "washU_text"), order = c("position", "score")[1])
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cd}{A \code{chicagoData} object.}
  \item{outfileprefix}{
  A character string that forms the prefix for each output file.
}
  \item{scoreCol}{
  The column of \code{cd@x} that contains the score.
}
  \item{cutoff}{
  The score cutoff.
}
  \item{b2bcutoff}{
  If desired, an alternative score cutoff for bait-to-bait interactions.
}
  \item{format}{
  The file format(s) to output. If a multiple formats are supplied as a vector, then all of these formats will be outputted.
}
  \item{order}{
  Should output be ordered by position or score?
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
NULL
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
  \dontrun{
    exportResults(cdUnitTest, outfileprefix = "unitTestOutput")
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
