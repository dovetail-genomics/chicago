\name{chicagoPipeline}
\alias{chicagoPipeline}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
chicagoPipeline(cd, outprefix, printMemory = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cd}{
%%     ~~Describe \code{cd} here~~
}
  \item{outprefix}{
%%     ~~Describe \code{outprefix} here~~
}
  \item{printMemory}{
%%     ~~Describe \code{printMemory} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
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
        print(gc(reset = T))
    }
    message("\n*** Running normaliseOtherEnds...\n")
    cd = normaliseOtherEnds(cd, outfile = paste0(outprefix, "_oeNorm.pdf"))
    if (printMemory) {
        print(gc(reset = T))
    }
    message("\n*** Running estimateTechnicalNoise...\n")
    cd = estimateTechnicalNoise(cd, outfile = paste0(outprefix, 
        "_techNoise.pdf"))
    if (printMemory) {
        print(gc(reset = T))
    }
    message("\n*** Running estimateDistFun...\n")
    cd = estimateDistFun(cd, outfile = paste0(outprefix, "_distFun.pdf"))
    if (printMemory) {
        print(gc(reset = T))
    }
    message("\n*** Running estimateBrownianNoise...\n")
    cd = estimateBrownianNoise(cd)
    if (printMemory) {
        print(gc(reset = T))
    }
    message("\n*** Running getPvals...\n")
    cd = getPvals(cd)
    if (printMemory) {
        print(gc(reset = T))
    }
    message("\n*** Running getScores...\n")
    cd = getScores(cd)
    if (printMemory) {
        print(gc(reset = T))
    }
    cd
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
