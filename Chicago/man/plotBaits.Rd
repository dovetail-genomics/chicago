\name{plotBaits}
\alias{plotBaits}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot Baits
}
\description{
Plot the read counts around baits.
}
\usage{
plotBaits(cd, pcol = "score", Ncol = "N", n = 16, baits = NULL,
      plotBaitNames = TRUE, plotBprof = FALSE, plevel1 = 5, plevel2 = 3,
      outfile = NULL, removeBait2bait = TRUE, width = 20, height = 20,
      maxD = NULL, bgCol = "black", lev2Col = "blue", lev1Col = "red",
      bgPch = 1, lev1Pch = 20, lev2Pch = 20, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cd}{A \code{chicagoData} object.}
  \item{pcol}{
The name of the column that contains the score.
}
  \item{Ncol}{
The name of the column that contains counts.
}
  \item{n}{
The number of baits to plot (ignored if \code{baits} is specified).
}
  \item{baits}{
The IDs of the baits to plot.
}
  \item{plotBaitNames}{
If \code{TRUE}, the names of the baits, rather than their IDs, will appear in the plot.
}
  \item{plotBprof}{
If \code{TRUE}, display a line representing the expected Brownian noise at each distance.
}
  \item{plevel1, plevel2}{
Thresholds used on the \code{pcol} column. \code{plevel1} should be the more stringent threshold.
}
  \item{outfile}{
If \code{NULL}, output to current plotting device. Otherwise, this specifies a pdf file to write to.
}
  \item{removeBait2bait}{
If \code{TRUE}, bait-to-bait interactions are not plotted.
}
  \item{width, height}{
Passed through to \code{\link{pdf}}
}
  \item{maxD}{
The maximum distance to plot.
}
  \item{bgCol, lev1Col, lev2Col}{
Colours to be used for background points, and for the two stringency levels defined by \code{plevel1} and \code{plevel2}, respectively.
}
  \item{bgPch, lev1Pch, lev2Pch}{
Plotting character for background points, and for points exceeding the two stringency levels defined by \code{plevel1} and \code{plevel2}, respectively. Specified as per \code{pch} in \code{\link{points}}.
}
  \item{\dots}{
Additional arguments passed to \code{\link{plot}}
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
Vector of the baitIDs plotted (useful if baitIDs were sampled randomly).
}
\author{
Mikhail Spivakov, Jonathan Cairns, Paula Freire Pritchett
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
data(cdUnitTest)

##modifications to cdUnitTest, ensuring it uses correct design directory
designDir <- file.path(system.file("extdata", package="Chicago"), "unitTestDesign")
new.settings <- paste0(file.path(designDir, "unitTest"), c(".rmap", ".baitmap", ".npb", ".nbpb", ".poe"))
new.settings <- as.list(new.settings)
names(new.settings) <- c("rmapfile", "baitmapfile", "nperbinfile", "nbaitsperbinfile", "proxOEfile")
cdUnitTest <- modifySettings(cd = cdUnitTest, new.settings)

plotBaits(cdUnitTest)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
