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
plotBaits(cd, pcol = "score", Ncol = "N", n = 16, baits = NULL, plotBaitNames = TRUE, plotBprof = FALSE, plevel1 = 5, plevel2 = 3, outfile = NULL, removeBait2bait = TRUE, width = 20, height = 20, maxD = NULL, bgCol = "black", lev2Col = "blue", lev1Col = "red", bgPch = 1, lev1Pch = 20, lev2Pch = 20, ...)
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
%%     ~~Describe \code{plevel1} here~~
}
%  \item{plevel2}{
%%     ~~Describe \code{plevel2} here~~
%}
  \item{outfile}{
%%     ~~Describe \code{outfile} here~~
}
  \item{removeBait2bait}{
%%     ~~Describe \code{removeBait2bait} here~~
}
  \item{width}{
%%     ~~Describe \code{width} here~~
}
  \item{height}{
%%     ~~Describe \code{height} here~~
}
  \item{maxD}{
%%     ~~Describe \code{maxD} here~~
}
  \item{bgCol}{
%%     ~~Describe \code{bgCol} here~~
}
  \item{lev2Col}{
%%     ~~Describe \code{lev2Col} here~~
}
  \item{lev1Col}{
%%     ~~Describe \code{lev1Col} here~~
}
  \item{bgPch}{
%%     ~~Describe \code{bgPch} here~~
}
  \item{lev1Pch}{
%%     ~~Describe \code{lev1Pch} here~~
}
  \item{lev2Pch}{
%%     ~~Describe \code{lev2Pch} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
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
function (cd, pcol = "score", Ncol = "N", n = 16, baits = NULL, 
    plotBaitNames = TRUE, plotBprof = FALSE, plevel1 = 5, plevel2 = 3, 
    outfile = NULL, removeBait2bait = TRUE, width = 20, height = 20, 
    maxD = NULL, bgCol = "black", lev2Col = "blue", lev1Col = "red", 
    bgPch = 1, lev1Pch = 20, lev2Pch = 20, ...) 
{
    if (plotBaitNames) {
        baitmap = fread(cd@settings$baitmapfile)
    }
    if (is.null(baits)) {
        baits = sample(unique(cd@x$baitID), n)
    }
    else {
        n = length(baits)
    }
    if (plotBprof) {
        disp = cd@params$dispersion
    }
    if (!is.null(outfile)) {
        pdf(outfile, width = width, height = height)
    }
    if (n >= 4) {
        par(mfrow = c(4, ceiling(n/4)))
    }
    else {
        par(mfrow = c(n, 1))
    }
    setkey(x@x, baitID)
    for (i in 1:n) {
        this = cd@x[baitID == baits[i]]
        this = this[is.na(distSign) == FALSE]
        if (!is.null(maxD)) {
            this = this[abs(distSign) <= maxD]
        }
        if (removeBait2bait) {
            this = this[isBait2bait == FALSE]
        }
        setDF(this)
        this = this[order(this$distSign), ]
        cols <- rep(bgCol, nrow(this))
        pchs <- rep(bgPch, nrow(this))
        sel1 <- this[, pcol] >= plevel1
        sel2 <- this[, pcol] >= plevel2
        cols[sel2] <- lev2Col
        cols[sel1] <- lev1Col
        pchs[sel2] <- lev2Pch
        pchs[sel1] <- lev2Pch
        title = paste(baits[i], sep = "")
        if (plotBaitNames) {
            baitName = baitmap[baitmap$V4 == baits[i]][, cd@settings$baitmapGeneIDcol, 
                with = F]
            if (length(grep(",", baitName))) {
                baitName = gsub("(\\S+,).+", "\\1", baitName)
                baitName = paste0(baitName, "...")
            }
            title = paste0(baitName, " (", title, ")")
        }
        plot(this$distSign, this[, Ncol], xlab = distcol, ylab = Ncol, 
            main = title, col = cols, pch = pchs, ...)
        abline(v = 0, col = "grey", lwd = 1)
        if (plotBprof) {
            lines(this[, distcol], this$Bmean, lwd = 1, col = "darkgrey")
            lines(this[, distcol], this$Bmean + 1.96 * sqrt(this$Bmean + 
                this$Bmean^2/disp), lwd = 1, lty = 2, col = "darkgrey")
        }
    }
    if (!is.null(outfile)) {
        dev.off()
    }
    baits
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
