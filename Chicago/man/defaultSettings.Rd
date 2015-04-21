\name{defaultSettings}
\alias{defaultSettings}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Default CHiCAGO settings
}
\description{
A function that gives the default settings used for a CHiCAGO experiment.
}
\usage{
defaultSettings(designDir = "")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{designDir}{A directory containing the digest information (see below).}
}
\details{

}
\value{
A list of the following settings:

\item{rmapfile}{Default: file.path(designDir, "Digest_Human_HindIII.bed")}
\item{baitmapfile}{Default: file.path(designDir, "Digest_Human_HindIII_baits_ID.bed")}
\item{nperbinfile}{Default: file.path(designDir, "Digest_Human_HindIII_NperBin.txt")}
\item{nbaitsperbinfile}{Default: file.path(designDir, "Digest_Human_HindIII_NbaitsPerBin.txt")}
\item{proxOEfile}{Default: file.path(designDir, "proxOE_out.txt")}
\item{Ncol}{Default: "N"}
\item{baitmapFragIDcol}{Default: 4}
\item{baitmapGeneIDcol}{Default: 5}
\item{maxLBrownEst}{Default: 1500000}
\item{minFragLen}{Default: 150}
\item{maxFragLen}{Default: 40000}
\item{minNPerBait}{Default: 250}
\item{binsize}{Default: 20000}
\item{removeAdjacent}{Default: TRUE}
\item{adjBait2bait}{Default: TRUE}
\item{tlb.filterTopPercent}{Default: 0.01}
\item{tlb.minProxOEPerBin}{Default: 1000}
\item{tlb.minProxB2BPerBin}{Default: 100}
\item{techNoise.minBaitsPerBin}{Default: 1000}
\item{brownianNoise.subset}{Default: 1000}
\item{brownianNoise.seed}{Default: NULL}
\item{baitIDcol}{Default: "baitID"}
\item{otherEndIDcol}{Default: "otherEndID"}
\item{otherEndLencol}{Default: "otherEndLen"}
\item{distcol}{Default: "distSign"}
\item{weightAlpha}{Default: 34.1157346557331}
\item{weightBeta}{Default: -2.58688050486759}
\item{weightGamma}{Default: -17.1347845819659}
\item{weightDelta}{Default: -7.07609245521541}
}
\author{
Mikhail Spivakov, Jonathan Cairns, Paula Freire Pritchett
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{setExperiment}}, \code{\link{modifySettings}}
}
\examples{
s <- defaultSettings()
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
