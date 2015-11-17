\name{defaultSettings}
\alias{defaultSettings}
\title{
Default CHiCAGO settings
}
\description{
A function that gives the default settings used for a CHiCAGO experiment.
}
\usage{
defaultSettings()
}
\value{
A list of the following settings:

\item{rmapfile}{Default: NA. The location of the restriction map file; see the vignette for a description of what this file should contain.}
\item{baitmapfile}{Default: NA. The location of the bait map file; see the vignette for a description of what this file should contain.}
\item{nperbinfile}{Default: NA. See vignette.}
\item{nbaitsperbinfile}{Default: NA. See vignette.}
\item{proxOEfile}{Default: NA. See vignette.}
\item{Ncol}{Default: "N". The column in intData(cd) that contains the number of reads.}
\item{baitmapFragIDcol}{Default: 4. In the bait map file, the number of the column that specifies the fragment ID of each bait.}
\item{baitmapGeneIDcol}{Default: 5. In the bait map file, the number of the column that specifies which gene(s) are on each fragment.}
\item{maxLBrownEst}{Default: 1500000. The maximum distance away from the bait for which Brownian noise is assumed to dominate.}
\item{minFragLen}{Default: 150. Fragments with size less than this are discarded.}
\item{maxFragLen}{Default: 40000. Fragments with size greater than this are discarded.}
\item{minNPerBait}{Default: 250. Baits must have at least this many associated reads. Otherwise, all such reads are discarded (except those that are bait-to-bait)}
\item{binsize}{Default: 20000. Bins of this size are used when estimating Brownian noise}
\item{removeAdjacent}{Default: TRUE}
\item{adjBait2bait}{Default: TRUE}
\item{tlb.filterTopPercent}{Default: 0.01}
\item{tlb.minProxOEPerBin}{Default: 1000}
\item{tlb.minProxB2BPerBin}{Default: 100}
\item{techNoise.minBaitsPerBin}{Default: 1000}
\item{brownianNoise.subset}{Default: 1000. Number of baits sampled from when estimating the Brownian noise dispersion. If set to NA, then all baits are used.}
\item{brownianNoise.seed}{Default: NA. If not NA, then \code{brownianNoise.seed} is used as the random number generator seed when subsampling baits. Set this to make your analysis reproducible.}
\item{baitIDcol}{Default: "baitID". The name of the baitID column in intData(cd).}
\item{otherEndIDcol}{Default: "otherEndID". The name of the otherEndID column in intData(cd).}
\item{otherEndLencol}{Default: "otherEndLen". The name of the column in intData(cd) that contains the lengths of the other end fragments.}
\item{distcol}{Default: "distSign". The name of the column in intData(cd) that contains the genomic distance that an interaction spans.}
\item{weightAlpha}{Default: 34.1157346557331. This, and the following parameters, are used in the p-value weighting procedure.}
\item{weightBeta}{Default: -2.58688050486759}
\item{weightGamma}{Default: -17.1347845819659}
\item{weightDelta}{Default: -7.07609245521541}
}
\author{
Mikhail Spivakov, Jonathan Cairns, Paula Freire Pritchett
}

\seealso{
\code{\link{setExperiment}}, \code{\link{modifySettings}}
}
\examples{
s <- defaultSettings()
print(s)
}
