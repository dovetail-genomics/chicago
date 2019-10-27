\name{defaultSettings}
\alias{defaultSettings}
\title{
Default CHiCAGO settings
}
\description{
A function that gives the default settings used for a CHiCAGO experiment.

IMPORTANT: from version 1.13, the following parameters are set based on the values in .npb file header and checked for consistency with the headers of .npbp and .poe files and custom-defined settings. They should therefore be provided to the makeDesignFiles.py script, which needs to be rerun if they need to be modified:

rmapfile (only the basename is checked; inconsistent baitmapfile will only generate a warning for compatibility with publicly released older designs), 
minFragLen, maxFragLen, binsize, removeAdjacent, adjBait2bait.

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
\item{maxLBrownEst}{Default: 1500000. The distance range to be used for estimating the Brownian component of the null model
The parameter setting should approximately reflect the maximum distance, at which the power-law distance dependence is still observable. }

\item{minFragLen}{Default: 150. (See maxFragLen.)}
\item{maxFragLen}{Default: 40000. minFragLen and maxFragLen correspond to the limits within which we observed no clear dependence between fragment length and the numbers of reads mapping to these fragments in HindIII PCHiC data.

These parameters need to be modified when using a restriction enzyme with a different cutting frequency (such as a 4-cutter) and can also be verified by users with their datasets in each individual case. However, we note that the fragment-level scaling factors (s_i and s_j) generally incorporate the effects of fragment size, so this filtering step only aims to remove the strongest bias.  }
\item{minNPerBait}{Default: 250. Minimum number of reads that a bait has to accumulate to be included in the analysis.

Reasonable numbers of per-bait reads are required for robust parameter estimation. If this value is too low, the confidence of interaction calling is reduced. If too high, too many baits may be unreasonably excluded from the analysis. If it is desirable to include baits below this threshold, we recommend decreasing this parameter and then visually examining the result bait profiles (for example, using plotBaits()).}
\item{binsize}{Default: 20000. The bin size (in bases) used when estimating the Brownian collision parameters.

The bin size should, on average, include several (~4-5) restriction fragments to increase the robustness of parameter estimation. However, using too large bins will reduce the precision of distance function estimation. Therefore, this value needs to be changed if using an enzyme with a different cutting frequency (such as a 4-cutter).}
\item{removeAdjacent}{Default: TRUE. Should fragments adjacent to baits be removed from analysis?

We remove fragments adjacent to baits by default, as the corresponding ligation products are indistinguishable from incomplete digestion. This setting however may be set to FALSE if the rmap and baitmap files represent bins over multiple fragments as opposed to fragment-level data (e.g., to address sparsity issues with low-coverage experiments).}
\item{adjBait2bait}{Default: TRUE. Should baited fragments be treated separately?
Baited fragments are treated separately from the rest in estimating other end-level scaling factors (si) and technical noise levels. It is a free parameter mainly for development purposes, and we do not recommend changing it.}

\item{tlb.filterTopPercent}{Default: 0.01. Top percent of fragments with respect to accumulated trans-counts to be filtered out in the binning procedure.

Other ends are pooled together when calculating their scaling factors and as part of technical noise estimation. Binning is performed by quantile, and for the most extreme outliers this approach is not going to be adequate. Increasing this value may potentially make the estimation for the highest-count bin more robust, but will exclude additional other ends from the analysis.}
\item{tlb.minProxOEPerBin}{Default: 50000. Minimum pool size (i.e. minimum number of other ends per pool), used when pooling other ends together based on trans-counts.

If this parameter is set too small, then estimates will be imprecise due to sparsity issues. If this parameter is set too large, then the model becomes inflexible and so the model fit is hindered. This parameter could be decreased in a dataset that has been sequenced to an extremely high depth. Alternatively, it may need to be decreased out of necessity, in a dataset with very few other ends - for example, the vignette decreases this setting to process the PCHiCdata package data (since these data sets span only a small subset of the genome, in each case).}
\item{tlb.minProxB2BPerBin}{Default: 2500. Minimum pool size, used when pooling other ends together (bait-to-bait interactions only). (See previous entry, tlb.minProxOEPerBin, for advice on setting parameter.)}
\item{techNoise.minBaitsPerBin}{Default: 1000. Minimum pool size, used when pooling baits together based on accumulated trans-counts. (See tlb.minProxOEPerBin for advice on setting parameter.)}
\item{brownianNoise.samples}{Default: 5. Number of times subsampling occurs when estimating the Brownian collision dispersion.

Dispersion estimation from a subset of baits has an error attached. Averaging over multiple subsamples allows us to decrease this error. Increasing this number improves the precision of dispersion estimation at the expense of greater runtime.}
\item{brownianNoise.subset}{Default: 1000. Number of baits sampled from when estimating the Brownian noise dispersion. If set to NA, then all baits are used.

Estimating dispersion from the entire dataset usually requires a prohibitively large amount of memory. A subset is chosen that is large enough to get a reasonably precise estimate of the dispersion, but small enough to stay in memory. A user with excess memory may wish to increase this number to further improve the estimate's precision.}
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
