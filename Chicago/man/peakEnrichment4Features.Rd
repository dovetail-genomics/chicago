\name{peakEnrichment4Features}
\alias{peakEnrichment4Features}
\title{Enrichment for Features}
\usage{
peakEnrichment4Features(x1 = NULL, score = 5, colname_score = "score",
  min_dist = 0, max_dist = NULL, folder = NULL, list_frag,
  sep = "\\t", filterB2B = TRUE, b2bcol = "isBait2bait", unique = TRUE,
  no_bins, sample_number, plot_name = NULL, position_otherEnd = NULL,
  colname_dist = NULL,trans=FALSE, plotPeakDensity=FALSE)
}
\arguments{
\item{x1}{a \code{\link{chicagoData}} object or a data table (\code{\link{data.table}}) containining other end IDs.}

\item{score}{the threshold above which interactions start being called as significant.}

\item{colname_score}{the name of the column which contains the score values which establish the level of significance of each interaction.}

\item{min_dist}{the minimum distance from bait required in the query. If this parameter is
set to \code{NULL} and \code{trans} is set to \code{TRUE}, cis interactions are disregarded from the analysis.
This parameter is also useful when the user only wants to look at cis distal interactions (very far from bait).}

\item{max_dist}{the maximum distance from bait required in the query.
This parameter is particularly useful when the user only wants to look at cis proximal interactions (interactions surrounding the bait).}

\item{folder}{the name of the folder where the files containing the features of interest are stored.}

\item{list_frag}{a list where each element is the name of a file containing a feature of interest (e. g. H3K4me1, CTCF, DHS etc.).
These files must have a bed format, with no header.
Each element of the list must be named.}

\item{sep}{the field separator character. Values are separated by this character on each line of the file containing
the coordinates of the restriction fragments (called by \code{position_otherEnd}).}

\item{filterB2B}{a logical value indicating whether bait-to-bait interactions should be removed from the analysis.}

\item{b2bcol}{the name of the column identifying bait-to-bait interactions in the \code{x1}.}

\item{unique}{a logical value indicating whether to removing duplicated other-ends from significant interactions and samples.}

\item{no_bins}{Number of bins to divide the range of \code{colname_dist} (after \code{colname_dist} has been trimmed
according to \code{min_dist} and \code{max_dist}). This will be important to determine how many interactions should be
sampled according to distance from bait.}

\item{sample_number}{Number of samples to be used in the permutation test.}

\item{plot_name}{the name of the file where to save the resulting plot. This parameter is only required if the user wants to save the plot.
Otherwise, the plot will be displayed on the screen, but not saved.}

\item{position_otherEnd}{the name of the file containing the coordinates of the restriction fragments
and the corresponding IDs. The coordinates should be "chromosome", "start" and "end", and the ID should be numeric.
\code{position_otherEnd} only needs to be specified if x is not a \code{\link{chicagoData}} object.}

\item{colname_dist}{the name of the column which contains the distances between bait and other end.
\code{colname_dist} only needs to be specified if x is not a \code{\link{chicagoData}} object.}

\item{trans}{a logical value indicating whether the enrichment is to be computed for trans interactions. 
If this parameter is set to \code{TRUE} and \code{min_dist} is set to \code{NULL}, cis interactions are disregarded from the analysis.}

\item{plotPeakDensity}{a logical value indicating whether to plot the density of interactions with distance. Setting this parameter to \code{TRUE} only applies to cis interactions.}

}
\value{
a data frame containing columns for the number of overlaps for each feature in our significant interactions,
the average number of overlaps for each feature in our samples, the corresponding standard deviations.
}
\description{
This function computes how many other-ends from a \code{\link{chicagoData}} object, that engage in significant interactions,
overlap with a set of genomic features.
In order to determine how those numbers compare to what would be expected if interaction significance had no effect on the overlaps,
this function samples different sets of interactions from the non-significant pool and assesses how they overlap with genomic features
(it computes the mean and confidence intervals).
Results are returned in a table and plotted in a barplot.
The difference between the results for the set of significant interactions and the random samples can be used as a measure of the
enrichment for genomic features.
Samples have the same size as the number of significant interactions called. Moreover, they follow the same distribution of distances between
bait and other-end. This is achieved by binning this distribution and drawing interactions per bin, according to the numbers observed in
the significant set.
}
\author{
Mikhail Spivakov, Jonathan Cairns, Paula Freire Pritchett
}
\examples{
FIXME <- function(FIXME)
}
