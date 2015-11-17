\name{mergeSamples}
\alias{mergeSamples}
\title{
Merge samples together.
}
\description{
Merge a number of \code{chicagoData} objects together, summarising their counts into a normalised value.
}
\usage{
mergeSamples(cdl, normalise = TRUE, NcolOut = "N",
      NcolNormPrefix = "NNorm", mergeMethod = c("weightedMean", "mean")[1], repNormCounts = (mergeMethod=="mean"))
}
\arguments{
  \item{cdl}{A \code{list} of \code{chicagoData} objects.}
  \item{normalise}{
If \code{TRUE}, use a normalisation procedure, specified by \code{mergeMethod}, to arrive at a normalised count. If \code{FALSE}, take the mean number of reads.
}
  \item{NcolOut}{
The column to store the normalised counts in.
}
  \item{NcolNormPrefix}{
Each sample gains a normalised count column, that begins with this prefix.
}
  \item{mergeMethod}{
If \code{mergeMethod == "weightedMean"}, then NcolOut is the weighted mean of the sample-wise counts adjusted by the samples' respective scaling factors \code{s_k}. If \code{mergeMethod == "mean"}, then sample-specific counts are first normalised by dividing by s_k, and NcolOut is computed as the mean of these normalised counts.
}

  \item{repNormCounts}{
Report normalised counts for each replicate (by dividing them by s_k) in the <NcolNormPrefix>.<sampleNo> column (by default, NNorm.1, NNorm.2, etc.). This option is on by default when \code{mergeMethod == "mean"}. However, it can also be used with \code{mergeMethod == "weightedMean"} (but the normalised counts will still be produced by dividing the raw counts for each replicate by s_k).
}

}
\note{
Raw per-sample counts will be stored in the N.<sampleNo> column (N.1, N.2, etc.)
}
\value{
An object of class \code{chicagoData}, with a \code{params(cd)$s_k} slot added representing the per-sample scaling factors used in normalisation.
}
\author{
Mikhail Spivakov, Jonathan Cairns, Paula Freire Pritchett
}

\seealso{
\code{\link{readAndMerge}}
}
\examples{
filesDir <- file.path(system.file("extdata", package="Chicago"), "unitTestData")
files <- file.path(filesDir, dir(filesDir))
print(files) ##we will read in and merge these files

designDir <- file.path(system.file("extdata", package="Chicago"), "unitTestDesign")

cdA <- setExperiment(designDir=designDir)
cdA <- readSample(files[1], cdA)

cdB <- setExperiment(designDir=designDir)
cdB <- readSample(files[2], cdB)

cdMerged <- mergeSamples(list(cdA, cdB))
}
