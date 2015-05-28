\name{overlapFragWithFeatures}
\alias{overlapFragWithFeatures}
\title{Overlap Other-Ends with Features}
\usage{
overlapFragWithFeatures(x = NULL, folder = NULL, list_frag, sep = "\\t",
  position_otherEnd = NULL)
}
\arguments{
\item{x}{a \code{\link{chicagoData}} object or a data table (\code{\link{data.table}}) containining other end IDs.}

\item{folder}{the name of the folder where the files containing the features of interest are stored.}

\item{list_frag}{a list where each element is the name of a file containing a feature of interest (e. g. H3K4me1, CTCF, DHS etc.).
These files must have a bed format, with no header.
Each element of the list must be named.}

\item{sep}{the field separator character. Values are separated by this character on each line of the file containing
the coordinates of the restriction fragments (called by \code{position_otherEnd}).}

\item{position_otherEnd}{the name of the file containing the coordinates of the restriction fragments
and the corresponding IDs. The coordinates should be "chromosome", "start" and "end", and the ID should be numeric.
\code{position_otherEnd} only needs to be specified if x is not a \code{\link{chicagoData}} object.}
}
\value{
a data table (\code{\link{data.table}}) built from \code{x}, where a column was added for each genomic feature present in \code{list_frag}.
The new columns contain logical values indicating whether there was an overlap or not.
}
\description{
This function checks which other-ends from a \code{\link{chicagoData}} object overlap with a set of genomic features.
}
\examples{
FIXME <- function(FIXME)
}
\author{
Mikhail Spivakov, Jonathan Cairns, Paula Freire Pritchett
}

