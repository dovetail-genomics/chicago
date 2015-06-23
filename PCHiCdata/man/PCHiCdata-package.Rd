\name{PCHiCdata-package}
\alias{PCHiCdata-package}
\alias{PCHiCdata}
\docType{package}
\title{
Promoter Capture Hi-C data
}
\description{
Promoter Capture data Hi-C conveniently packaged for consumption by R users. Data includes interactions detected for chromosomes 20 and 21 in GM12878 cells and for chromosomes 18 and 19 in mESC.
}
\details{
\tabular{ll}{
Package: \tab PCHiCdata\cr
Type: \tab Package\cr
Version: \tab 0.1\cr
Date: \tab 2015-05-22\cr
License: \tab Artistic-2.0\cr
Depends: \tab R (>= 3.1.2), Chicago\cr
LazyData: \tab true\cr
}

Index:
\preformatted{
sGM12878                Promoter Capture data for chromosomes 20 and 21
                        in GM12878
smESC                   Promoter Capture data for chromosomes 18 and 19
                        in smESC
}
}
\author{
Mikhail Spivakov, Jonathan Cairns, Paula Freire Pritchett, Steven Wingett

Maintainer: Paula Freire-Pritchett <paula.freire-pritchett@babraham.ac.uk>
}
\references{
Mifsud, B. et al. "Mapping long-range promoter contacts in human cells with high-resolution capture Hi-C." Nature Genetics (2015) doi:10.1038/ng.3286

Schoenfelder, S. et al. "The pluripotent regulatory circuitry connecting promoters to their long-range interacting elements." Genome research 25.4 (2015): 582-597.
}
\examples{
data(sGM12878)
data(smESC)
}
