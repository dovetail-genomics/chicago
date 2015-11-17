\docType{data}
\name{sGM12878}
\alias{sGM12878}
\title{Promoter Capture data for chromosomes 20 and 21 in GM12878}
\format{A ChicagoData object.}
\source{
Mifsud, B. et al. "Mapping long-range promoter contacts in human cells with high-resolution capture Hi-C." Nature Genetics (2015) doi:10.1038/ng.3286
}
\usage{
sGM12878
}
\description{
Promoter Capture data to be used as a toy example to run all steps of Chicago. This data only incorporates read pairs
including both chromosomes 20 and 21 in order to minimize processing time and memory storage. Thus, this data set
includes all cis read pairs for these two chromosomes and all trans read pairs between them.

The hg19 genome was used.

The package includes ChIP-seq data from the ENCODE consortium, also restricted to chr20 and chr21. (The ENCODE Project Consortium, 2012 "An Integrated Encyclopedia of DNA Elements in the Human Genome." Nature 489 (September): 57-74. Data accession numbers: Bernstein lab GSM733752, GSM733772, GSM733708, GSM733664, GSM733771, GSM733758)
}
\value{
A \code{\link{chicagoData}} object.
}
\examples{
data(sGM12878)
##modifications to sGM12878, ensuring it uses correct design directory
designDir <- file.path(system.file("extdata", package="PCHiCdata"), "hg19TestDesign")
sGM12878 <- modifySettings(cd = sGM12878, designDir=designDir)
}
\seealso{
\code{\link{chicagoData}}
}
\keyword{datasets}

