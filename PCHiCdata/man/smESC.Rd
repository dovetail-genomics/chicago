\docType{data}
\name{smESC}
\alias{smESC}
\title{Promoter Capture data for chromosomes 18 and 19 in smESC}
\format{A ChicagoData object.}
\source{
Schoenfelder, S. et al. "The pluripotent regulatory circuitry connecting promoters to their long-range interacting elements." Genome research 25.4 (2015): 582-597.
}
\usage{
smESC
}
\description{
Promoter Capture data to be used as a toy example to run all steps of Chicago. This data only incorporates read pairs
including both chromosomes 18 and 19 in order to minimize processing time and memory storage. Thus, this data set
includes all cis read pairs for these two chromosomes and all trans read pairs between them.

The mm9 genome was used.

The package includes ChIP-seq data from the ENCODE consortium, also restricted to chr18 and chr19. (The ENCODE Project Consortium, 2012 "An Integrated Encyclopedia of DNA Elements in the Human Genome." Nature 489 (September): 57-74. Data accession numbers: Snyder lab GSM1003750, GSM1003756, GSM1003751; Ren lab GSM1000126, Hardison lab ENCFF001ZJO, ENCFF001ZQY)
}
\value{
A \code{\link{chicagoData}} object.
}
\examples{
data(smESC)
##modifications to smESC, ensuring it uses correct design directory
designDir <- file.path(system.file("extdata", package="PCHiCdata"), "mm9TestDesign")
smESC <- modifySettings(cd = smESC, designDir=designDir)
}
\seealso{
\code{\link{chicagoData}}
}
\keyword{datasets}

