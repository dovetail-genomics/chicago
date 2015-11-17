\name{exportResults}
\alias{exportResults}
\alias{exportToGI}

\title{
Export Results
}
\description{
Export the results from a \code{chicagoData} object to disk, or to a \code{\link[GenomicInteractions:GenomicInteractions]{GenomicInteractions}} object.
}
\usage{
exportResults(cd, outfileprefix, scoreCol = "score", cutoff = 5, b2bcutoff = NULL,
        format = c("seqMonk", "interBed", "washU_text"),
        order = c("position", "score")[1], removeMT=TRUE)
exportToGI(cd, scoreCol="score", cutoff=5, b2bcutoff=NULL,
                       order=c("position", "score")[1], removeMT=TRUE)
}
\arguments{
  \item{cd}{A \code{chicagoData} object.}
  \item{outfileprefix}{
  A character string that forms the prefix for each output file.
}
  \item{scoreCol}{
  The column of \code{intData(cd)} that contains the score.
}
  \item{cutoff}{
  The score cutoff.
}
  \item{b2bcutoff}{
  If desired, an alternative score cutoff for bait-to-bait interactions.
}
  \item{format}{
  The file format(s) to output. If a multiple formats are supplied as a vector, then all of these formats will be outputted.

  Supported formats are: "seqMonk", "interBed", "washU_text" and, for advanced users, "washU_track".
}
  \item{order}{
  Should output be ordered by position or score?
}
  \item{removeMT}{
  Logical. If \code{TRUE}, remove any interactions involving mitochondrial DNA from the output.
}
}
\details{
Important notes on the washU formats: Most users will prefer "washU_text" output to "washU_track" output. The "washU_text" output can be uploaded to the washU browser directly. To do this, open the browser, select "Add custom tracks", and use the "Got text files instead? Upload them from your computer" link near the bottom of the page.

The "washU_track" output needs to be hosted elsewhere. You can then link the browser to the data via the "Interaction - pairwise interaction" button on the "Add custom tracks" page.

If you get the warning "WashU Browser track format could not be finalized due to absence of bgzip or tabix", this could be because you have not installed SAMtools and htslib. You can check with system2("tabix") and system2("bgzip"). Sometimes RStudio has issues with reading \code{$PATH} - you can check this with \code{system2("echo", "$PATH")}. Consider running the command in R, outside of RStudio, to fix this problem.

If all else fails, and you need "washU_track" output, then you can manually perform the final steps yourself by running: \code{bgzip <outfileprefix>_washU_track.txt} and \code{tabix -p bed <outfileprefix>.txt.gz}.
}
\value{
\code{exportResults()}: NULL.

\code{exportToGI()}: a \code{\link[GenomicInteractions:GenomicInteractions]{GenomicInteractions}} object. Anchor one is the bait, anchor two is the other end.
}
\author{
Mikhail Spivakov, Jonathan Cairns, Paula Freire Pritchett
}

\seealso{
\code{\link{chicagoPipeline}}
}
\examples{
data(cdUnitTest)

##modifications to cdUnitTest, ensuring it uses correct design directory
designDir <- file.path(system.file("extdata", package="Chicago"), "unitTestDesign")
cdUnitTest <- modifySettings(cd=cdUnitTest, designDir=designDir)

##create a temporary directory, export output there
tempDirectory <- tempdir()
print(tempDirectory)
exportResults(cdUnitTest, outfileprefix = file.path(tempDirectory, "unitTestOutput"))

GI <- exportToGI(cdUnitTest)
}