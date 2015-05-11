\name{normaliseBaits}
\alias{normaliseBaits}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Normalise Baits
}
\description{
Calculate normalisation factors \code{s_j} for each bait.
}
\usage{
normaliseBaits(cd, normNcol = "NNb", shrink = FALSE,
      plot = TRUE, outfile = NULL, debug = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cd}{A \code{chicagoData} object.}
  \item{normNcol}{
The name of the column in \code{cd} that contains normalised counts.
}
  \item{shrink}{
Deprecated.
}
  \item{plot}{
If \code{TRUE}, output a diagnostic plot.
}
  \item{outfile}{\code{NULL}, or a character string. If \code{NULL}, the diagnostic plot is outputted to the current plotting device. If a character string, e.g. \code{outfile="tech.pdf"}, then the plot will be outputted to that file.}
  \item{debug}{
Deprecated.
}
}
\section{Warning}{
The object \code{cd@x} is updated by reference. Thus, \code{cd@x} will be altered. See vignette for further information.
}
\details{
A \code{chicagoData} object: \code{cd@x} gains a new column \code{s_j}, and normalised output \code{NNb} (unless the \code{normNcol} parameter is altered).
}
\value{
An object of class \code{chicagoData}.
}
\author{
Mikhail Spivakov, Jonathan Cairns, Paula Freire Pritchett
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
data(cdUnitTest)

##modifications to cdUnitTest, ensuring it uses correct design directory
designDir <- file.path(system.file("extdata", package="Chicago"), "unitTestDesign")
new.settings <- paste0(file.path(designDir, "unitTest"), c(".rmap", ".baitmap", ".npb", ".nbpb", ".poe"))
new.settings <- as.list(new.settings)
names(new.settings) <- c("rmapfile", "baitmapfile", "nperbinfile", "nbaitsperbinfile", "proxOEfile")
cdUnitTest <- modifySettings(cd = cdUnitTest, new.settings)

cdUnitTest <- normaliseBaits(cdUnitTest)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
