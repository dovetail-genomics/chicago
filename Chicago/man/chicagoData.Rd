\name{chicagoData}
\alias{chicagoData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The \code{chicagoData} class.
}
\description{
Constructor for the \code{chicagoData} class.
}
\usage{
chicagoData(...)
}
\arguments{
  \item{\dots}{
Arguments passed to \code{new()}.
}
}
\details{
While this function can be used to create a \code{chicagoData} object, most users will use the \code{\link{setExperiment}} function instead.
}
\value{
A \code{chicagoData} object has three slots:

* \code{cd@x} is a data.table (note: not a data.frame) that contains information about fragment pairs.
* \code{cd@settings} is a list of settings, usually set with the setExperiment() function. For more information about valid settings, please see \code{\link{defaultSettings}}.
* \code{cd@params} is a list of parameters. CHiCAGO estimates these automatically, as part of the pipeline.

}
\author{
Mikhail Spivakov, Jonathan Cairns, Paula Freire Pritchett
}

\seealso{
\code{\link{setExperiment}}, \code{\link{defaultSettings}}
}
\examples{
  cd <- chicagoData()
}
