\name{getPvals}
\alias{getPvals}

\title{
Get P-values
}
\description{
Based on a Delaporte model, calculate the P-value associated with each observation.
}
\usage{
getPvals(cd)
}
\arguments{
  \item{cd}{A \code{chicagoData} object.}
}
\section{Warning}{
The object \code{intData(cd)} is updated by reference. Thus, \code{intData(cd)} will be altered. See vignette for further information.
}
\details{
The parameters for the Delaporte distribution are obtained as follows: the NB mean from the column \code{intData(cd)$Bmean}, the Poisson mean from the column \code{intData(cd)$Tmean}, and the dispersion from \code{params(cd)$dispersion}.
}
\value{
An object of class \code{chicagoData}, with new column \code{log.p}.
}
\author{
Mikhail Spivakov, Jonathan Cairns, Paula Freire Pritchett
}

\seealso{
\code{\link{chicagoPipeline}}
}
\examples{
data(cdUnitTest)
cdUnitTest <- getPvals(cdUnitTest)
}
