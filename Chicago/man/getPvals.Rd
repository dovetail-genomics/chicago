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
The object \code{cd@x} is updated by reference. Thus, \code{cd@x} will be altered. See vignette for further information.
}
\details{
The parameters for the Delaporte distribution are obtained as follows: the NB mean from the column \code{cd@x$Bmean}, the Poisson mean from the column \code{cd@x$Tmean}, and the dispersion from \code{cd@params$dispersion}.
}
\value{
An object of class \code{chicagoData}, with new column \code{log.p}.
}
\author{
Mikhail Spivakov, Jonathan Cairns, Paula Freire Pritchett
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{chicagoPipeline}}
}
\examples{
data(cdUnitTest)
cdUnitTest <- getPvals(cdUnitTest)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
