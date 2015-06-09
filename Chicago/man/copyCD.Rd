\name{copyCD}
\alias{copyCD}

\title{
Copy chicagoData object
}
\description{
Copies a chicagoData object. (Failing to use this function may mean that an object is updated by reference when its 'copy' is altered.)
}
\usage{
copyCD(cd)
}

\arguments{
  \item{cd}{
\code{\link{chicagoData}} object.
}
}
\value{
\code{\link{chicagoData}} object.
}
\author{
Jonathan Cairns
}

\examples{
data(cdUnitTest)
x <- copyCD(cdUnitTest)
}

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
