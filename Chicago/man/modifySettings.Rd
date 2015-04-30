\name{modifySettings}
\alias{modifySettings}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Modify Settings
}
\description{
Modify the settings in a \code{chicagoData} object.
}
\usage{
modifySettings(cd, settings)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cd}{A \code{chicagoData} object.}
  \item{settings}{A named list.}
}
\details{
\code{cd}'s settings are updated. Setting \code{names(settings)[1]} is set to (settings)[[1]], and so on. For a list of available settings, see \code{\link{defaultSettings}}.
}
\value{
An object of class \code{chicagoData}.
}
\author{
Mikhail Spivakov, Jonathan Cairns, Paula Freire Pritchett
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{defaultSettings}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (cd, settings) 
{
    message("Warning: settings are not checked for consistency with the previous ones.")
    for (s in names(settings)) {
        cd@settings[[s]] = settings[[s]]
    }
    cd
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
