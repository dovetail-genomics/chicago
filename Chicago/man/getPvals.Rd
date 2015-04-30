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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (cd) 
{
    alpha = cd@params$dispersion
    x = cd@x
    if (is.null(alpha)) {
        stop("getPvals: 'dispersion' parameter of x not found.")
    }
    message("Calculating p-values...")
    x[, `:=`(log.p, ifelse(Bmean < .Machine$double.eps, ppois(N - 
        1L, lambda = Tmean, lower.tail = FALSE, log.p = TRUE), 
        pdelap(N - 1L, alpha, beta = Bmean/alpha, lambda = Tmean, 
            lower.tail = FALSE, log.p = TRUE)))]
    sel <- which(is.infinite(x$log.p) | is.nan(x$log.p))
    if (length(sel) > 0) {
        message("Approximating ", length(sel), " very small p-values.")
        gamma <- x[sel, alpha * (1 + Tmean/Bmean)^2]
        gamma <- pmin(gamma, 1e+10)
        x[sel, `:=`(log.p, pnbinom(N - 1L, size = gamma, mu = Bmean + 
            Tmean, lower.tail = FALSE, log.p = TRUE))]
        if (any(is.infinite(x[sel, log.p]))) {
            warning("Some log-p-values were infinite.")
        }
        if (any(is.nan(x[sel, log.p]))) {
            warning("Some log-p-values were NaNs.")
        }
    }
    if (any(is.na(x$log.p))) {
        warning("Some log-p-values were NA.")
    }
    cd
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
