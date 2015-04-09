\name{estimateDistFun}
\alias{estimateDistFun}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
estimateDistFun(cd, method = "cubic", n.obs.head = 10, n.obs.tail = 25, logScale = FALSE, outfile = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cd}{
%%     ~~Describe \code{cd} here~~
}
  \item{method}{
%%     ~~Describe \code{method} here~~
}
  \item{n.obs.head}{
%%     ~~Describe \code{n.obs.head} here~~
}
  \item{n.obs.tail}{
%%     ~~Describe \code{n.obs.tail} here~~
}
  \item{logScale}{
%%     ~~Describe \code{logScale} here~~
}
  \item{outfile}{
%%     ~~Describe \code{outfile} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (cd, method = "cubic", n.obs.head = 10, n.obs.tail = 25, 
    logScale = FALSE, outfile = NULL) 
{
    if (!method \%in\% c("lm", "cubic")) {
        stop("Unknown method.\n")
    }
    setkey(cd@x, distbin, refBinMean)
    f.d <- unique(cd@x)[is.na(refBinMean) == FALSE][, c("distbin", 
        "refBinMean"), with = F]
    f.d <- f.d[order(refBinMean, decreasing = TRUE)]
    setDF(f.d)
    f.d$midpoint <- seq(from = round(cd@settings$binsize/2), 
        by = cd@settings$binsize, length.out = nrow(f.d))
    obs.min <- log(min(f.d$midpoint))
    obs.max <- log(max(f.d$midpoint))
    if (method == "lm") {
        log.f.obs <- approxfun(log(f.d$midpoint), log(f.d$refBinMean), 
            rule = c(1, 1))
        head.coef <- coefficients(lm(log(refBinMean) ~ log(midpoint), 
            data = head(f.d, n.obs.head)))
        tail.coef <- coefficients(lm(log(refBinMean) ~ log(midpoint), 
            data = tail(f.d, n.obs.tail)))
        log.f.head <- function(x, head.coef) {
            head.coef[1] + x * head.coef[2]
        }
        log.f.tail <- function(x, tail.coef) {
            tail.coef[1] + x * tail.coef[2]
        }
    }
    if (method == "cubic") {
        f.d.cubic <- lm(log(refBinMean) ~ log(midpoint) + I(log(midpoint)^2) + 
            I(log(midpoint)^3), data = f.d)
        fit <- f.d.cubic$coefficients
        log.f.obs <- function(x, fit. = fit) {
            fit.[1] + fit.[2] * x + fit.[3] * (x^2) + fit.[4] * 
                (x^3)
        }
        obs.min <- log(min(f.d$midpoint))
        obs.max <- log(max(f.d$midpoint))
        beta <- fit[2] + 2 * fit[3] * c(obs.min, obs.max) + 3 * 
            fit[4] * (c(obs.min, obs.max)^2)
        alpha <- fit[1] + (fit[2] - beta) * c(obs.min, obs.max) + 
            fit[3] * c(obs.min, obs.max)^2 + fit[4] * c(obs.min, 
            obs.max)^3
        head.coef <- c(alpha[1], beta[1])
        tail.coef <- c(alpha[2], beta[2])
        log.f.head <- function(x, head.coef. = head.coef) {
            head.coef.[1] + x * head.coef.[2]
        }
        log.f.tail <- function(x, tail.coef. = tail.coef) {
            tail.coef.[1] + x * tail.coef.[2]
        }
    }
    log.f <- function(x, head.coef. = head.coef, tail.coef. = tail.coef, 
        obs.min. = obs.min, obs.max. = obs.max) {
        ifelse(x > obs.max., log.f.tail(x, tail.coef.), ifelse(x < 
            obs.min., log.f.head(x, head.coef.), log.f.obs(x)))
    }
    if (logScale) {
        f <- log.f
    }
    else {
        f <- function(x) exp(log.f(log(x)))
    }
    if (!is.null(outfile)) {
        pdf(outfile)
    }
    curve(log.f.obs, obs.min, obs.max, main = paste0("Distance function (points = obs, line = ", 
        method, " fit)"), xlab = "log(distance)", ylab = "log(f(d))")
    with(f.d, points(log(midpoint), log(refBinMean)))
    if (!is.null(outfile)) {
        dev.off()
    }
    cd@params$f = f
    cd
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
