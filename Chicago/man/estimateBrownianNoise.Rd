\name{estimateBrownianNoise}
\alias{estimateBrownianNoise}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
estimateBrownianNoise(cd, distFun, Ncol = "N", reEstimateMean = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cd}{
%%     ~~Describe \code{cd} here~~
}
  \item{distFun}{
%%     ~~Describe \code{distFun} here~~
}
  \item{Ncol}{
%%     ~~Describe \code{Ncol} here~~
}
  \item{reEstimateMean}{
%%     ~~Describe \code{reEstimateMean} here~~
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
function (cd, distFun, Ncol = "N", reEstimateMean = FALSE) 
{
    s = cd@settings
    adjBait2bait = s$adjBait2bait
    subset = s$brownianNoise.subset
    seed = s$brownianNoise.seed
    maxLBrownEst = s$maxLBrownEst
    if (!is.null(seed)) {
        set.seed(seed)
    }
    if (reEstimateMean) {
        stop("reEstimateMean=TRUE not implemented yet.")
    }
    siPresent <- "s_i" \%in\% colnames(cd@x)
    if (siPresent) {
        message("s_i factors found - estimating Brownian noise...")
    }
    else {
        message("s_i factors NOT found - variance will increase, estimating Brownian noise anyway...")
    }
    if (!is.null(subset)) {
        if (!class(subset) \%in\% c("numeric", "integer")) {
            stop("'subset' must be an integer.")
        }
        setkey(cd@x, baitID)
        if (nrow(cd@x[, .I[1], by = baitID]) > subset) {
            sel.sub <- sort(sample(unique(cd@x$baitID), subset))
            x <- cd@x[J(sel.sub)]
        }
        else {
            x <- cd@x
            subset = NULL
            warning("subset > number of baits in data, so use the full dataset.\n")
        }
    }
    else {
        x <- cd@x
    }
    setkey(x, distSign)
    x = x[abs(distSign) < maxLBrownEst & is.na(distSign) == F, 
        ]
    if (adjBait2bait) {
        if (!"isBait2bait" \%in\% names(x)) {
            x[, `:=`(isBait2bait, FALSE)]
            x[wb2b(otherEndID), `:=`(isBait2bait, TRUE)]
        }
        x = x[isBait2bait == FALSE]
    }
    proxOE <- .readProxOEfile(s)
    if (!is.null(subset)) {
        sel.baits <- sel.sub
    }
    else {
        sel.baits <- unique(x$baitID)
    }
    setkey(proxOE, baitID)
    proxOE <- proxOE[J(sel.baits), ]
    setkey(x, baitID, otherEndID)
    sjLookup <- unique(x[, c("baitID", "s_j"), with = FALSE])
    setkey(sjLookup, baitID)
    if (siPresent) {
        siLookup <- unique(x[, c("otherEndID", "s_i"), with = FALSE])
        setkey(siLookup, otherEndID)
    }
    setkey(proxOE, baitID, otherEndID)
    if (siPresent) {
        x <- merge(x, proxOE, all.y = TRUE)[, c("baitID", "otherEndID", 
            "s_i", "s_j", "N", "distSign", "dist"), with = FALSE]
    }
    else {
        x <- merge(x, proxOE, all.y = TRUE)[, c("baitID", "s_j", 
            "N", "distSign", "dist"), with = FALSE]
    }
    x[is.na(N), `:=`(N, 0)]
    x[, `:=`(s_j, sjLookup[J(x$baitID)]$s_j)]
    if (siPresent) {
        x[, `:=`(s_i, siLookup[J(x$otherEndID)]$s_i)]
        if (any(is.na(x$s_i))) {
            warning("Some other ends did not have s_i factors. Assuming s_i = 1 for these.")
            x[, `:=`(s_i, ifelse(is.na(s_i), 1, s_i))]
        }
    }
    if (any(removeNAs(abs(x$dist - abs(x$distSign))) > 1)) {
        warning("estimateBrownianNoise: Distances in precomputed ProxOE file did not match distances supplied.")
    }
    x[, `:=`(distSign, dist)]
    x <- .estimateBMean(x, distFun = cd@params$f)
    message("Calculating dispersion...")
    model <- glm.nb(formula = x$N ~ offset(log(x$Bmean)) + 0)
    cd@params$dispersion <- model$theta
    if (reEstimateMean) {
        stop("Not implemented yet")
    }
    else {
        cd@x <- .estimateBMean(cd@x, cd@params$f)
    }
    cd
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
