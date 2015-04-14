\name{readSample}
\alias{readSample}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Read Sample
}
\description{
This function reads input data from a file, into a \code{chicagoData} object.
}
\usage{
readSample(file, cd)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file}{
The location of an input file FIXME more details!
}
  \item{cd}{A \code{chicagoData} object.}
}
\section{Warning}{
The object \code{cd} is updated by reference. Thus, \code{cd} will be altered. See vignette for further information.
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
An object of class \code{chicagoData}.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Mikhail Spivakov, Jonathan Cairns, Paula Freire Pritchett
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
function (file, cd) 
{
    message(paste("Reading", file))
    x = fread(file)
    message("Processing input...")
    s = cd@settings
    if ((!s$baitIDcol \%in\% names(x)) | (!s$otherEndIDcol \%in\% 
        names(x)) | (!s$Ncol \%in\% names(x)) | (!s$otherEndLencol \%in\% 
        names(x)) | (!s$distcol \%in\% names(x))) {
        stop("Named columns baitIDcol = ", s$baitIDcol, ", otherEndIDcol = ", 
            s$otherEndIDcol, ", Ncol = ", s$Ncol, ", otherEndLencol = ", 
            s$otherEndLencol, " and distcol = ", s$distcol, " must be present in the input file. Change these global parameters if names do not macth\n")
    }
    setnames(x, s$baitIDcol, "baitID")
    setnames(x, s$otherEndIDcol, "otherEndID")
    setnames(x, s$Ncol, "N")
    setnames(x, s$distcol, "distSign")
    setnames(x, s$otherEndLencol, "otherEndLen")
    xlen = nrow(x)
    x = x[otherEndLen \%between\% c(s$minFragLen, s$maxFragLen)]
    message("minFragLen = ", s$minFragLen, " maxFragLen = ", 
        s$maxFragLen)
    message("Filtered out ", xlen - nrow(x), " interactions involving other ends < minFragLen or > maxFragLen.")
    setkey(x, baitID)
    baitlen = length(unique(x$baitID))
    x = x[, `:=`(nperbait, sum(N)), by = baitID]
    x = x[nperbait >= s$minNPerBait]
    message("minNPerBait = ", s$minNPerBait)
    message("Filtered out ", baitlen - length(unique(x$baitID)), 
        " baits with < minNPerBait reads.\n")
    set(x, NULL, "nperbait", NULL)
    if (s$removeAdjacent) {
        x[, `:=`(isAdjacent, abs(baitID - otherEndID) == 1), 
            by = baitID]
        x = x[isAdjacent == FALSE]
        set(x, NULL, "isAdjacent", NULL)
        message("Removed interactions with fragments adjacent to baits.")
    }
    baitlen = length(unique(x$baitID))
    x[, `:=`(isBait2bait, FALSE)]
    x[wb2b(otherEndID, s), `:=`(isBait2bait, TRUE)]
    x[, `:=`(isAllB2BProx, {
        prox = abs(distSign) < s$maxLBrownEst & !is.na(distSign)
        if (!length(prox)) {
            TRUE
        }
        else {
            all(isBait2bait[prox])
        }
    }), by = baitID]
    x = x[isAllB2BProx == FALSE]
    set(x, NULL, "isAllB2BProx", NULL)
    message("Filtered out ", baitlen - length(unique(x$baitID)), 
        " baits without proximal non-Bait2bait interactions\n")
    chicagoData(x = x, settings = cd@settings, params = list())
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
