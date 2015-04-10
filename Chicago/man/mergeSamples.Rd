\name{mergeSamples}
\alias{mergeSamples}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
mergeSamples(cdl, normalise = TRUE, NcolOut = "N", NcolNormPrefix = "NNorm", mergeMethod = c("weightedMean", "mean")[1])
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cdl}{A \code{list} of \code{chicagoData} objects.}
  \item{normalise}{
%%     ~~Describe \code{normalise} here~~
}
  \item{NcolOut}{
%%     ~~Describe \code{NcolOut} here~~
}
  \item{NcolNormPrefix}{
%%     ~~Describe \code{NcolNormPrefix} here~~
}
  \item{mergeMethod}{
%%     ~~Describe \code{mergeMethod} here~~
}
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
function (cdl, normalise = TRUE, NcolOut = "N", NcolNormPrefix = "NNorm", 
    mergeMethod = c("weightedMean", "mean")[1]) 
{
    if (!mergeMethod \%in\% c("weightedMean", "mean")) {
        stop("Unknown mergeMethod.\n")
    }
    attr = vector("list")
    x = as.dataTableList(cdl)
    for (i in 1:length(x)) {
        if (i > 1) {
            if (!identical(cdl[[i]]@settings, cdl[[i - 1]]@settings)) {
                stop("All samples to merge should have identical experiment settings")
            }
        }
        Ncol = grep("^N$", names(x[[i]]), value = T)
        if (!length(Ncol)) {
            Ncol = grep("^N\\.", names(x[[i]]), value = T)
            if (length(Ncol) == 1) {
                warning("Could not find column name \"N\" in element ", 
                  i, " of the input list. Using \"", Ncol, "\" instead\n")
            }
            else {
                stop(paste("Could not find column name \"N\" in element ", 
                  i, " of the input list.\n"))
            }
        }
        setnames(x[[i]], Ncol, paste0("N.", i))
        attr[[i]] = copy(x[[i]])
        set(attr[[i]], NULL, paste0("N.", i), NULL)
        remCols = names(x[[i]])[!names(x[[i]]) \%in\% c("baitID", 
            "otherEndID", paste0("N.", i))]
        for (rc in remCols) {
            set(x[[i]], NULL, rc, NULL)
        }
        setkey(x[[i]], baitID, otherEndID)
    }
    attr = rbindlist(attr)
    setkey(attr, baitID, otherEndID)
    attr = unique(attr)
    message("Merging samples...")
    xmerge = Reduce(function(...) merge(..., by = c("baitID", 
        "otherEndID"), all = T), x)
    setkey(xmerge, baitID, otherEndID)
    xmerge = attr[xmerge]
    for (i in 1:length(x)) {
        iNcol = paste0("N.", i)
        set(xmerge, which(is.na(xmerge[[iNcol]])), iNcol, 0)
    }
    s = cdl[[1]]@settings
    message("Computing merged scores...\n")
    if (normalise) {
        Nnames = names(xmerge)[names(xmerge) \%in\% paste0("N.", 
            1:length(x))]
        s_ks = .getSampleScalingFactors(xmerge, s)
        if (mergeMethod == "weightedMean") {
            xmerge[, eval(parse(text = paste0(NcolOut, ":=round((", 
                paste0(Nnames, "*s_ks[\"", Nnames, "\"]", collapse = "+"), 
                ")/sum(s_ks))")))]
        }
        else {
            xmerge[, eval(parse(text = paste0(NcolOut, ":=round((", 
                paste0(Nnames, "/s_ks[\"", Nnames, "\"]", collapse = "+"), 
                ")/length(Nnames))")))]
            for (i in 1:length(x)) {
                xmerge[, eval(parse(text = paste0(NcolNormPrefix, 
                  ".", i, ":=round(", Nnames[i], "/s_ks[\"", 
                  Nnames[i], "\"])")))]
            }
        }
    }
    else {
        xmerge[, eval(parse(text = paste0(NcolOut, ":=round((", 
            paste0(Nnames, collapse = "+"), ")/length(Nnames))")))]
    }
    set(xmerge, which(!xmerge[[NcolOut]]), NcolOut, 1)
    chicagoData(x = xmerge, params = list(s_k = s_ks), settings = s)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
