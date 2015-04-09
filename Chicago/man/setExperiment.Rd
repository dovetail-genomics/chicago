\name{setExperiment}
\alias{setExperiment}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
setExperiment(designDir = "", settings = list(), settingsFile = NULL, def.settings = list(rmapfile = file.path(designDir, "Digest_Human_HindIII.bed"), baitmapfile = file.path(designDir, "Digest_Human_HindIII_baits_ID.bed"), nperbinfile = file.path(designDir, "Digest_Human_HindIII_NperBin.txt"), nbaitsperbinfile = file.path(designDir, "Digest_Human_HindIII_NbaitsPerBin.txt"), proxOEfile = file.path(designDir, "proxOE_out.txt"), Ncol = "N", baitmapFragIDcol = 4, baitmapGeneIDcol = 5, maxLBrownEst = 1500000, minFragLen = 150, maxFragLen = 40000, minNPerBait = 250, binsize = 20000, removeAdjacent = TRUE, 
    adjBait2bait = TRUE, tlb.filterTopPercent = 0.01, tlb.minProxOEPerBin = 1000, tlb.minProxB2BPerBin = 100, techNoise.minBaitsPerBin = 1000, brownianNoise.subset = 1000, brownianNoise.seed = NULL, baitIDcol = "baitID", otherEndIDcol = "otherEndID", otherEndLencol = "otherEndLen", distcol = "distSign", weightAlpha = 34.1157346557331, weightBeta = -2.58688050486759, weightGamma = -17.1347845819659, weightDelta = -7.07609245521541))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{designDir}{
%%     ~~Describe \code{designDir} here~~
}
  \item{settings}{
%%     ~~Describe \code{settings} here~~
}
  \item{settingsFile}{
%%     ~~Describe \code{settingsFile} here~~
}
  \item{def.settings}{
%%     ~~Describe \code{def.settings} here~~
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
function (designDir = "", settings = list(), settingsFile = NULL, 
    def.settings = list(rmapfile = file.path(designDir, "Digest_Human_HindIII.bed"), 
        baitmapfile = file.path(designDir, "Digest_Human_HindIII_baits_ID.bed"), 
        nperbinfile = file.path(designDir, "Digest_Human_HindIII_NperBin.txt"), 
        nbaitsperbinfile = file.path(designDir, "Digest_Human_HindIII_NbaitsPerBin.txt"), 
        proxOEfile = file.path(designDir, "proxOE_out.txt"), 
        Ncol = "N", baitmapFragIDcol = 4, baitmapGeneIDcol = 5, 
        maxLBrownEst = 1500000, minFragLen = 150, maxFragLen = 40000, 
        minNPerBait = 250, binsize = 20000, removeAdjacent = TRUE, 
        adjBait2bait = TRUE, tlb.filterTopPercent = 0.01, tlb.minProxOEPerBin = 1000, 
        tlb.minProxB2BPerBin = 100, techNoise.minBaitsPerBin = 1000, 
        brownianNoise.subset = 1000, brownianNoise.seed = NULL, 
        baitIDcol = "baitID", otherEndIDcol = "otherEndID", otherEndLencol = "otherEndLen", 
        distcol = "distSign", weightAlpha = 34.1157346557331, 
        weightBeta = -2.58688050486759, weightGamma = -17.1347845819659, 
        weightDelta = -7.07609245521541)) 
{
    modSettings = vector("list")
    if (!is.null(settingsFile)) {
        message(paste0("Reading custom experimental settings from ", 
            settingsFile, "..."))
        sf <- scan(settingsFile, what = "", sep = "\n")
        modSettings <- strsplit(sf, "[[:space:]]+")
        names(modSettings) <- sapply(modSettings, `[[`, 1)
        modSettings <- lapply(modSettings, `[`, -1)
        suppressWarnings({
            modSettings <- lapply(modSettings, function(s) {
                num_s = as.numeric(s)
                if (!is.na(num_s)) {
                  return(num_s)
                }
                bool_s = as.logical(s)
                if (!is.na(bool_s)) {
                  return(bool_s)
                }
                s
            })
        })
    }
    for (s in names(settings)) {
        modSettings[[s]] = settings[[s]]
    }
    for (s in names(modSettings)) {
        def.settings[[s]] = modSettings[[s]]
    }
    cd = chicagoData(x = data.table(), params = list(), settings = def.settings)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
