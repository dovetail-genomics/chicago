## ----style, echo = FALSE, results = 'asis'-------------------------------
BiocStyle::markdown()

## ------------------------------------------------------------------------
library(PCHiCdata)

## ------------------------------------------------------------------------
dataPath <- system.file("extdata", package="PCHiCdata")
testDesignDir <- file.path(dataPath, "hg19TestDesign")
dir(testDesignDir)

## ------------------------------------------------------------------------
testDataPath <- file.path(dataPath, "GMchinputFiles")
dir(testDataPath)

files <- c(
    file.path(testDataPath, "GM_rep1.chinput"),
    file.path(testDataPath, "GM_rep2.chinput"),
    file.path(testDataPath, "GM_rep3.chinput")
  )

## ------------------------------------------------------------------------
settingsFile <- file.path(system.file("extdata", package="PCHiCdata"),
                          "sGM12878Settings", "sGM12878.settingsFile")

## ---- message=FALSE------------------------------------------------------
library(Chicago)

cd <- setExperiment(designDir = testDesignDir, settingsFile = settingsFile)

## ---- message=FALSE------------------------------------------------------
cd <- readAndMerge(files=files, cd=cd)

## ---- eval=FALSE---------------------------------------------------------
#  cd <- chicagoPipeline(cd)

## ---- echo=FALSE, message=FALSE------------------------------------------
cd <- chicagoPipeline(cd)

## ------------------------------------------------------------------------
outputDirectory <- tempdir()
exportResults(cd, file.path(outputDirectory, "vignetteOutput"))

## ------------------------------------------------------------------------
featuresFolder <- file.path(dataPath, "GMfeatures")
dir(featuresFolder)

featuresFile <- file.path(featuresFolder, "featuresGM.txt")
featuresTable <- read.delim(featuresFile, header=FALSE, as.is=TRUE)
featuresList <- as.list(featuresTable$V2)
names(featuresList) <- featuresTable$V1
featuresList

## ---- message=FALSE, fig.width=12, fig.height=7--------------------------
no_bins <- ceiling(max(abs(cd@x$distSign), na.rm = T)/1e4)

enrichmentResults <- peakEnrichment4Features(cd, folder=featuresFolder, list_frag=featuresList,
                                             no_bins=no_bins, sample_number=100)

## ------------------------------------------------------------------------
enrichmentResults

## ------------------------------------------------------------------------
head(cd@x, 2)

## ------------------------------------------------------------------------
newCd = copyCD(cd)

## ------------------------------------------------------------------------
sessionInfo()

