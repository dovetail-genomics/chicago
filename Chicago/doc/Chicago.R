## ----style, echo = FALSE, results = 'asis'-------------------------------
BiocStyle::markdown()

## ---- message=FALSE------------------------------------------------------
library(Chicago)
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

## ----echo=FALSE----------------------------------------------------------

a <- read.table(file.path(outputDirectory, "vignetteOutput.ibed"), header=TRUE)
head(a)


## ----echo=FALSE----------------------------------------------------------

a <- read.table(file.path(outputDirectory, "vignetteOutput_seqmonk.txt"), header=FALSE)
head(a)


## ----echo=FALSE----------------------------------------------------------

a <- read.table(file.path(outputDirectory, "vignetteOutput_washU_text.txt"), header=FALSE)
head(a)


## ----eval=FALSE----------------------------------------------------------
#  
#  plottedBaitIDs <- plotBaits(cd, n=6)
#  

## ----echo=FALSE, fig.height=10-------------------------------------------

invisible(plotBaits(cd, baits=c(415839, 404491, 425847, 417632, 409335, 414114)))


## ------------------------------------------------------------------------
featuresFolder <- file.path(dataPath, "GMfeatures")
dir(featuresFolder)

featuresFile <- file.path(featuresFolder, "featuresGM.txt")
featuresTable <- read.delim(featuresFile, header=FALSE, as.is=TRUE)
featuresList <- as.list(featuresTable$V2)
names(featuresList) <- featuresTable$V1
featuresList

## ---- message=FALSE, fig.width=12, fig.height=7--------------------------
no_bins <- ceiling(max(abs(intData(cd)$distSign), na.rm = TRUE)/1e4)

enrichmentResults <- peakEnrichment4Features(cd, folder=featuresFolder,
              list_frag=featuresList, no_bins=no_bins, sample_number=100)

## ------------------------------------------------------------------------
enrichmentResults

## ---- message=FALSE------------------------------------------------------
library(GenomicInteractions)
library(GenomicRanges)
gi <- exportToGI(cd)

## ---- message=FALSE------------------------------------------------------
library(AnnotationHub)
ah <- AnnotationHub()
hs <- query(ah, c("GRanges", "EncodeDCC", "Homo sapiens", "H3k4me1"))
enhancerTrack <- hs[["AH23254"]]

## ------------------------------------------------------------------------
otherEnds <- anchorTwo(gi)
otherEnds <- renameSeqlevels(otherEnds, c("chr20","chr21"))

## ------------------------------------------------------------------------
findOverlaps(otherEnds, enhancerTrack)

## ------------------------------------------------------------------------
hs["AH23254"]$genome

## ------------------------------------------------------------------------
head(intData(cd), 2)

## ------------------------------------------------------------------------
newCd = copyCD(cd)

## ------------------------------------------------------------------------
weightsPath <- file.path(system.file("extdata", package="Chicago"),
                         "weights")
dir(weightsPath)

## ---- message=FALSE------------------------------------------------------
weightSettings <- file.path(weightsPath, "GM12878-2reps.settings")
cd <- setExperiment(designDir = testDesignDir, settingsFile = weightSettings)

## ------------------------------------------------------------------------
sessionInfo()

