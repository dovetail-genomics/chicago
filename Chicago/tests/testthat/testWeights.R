library(Chicago)

data("cdUnitTest")

##modifications to cdUnitTest, ensuring it uses correct design directory
designDir <- file.path(system.file("extdata", package="Chicago"), "unitTestDesign")
cdUnitTest <- modifySettings(cd=cdUnitTest, designDir=designDir,
                             settings = list(brownianNoise.samples=1, brownianNoise.subset=NA))

cd <- copyCD(cdUnitTest) ##NB not cd <- cdUnitTest. this line is crucial!

weightsPath <- file.path(system.file("extdata", package="Chicago"), "weights")

context("Weights")

##move weights

##FIXME!

setnames(cd@x, old = "log.w", new="logW.cdUnitTest")
setnames(cd@x, old = "log.q", new="logQ.cdUnitTest")

GMweightSettings <- file.path(weightsPath, "GM12878-2reps.settings")
cd <- modifySettings(cd, settingsFile = GMweightSettings)
cd <- getScores(cd)
setnames(cd@x, old = "log.w", new="logW.GM")
setnames(cd@x, old = "log.q", new="logQ.GM")

hMweightSettings <- file.path(weightsPath, "humanMacrophage-7reps.settings")
cd <- modifySettings(cd, settingsFile = hMweightSettings)
cd <- getScores(cd)
setnames(cd@x, old = "log.w", new="logW.hM")
setnames(cd@x, old = "log.q", new="logQ.hM")

mESCweightSettings <- file.path(weightsPath, "mESC-2reps.settings")
cd <- modifySettings(cd, settingsFile = mESCweightSettings)
cd <- getScores(cd)
setnames(cd@x, old = "log.w", new="logW.mESC")
setnames(cd@x, old = "log.q", new="logQ.mESC")

test_that("using alternative weight settings results in differing weights", {
  expect_true(all(cd@x[, logW.GM != logW.hM]))
  expect_true(all(cd@x[, logW.GM != logW.mESC]))
  expect_true(all(cd@x[, logW.hM != logW.mESC]))
})