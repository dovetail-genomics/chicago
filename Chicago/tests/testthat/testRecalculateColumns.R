library(Chicago)

data("cdUnitTest")

##modifications to cdUnitTest, ensuring it uses correct design directory
designDir <- file.path(system.file("extdata", package="Chicago"), "unitTestDesign")
cdUnitTest <- modifySettings(cd=cdUnitTest, designDir=designDir,
                             settings = list(brownianNoise.samples=1, brownianNoise.subset=NA))
setkey(cdUnitTest@x, baitID, otherEndID)

context("Recalculate cdUnitTest columns")


# normaliseBaits ##Not expected to be recapitulated
# normaliseOtherEnds ##Not expected to be recapitulated
# estimateTechnicalNoise ##Not expected to be recapitulated
# estimateDistFun ## Not expected to be recapitulated
# estimateBrownianNoise
test_that("estimateBrownianNoise reproduces cdUnitTest data", {
  cdRecalc <- copyCD(cdUnitTest)
  cdRecalc <- estimateBrownianNoise(cdRecalc)
  setkey(cdRecalc@x, baitID, otherEndID)
  expect_equal(cdUnitTest@x$Bmean, cdRecalc@x$Bmean)
})

# getPvals
test_that("getPvals reproduces cdUnitTest data", {
  cdRecalc <- copyCD(cdUnitTest)
  cdRecalc <- getPvals(cdRecalc)
  setkey(cdRecalc@x, baitID, otherEndID)
  expect_equal(cdUnitTest@x$log.p, cdRecalc@x$log.p)
})

# getScores
test_that("getScores reproduces cdUnitTest data", {
  cdRecalc <- copyCD(cdUnitTest)
  cdRecalc <- getScores(cdRecalc)
  setkey(cdRecalc@x, baitID, otherEndID)
  expect_equal(cdUnitTest@x$score, cdRecalc@x$score)
})