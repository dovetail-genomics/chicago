library(Chicago)

data("cdUnitTest")

##modifications to cdUnitTest, ensuring it uses correct design directory
designDir <- file.path(system.file("extdata", package="Chicago"), "unitTestDesign")
cdUnitTest <- modifySettings(cd=cdUnitTest, designDir=designDir,
                             settings = list(brownianNoise.samples=1, brownianNoise.subset=NA))

context("f(x) = f(f(x)) for key pipeline functions (idempotence)")


# normaliseBaits
test_that("normaliseBaits is idempotent", {
  cd1 <- copyCD(cdUnitTest)
  set(cd1@x, NULL, "s_j", NULL)
  set(cd1@x, NULL, "refBinMean", NULL)
  cd1 <- normaliseBaits(cd1)
  cd2 <- copyCD(cd1)
  set(cd2@x, NULL, "s_j", NULL)
  set(cd2@x, NULL, "refBinMean", NULL)
  cd2 <- normaliseBaits(cd2)
  expect_identical(cd1, cd2)
})

# normaliseOtherEnds

##Not enough data to test here

# estimateTechnicalNoise
test_that("estimateTechnicalNoise is idempotent", {
  cd1 <- copyCD(cdUnitTest)
  set(cd1@x, NULL, "tblb", NULL)
  set(cd1@x, NULL, "Tmean", NULL)
  cd1 <- estimateTechnicalNoise(cd1)
  cd2 <- copyCD(cd1)
  set(cd2@x, NULL, "tblb", NULL)
  set(cd2@x, NULL, "Tmean", NULL)
  cd2 <- estimateTechnicalNoise(cd2)
  expect_identical(cd1, cd2)
})

# estimateDistFun
 test_that("estimateDistFun is idempotent", {
   cd1 <- copyCD(cdUnitTest)
   cd1 <- estimateDistFun(cd1)
   cd2 <- copyCD(cd1)
   cd2 <- estimateDistFun(cd2)
   expect_identical(cd1, cd2)
 })

# getPvals
 test_that("getPvals is idempotent", {
   cd1 <- copyCD(cdUnitTest)
   cd1 <- getPvals(cd1)
   cd2 <- copyCD(cd1)
   cd2 <- getPvals(cd2)
   expect_identical(cd1, cd2)
 })

# getScores
 test_that("getScores is idempotent", {
   cd1 <- copyCD(cdUnitTest)
   cd1 <- getScores(cd1)
   cd2 <- copyCD(cd1)
   cd2 <- getScores(cd2)
   expect_identical(cd1, cd2)
 })
