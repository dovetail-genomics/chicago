library(Chicago)

data("cdUnitTest")

##modifications to cdUnitTest, ensuring it uses correct design directory
designDir <- file.path(system.file("extdata", package="Chicago"), "unitTestDesign")
cdUnitTest <- modifySettings(cd=cdUnitTest, designDir=designDir)

context("f(x) = f(f(x)) for key pipeline functions (idempotence)")


# normaliseBaits
test_that("normaliseBaits is idempotent", {
  cd1 <- copy(cdUnitTest)
  cd1 <- normaliseBaits(cd1)
  cd2 <- copy(cd1)
  cd2 <- normaliseBaits(cd2)
  expect_identical(cd1, cd2)
})

# normaliseOtherEnds

##Not enough data to test here

# estimateTechnicalNoise
 test_that("estimateTechnicalNoise is idempotent", {
  cd1 <- copy(cdUnitTest)
  cd1 <- estimateTechnicalNoise(cd1)
  cd2 <- copy(cd1)
  cd2 <- estimateTechnicalNoise(cd2)
   expect_identical(cd1, cd2)
 })

# estimateDistFun
 test_that("estimateDistFun is idempotent", {
   cd1 <- copy(cdUnitTest)
   cd1 <- estimateDistFun(cd1)
   cd2 <- copy(cd1)
   cd2 <- estimateDistFun(cd2)
   expect_identical(cd1, cd2)
 })

# getPvals
 test_that("getPvals is idempotent", {
   cd1 <- copy(cdUnitTest)
   cd1 <- getPvals(cd1)
   cd2 <- copy(cd1)
   cd2 <- getPvals(cd2)
   expect_identical(cd1, cd2)
 })

# getScores
 test_that("getScores is idempotent", {
   cd1 <- copy(cdUnitTest)
   cd1 <- getScores(cd1)
   cd2 <- copy(cd1)
   cd2 <- getScores(cd2)
   expect_identical(cd1, cd2)
 })
