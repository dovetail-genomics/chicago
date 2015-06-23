library(PCHiCdata)

data(sGM12878)
data(smESC)

# check that Chicago recognizes these as valid objects
 test_that("sGM12878 is a valid object", {
   expect_true(validObject(sGM12878))
 })

 test_that("smESC is a valid object", {
   expect_true(validObject(smESC))
 })
