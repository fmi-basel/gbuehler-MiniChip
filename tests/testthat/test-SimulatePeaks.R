context("SimulatePeaks")

test_that("SimulatePeaks works as expected", {
  expect_error(SimulatePeaks())
  expect_error(SimulatePeaks(1000,100,100))
  simPeaks <- SimulatePeaks(1000,rep(100,1000),chromosomeSizes=system.file("extdata", "chrNameLength_mm10_chr11.txt", package = "MiniChip"))
  expect_that(length(simPeaks), equals(1000))
})
