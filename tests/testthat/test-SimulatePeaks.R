context("SimulatePeaks")

test_that("SimulatePeaks works as expected", {
  expect_error(SimulatePeaks())
  expect_error(SimulatePeaks(1000,100,100))
})
