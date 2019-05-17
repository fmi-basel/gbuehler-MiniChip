context("SimulatePeaks")

test_that("SimulatePeaks works as expected", {
  expect_error(SimulatePeaks())
  expect_error(SimulatePeaks(1000,100,100))
  expect_error(SimulatePeaks(1000,100,chromosomeSizes="work/gbioinfo/DB/genomes/mm10/starIndex_BSgenome.Mmusculus.mm10/chrNameLength.txt"))
})
