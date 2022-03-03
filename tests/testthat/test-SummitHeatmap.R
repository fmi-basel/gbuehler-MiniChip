context("SummitHeatmap")

test_that("SummitHeatmap crashes with wrong options and gives a warning if all windows are 0", {
  expect_error(SummitHeatmap(plotHM=FALSE))

  peaks <- GenomicRanges::GRanges(
  seqnames = Rle(c("chr1", "chr2", "chr1", "chr3"), c(1, 3, 2, 4)),
  ranges = IRanges(50101:50110, end = 51111:51120),
  strand = Rle(strand(c("-", "+", "*", "+", "-")), c(1, 2, 2, 3, 2)),
  name = head(letters, 10), score = 1:10)
  names(peaks) <- NULL
  bamFiles <- list.files(system.file("extdata", package = "MiniChip"), full.names=TRUE,pattern="*bam$")
  expect_warning(SummitHeatmap(peaks=peaks,bamFiles=bamFiles))

})

test_that("SummitHeatmap returns list of count matrices as expected", {
  peaks <- SimulatePeaks(5000,rep(100,5000),chromosomeSizes=
                            system.file("extdata", "chrNameLength_mm10_chr11.txt", package = "MiniChip"))
  bamFiles <- list.files(system.file("extdata", package = "MiniChip"), full.names=TRUE,pattern="*bam$")
  res <- SummitHeatmap(peaks=peaks,bamFiles=bamFiles)
  expect_that(res, is_a("list"))
})
