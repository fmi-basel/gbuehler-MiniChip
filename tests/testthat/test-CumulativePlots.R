context("CumulativePlots")

#test_that("CumulativePlots works as expected", {
#  expect_error(CumulativePlots(counts=counts))
#})

test_that("CumulativePlots uses a list of count matrices and generates plots as expected", {
  peaks <- GRanges(
    seqnames = Rle(c("chr1", "chr2", "chr1", "chr3"), c(1, 3, 2, 4)),
    ranges = IRanges(40101:40110, end = 51111:51120),
    strand = Rle(strand(c("-", "+", "*", "+", "-")), c(1, 2, 2, 3, 2)),
    name = head(letters, 10), summit = 1:10)
  names(peaks) <- elementMetadata(peaks)[,"name"]
  bamFiles <- c("/work2/gbuehler/deepSeqRepos/bam//HP1a_wt_ChIP_r1_818F1_multi.bam", "/work2/gbuehler/deepSeqRepos/bam//HP1a_wt_ChIP_r2_818F3_multi.bam")
  res <- SummitHeatmap(peaks=peaks,bamFiles=bamFiles,bamNames=c("reads1","reads2"))

  mean.plots <- CumulativePlots(res,bamNames = names(res),
                               summarizing = "mean",overlapNames = names(peaks))
  expect_that(is.numeric(mean.plots$overlap1)=="TRUE")
#  expect_error(CumulativePlots(res, bamNames = names(res), summarizing = "mean"))
#  expect_error(CumulativePlots(res, bamNames = names(res), summarizing = "average",overlapNames = names(peaks)))

})
