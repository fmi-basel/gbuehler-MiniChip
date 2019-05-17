context("SummitHeatmap")

test_that("SummitHeatmap works as expected", {
  expect_error(SummitHeatmap())
  library(GenomicRanges)
  peaks <- GRanges(
  seqnames = Rle(c("chr1", "chr2", "chr1", "chr3"), c(1, 3, 2, 4)),
  ranges = IRanges(101:110, end = 111:120),
  strand = Rle(strand(c("-", "+", "*", "+", "-")), c(1, 2, 2, 3, 2)),
  name = head(letters, 10), score = 1:10)
  bamFiles <- c("/work2/gbuehler/deepSeqRepos/bam//HP1a_wt_ChIP_r1_818F1_multi.bam", "/work2/gbuehler/deepSeqRepos/bam//HP1a_wt_ChIP_r2_818F3_multi.bam")
  expect_error(SummitHeatmap(peaks,bamFiles,plotHM=FALSE))

  peaks <- GRanges(
    seqnames = Rle(c("chr1", "chr2", "chr1", "chr3"), c(1, 3, 2, 4)),
    ranges = IRanges(101:110, end = 111:120),
    strand = Rle(strand(c("-", "+", "*", "+", "-")), c(1, 2, 2, 3, 2)),
    name = head(letters, 10), summit = 1:10)
  bamFiles <- c("work2/gbuehler/deepSeqRepos/bam//HP1a_wt_ChIP_r1_818F1_multi.bam", "work2/gbuehler/deepSeqRepos/bam//HP1a_wt_ChIP_r2_818F3_multi.bam")
  expect_error(SummitHeatmap(peaks,bamFiles,plotHM=FALSE))
})
