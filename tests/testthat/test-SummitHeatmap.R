context("SummitHeatmap")

test_that("SummitHeatmap generates window counts as expected", {
  expect_error(SummitHeatmap(plotHM=FALSE))

  #library(GenomicRanges)
  peaks <- GenomicRanges::GRanges(
  seqnames = Rle(c("chr1", "chr2", "chr1", "chr3"), c(1, 3, 2, 4)),
  ranges = IRanges(50101:50110, end = 51111:51120),
  strand = Rle(strand(c("-", "+", "*", "+", "-")), c(1, 2, 2, 3, 2)),
  name = head(letters, 10), score = 1:10)
  bamFiles <- c("/work2/gbuehler/deepSeqRepos/bam//HP1a_wt_ChIP_r1_818F1_multi.bam", "/work2/gbuehler/deepSeqRepos/bam//HP1a_wt_ChIP_r2_818F3_multi.bam")
  expect_error(SummitHeatmap(peaks=peaks,bamFiles=bamFiles,plotHM=FALSE))

  peaks <- GenomicRanges::GRanges(
    seqnames = Rle(c("chr1", "chr2", "chr1", "chr3"), c(1, 3, 2, 4)),
    ranges = IRanges(50101:50110, end = 51111:51120),
    strand = Rle(strand(c("-", "+", "*", "+", "-")), c(1, 2, 2, 3, 2)),
    name = head(letters, 10), summit = 1:10)
  bamFiles <- c("work2/gbuehler/deepSeqRepos/bam//HP1a_wt_ChIP_r1_818F1_multi.bam", "work2/gbuehler/deepSeqRepos/bam//HP1a_wt_ChIP_r2_818F3_multi.bam")
  expect_error(SummitHeatmap(peaks=peaks,bamFiles=bamFiles,plotHM=FALSE))
})

test_that("SummitHeatmap returns list of count matrices and plots heatmap as expected", {
  #library(GenomicRanges)
  peaks <- GRanges(
    seqnames = Rle(c("chr1", "chr2", "chr1", "chr3"), c(1, 3, 2, 4)),
    ranges = IRanges(40101:40110, end = 51111:51120),
    strand = Rle(strand(c("-", "+", "*", "+", "-")), c(1, 2, 2, 3, 2)),
    name = head(letters, 10), summit = 1:10)
  bamFiles <- c("/work2/gbuehler/deepSeqRepos/bam//HP1a_wt_ChIP_r1_818F1_multi.bam", "/work2/gbuehler/deepSeqRepos/bam//HP1a_wt_ChIP_r2_818F3_multi.bam")
  res <- SummitHeatmap(peaks=peaks,bamFiles=bamFiles,plotHM=TRUE)
  expect_that(res, is_a("list"))
  unlink("heatmaps-mypeaks",recursive=TRUE)
})
