context("DrawSummitHeatmaps")

test_that("DrawSummitHeatmaps works as expected", {
  expect_error(DrawSummitHeatmaps(counts=counts))
})

test_that("DrawSummitHeatmaps uses a list of count matrices and plots heatmap as expected", {
  peaks <- GRanges(
    seqnames = Rle(c("chr1", "chr2", "chr1", "chr3"), c(1, 3, 2, 4)),
    ranges = IRanges(40101:40110, end = 51111:51120),
    strand = Rle(strand(c("-", "+", "*", "+", "-")), c(1, 2, 2, 3, 2)),
    name = head(letters, 10), summit = 1:10)
  names(peaks) <- elementMetadata(peaks)[,"name"]
  bamFiles <- c("/work2/gbuehler/deepSeqRepos/bam//HP1a_wt_ChIP_r1_818F1_multi.bam", "/work2/gbuehler/deepSeqRepos/bam//HP1a_wt_ChIP_r2_818F3_multi.bam")
  res <- SummitHeatmap(peaks=peaks,bamFiles=bamFiles,bamNames=c("reads1","reads2"))

  heatmap_list <- DrawSummitHeatmaps(res, names(res), orderSample = 1, summarizing = "mean", bottomCpm = 0, topCpm=10,
                                     use.log=TRUE,orderWindows=2, TargetHeight=500)
  draw(heatmap_list, padding = unit(c(3, 8, 8, 2), "mm"),show_heatmap_legend=FALSE)
  expect_that(heatmap_list, is_a("HeatmapList"))
  expect_error(DrawSummitHeatmaps(res, names(res), orderSample = 5, summarizing = "mean"))
  expect_error(DrawSummitHeatmaps(res, names(res), orderSample = 1, summarizing = "average"))

})
