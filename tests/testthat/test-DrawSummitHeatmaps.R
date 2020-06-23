context("DrawSummitHeatmaps")

test_that("DrawSummitHeatmaps works as expected", {
  expect_error(DrawSummitHeatmaps(counts=counts))
})

test_that("DrawSummitHeatmaps uses a list of count matrices and plots heatmap as expected", {
  peaks1.d <- read.table(system.file("extdata", "Adnp_rep1_chr11_peaks.narrowPeak", package = "MiniChip"),header=FALSE)
  names(peaks1.d) <- c("chr","start","end","name","score","empty","foldchange","pvalue","qvalue","summit")
  peaks <- makeGRangesFromDataFrame(peaks1.d,
                                    keep.extra.columns=TRUE,
                                    ignore.strand=TRUE,
                                    seqinfo=NULL,
                                    seqnames.field=c("chr"),
                                    start.field=c("start"),
                                    end.field=c("end"),
                                    starts.in.df.are.0based=FALSE)
  bamFiles <- list.files(system.file("extdata", package = "MiniChip"), full.names=TRUE,pattern="*bam$")
  bamNames <- gsub(paste(system.file("extdata", package = "MiniChip"),"/",sep=""),"",bamFiles)
  bamNames <- gsub("_chr11.bam","",bamNames)
  res <- SummitHeatmap(peaks=peaks,bamFiles=bamFiles,bamNames=bamNames)

  heatmap_list <- DrawSummitHeatmaps(res, names(res), orderSample = 1, summarizing = "mean",
                                     use.log=TRUE,orderWindows=2)
  #draw(heatmap_list, padding = unit(c(3, 8, 8, 2), "mm"),show_heatmap_legend=FALSE)
  expect_that(heatmap_list, is_a("HeatmapList"))
  expect_error(DrawSummitHeatmaps(res, names(res), orderSample = 10, summarizing = "mean"))
  expect_error(DrawSummitHeatmaps(res, names(res), orderSample = 1, summarizing = "average"))

})
