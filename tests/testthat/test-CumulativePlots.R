context("CumulativePlots")

#test_that("CumulativePlots works as expected", {
#  expect_error(CumulativePlots(counts=counts))
#})

test_that("CumulativePlots uses a list of count matrices and generates plots as expected", {
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

  mean.plots <- CumulativePlots(res,bamNames = names(res),
                               summarizing = "mean",overlapNames = "NA",plot=FALSE)
 # expect_that(is.numeric(mean.plots$overlap2)=="TRUE")
  expect_error(CumulativePlots(res, bamNames = names(res), summarizing = "mean",span=50,step=2000))
  expect_error(CumulativePlots(res, bamNames = names(res), summarizing = "average",overlapNames = names(peaks)))

})
