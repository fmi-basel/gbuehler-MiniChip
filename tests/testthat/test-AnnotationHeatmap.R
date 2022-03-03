context("AnnotationHeatmap")

test_that("AnnotationHeatmap does not work without the correct input", {
  expect_error(AnnotationHeatmap(plotHM=FALSE))

  peaks <- GenomicRanges::GRanges(
  seqnames = Rle(c("chr1", "chr2", "chr1", "chr3"), c(1, 3, 2, 4)),
  ranges = IRanges(50101:50110, end = 51111:51120),
  strand = Rle(strand(c("-", "+", "*", "+", "-")), c(1, 2, 2, 3, 2)),
  name = head(letters, 10), score = 1:10)
  names(peaks) <- peaks$name

  annotation <- GenomicRanges::GRanges(
      seqnames = Rle(c("chr1", "chr4", "chr1", "chr3"), c(1, 3, 2, 4)),
      ranges = IRanges(50101:50110, end = 51111:51120),
      strand = Rle(strand(c("-", "+", "*", "-", "-")), c(1, 2, 2, 3, 2)),
      name = head(letters, 10))
  names(peaks) <- peaks$name
  expect_error(SummitHeatmap(peaks=peaks,annotation=annotation))

  peaks <- GenomicRanges::GRanges(
    seqnames = Rle(c("chr1", "chr2", "chr1", "chr3"), c(1, 3, 2, 4)),
    ranges = IRanges(50101:50110, end = 51111:51120),
    strand = Rle(strand(c("-", "+", "*", "+", "-")), c(1, 2, 2, 3, 2)),
    name = head(letters, 10), summit = 1:10)
  names(peaks) <- peaks$name
  annotation <- data.frame(chr=c("chr1", "chr4", "chr1", "chr3"),
                         start=c(50101,100101,500101,5000101),
                         end=c(50501,100501,500501,5000501))
  expect_error(SummitHeatmap(peaks=peaks,annotation=annotation))
})

test_that("AnnotationHeatmap returns a  matrix and plots heatmap as expected", {
  peaks <- GRanges(
    seqnames = Rle(c("chr1", "chr2", "chr1", "chr3"), c(1, 3, 2, 4)),
    ranges = IRanges(40101:40110, end = 51111:51120),
    strand = Rle(strand(c("-", "+", "*", "+", "-")), c(1, 2, 2, 3, 2)),
    name = head(letters, 10), summit = 1:10)
  names(peaks) <- peaks$name
  annotation <- GenomicRanges::GRanges(
    seqnames = Rle(c("chr1", "chr4", "chr1", "chr3"), c(1, 3, 2, 4)),
    ranges = IRanges(50101:50110, end = 51111:51120),
    strand = Rle(strand(c("-", "+", "*", "-", "-")), c(1, 2, 2, 3, 2)),
    name = head(letters, 10))
  names(peaks) <- peaks$name
  res <- AnnotationHeatmap(peaks=peaks,annotation=annotation)
  expect_that(res, is_a("matrix"))
 # unlink("heatmaps-mypeaks",recursive=TRUE)
})
