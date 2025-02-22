#' @title AnnotationHeatmap
#'
#' @description Generate heatmaps of Genomic Ranges around positions of interest in a genome.
#'
#' @details This function generates a heatmap of annotations around positions of interest in a genome,
#' for example, the summits of ChIP peaks. An annotation GRanges object must be provided.
#' This can contain regions representing anything of interest, for example promoters, genes, repeats,
#' or transcription factor motif occurences in the genome. The center of your GRanges object defines the coordinates
#' of these positions of interest (eg peak summits, TSSs).  Features that lie on the minus strand will be reversed in the final output. 
#' If at least one region in the annotation GRanges object overlaps a window
#' in the heatmap by 1/2 of the step size provided, this window will be called overlapping (= 1), otherwise it will be called non-overlapping (= 0).
#'
#' @param peaks A Granges object containing your positions of interest in a genome (eg ChIP peak summits). Must include seqnames (chromosomes), start, end, strand, and name.
#' @param annotation A GRanges object containing the annotation ranges you want to plot around peak summits.For example, promoter regions.
#' @param annoname A character scalar to describe the annotation that was provided (for example: "promoters"). Defaults to "annotation".
#' @param span Integer scalar specifying the distance from the peak center to the left and right that you want your heatmap to extend to. Default is 2025.
#' @param step Integer scalar specifying the window size in which annotation overlaps are counted. Default is 50.
#' @param ignoreStrand Logical scalar indicating if hould an overlap should be counted only if on the same strand (ignore.strand=FALSE), or on any strand (ignore.strand=TRUE, default).
#' @param minoverlap Integer scalar indicating the desired minimum overlap of the annotation region with the window of the heatmap 
#' for a window to be counted as overlapping the annotation. Default is 1/2 * \code{step} (half the step size).
#'
#' @return
#' A matrix that contains the overlap with the annotation in each window (column headers=middle of the window)
#' for each peak (rownames=peaknames). 1 = at least half of the window overlaps with a region in the annotation
#' object, 0 = less than half of the window overlaps with the annotation object.
#'
#' @examples
#' peaks <- GenomicRanges::GRanges(
#' seqnames = Rle(c("chr1", "chr2", "chr1", "chr3"), c(1, 3, 2, 4)),
#' ranges = IRanges(50101:50110, end = 51111:51120),
#' strand = Rle(strand(c("-", "+", "*", "+", "-")), c(1, 2, 2, 3, 2)),
#' summit = 1:10, name = head(letters, 10))
#' annotation <- GenomicRanges::GRanges(
#' seqnames = Rle(c("chr1", "chr2", "chr2", "chr3"), c(1, 3, 2, 4)),
#' ranges = IRanges(50101:50110, end = 51111:51120),
#' strand = Rle(strand(c("-", "+", "*", "-", "-")), c(1, 2, 2, 3, 2)),
#' name = head(letters, 10))
#' AnnotationHeatmap(peaks=peaks,annotation=annotation)
#'
#' @importFrom GenomicRanges start
#' @importFrom GenomicRanges resize
#' @importFrom GenomicRanges width
#' @importFrom GenomicRanges strand
#' @importFrom GenomicRanges strand<-
#' @importFrom GenomicRanges slidingWindows
#' @importFrom GenomicRanges seqnames
#' @importFrom GenomicRanges GRanges
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom IRanges overlapsAny
#'
#' @export
AnnotationHeatmap <- function(peaks,annotation,annoname = "annotation", span=2025,step=50,ignoreStrand=TRUE,
                              minoverlap = ceiling(step/2)){

  #take the middle of the GRanges region, then define whole heatmap region
  nwindows <- ceiling((span*2)/step)
  if(nwindows %% 2 == 0){
    nwindows <- nwindows + 1
  }
  regionwidth <- step * nwindows
  peaks <- resize(peaks,width=regionwidth,fix="center")

   if(class(names(peaks)) == "NULL" & class(peaks$name) != "NULL"){
    names(peaks) <- peaks$name
   }

  if(class(strand(peaks)) == "NULL" & class(peaks$strand) != "NULL"){
    strand(peaks) <- peaks$strand
  }

  #remove peaks with negative start values
  peaks <- peaks[start(peaks) >= 0 & width(peaks)== regionwidth]

  #generate window starts and ends across span
  windows <- seq(from=0,to=regionwidth-step,by=step)
  binmids <- windows - regionwidth/2 + step/2

  peaks2plus <- unlist(GenomicRanges::slidingWindows(peaks[strand(peaks)!="-"], width = step, step = step),use.names=FALSE)
  peaks2minus <- rev(unlist(GenomicRanges::slidingWindows(peaks[strand(peaks)=="-"], width = step, step = step),use.names=FALSE))
  peaks2 <- c(peaks2plus,peaks2minus)

  #calculate the overlap in each window with annotation
  overlap <- overlapsAny(peaks2, annotation,minoverlap=minoverlap,ignore.strand=ignoreStrand)

  overlap <- ifelse(overlap ==TRUE,1,0)
  overlap <- matrix(overlap,nrow=length(peaks),ncol=length(windows), byrow=TRUE)
  colnames(overlap) <- binmids
  rownames(overlap) <- matrix(names(peaks2),nrow=length(peaks),ncol=length(windows), byrow=TRUE)[,1]

  #sort the rows by the original peak GRanges object order
  overlap <- overlap[names(peaks),]
  
  #return the results
  return(overlap)
}
