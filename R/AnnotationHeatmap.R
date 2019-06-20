#' @title AnnotationHeatmap
#'
#' @description This function generates a heatmap of Genomic Ranges around the summits of ChIP peaks.
#'
#' @details This function generates a heatmap of annotations around the summits of ChIP peaks. An annotation GRanges object must be provided.
#' This can contain regions representing anything of interest, for example promoters, genes, repeats,
#' or transcription factor motif occurences in the genome. If at least oen region in the annotation GRanges object overlaps a window
#' in the heatmap by 1/2 of the step size provided, this window will be called overlapping (= 1), otherwise it will be called non-overlapping (= 0).
#'
#' @param peaks A Granges object containing your ChIP peaks. Must include seqnames (chromosomes), start, end, strand, summit, and name.
#' @param peaknames A character string that describes your peaks. Default is "mypeaks".
#' @param annotation A GRanges object containing the annotation ranges you want to plot around peak summits.For example, promoter regions.
#' @param annoname A name to describe the annotationthat was provided (for example: "promoters"). Defaults to "annotation".
#' @param span The distance from the peak summit to the left and right that you want your heatmap to extand to. Default is 2025.
#' @param step The window size in which reads are counted/plotted in the heatmap. Default is 50.
#' @param plotcols The colors for the heatmap. Default is circlize::colorRamp2(c(0, 5), c("white", "darkblue")).
#' @param plotHM Also plot the heatmap or only generate overlap counts in the windows? Default is TRUE.
#'
#'
#' @return
#' A matrix that contains the overlap with the annotation in each window (column headers=middle of the window)
#' for each peak (rownames=peaknames). 1 = at least half of the window overlaps with a region in the annotation
#' object, 0 = less than half of the window overlaps with the annotation object. If plotHM is TRUE, a
#' heatmap is plotted as well and saved in a directory (in your current working directory. if it already exists,
#' it is not overwritten) that is named after the peaknames object. This heatmaps is clustered.
#'
#' @examples
#' peaks <- GenomicRanges::GRanges(
#' seqnames = Rle(c("chr1", "chr2", "chr1", "chr3"), c(1, 3, 2, 4)),
#' ranges = IRanges(50101:50110, end = 51111:51120),
#' strand = Rle(strand(c("-", "+", "*", "+", "-")), c(1, 2, 2, 3, 2)),
#' summit = 1:10, name = head(letters, 10))
#' annotation <- GenomicRanges::GRanges(
#' seqnames = Rle(c("chr1", "chr4", "chr1", "chr3"), c(1, 3, 2, 4)),
#' ranges = IRanges(50101:50110, end = 51111:51120),
#' strand = Rle(strand(c("-", "+", "*", "-", "-")), c(1, 2, 2, 3, 2)),
#' name = head(letters, 10))
#' AnnotationHeatmap(peaks=peaks,annotation=annotation,plotHM=FALSE)
#'
#' @importFrom GenomicRanges start
#' @importFrom GenomicRanges start<-
#' @importFrom GenomicRanges end
#' @importFrom GenomicRanges end<-
#' @importFrom GenomicRanges strand
#' @importFrom GenomicRanges elementMetadata
#' @importFrom GenomicRanges seqnames
#' @importFrom GenomicRanges GRanges
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom IRanges overlapsAny
#' @importFrom ComplexHeatmap Heatmap
#' @importFrom ComplexHeatmap HeatmapAnnotation
#' @importFrom ComplexHeatmap anno_text
#' @importFrom ComplexHeatmap max_text_width
#' @importFrom ComplexHeatmap draw
#' @importFrom grid gpar
#' @importFrom grid unit
#' @importFrom circlize colorRamp2
#' @importFrom grDevices pdf
#' @importFrom grDevices dev.off
#'
#' @export
AnnotationHeatmap <- function(peaks,peaknames="mypeaks",annotation,annoname = "annotation", span=2025,step=50,
                        plotcols=circlize::colorRamp2(c(0, 5), c("white", "darkblue")),plotHM=TRUE){
  #test ifGranges object has all required columns
  if(is.null(peaks$summit)) { stop("Please make sure your peaks GRanges object has a summit metadata culumn")}
  #make the summit of the peak the start and end, then define the whole heatmap region
  start(peaks) <- start(peaks) + peaks$summit - span
  end(peaks) <- start(peaks) + span*2
  #remove peaks with negative start values
  peaks <- peaks[start(peaks) >= 0]

  #generate window starts and ends across span
  windows <- seq(from=0,to=span*2-step,by=step)
  binmids <- windows - span + step/2

  #generate empty saf format data frame
  saf <- data.frame(GeneID= NA,Chr=NA, Start=NA,End=NA,Strand=NA)
  saf <- saf[FALSE, ]

  #generate genomic coordinates for each window across all peaks, concatenate these windows into one big saf format file
  for (w in seq_along(windows)){
    saf.w <- data.frame(GeneID= elementMetadata(peaks)[,"name"], Chr=seqnames(peaks),
                        Start=start(peaks) + windows[[w]], End=start(peaks) + windows[[w]] + step,Strand=strand(peaks))
    saf <- rbind(saf,saf.w)
  }

  #convert back to GRanges object
  saf.gr <- makeGRangesFromDataFrame(saf,
                                     keep.extra.columns=TRUE,
                                     ignore.strand=TRUE,
                                     seqinfo=NULL,
                                     seqnames.field=c("Chr"),
                                     start.field=c("Start"),
                                     end.field=c("End"),
                                     starts.in.df.are.0based=FALSE)


  #calculate the overlap in each window with annotation
  overlap <- overlapsAny(saf.gr, annotation,minoverlap=round(step/2))
  overlap <- ifelse(overlap ==TRUE,1,0)
  overlap <- matrix(overlap,nrow=length(peaks),ncol=length(windows), byrow=FALSE)
  colnames(overlap) <- binmids
  rownames(overlap) <- peaks$name

  if (plotHM == TRUE){
    #create directory for heatmaps
    ifelse(!dir.exists(sprintf("heatmaps-%s",peaknames)), dir.create(sprintf("heatmaps-%s",peaknames)), FALSE)
    #define colors and window annotations for heatmaps
    #col_fun = circlize::colorRamp2(c(0, 3), c("white", "darkblue"))
    binnames <- ifelse((binmids/1000)%%1==0,binmids,"")

    #plot and save heatmaps
    ht <- ComplexHeatmap::Heatmap(overlap,name = annoname, cluster_rows = TRUE, cluster_columns=FALSE,
                                  column_order = colnames(overlap), col = plotcols,
                                  column_title = annoname, column_title_gp = gpar(fontsize = 10, fontface = "bold"),
                                  show_row_names = FALSE, show_column_names = FALSE,show_row_dend = FALSE,
                                  bottom_annotation = HeatmapAnnotation(
                                    text = anno_text(binnames, rot = 45, offset = unit(1, "npc"), just = "right",gp=gpar(fontsize=10,fontface = "italic")),
                                    annotation_height = max_text_width(binnames)),use_raster = TRUE, raster_device = "png",raster_quality=1
    )

    pdf(sprintf("heatmaps-%s/heatmap-%s-%s.pdf",peaknames,peaknames,annoname),height=20,width=10)
    ComplexHeatmap::draw(ht)
    dev.off()

  } else {
    next
  }

  #return the results
  return(overlap)
}
