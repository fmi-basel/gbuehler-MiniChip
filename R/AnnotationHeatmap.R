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
AnnotationHeatmap <- function(peaks, annotation, annoname = "annotation", span = 2025, step = 50, 
                              ignoreStrand = TRUE, minoverlap = ceiling(step/2)) {
  
  # Check required packages
  if (!requireNamespace("GenomicRanges", quietly = TRUE)) {
    stop("Package 'GenomicRanges' is required but not installed. Please install it using BiocManager::install('GenomicRanges').")
  }
  if (!requireNamespace("IRanges", quietly = TRUE)) {
    stop("Package 'IRanges' is required but not installed. Please install it using BiocManager::install('IRanges').")
  }
  
  # Check input parameters
  if (missing(peaks)) {
    stop("'peaks' parameter is required but missing.")
  }
  if (missing(annotation)) {
    stop("'annotation' parameter is required but missing.")
  }
  
  # Check class of input objects
  if (!inherits(peaks, "GRanges")) {
    stop("'peaks' must be a GRanges object. Please provide a valid GRanges object.")
  }
  if (!inherits(annotation, "GRanges")) {
    stop("'annotation' must be a GRanges object. Please provide a valid GRanges object.")
  }
  
  # Check if peaks is empty
  if (length(peaks) == 0) {
    stop("'peaks' GRanges object is empty. Please provide a non-empty GRanges object.")
  }
  
  # Check if annotation is empty
  if (length(annotation) == 0) {
    stop("'annotation' GRanges object is empty. Please provide a non-empty GRanges object.")
  }
  
  # Check scalar parameters
  if (!is.character(annoname) || length(annoname) != 1) {
    stop("'annoname' must be a character scalar.")
  }
  if (!is.numeric(span) || length(span) != 1 || span <= 0) {
    stop("'span' must be a positive numeric scalar.")
  }
  if (!is.numeric(step) || length(step) != 1 || step <= 0) {
    stop("'step' must be a positive numeric scalar.")
  }
  if (step > span) {
    warning("'step' is larger than 'span'. This may result in poor resolution.")
  }
  if (!is.logical(ignoreStrand) || length(ignoreStrand) != 1) {
    stop("'ignoreStrand' must be a logical scalar (TRUE or FALSE).")
  }
  if (!is.numeric(minoverlap) || length(minoverlap) != 1 || minoverlap < 0) {
    stop("'minoverlap' must be a non-negative numeric scalar.")
  }
  if (minoverlap > step) {
    warning("'minoverlap' is larger than 'step'. This may result in no overlaps being found.")
  }
  
  # Check if peaks have names
  if (is.null(names(peaks)) && is.null(peaks$name)) {
    stop("'peaks' must have either 'names' attribute or a 'name' metadata column.")
  }
  
  # Check if peaks have strand information
  if (is.null(strand(peaks)) && is.null(peaks$strand)) {
    warning("No strand information found in 'peaks'. All peaks will be treated as being on the positive strand.")
    strand(peaks) <- "+"
  }
  
  # take the middle of the GRanges region, then define whole heatmap region
  nwindows <- ceiling((span*2)/step)
  if(nwindows %% 2 == 0){
    nwindows <- nwindows + 1
  }
  regionwidth <- step * nwindows
  
  # Try-catch for resizing peaks
  tryCatch({
    peaks <- resize(peaks, width=regionwidth, fix="center")
  }, error = function(e) {
    stop("Error resizing peaks: ", e$message, 
         "\nMake sure your peaks have valid coordinates and are properly formatted.")
  })
  
  # Assign names if needed
  if(is.null(names(peaks)) && !is.null(peaks$name)){
    names(peaks) <- peaks$name
  }
  
  # Assign strand if needed
  if(is.null(strand(peaks)) && !is.null(peaks$strand)){
    strand(peaks) <- peaks$strand
  }
  
  # Remove peaks with negative start values and check width
  originalLength <- length(peaks)
  peaks <- peaks[start(peaks) >= 0 & width(peaks) == regionwidth]
  if (length(peaks) == 0) {
    stop("After filtering out peaks with negative start values or incorrect width, no peaks remain. Check your input data.")
  }
  if (length(peaks) < originalLength) {
    warning(paste0(originalLength - length(peaks), " peaks were removed due to negative start values or incorrect width."))
  }
  
  # Generate window starts and ends across span
  windows <- seq(from=0, to=regionwidth-step, by=step)
  binmids <- windows - regionwidth/2 + step/2
  
  # Try-catch for slidingWindows operation
  tryCatch({
    # Process peaks on plus strand
    if (sum(strand(peaks) != "-") > 0) {
      peaks2plus <- unlist(GenomicRanges::slidingWindows(peaks[strand(peaks) != "-"], 
                                                         width = step, step = step), 
                           use.names = FALSE)
    } else {
      peaks2plus <- GenomicRanges::GRanges()
    }
    
    # Process peaks on minus strand
    if (sum(strand(peaks) == "-") > 0) {
      peaks2minus <- rev(unlist(GenomicRanges::slidingWindows(peaks[strand(peaks) == "-"], 
                                                              width = step, step = step), 
                                use.names = FALSE))
    } else {
      peaks2minus <- GenomicRanges::GRanges()
    }
    
    peaks2 <- c(peaks2plus, peaks2minus)
    
    if (length(peaks2) == 0) {
      stop("Error generating sliding windows. Check your input data.")
    }
    
  }, error = function(e) {
    stop("Error generating sliding windows: ", e$message, 
         "\nMake sure your peaks have valid coordinates and are properly formatted.")
  })
  
  # Calculate the overlap in each window with annotation
  tryCatch({
    overlap <- overlapsAny(peaks2, annotation, minoverlap = minoverlap, ignore.strand = ignoreStrand)
  }, error = function(e) {
    stop("Error calculating overlaps: ", e$message, 
         "\nCheck that your annotation GRanges object is valid and compatible with your peaks.")
  })
  
  overlap <- ifelse(overlap == TRUE, 1, 0)
  
  # Create and populate the matrix
  tryCatch({
    overlap <- matrix(overlap, nrow = length(peaks), ncol = length(windows), byrow = TRUE)
    colnames(overlap) <- binmids
    rownames(overlap) <- matrix(names(peaks2), nrow = length(peaks), ncol = length(windows), byrow = TRUE)[,1]
    
    # Sort the rows by the original peak GRanges object order
    if (all(names(peaks) %in% rownames(overlap))) {
      overlap <- overlap[names(peaks),]
    } else {
      warning("Some peak names not found in result matrix. Check your input data for duplicated or missing names.")
    }
  }, error = function(e) {
    stop("Error creating result matrix: ", e$message, 
         "\nThis may be due to inconsistent data dimensions or naming issues.")
  })
  
  # Return the results
  return(overlap)
}