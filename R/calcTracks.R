#' @title calcTracks
#'
#' @description Prepare ChIP, repeat, and gene annotation data for plotting tracks 
#' for a given region of interest (ROI).
#'
#' @details calcTracks is used to prepare ChIP, repeat, and gene annotation data for plotting tracks 
#' for a given region of interest (ROI).
#'
#' @param bigwigFiles Character vector containing the filenames (including the full path) of read alignment files in bigwig format.
#' @param bigwigNames Character vector containing the names to describe the \code{bigwigFiles} you are using (for example: "H3K9me3_replicate1"). 
#' @param bigwigGroupNames Character vector containing the group names to describe the \code{bigwigFiles} you are using (for example: "H3K9me3"). 
#' @param ROI A GRanges object containing your single region of interest that should be plotted. Must include seqnames (chromosomes), start, end.
#' @param CoverageSmoothing Logical scalar, indicating whether the coverage tracks should be smoothed. Default: TRUE.
#' @param smoothing_window Numeric scalar providing the window size for smoothing. Default: 20.
#' @param reps A GRanges object containing all repeats in the genome, or at least the ones overlapping the \code{ROI}.
#' @param genes A GRanges object containing all genes in the genome, or at least the ones overlapping the \code{ROI}.
#' @param exons A GRanges object containing all exons in the genome, or at least the ones overlapping the \code{ROI}.
#'
#' @return A list of data frames containing all necessary information to plot
#'  the screenshots using the \code{plotTracks} function.
#'
#' @examples
#' #see vignette
#'
#' @importFrom rtracklayer import.bw
#' @importFrom dplyr group_by mutate ungroup summarise
#' @importFrom zoo rollmean
#' @importFrom GenomicRanges GRanges start end strand seqnames promoters terminators start<- end<-
#' @importFrom IRanges IRanges subsetByOverlaps
#' @importFrom methods is
#'
#' @export
calcTracks <- function(bigwigFiles, bigwigNames, bigwigGroupNames, ROI, CoverageSmoothing=TRUE, smoothing_window=20,
                       reps, genes, exons){ 
  
  # Input validation
  
  # Check if required packages are installed
  required_packages <- c("rtracklayer", "GenomicRanges", "IRanges", "dplyr", "zoo")
  missing_packages <- required_packages[!sapply(required_packages, function(p) requireNamespace(p, quietly = TRUE))]
  if (length(missing_packages) > 0) {
    stop(paste("Missing required packages:", paste(missing_packages, collapse=", "), 
               "\nPlease install them using BiocManager::install() or install.packages()"))
  }
  
  # Check bigwigFiles parameter
  if (missing(bigwigFiles)) {
    stop("Error: 'bigwigFiles' argument is missing. Please provide paths to your bigwig files.")
  }
  if (!is.character(bigwigFiles)) {
    stop("Error: 'bigwigFiles' must be a character vector containing file paths.")
  }
  
  # Check file existence
  non_existent_files <- bigwigFiles[!file.exists(bigwigFiles)]
  if (length(non_existent_files) > 0) {
    stop(paste("Error: The following bigwig files do not exist:", 
               paste(non_existent_files, collapse=", ")))
  }
  
  # Check bigwigNames parameter
  if (missing(bigwigNames)) {
    stop("Error: 'bigwigNames' argument is missing. Please provide names for your bigwig files.")
  }
  if (!is.character(bigwigNames)) {
    stop("Error: 'bigwigNames' must be a character vector.")
  }
  if (length(bigwigNames) != length(bigwigFiles)) {
    stop(paste("Error: The length of 'bigwigNames' (", length(bigwigNames), 
               ") does not match the length of 'bigwigFiles' (", length(bigwigFiles), ")."))
  }
  
  # Check bigwigGroupNames parameter
  if (missing(bigwigGroupNames)) {
    stop("Error: 'bigwigGroupNames' argument is missing. Please provide group names for your bigwig files.")
  }
  if (!is.character(bigwigGroupNames)) {
    stop("Error: 'bigwigGroupNames' must be a character vector.")
  }
  if (length(bigwigGroupNames) != length(bigwigFiles)) {
    stop(paste("Error: The length of 'bigwigGroupNames' (", length(bigwigGroupNames), 
               ") does not match the length of 'bigwigFiles' (", length(bigwigFiles), ")."))
  }
  
  # Check ROI parameter
  if (missing(ROI)) {
    stop("Error: 'ROI' argument is missing. Please provide a GRanges object with your region of interest.")
  }
  if (!methods::is(ROI, "GRanges")) {
    stop("Error: 'ROI' must be a GRanges object.")
  }
  if (length(ROI) != 1) {
    stop(paste("Error: 'ROI' must contain exactly one region, but contains", length(ROI), "regions."))
  }
  if (any(is.na(GenomicRanges::start(ROI))) || any(is.na(GenomicRanges::end(ROI)))) {
    stop("Error: 'ROI' contains NA values in start or end positions.")
  }
  
  # Check CoverageSmoothing parameter
  if (!is.logical(CoverageSmoothing)) {
    stop("Error: 'CoverageSmoothing' must be a logical value (TRUE or FALSE).")
  }
  
  # Check smoothing_window parameter
  if (!is.numeric(smoothing_window) || length(smoothing_window) != 1) {
    stop("Error: 'smoothing_window' must be a single numeric value.")
  }
  if (smoothing_window < 1) {
    stop("Error: 'smoothing_window' must be at least 1.")
  }
  if (smoothing_window %% 1 != 0) {
    warning("Warning: 'smoothing_window' is not an integer and will be rounded.")
    smoothing_window <- round(smoothing_window)
  }
  
  # Check reps parameter
  if (missing(reps)) {
    warning("Warning: 'reps' argument is missing. No repeat annotations will be included.")
    reps <- GenomicRanges::GRanges()
  } else if (!methods::is(reps, "GRanges")) {
    stop("Error: 'reps' must be a GRanges object.")
  } else {
    # Check for required metadata columns in reps
    if (!"repeat_name" %in% colnames(GenomicRanges::mcols(reps))) {
      stop("Error: 'reps' must contain a metadata column named 'repeat_name'.")
    }
    if (any(is.na(GenomicRanges::mcols(reps)$repeat_name))) {
      warning("Warning: NA values found in repeat_name column of 'reps'.")
    }
  }
  
  # Check genes parameter
  if (missing(genes)) {
    warning("Warning: 'genes' argument is missing. No gene annotations will be included.")
    genes <- GenomicRanges::GRanges()
  } else if (!methods::is(genes, "GRanges")) {
    stop("Error: 'genes' must be a GRanges object.")
  } else {
    # Check for required metadata columns in genes
    if (!"gene_name" %in% colnames(GenomicRanges::mcols(genes))) {
      stop("Error: 'genes' must contain a metadata column named 'gene_name'.")
    }
    if (any(is.na(GenomicRanges::mcols(genes)$gene_name))) {
      warning("Warning: NA values found in gene_name column of 'genes'.")
    }
  }
  
  # Check exons parameter
  if (missing(exons)) {
    warning("Warning: 'exons' argument is missing. No exon annotations will be included.")
    exons <- GenomicRanges::GRanges()
  } else if (!methods::is(exons, "GRanges")) {
    stop("Error: 'exons' must be a GRanges object.")
  }
  
  # Initialize a list to store data for all BigWig files
  coverage_data_list <- list()
  
  # Loop through each BigWig file
  for (i in seq_along(bigwigFiles)) {
    tryCatch({
      # Import the BigWig file for the given region
      bwf <- rtracklayer::import.bw(bigwigFiles[i], which = ROI)
      
      if (length(bwf) == 0) {
        warning(paste("Warning: No data found in", bigwigFiles[i], "for the specified ROI."))
        next
      }
      
      # Convert to data frame for processing
      for_plotting <- as.data.frame(bwf)
      
      # Expand scores and indices
      expanded_scores <- unlist(mapply(rep, for_plotting$score, for_plotting$width))
      expanded_seqnames <- unlist(mapply(rep, for_plotting$seqnames, for_plotting$width))
      expanded_index <- unlist(mapply(function(start, width) seq(start, by = 1, length.out = width), 
                                      for_plotting$start, for_plotting$width))
      
      # Create a new data frame for coverage
      coverage_table <- data.frame(
        seqnames = expanded_seqnames,
        index = expanded_index,
        score = expanded_scores,
        file_name = bigwigNames[i],  # Add a column for the BigWig file name
        group = bigwigGroupNames[i] # Add a column for the sample group name
      )
      
      # Store the result in the list
      coverage_data_list[[i]] <- coverage_table
      
    }, error = function(e) {
      stop(paste("Error processing", bigwigFiles[i], ":", e$message))
    })
  }
  
  # Check if any coverage data was successfully processed
  if (length(coverage_data_list) == 0) {
    stop("Error: No coverage data could be extracted from any of the bigwig files. Please check your files and ROI.")
  }
  
  # Combine all the coverage data into a single data frame
  combined_coverage <- do.call(rbind, coverage_data_list)
  
  if (nrow(combined_coverage) == 0) {
    stop("Error: No coverage data found. Please check if your ROI overlaps with data in the bigwig files.")
  }
  
  # Define a smoothing function
  smooth_coverage <- function(data, smoothing_window = 20) {
    data %>%
      group_by(file_name) %>%
      mutate(score = zoo::rollmean(score, k = smoothing_window, fill = NA, align = "center")) %>%
      ungroup()
  }
  
  if (CoverageSmoothing == TRUE) {
    # Apply smoothing to the combined coverage data
    tryCatch({
      combined_coverage <- smooth_coverage(combined_coverage, smoothing_window = smoothing_window)
    }, error = function(e) {
      warning(paste("Warning: Error during smoothing:", e$message, 
                    "Continuing with unsmoothed data."))
    })
  }
  
  # Filter annotations to the combined index range
  
  ### generate window GRanges
  score_window <- GenomicRanges::GRanges(
    seqnames = unique(combined_coverage$seqnames),
    ranges = IRanges::IRanges(start = min(combined_coverage$index), end = max(combined_coverage$index))
  )
  
  # Initialize empty data frames for return values
  reps_df <- data.frame(
    start = numeric(0),
    end = numeric(0),
    repeat_name = character(0),
    ypos = numeric(0)
  )
  
  genes_df <- data.frame(
    start = numeric(0),
    end = numeric(0),
    gene_name = character(0),
    ypos = numeric(0)
  )
  
  exons_df <- data.frame(
    start = numeric(0),
    end = numeric(0),
    ypos = numeric(0)
  )
  
  # Process repeats if they exist
  if (length(reps) > 0) {
    #filter repeats
    reps_filtered <- IRanges::subsetByOverlaps(reps, score_window, ignore.strand=TRUE)
    
    if (length(reps_filtered) > 0) {
      #truncate start and end coordinates to fit into plotting window
      GenomicRanges::start(reps_filtered) <- ifelse(GenomicRanges::start(reps_filtered) >= GenomicRanges::start(score_window),
                                                    GenomicRanges::start(reps_filtered),
                                                    GenomicRanges::start(score_window))
      GenomicRanges::end(reps_filtered) <- ifelse(GenomicRanges::end(reps_filtered) <= GenomicRanges::end(score_window),
                                                  GenomicRanges::end(reps_filtered),
                                                  GenomicRanges::end(score_window))
      
      # Convert filtered reps to a data frame for plotting
      reps_df <- data.frame(
        start = GenomicRanges::start(GenomicRanges::promoters(reps_filtered, upstream=0, downstream=1)),
        end = GenomicRanges::end(GenomicRanges::terminators(reps_filtered, upstream=0, downstream=1)),
        repeat_name = GenomicRanges::mcols(reps_filtered)$repeat_name,
        ypos = as.numeric(paste0(GenomicRanges::strand(reps_filtered), 5))
      )
    }
  }
  
  # Process genes if they exist
  if (length(genes) > 0) {
    #filter genes
    genes_filtered <- IRanges::subsetByOverlaps(genes, score_window, ignore.strand=TRUE)
    
    if (length(genes_filtered) > 0) {
      #truncate start and end coordinates to fit into plotting window
      GenomicRanges::start(genes_filtered) <- ifelse(GenomicRanges::start(genes_filtered) >= GenomicRanges::start(score_window),
                                                     GenomicRanges::start(genes_filtered),
                                                     GenomicRanges::start(score_window))
      GenomicRanges::end(genes_filtered) <- ifelse(GenomicRanges::end(genes_filtered) <= GenomicRanges::end(score_window),
                                                   GenomicRanges::end(genes_filtered),
                                                   GenomicRanges::end(score_window))
      
      # Convert filtered genes to a data frame for plotting
      genes_df <- data.frame(
        start = GenomicRanges::start(GenomicRanges::promoters(genes_filtered, upstream=0, downstream=1)),
        end = GenomicRanges::end(GenomicRanges::terminators(genes_filtered, upstream=0, downstream=1)),
        gene_name = GenomicRanges::mcols(genes_filtered)$gene_name,
        ypos = as.numeric(paste0(GenomicRanges::strand(genes_filtered), 5))
      )
    }
  }
  
  # Process exons if they exist
  if (length(exons) > 0) {
    #filter exons
    exons_filtered <- IRanges::subsetByOverlaps(exons, score_window, ignore.strand=TRUE)
    
    if (length(exons_filtered) > 0) {
      #truncate start and end coordinates to fit into plotting window
      GenomicRanges::start(exons_filtered) <- ifelse(GenomicRanges::start(exons_filtered) >= GenomicRanges::start(score_window),
                                                     GenomicRanges::start(exons_filtered),
                                                     GenomicRanges::start(score_window))
      GenomicRanges::end(exons_filtered) <- ifelse(GenomicRanges::end(exons_filtered) <= GenomicRanges::end(score_window),
                                                   GenomicRanges::end(exons_filtered),
                                                   GenomicRanges::end(score_window))
      
      # Convert filtered exons to a data frame for plotting
      exons_df <- data.frame(
        start = GenomicRanges::start(GenomicRanges::promoters(exons_filtered, upstream=0, downstream=1)),
        end = GenomicRanges::end(GenomicRanges::terminators(exons_filtered, upstream=0, downstream=1)),
        ypos = as.numeric(paste0(GenomicRanges::strand(exons_filtered), 5))
      )
    }
  }
  
  # Final check that all components have the correct structure
  if (!all(c("seqnames", "index", "score", "file_name", "group") %in% colnames(combined_coverage))) {
    stop("Error: combined_coverage is missing required columns. This may be due to an internal error.")
  }
  
  # Return the combined data
  return(list(combined_coverage = combined_coverage, 
              reps_df = reps_df, 
              genes_df = genes_df, 
              exons_df = exons_df))
}