#' @title calcTracks
#'
#' @description Prepare ChIP, repeat, and gene annotation data for plotting tracks 
#' for a given region of interest (ROI).
#'
#' @details calcTracks is used to prepare ChIP, repeat, and gene annotation data for plotting tracks 
#' for a given region of interest (ROI).
#'
#' @param bigwigFiles Character vector containing the filenames (including the full path) of read alignment files in bigwig format.
#' @param bigwigNames  Character vector containing the names to describe the \code{bigwigFiles} you are using (for example: "H3K9me3_replicate1"). 
#' @param bigwigGroupNames  Character vector containing the group names to describe the \code{bigwigFiles} you are using (for example: "H3K9me3"). 
#' @param ROI A GRanges object containing your single region of interest that should be plotted. Must include seqnames (chromosomes), start, end.
#' @param CoverageSmoothing Logical scalar, indicating whether the coverage tracks should be smoothed.
#' @param smoothing_window Numeric scalar providing the window size for smoothing.
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
#'
#' @export
calcTracks <- function(bigwigFiles, bigwigNames, bigwigGroupNames, ROI, CoverageSmoothing=TRUE, smoothing_window=20,
                           reps,genes,exons){ 
  
  # Initialize a list to store data for all BigWig files
  coverage_data_list <- list()
  
  # Loop through each BigWig file
  for (i in seq_along(bigwigFiles)) {
    
    # Import the BigWig file for the given region (reps2_sub1 assumed to be defined elsewhere)
    bwf <- rtracklayer::import.bw(bigwigFiles[i], which = ROI)
    
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
  }
  
  
  # Combine all the coverage data into a single data frame
  combined_coverage <- do.call(rbind, coverage_data_list)
  
  
  # Define a smoothing function
  smooth_coverage <- function(data, smoothing_window = 20) {
    data %>%
      group_by(file_name) %>%
      mutate(score = zoo::rollmean(score, k = smoothing_window, fill = NA, align = "center")) %>%
      ungroup()
  }
  
  if (CoverageSmoothing==TRUE){
    # Apply smoothing to the combined coverage data
    combined_coverage <- smooth_coverage(combined_coverage,smoothing_window=smoothing_window)
  }
  

  # Filter annotations to the combined index range
  
  ### generate window GRanges
  score_window <- GRanges(
    seqnames = unique(combined_coverage$seqnames),
    ranges = IRanges(start = min(combined_coverage$index), end = max(combined_coverage$index))
  )
  
  #filter repeats
  reps_filtered <- IRanges::subsetByOverlaps(reps, score_window, ignore.strand=TRUE)
  
  #truncate start and end coordinates to fit into plotting window
  start(reps_filtered) <- ifelse(start(reps_filtered) >=  start(score_window),start(reps_filtered),start(score_window))
  end(reps_filtered) <- ifelse(end(reps_filtered) <=  end(score_window),end(reps_filtered),end(score_window))
  
  # Convert filtered reps2 to a data frame for plotting
  if(length(reps_filtered) > 0){
    
    reps_df <- data.frame(
      start = start(GenomicRanges::promoters(reps_filtered,upstream=0,downstream=1)),
      end = end(GenomicRanges::terminators(reps_filtered,upstream=0,downstream=1)),
      repeat_name = mcols(reps_filtered)$repeat_name,
      ypos= as.numeric(paste0(strand(reps_filtered),5))
    )
  } else {
    reps_df <- data.frame(
      start = 0,
      end = 0,
      repeat_name = "",
      ypos= 5)
  }
  
  #filter genes
  genes_filtered <- IRanges::subsetByOverlaps(genes, score_window, ignore.strand=TRUE)
  exons_filtered <- IRanges::subsetByOverlaps(exons, score_window, ignore.strand=TRUE)
  
  #truncate start and end coordinates to fit into plotting window
  start(genes_filtered) <- ifelse(start(genes_filtered) >=  start(score_window),start(genes_filtered),start(score_window))
  end(genes_filtered) <- ifelse(end(genes_filtered) <=  end(score_window),end(genes_filtered),end(score_window))
  start(exons_filtered) <- ifelse(start(exons_filtered) >=  start(score_window),start(exons_filtered),start(score_window))
  end(exons_filtered) <- ifelse(end(exons_filtered) <=  end(score_window),end(exons_filtered),end(score_window))
  
  # Convert filtered genes to a data frame for plotting
  if(length(genes_filtered) > 0){
    genes_df <- data.frame(
      start = start(GenomicRanges::promoters(genes_filtered,upstream=0,downstream=1)),
      end = end(GenomicRanges::terminators(genes_filtered,upstream=0,downstream=1)),
      gene_name = mcols(genes_filtered)$gene_name,
      ypos= as.numeric(paste0(strand(genes_filtered),5))
    )
  } else {
    genes_df <- data.frame(
      start = 0,
      end = 0,
      gene_name = "",
      ypos= 5)
  }
  
  if(length(exons_filtered) > 0){
    exons_df <- data.frame(
      start = start(GenomicRanges::promoters(exons_filtered,upstream=0,downstream=1)),
      end = end(GenomicRanges::terminators(exons_filtered,upstream=0,downstream=1)),
      ypos= as.numeric(paste0(strand(exons_filtered),5))
    )
  } else {
    exons_df <- data.frame(
      start = 0,
      end = 0,
      ypos= 5)
    
  }
  
  return(list(combined_coverage=combined_coverage,reps_df=reps_df,genes_df=genes_df,exons_df=exons_df))
}