#' @title GCbias
#'
#' @description Plot GC content versus read counts.
#'
#' @details This function generates a scatter plot of the number of Gs and Cs on the x-axis and the read count (cpm) on the y-axis in windows
#' of size \code{winWidth} bp across the genome. A seperate plot is generated for each read alignment file in \code{bamFiles}. 
#' Supports both single and paired-end experiments. These plots allow the user to check if there is a potential GCbias in the (ChIPseq) data. 
#'
#' @param bamFiles Character vector containing the filenames filenames (including the full path) of read alignment files in bam format.
#' @param bamNames Character vector containing the names to describe the \code{bamFiles} you are using (for example: "H3K9me3_reads"). If no names are supplied, the full \code{bamFiles} names are used.
#' @param minMQS Integer scalar, specifying the minimum mapping quality that a read must have to be included. Default is 255, which eliminates multimapping reads in case the STAR aligner was used to generate the \code{bamFiles}.
#' @param maxFrag Integer scalar, specifying the maximum fragment length corresponding to a read pair. Defaults to 500 base pairs.
#' @param pe Character scalar indicating whether paired-end data is present; set to "none" (the default), "both", "first" or "second".
#' @param restrict Character vector containing the names of allowable chromosomes from which reads will be extracted. Default is "chr11".
#' @param winWidth Integer scalar specifying the width of the window, in which reads are counted and GC content calculated. Default is 5000 base pairs.
#' @param col Color scheme for the smooth scatter plots. If not provided, viridis::inferno is used. 
#' @param genome BSGenome object. Required parameter. For example,use BSgenome.Mmusculus.UCSC.mm10 for mouse.
#' @param GCprob Logical scalar, indicating whether the GC content should be displayed as absolute counts (GCprob=FALSE) or as fraction of GCs (GCprob=TRUE,default).
#' @param span Numeric scalar specifying the span that is used for loess trendline. Default= 0.1
#' @param plot If TRUE, the output will be plotted, otherwise the matrix to generate the plots will be returned.
#' @param logCPM Logical, should the cpm be reported on log scale? 
#' @param priorCount Prior Count for calculating cpm.
#'
#' @return This function generates a scatter plot of the number of Gs and Cs on the x-axis and the read count (cpm) on the y-axis in windows
#' of size \code{winWidth} bp across the genome. A loess trendline is added to allow the user to see a potential GCbias trend in the data provided.  
#' 
#'
#' @examples
#' library(BSgenome.Mmusculus.UCSC.mm10)
#' bamFiles <- list.files(system.file("extdata", package = "MiniChip"),
#'  full.names=TRUE,pattern="*bam$")[1:2]
#' bamNames <- gsub(paste(system.file("extdata", package = "MiniChip"),
#' "/",sep=""),"",bamFiles)
#' bamNames <- gsub("_chr11.bam","",bamNames)
#' GCbias(bamFiles=bamFiles,bamNames=bamNames,
#' genome=BSgenome.Mmusculus.UCSC.mm10)
#'
#' @importFrom csaw readParam calculateCPM windowCounts calculateCPM
#' @importFrom SummarizedExperiment assay rowRanges
#' @importFrom BSgenome getSeq
#' @importFrom Biostrings letterFrequency
#' @importFrom viridis inferno
#' @importFrom graphics smoothScatter lines par
#' @importFrom stats loess predict
#' @importFrom grDevices colorRampPalette
#'
#' @export
GCbias <- function(bamFiles, bamNames=bamFiles, minMQS=255, maxFrag=500, pe="none", restrict="chr11",
                   winWidth=5000, col=inferno, genome, GCprob=TRUE, span=0.1, plot=TRUE, 
                   logCPM=TRUE, priorCount=1) {
  
  # Check required packages
  required_packages <- c("csaw", "SummarizedExperiment", "BSgenome", "Biostrings",
                         "viridis", "graphics", "stats")
  missing_packages <- character(0)
  
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      missing_packages <- c(missing_packages, pkg)
    }
  }
  
  if (length(missing_packages) > 0) {
    pkg_list <- paste(missing_packages, collapse = ", ")
    bioc_pkgs <- intersect(missing_packages, c("csaw", "SummarizedExperiment", "BSgenome", "Biostrings"))
    if (length(bioc_pkgs) > 0) {
      stop("The following Bioconductor packages are required but not installed: ", 
           paste(bioc_pkgs, collapse = ", "), 
           ". Please install them using BiocManager::install(c('", 
           paste(bioc_pkgs, collapse = "', '"), "')).")
    }
    cran_pkgs <- setdiff(missing_packages, bioc_pkgs)
    if (length(cran_pkgs) > 0) {
      stop("The following CRAN packages are required but not installed: ", 
           paste(cran_pkgs, collapse = ", "), 
           ". Please install them using install.packages(c('", 
           paste(cran_pkgs, collapse = "', '"), "')).")
    }
  }
  
  # Check for required parameters
  if (missing(bamFiles)) {
    stop("'bamFiles' parameter is required but missing.")
  }
  if (missing(genome)) {
    stop("'genome' parameter is required but missing.")
  }
  
  # Validate bamFiles parameter
  if (!is.character(bamFiles)) {
    stop("'bamFiles' must be a character vector containing file paths.")
  }
  if (length(bamFiles) == 0) {
    stop("'bamFiles' cannot be empty. Please provide at least one BAM file path.")
  }
  
  # Check if BAM files exist
  missing_files <- bamFiles[!file.exists(bamFiles)]
  if (length(missing_files) > 0) {
    stop("The following BAM files do not exist: ", paste(missing_files, collapse = ", "))
  }
  
  # Check if BAM files are readable
  unreadable_files <- character(0)
  for (file in bamFiles) {
    if (!file.access(file, mode = 4) == 0) {
      unreadable_files <- c(unreadable_files, file)
    }
  }
  if (length(unreadable_files) > 0) {
    stop("The following BAM files are not readable: ", paste(unreadable_files, collapse = ", "))
  }
  
  # Validate bamNames parameter
  if (!is.character(bamNames)) {
    stop("'bamNames' must be a character vector.")
  }
  if (length(bamNames) != length(bamFiles)) {
    warning("Length of 'bamNames' (", length(bamNames), 
            ") differs from length of 'bamFiles' (", length(bamFiles), 
            "). Using 'bamFiles' as names instead.")
    bamNames <- bamFiles
  }
  
  # Check for duplicate names
  if (length(unique(bamNames)) != length(bamNames)) {
    warning("Duplicate names detected in 'bamNames'. This may cause issues with labeling.")
  }
  
  # Validate numeric parameters
  if (!is.numeric(minMQS) || length(minMQS) != 1 || minMQS < 0) {
    stop("'minMQS' must be a non-negative numeric scalar.")
  }
  
  if (!is.numeric(maxFrag) || length(maxFrag) != 1 || maxFrag <= 0) {
    stop("'maxFrag' must be a positive numeric scalar.")
  }
  
  if (!is.numeric(winWidth) || length(winWidth) != 1 || winWidth <= 0) {
    stop("'winWidth' must be a positive numeric scalar.")
  }
  
  if (!is.numeric(span) || length(span) != 1 || span <= 0 || span > 1) {
    stop("'span' must be a numeric scalar between 0 and 1.")
  }
  
  if (!is.numeric(priorCount) || length(priorCount) != 1 || priorCount < 0) {
    stop("'priorCount' must be a non-negative numeric scalar.")
  }
  
  # Validate logical parameters
  if (!is.logical(GCprob) || length(GCprob) != 1) {
    stop("'GCprob' must be a logical scalar (TRUE or FALSE).")
  }
  
  if (!is.logical(plot) || length(plot) != 1) {
    stop("'plot' must be a logical scalar (TRUE or FALSE).")
  }
  
  if (!is.logical(logCPM) || length(logCPM) != 1) {
    stop("'logCPM' must be a logical scalar (TRUE or FALSE).")
  }
  
  # Validate pe parameter
  valid_pe_values <- c("none", "both", "first", "second")
  if (!is.character(pe) || length(pe) != 1 || !(pe %in% valid_pe_values)) {
    stop("'pe' must be one of the following: ", paste(valid_pe_values, collapse = ", "))
  }
  
  # Validate restrict parameter
  if (!is.character(restrict)) {
    stop("'restrict' must be a character vector.")
  }
  
  # Validate genome parameter
  if (!inherits(genome, "BSgenome")) {
    stop("'genome' must be a BSgenome object. Please provide a valid BSgenome object.")
  }
  
  # Validate col parameter
  if (!is.function(col)) {
    if (exists("inferno", mode = "function")) {
      col <- inferno
    } else {
      # If viridis is installed but inferno is not accessible directly
      if (requireNamespace("viridis", quietly = TRUE)) {
        col <- viridis::inferno
      } else {
        # Fallback to default R color palette
        col <- colorRampPalette(c("blue", "purple", "red", "orange", "yellow"))
        warning("'col' should be a color ramp function. Using default color palette instead.")
      }
    }
  }
  
  # Start the analysis with error handling
  
  # Define parameters for read extraction
  tryCatch({
    paramPE <- csaw::readParam(minq = minMQS, max.frag = maxFrag, pe = pe, restrict = restrict)
  }, error = function(e) {
    stop("Error setting read parameters: ", e$message, 
         "\nCheck that your parameter values are valid.")
  })
  
  # Count number of reads in bins
  tryCatch({
    chip_bins <- csaw::windowCounts(bamFiles, width = winWidth, param = paramPE, bin = TRUE, filter = 0)
    
    # Check if any data was read
    if (nrow(SummarizedExperiment::assay(chip_bins)) == 0) {
      stop("No reads were counted in the specified windows. Check your BAM files and parameters.")
    }
    
  }, error = function(e) {
    stop("Error counting reads: ", e$message, 
         "\nThis could be due to issues with your BAM files or invalid parameter values.")
  })
  
  # Extract bin data
  tryCatch({
    bins_data <- data.frame(SummarizedExperiment::assay(chip_bins))
    colnames(bins_data) <- bamNames
    
    # Calculate CPM
    bins_cpm <- data.frame(csaw::calculateCPM(chip_bins, use.norm.factors = TRUE, use.offsets = FALSE,
                                              log = logCPM, prior.count = priorCount, assay.id = "counts"))
    colnames(bins_cpm) <- bamNames
    
  }, error = function(e) {
    stop("Error processing bin data: ", e$message, 
         "\nCheck that your input data is valid.")
  })
  
  # Extract window coordinates and sequences
  tryCatch({
    bins_ranges <- SummarizedExperiment::rowRanges(chip_bins)
    
    # Check if all chromosomes in bins_ranges exist in the genome
    missing_chromosomes <- setdiff(unique(as.character(seqnames(bins_ranges))), 
                                   seqnames(genome))
    
    if (length(missing_chromosomes) > 0) {
      warning("The following chromosomes from your BAM files are not present in the provided genome: ", 
              paste(missing_chromosomes, collapse = ", "), 
              ". These will be excluded from the analysis.")
      
      # Remove bins from missing chromosomes
      bins_ranges <- bins_ranges[!as.character(seqnames(bins_ranges)) %in% missing_chromosomes]
      
      if (length(bins_ranges) == 0) {
        stop("No valid chromosomes remaining after filtering. Ensure your BAM files match the provided genome.")
      }
    }
    
    # Extract sequences
    bins_seqs <- BSgenome::getSeq(genome, bins_ranges)
    
    # Calculate GC content
    bins_cpm$bins_gc <- Biostrings::letterFrequency(bins_seqs, "GC", as.prob = GCprob)
    
  }, error = function(e) {
    stop("Error calculating GC content: ", e$message, 
         "\nEnsure that your genome object is valid and contains the chromosomes in your BAM files.")
  })
  
  # Set x-axis label based on GCprob
  if (GCprob == FALSE) {
    xaxislab <- "number of Gs and Cs"
  } else {
    xaxislab <- "fraction GC"
  }
  
  # Order table by GC content
  bins_cpm <- bins_cpm[order(bins_cpm$bins_gc), ]
  
  # Return data if plot=FALSE
  if (plot == FALSE) {
    return(bins_cpm)
  } else {
    # Plot data with error handling
    tryCatch({
      # Set up plotting area
      opar <- par(no.readonly = TRUE)  # Save original par settings
      on.exit(par(opar))  # Restore original par settings when function exits
      
      par(mfrow = c(ceiling(length(bamNames) / 4), min(4, length(bamNames))))
      
      # Create plots
      for (i in seq_along(bamNames)) {
        smoothScatter(bins_cpm$bins_gc, bins_cpm[, i], colramp = col, 
                      xlab = xaxislab, ylab = if(logCPM) "log CPM" else "CPM", 
                      main = bamNames[i])
        
        # Add trendline
        if (nrow(bins_cpm) >= 4) {  # Need at least 4 points for loess
          loessGC <- loess(bins_cpm[, i] ~ bins_cpm$bins_gc, span = span)
          lines(bins_cpm$bins_gc, predict(loessGC), lwd = 2, col = "white")
        } else {
          warning("Too few points to fit loess curve for ", bamNames[i])
        }
      }
      
      invisible(bins_cpm)  # Return the data invisibly
      
    }, error = function(e) {
      # Restore original par settings in case of error
      par(opar)
      stop("Error generating plots: ", e$message)
    })
  }
}