#' @title SummitHeatmap
#'
#' @description Generate heatmaps of read counts around positions of interest in a genome.
#'
#' @details This function generates heatmaps of different experiments (read counts) around positions of interest in the genome,
#' for example the summits of ChIP peaks, or transcription start sites. The center of your GRanges object defines the coordinates
#' of these positions of interest (eg peak summits, TSSs). Features that lie on the minus strand will be reversed in the final output.
#'
#' @param peaks A GRanges object containing your regions of interest. Must include seqnames (chromosomes), start, end, strand, and name.
#' @param bamFiles Character vector containing the filenames (including the full path) of read alignment files in bam format.
#' @param bamNames Character vector containing the names to describe the \code{bamFiles} you are using (for example: "H3K9me3_reads"). 
#' If no names are supplied, the full \code{bamFiles} names are used.
#' @param span Integer scalar specifyig the distance from the peak center to the left and right that you want your heatmap to extend to. Default is 2025.
#' @param step Integer scalar specifyig the window size in which reads are counted. Default is 50.
#' @param useCPM Logical scalar indicationg wether to normalize the number of reads per window to the total number of reads in the sample (bam file). Default is TRUE.
#' @param PairedEnd Logical scalar, indicating wether reads in bam files were generated with paired-end or single-end sequencing. Default is FALSE (=single-end).
#' @param minMQS Integer scalar, specifying the minimum mapping quality that a read must have to be included. Default is 255, which eliminates multimapping 
#' reads in case the STAR aligner was used to generate the \code{bamFiles}.
#' @param strand Integer vector indicating if strand-specific read counting should be performed. Length of the vector should be either 1 
#' (meaning that the value is applied to all input files), or equal to the total number of input files provided. Each 
#' vector element should have one of the following three values: 0 (unstranded), 1 (stranded) and 2 (reversely stranded). 
#' Default value of this parameter is 0 (ie. unstranded read counting is performed for all input files).
#' @param splitOnly Logical scalar indicating whether only split alignments (their CIGAR strings contain letter 'N') should be included. FALSE by default.
#' @param nonSplitOnly Logical scalar indicating whether only non-split alignments (their CIGAR strings do not contain letter 'N') should be included. FALSE by default.
#' @param minOverlap Integer scalar giving the minimum number of overlapping bases required for assigning a read to a heatmap window. 
#' For assignment of read pairs, number of overlapping bases from each read in the same pair will be summed. 
#' If a negative value is provided, then a gap of up to specified size will be allowed between read and the feature that the read is assigned to. 1 by default.
#' @param readExtension3 Integer scalar giving the number of bases extended downstream from 3' end of each read. 0 by default. Negative value is not allowed.
#' @param readShiftSize Integer scalar specifying the number of bases the reads will be shifted downstream by. 0 by default. Negative value is not allowed. 
#' In case of mode = "Q" and PairedEnd = TRUE, it should be set to 'halfInsert' to shift each read by half the fragment size. 
#' @param requireBothEndsMapped Logical scalar indicating if both ends from the same fragment are required to be successfully aligned before the fragment can be assigned to a feature or meta-feature. 
#' This parameter is only appliable when PairedEnd is TRUE.
#' @param read2pos Specifying whether each read should be reduced to its 5' most base or 3' most base. It has three possible values: NULL, 5 (denoting 5' most base) and 3 (denoting 3' most base).
#' Default value is 5, ie. only the 5' end of the read will be counted. If a read is reduced to a single base, only that base will be considered for the read assignment.
#' Read reduction is performed after read shifting and extension.
#' @param mode Specify if you want to use QuasR (mode = "Q") or FeatureCounts (mode = "F", default) to count the number of rads per window.
#' @param genome The reference genome, either a string specifying a BSgenome object or the name of a genome fasta file. Default is "BSgenome.Mmusculus.UCSC.mm10".
#' 
#' @return A list of matrices of the same length as the number of samples supplied in \code{bamFiles}.
#' Each matrix contains the number of reads (or cpm) in each window (column headers indicate the middle of the window)
#' around the center of each region provided in \code{peaks}.
#'
#' @examples
#' peaks <- SimulatePeaks(1000,rep(100,1000),chromosomeSizes=
#' system.file("extdata", "chrNameLength_mm10_chr11.txt", package = "MiniChip"))
#' #peaks <- GenomicRanges::GRanges(
#' #seqnames = Rle(c("chr1", "chr2", "chr1", "chr3"), c(1, 3, 2, 4)),
#' #ranges = IRanges(50101:50110, end = 51111:51120),
#' #strand = Rle(strand(c("-", "+", "*", "+", "-")), c(1, 2, 2, 3, 2)),
#' #summit = 1:10, name = head(letters, 10))
#' bamFiles <- list.files(system.file("extdata", package = "MiniChip"),
#'  full.names=TRUE,pattern="*bam$")
#' SummitHeatmap(peaks=peaks,bamFiles=bamFiles)
#'
#' @importFrom GenomicRanges start
#' @importFrom GenomicRanges end
#' @importFrom GenomicRanges resize
#' @importFrom GenomicRanges width
#' @importFrom GenomicRanges slidingWindows
#' @importFrom GenomicRanges strand
#' @importFrom GenomicRanges strand<-
#' @importFrom GenomicRanges elementMetadata
#' @importFrom GenomicRanges seqnames
#' @importFrom GenomicRanges GRanges
#' @importFrom Rsubread featureCounts
#' @importFrom parallel makeCluster
#' @importFrom QuasR qAlign
#' @importFrom QuasR qProfile
#' @importFrom QuasR alignmentStats
#'
#' @export
SummitHeatmap <- function(peaks, bamFiles, bamNames="myreads", span=2025, step=50, minOverlap = 1,
                          useCPM=TRUE, PairedEnd=FALSE, minMQS=255, strand=0, splitOnly=FALSE, nonSplitOnly=FALSE,
                          readExtension3=0, readShiftSize=0, requireBothEndsMapped=FALSE, read2pos=5, mode="F",
                          genome="BSgenome.Mmusculus.UCSC.mm10"){
  
  # Input validation
  
  # Check if peaks is a GRanges object
  if (!inherits(peaks, "GRanges")) {
    stop("'peaks' must be a GRanges object")
  }
  
  # Check if bamFiles is provided and is a character vector
  if (missing(bamFiles)) {
    stop("'bamFiles' parameter is required")
  }
  if (!is.character(bamFiles)) {
    stop("'bamFiles' must be a character vector containing file paths")
  }
  
  # Check if all BAM files exist
  non_existing_files <- bamFiles[!file.exists(bamFiles)]
  if (length(non_existing_files) > 0) {
    stop(paste0("The following BAM files don't exist: ", 
                paste(non_existing_files, collapse=", ")))
  }
  
  # Check if bamNames has the same length as bamFiles
  if (length(bamNames) != 1 && length(bamNames) != length(bamFiles)) {
    stop(paste0("'bamNames' must have the same length as 'bamFiles' or be a single value. ",
                "Length of bamFiles: ", length(bamFiles), 
                ", length of bamNames: ", length(bamNames)))
  }
  
  # If bamNames is a single value and bamFiles has multiple files, replicate bamNames
  if (length(bamNames) == 1 && length(bamFiles) > 1) {
    bamNames <- rep(bamNames, length(bamFiles))
  }
  
  # Check numeric parameters
  if (!is.numeric(span) || span <= 0) {
    stop("'span' must be a positive numeric value")
  }
  
  if (!is.numeric(step) || step <= 0) {
    stop("'step' must be a positive numeric value")
  }
  
  if (!is.numeric(minOverlap)) {
    stop("'minOverlap' must be a numeric value")
  }
  
  if (!is.numeric(minMQS) || minMQS < 0) {
    stop("'minMQS' must be a non-negative numeric value")
  }
  
  if (!is.numeric(readExtension3) || readExtension3 < 0) {
    stop("'readExtension3' must be a non-negative numeric value")
  }
  
  if (!is.numeric(readShiftSize) || readShiftSize < 0) {
    if (!(mode == "Q" && PairedEnd && readShiftSize == "halfInsert")) {
      stop("'readShiftSize' must be a non-negative numeric value or 'halfInsert' (only when mode='Q' and PairedEnd=TRUE)")
    }
  }
  
  # Check boolean parameters
  if (!is.logical(useCPM)) {
    stop("'useCPM' must be a logical value (TRUE or FALSE)")
  }
  
  if (!is.logical(PairedEnd)) {
    stop("'PairedEnd' must be a logical value (TRUE or FALSE)")
  }
  
  if (!is.logical(splitOnly)) {
    stop("'splitOnly' must be a logical value (TRUE or FALSE)")
  }
  
  if (!is.logical(nonSplitOnly)) {
    stop("'nonSplitOnly' must be a logical value (TRUE or FALSE)")
  }
  
  if (!is.logical(requireBothEndsMapped)) {
    stop("'requireBothEndsMapped' must be a logical value (TRUE or FALSE)")
  }
  
  # Check if splitOnly and nonSplitOnly are not both TRUE
  if (splitOnly && nonSplitOnly) {
    stop("'splitOnly' and 'nonSplitOnly' cannot both be TRUE")
  }
  
  # Check requireBothEndsMapped only applicable when PairedEnd is TRUE
  if (requireBothEndsMapped && !PairedEnd) {
    warning("'requireBothEndsMapped' is only applicable when PairedEnd=TRUE. Setting will be ignored.")
  }
  
  # Check strand parameter
  if (!is.numeric(strand) || !all(strand %in% c(0, 1, 2))) {
    stop("'strand' must be one of the following values: 0 (unstranded), 1 (stranded), or 2 (reversely stranded)")
  }
  
  if (length(strand) != 1 && length(strand) != length(bamFiles)) {
    stop(paste0("Length of 'strand' must be either 1 or equal to the number of BAM files. ",
                "Current length: ", length(strand), ", number of BAM files: ", length(bamFiles)))
  }
  
  # Check read2pos parameter
  if (!is.null(read2pos) && !(read2pos %in% c(5, 3))) {
    stop("'read2pos' must be NULL, 5, or 3")
  }
  
  # Check mode parameter
  if (!(mode %in% c("F", "Q"))) {
    stop("'mode' must be either 'F' (FeatureCounts) or 'Q' (QuasR)")
  }
  
  # Check if QuasR is installed if mode is "Q"
  if (mode == "Q" && !requireNamespace("QuasR", quietly = TRUE)) {
    stop("'QuasR' package is required when mode='Q'. Please install it with BiocManager::install('QuasR')")
  }
  
  # Check if Rsubread is installed if mode is "F"
  if (mode == "F" && !requireNamespace("Rsubread", quietly = TRUE)) {
    stop("'Rsubread' package is required when mode='F'. Please install it with BiocManager::install('Rsubread')")
  }
  
  # Now proceed with the original function logic
  if (is.null(names(peaks)) && !is.null(peaks$name)) {
    names(peaks) <- peaks$name
  }
  
  if (is.null(strand(peaks)) && !is.null(peaks$strand)) {
    strand(peaks) <- peaks$strand
  }
  
  # Check if names are available for peaks
  if (is.null(names(peaks))) {
    warning("No names found for peaks. Using automatically generated names.")
    names(peaks) <- paste0("peak_", seq_along(peaks))
  }
  
  # Check if strand information is available
  if (all(is.na(strand(peaks)))) {
    warning("No strand information found for peaks. Setting all strands to '*'.")
    strand(peaks) <- "*"
  }
  
  if (mode == "Q") { # use QUASR version 
    #write a table to read in samples for QUASR
    write.table(data.frame(FileName=bamFiles, SampleName=bamNames), file="QUASR.txt", sep="\t", col.names=TRUE, row.names=FALSE, append=FALSE, quote=FALSE)
    
    #translate options
    cl <- tryCatch({
      makeCluster(20)
    }, error = function(e) {
      warning("Could not create cluster with 20 cores. Using 1 core instead.")
      makeCluster(1)
    })
    
    Qpaired <- ifelse(PairedEnd == TRUE, "fr", "no")
    selectReadPosition <- ifelse(read2pos == 3, "end", "start")
    orientation <- ifelse(strand == 1, "same", ifelse(strand == 2, "opposite", "any"))
    
    #generate project
    tryCatch({
      proj <- qAlign("QUASR.txt", genome, paired = Qpaired, clObj = cl)
    }, error = function(e) {
      stop(paste0("Error in qAlign: ", conditionMessage(e), 
                  "\nMake sure 'genome' parameter is valid and QuasR is properly configured."))
    })
    
    #generate counts matrices
    tryCatch({
      prf <- qProfile(proj, resize(peaks, width = 1L, fix = "center"), upstream=span, binSize = step,
                      selectReadPosition= selectReadPosition, orientation = orientation, shift = readShiftSize, 
                      useRead="any", clObj = cl, mapqMin = minMQS)[-1]
    }, error = function(e) {
      stop(paste0("Error in qProfile: ", conditionMessage(e)))
    })
    
    if (useCPM == TRUE) { 
      #normalize to million reads (cpm)
      normfactor <- alignmentStats(proj)[,"mapped"]/1000000
      prfcpm <- prf
      for (i in seq_along(prf)) {
        prfcpm[[i]] <- prf[[i]]/normfactor[i]
      }
      unlink("QUASR.txt")
      return(prfcpm)
    } else {
      unlink("QUASR.txt")
      return(prf)
    }
  }
  
  if (mode == "F") { #use Feature Counts version 
    
    #take the middle of the GRanges region, then define whole heatmap region
    nwindows <- ceiling((span*2)/step)
    if (nwindows %% 2 == 0) {
      nwindows <- nwindows + 1
    }
    regionwidth <- step * nwindows
    
    # Try to resize peaks and catch potential errors
    tryCatch({
      peaks <- resize(peaks, width=regionwidth, fix="center")
    }, error = function(e) {
      stop(paste0("Error resizing peaks: ", conditionMessage(e), 
                  "\nMake sure your peaks GRanges object is properly constructed."))
    })
    
    #remove peaks with negative start values
    original_peak_count <- length(peaks)
    peaks <- peaks[start(peaks) >= 0 & width(peaks) == regionwidth]
    if (length(peaks) < original_peak_count) {
      warning(paste0(original_peak_count - length(peaks), 
                     " peaks were removed because they had negative start positions or incorrect width."))
    }
    
    if (length(peaks) == 0) {
      stop("No valid peaks left after filtering. Please check your input peaks.")
    }
    
    #generate window starts and ends across span
    windows <- seq(from=0, to=regionwidth-step, by=step)
    binmids <- windows - regionwidth/2 + step/2
    
    #generate sliding windows across peaks 
    tryCatch({
      peaks2plus <- unlist(GenomicRanges::slidingWindows(peaks[strand(peaks) != "-"], width = step, step = step), use.names=FALSE)
      peaks2minus <- rev(unlist(GenomicRanges::slidingWindows(peaks[strand(peaks) == "-"], width = step, step = step), use.names=FALSE))
      peaks2 <- c(peaks2plus, peaks2minus)
    }, error = function(e) {
      stop(paste0("Error generating sliding windows: ", conditionMessage(e)))
    })
    
    #generate saf format data frame
    saf <- data.frame(GeneID = names(peaks2), Chr = seqnames(peaks2),
                      Start = start(peaks2), End = end(peaks2), Strand = strand(peaks2))
    
    #calculate the number of reads in each window for all bam files specified
    tryCatch({
      f_counts <- featureCounts(bamFiles, annot.ext=saf, useMetaFeatures=FALSE, allowMultiOverlap=TRUE, read2pos=read2pos,
                                minOverlap=minOverlap, readExtension3=readExtension3, countMultiMappingReads=FALSE, fraction=FALSE,
                                minMQS=minMQS, strandSpecific=strand, nthreads=20, verbose=FALSE, isPairedEnd=PairedEnd,
                                splitOnly=splitOnly, nonSplitOnly=nonSplitOnly, readShiftType="downstream", readShiftSize=readShiftSize, 
                                requireBothEndsMapped=requireBothEndsMapped)
    }, error = function(e) {
      stop(paste0("Error in featureCounts: ", conditionMessage(e), 
                  "\nMake sure your BAM files are properly formatted and accessible."))
    })
    
    #extract gene annotation (for peak names)
    GeneID <- matrix(names(peaks2), nrow=length(peaks), ncol=length(windows), byrow=TRUE)[,1]
    
    #extract number of mapped reads for all bam samples
    if (length(bamNames) > 1) {
      mapped.reads <- apply(f_counts$stat[c(1,11,12), -1], 2, sum)
    } else {
      mapped.reads <- sum(f_counts$stat[c(1,11,12), -1])
    }
    
    #extract counts for all samples, name columns and rows, calculate log CPM
    all.counts <- vector("list", length(bamNames))
    
    for (bam.sample in seq_along(bamNames)) {
      counts <- matrix(f_counts$counts[, bam.sample], nrow=length(peaks), ncol=length(windows), byrow=TRUE)
      
      colnames(counts) <- binmids
      rownames(counts) <- GeneID
      
      #prepare counts or cpm for saving as list of matrices
      if (useCPM == TRUE) {
        all.counts[[bam.sample]] <- ((counts)/mapped.reads[bam.sample])*1000000
      } else {
        all.counts[[bam.sample]] <- counts
      }
      
      # Check if there are names to sort by
      if (!is.null(names(peaks)) && length(names(peaks)) == nrow(counts)) {
        all.counts[[bam.sample]] <- all.counts[[bam.sample]][names(peaks), ]
      }
    }
    
    # give warning if all counts are 0
    if (max(all.counts[[1]]) == 0) {
      warning("All windows have 0 counts, make sure the chromosome names (seqnames) match between your GRanges object and bam files!")
    }
    
    #return the results
    names(all.counts) <- bamNames
    return(all.counts)
    
  } else {
    stop("mode must be one of Q (QuasR) or F (FeatureCounts)!")
  }
}