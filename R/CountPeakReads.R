#' @title CountPeakReads
#'
#' @description Generate normalized counts in genomic regions of interest (for example, ChIPseq peaks).
#'
#' @details Takes regions of interest in the genome (\code{peaks}), for example ChIP peaks, or transcription start sites. 
#' Then the number of reads in each one of the prvovided \code{bamFiles} is counted. The read counts in a ChIP sample (provided in \code{chips}) 
#' are normalized to the read counts in the corresponding Input sample (provided in the same position in \code{inputs} as the corresponding 
#' ChIP sample in \code{chips}), if provided, taking the total library sizes into account. 
#' If no Input sample names are provided, the counts per milliion reads (cpm) for each ChIP sample is returned. 
#'
#' @param peaks A GRanges object containing your regions of interest. Must include seqnames (chromosomes), start, end, strand, and name.
#' @param bamFiles Character vector containing the filenames (including the full path) of read alignment files in bam format.
#' @param bamNames Character vector containing the names to describe the \code{bamFiles} you are using (for example: "H3K9me3_reads"). 
#' If no names are supplied, the full \code{bamFiles} names are used.
#' @param chips Character vector containing the names of the ChIP samples, corresponding exactly to the ChIP sample names in \code{bamNames}.
#' If no names are supplied, the full \code{bamNames} names and all samples listed in \code{bamNames} are used
#' as ChIP experiments without Input normalization.
#' @param inputs Character vector containing the names of the Input samples corresponding exactly to the Input sample names in \code{bamNames}, 
#' provided in an order crresponding to the ChIP samples (This is because the function pairs the ChIP and Input samples based on the order in
#' the \code{chips} and \code{inputs} vectors.) If missing, no Input normalization is performed.
#' @param width Integer scalar providing the width around the center of the peaks GRanges object provided in which the reads will be counted,default = 0. 
#' If width=0, the provided regions will remain the original size. 
#' @param pseudocount Numeric scalar providing the pseudocount that is added to the read counts of Chip and Input before normalization. Default is 8.
#' @param PairedEnd Logical scalar, indicating wether reads in bam files were generated with paired-end or single-end sequencing. Default is FALSE (=single-end).
#' @param minMQS Integer scalar, specifying the minimum mapping quality that a read must have to be included. Default is 255, which eliminates multimapping 
#' reads in case the STAR aligner was used to generate the \code{bamFiles}.
#' @param strand Integer vector indicating if strand-specific read counting should be performed. Length of the vector should be either 1 
#' (meaning that the value is applied to all input files), or equal to the total number of input files provided. Each 
#' vector element should have one of the following three values: 0 (unstranded), 1 (stranded) and 2 (reversely stranded). 
#' Default value of this parameter is 0 (ie. unstranded read counting is performed for all input files).
#' @param splitOnly Logical scalar indicating whether only split alignments (their CIGAR strings contain letter 'N') should be included. FALSE by default.
#' @param nonSplitOnly Logical scalar indicating whether only non-split alignments (their CIGAR strings do not contain letter 'N') should be included. FALSE by default.
#' @param minOverlap Integer scalar giving the minimum number of overlapping bases required for assigning a read to a genomic region (provided in \code{peaks}). 
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
#' @param mode Specify if you want to use QuasR (mode = "Q") or FeatureCounts (mode = "F", default) to count the number of reads per region
#' @param genome The reference genome, either a string specifying a BSgenome object or the name of a genome fasta file. Default is "BSgenome.Mmusculus.UCSC.mm10".

#'
#' @return A list of 2 matrices. Both matrices contain a column for each chip sample provided, and the same number of rows as regions provided in the GRanges object.
#' The first matrix contains the (Input) normalized number of reads in each region (log2 of Chip + \code{pseudocount}/Input + \code{pseudocount}). 
#' If no Input samples are provided, then the first matrix contains the counts per million (cpm) values for the ChIP samples.
#' The second matrix contains the total number of reads in the the chip sample, or in the chip plus the input sample, if an input sample is provided. 
#' 
#'
#' @examples
#' peaks <- SimulatePeaks(1000,rep(100,1000),chromosomeSizes=
#' system.file("extdata", "chrNameLength_mm10_chr11.txt", package = "MiniChip"))
#' bamFiles <- list.files(system.file("extdata", 
#' package = "MiniChip"), full.names=TRUE,pattern="*bam$")
#' bamNames <- gsub(paste(system.file("extdata", 
#' package = "MiniChip"),"/",sep=""),"",bamFiles)
#' bamNames <- gsub("_chr11.bam","",bamNames)
#' CountPeakReads(peaks=peaks,bamFiles=bamFiles,
#' bamNames=bamNames,chips=bamNames[1:3],inputs=bamNames[4:6])
#'
#' @importFrom GenomicRanges start
#' @importFrom GenomicRanges end
#' @importFrom GenomicRanges resize
#' @importFrom GenomicRanges strand
#' @importFrom GenomicRanges strand<-
#' @importFrom GenomicRanges elementMetadata
#' @importFrom GenomicRanges seqnames
#' @importFrom GenomicRanges GRanges
#' @importFrom Rsubread featureCounts
#' @importFrom parallel makeCluster
#' @importFrom QuasR qAlign
#' @importFrom QuasR qCount
#' @importFrom QuasR alignmentStats
#'
#' @export
CountPeakReads <- function(peaks,bamFiles,bamNames=bamFiles,chips=bamNames,inputs,width=500,minOverlap = 1,
                           PairedEnd=FALSE, minMQS=255,strand=0, splitOnly=FALSE, nonSplitOnly=FALSE,
                           readExtension3=0,readShiftSize=0,requireBothEndsMapped=FALSE,read2pos=5,pseudocount=8,mode="F",
                           genome="BSgenome.Mmusculus.UCSC.mm10"){
  
  # Input validation section
  #---------------------------------------------------------------------------------------
  # Check peaks
  if(missing(peaks)) {
    stop("'peaks' argument is missing. Please provide a GRanges object with your regions of interest.")
  }
  if(!inherits(peaks, "GRanges")) {
    stop("'peaks' must be a GRanges object. Current class: ", class(peaks))
  }
  if(length(peaks) == 0) {
    stop("'peaks' GRanges object is empty. Please provide regions of interest.")
  }
  
  # Check bamFiles
  if(missing(bamFiles)) {
    stop("'bamFiles' argument is missing. Please provide paths to your BAM files.")
  }
  if(!is.character(bamFiles)) {
    stop("'bamFiles' must be a character vector. Current class: ", class(bamFiles))
  }
  if(length(bamFiles) == 0) {
    stop("'bamFiles' is empty. Please provide at least one BAM file path.")
  }
  for(i in seq_along(bamFiles)) {
    if(!file.exists(bamFiles[i])) {
      stop("BAM file does not exist: ", bamFiles[i])
    }
    if(!grepl("\\.bam$", bamFiles[i], ignore.case = TRUE)) {
      warning("File may not be in BAM format: ", bamFiles[i])
    }
    # Check if BAM index exists
    if(!file.exists(paste0(bamFiles[i], ".bai")) && !file.exists(gsub("\\.bam$", ".bai", bamFiles[i]))) {
      warning("BAM index (.bai) file not found for: ", bamFiles[i], 
              ". This may cause failures. Create an index with 'samtools index'.")
    }
  }
  
  # Check bamNames 
  if(length(bamNames) != length(bamFiles)) {
    stop("'bamNames' must have the same length as 'bamFiles'. ",
         "bamFiles length: ", length(bamFiles), ", bamNames length: ", length(bamNames))
  }
  
  # Check for duplicate bamNames
  if(length(unique(bamNames)) != length(bamNames)) {
    duplicates <- bamNames[duplicated(bamNames)]
    stop("Duplicate names found in 'bamNames': ", paste(duplicates, collapse=", "), 
         ". Please provide unique names.")
  }
  
  # Check chips
  if(!all(chips %in% bamNames)) {
    missing_chips <- chips[!chips %in% bamNames]
    stop("Some ChIP sample names are not in 'bamNames': ", paste(missing_chips, collapse=", "))
  }
  
  # Check inputs if provided
  if(!missing(inputs)) {
    if(!all(inputs %in% bamNames)) {
      missing_inputs <- inputs[!inputs %in% bamNames]
      stop("Some input sample names are not in 'bamNames': ", paste(missing_inputs, collapse=", "))
    }
    if(length(inputs) != length(chips)) {
      stop("'inputs' must have the same length as 'chips'. ",
           "chips length: ", length(chips), ", inputs length: ", length(inputs))
    }
    if(any(inputs %in% chips)) {
      overlap <- inputs[inputs %in% chips]
      warning("Some samples are used as both ChIP and input: ", paste(overlap, collapse=", "), 
              ". This is unusual but allowed.")
    }
  }
  
  # Check width
  if(!is.numeric(width) || length(width) != 1) {
    stop("'width' must be a single numeric value. Current value: ", paste(width, collapse=", "))
  }
  if(width < 0) {
    stop("'width' must be non-negative. Current value: ", width)
  }
  
  # Check pseudocount
  if(!is.numeric(pseudocount) || length(pseudocount) != 1) {
    stop("'pseudocount' must be a single numeric value. Current value: ", 
         paste(pseudocount, collapse=", "))
  }
  if(pseudocount <= 0) {
    warning("'pseudocount' should generally be positive to avoid division by zero or log of zero. Current value: ", pseudocount)
  }
  
  # Check PairedEnd
  if(!is.logical(PairedEnd) || length(PairedEnd) != 1) {
    stop("'PairedEnd' must be a single logical value (TRUE/FALSE). Current value: ", 
         paste(PairedEnd, collapse=", "))
  }
  
  # Check requireBothEndsMapped consistency with PairedEnd
  if(requireBothEndsMapped && !PairedEnd) {
    warning("'requireBothEndsMapped' is TRUE but 'PairedEnd' is FALSE. Setting 'requireBothEndsMapped' to FALSE.")
    requireBothEndsMapped <- FALSE
  }
  
  # Check minMQS
  if(!is.numeric(minMQS) || length(minMQS) != 1 || minMQS < 0) {
    stop("'minMQS' must be a single non-negative numeric value. Current value: ", 
         paste(minMQS, collapse=", "))
  }
  
  # Check strand
  if(!is.numeric(strand) || any(!(strand %in% c(0, 1, 2)))) {
    stop("'strand' must contain only values 0 (unstranded), 1 (stranded), or 2 (reversely stranded). ",
         "Current values: ", paste(strand, collapse=", "))
  }
  if(length(strand) != 1 && length(strand) != length(bamFiles)) {
    stop("'strand' must have length 1 or the same length as 'bamFiles'. ",
         "Current length: ", length(strand), ", bamFiles length: ", length(bamFiles))
  }
  
  # Check splitOnly and nonSplitOnly
  if(splitOnly && nonSplitOnly) {
    stop("'splitOnly' and 'nonSplitOnly' cannot both be TRUE.")
  }
  
  # Check minOverlap
  if(!is.numeric(minOverlap) || length(minOverlap) != 1) {
    stop("'minOverlap' must be a single numeric value. Current value: ", 
         paste(minOverlap, collapse=", "))
  }
  
  # Check readExtension3
  if(!is.numeric(readExtension3) || length(readExtension3) != 1 || readExtension3 < 0) {
    stop("'readExtension3' must be a single non-negative numeric value. Current value: ", 
         paste(readExtension3, collapse=", "))
  }
  
  if (!is.numeric(readShiftSize) || readShiftSize < 0) {
    if (!(mode == "Q" && PairedEnd && readShiftSize == "halfInsert")) {
      stop("'readShiftSize' must be a non-negative numeric value or 'halfInsert' (only when mode='Q' and PairedEnd=TRUE)")
    }
  }
  
  # Check read2pos
  if(!is.null(read2pos) && !(read2pos %in% c(3, 5))) {
    stop("'read2pos' must be NULL, 3, or 5. Current value: ", read2pos)
  }
  
  # Check mode
  if(!(mode %in% c("F", "Q"))) {
    stop("'mode' must be either 'F' (FeatureCounts) or 'Q' (QuasR). Current value: ", mode)
  }
  if(mode == "Q") {
    # Check QuasR availability
    if(!requireNamespace("QuasR", quietly = TRUE)) {
      stop("The 'QuasR' package is needed when mode='Q'. Please install it with BiocManager::install('QuasR').")
    }
  } else {
    # Check Rsubread availability
    if(!requireNamespace("Rsubread", quietly = TRUE)) {
      stop("The 'Rsubread' package is needed when mode='F'. Please install it with BiocManager::install('Rsubread').")
    }
  }
  
  # Check GRanges structure - are required fields present?
  if(is.null(names(peaks)) && is.null(peaks$name)) {
    warning("Neither 'name' attribute nor 'names' found in 'peaks' GRanges. Using row numbers as identifiers.")
    names(peaks) <- paste0("peak_", seq_along(peaks))
  }
  
  if(is.null(strand(peaks)) && is.null(peaks$strand)) {
    warning("No strand information found in 'peaks' GRanges. Assuming all features are on the '+' strand.")
    strand(peaks) <- "+"
  }
  
  # Original function execution
  #---------------------------------------------------------------------------------------
  
  # Resize peaks if width > 0
  if (width==0){
    peaks <- peaks
  } else {
    peaks <- resize(peaks, width=width, fix="center")
    peaks <- peaks[start(peaks) >= 0 & width(peaks) == width]
    
    if(length(peaks) == 0) {
      stop("After resizing, no valid peaks remain. Check if your width parameter is appropriate.")
    }
  }
  
  # Ensure peaks has names
  if(is.null(names(peaks)) && !is.null(peaks$name)){
    names(peaks) <- peaks$name
  }
  
  # Ensure peaks has strand information
  if(is.null(strand(peaks)) && !is.null(peaks$strand)){
    strand(peaks) <- peaks$strand
  }
  
  # QuasR mode
  if (mode == "Q"){ 
    tryCatch({
      # write a table to read in samples for QUASR
      quasr_file <- tempfile(pattern = "QUASR", fileext = ".txt")
      write.table(data.frame(FileName=bamFiles, SampleName=bamNames), 
                  file=quasr_file, sep="\t", col.names=TRUE, 
                  row.names=FALSE, append=FALSE, quote=FALSE)
      
      # translate options
      cl <- makeCluster(min(20, parallel::detectCores()))
      on.exit(parallel::stopCluster(cl), add = TRUE)
      
      Qpaired <- ifelse(PairedEnd == TRUE, "fr", "no")
      selectReadPosition <- ifelse(read2pos == 3, "end", "start")
      orientation <- ifelse(strand == 1, "same", ifelse(strand == 2, "opposite", "any"))
      
      # generate project
      message("Creating QuasR project...")
      proj <- qAlign(quasr_file, genome, paired = Qpaired, clObj = cl)
      
      # generate counts matrices
      message("Counting reads with QuasR...")
      counts <- qCount(proj, peaks,
                       selectReadPosition= selectReadPosition, orientation = orientation, 
                       shift = readShiftSize, useRead="any", clObj = cl, mapqMin = minMQS)[,-1]
      
      mapped.reads <- alignmentStats(proj)[,"mapped"]
      unlink(quasr_file)
    }, error = function(e) {
      stop("Error in QuasR processing: ", e$message, 
           "\nConsider switching to mode='F' if the problem persists.")
    })
  }  
  
  # FeatureCounts mode
  if (mode == "F"){
    tryCatch({
      saf <- data.frame(GeneID= names(peaks), Chr=seqnames(peaks),
                        Start=start(peaks), End=end(peaks), Strand=strand(peaks))
      
      message("Counting reads with featureCounts...")
      
      # Calculate the number of reads in each peak for all bam files specified
      f_counts <- featureCounts(bamFiles, annot.ext=saf, useMetaFeatures=FALSE, 
                                allowMultiOverlap=TRUE, minOverlap=minOverlap, 
                                countMultiMappingReads=FALSE, fraction=FALSE,
                                minMQS=minMQS, strandSpecific=strand, 
                                nthreads=min(10, parallel::detectCores()), verbose=FALSE,
                                isPairedEnd=PairedEnd, requireBothEndsMapped=requireBothEndsMapped,
                                readShiftType="downstream", readShiftSize=readShiftSize, read2pos=read2pos,
                                splitOnly=splitOnly, nonSplitOnly=nonSplitOnly, readExtension3=readExtension3)
      
      counts <- f_counts$counts
      colnames(counts) <- bamNames
      
      # Get the number of mapped reads for each bam sample
      if (length(bamFiles) == 1){
        mapped.reads <- sum(f_counts$stat[c(1, 12), -1])
      } else {
        mapped.reads <- apply(f_counts$stat[c(1, 12), -1], 2, sum)   
      }
      
      # Check if mapping statistics look reasonable
      if(any(mapped.reads == 0)) {
        zero_mapped <- names(mapped.reads)[mapped.reads == 0]
        warning("No mapped reads found for: ", paste(zero_mapped, collapse=", "), 
                ". Check that your BAM files are properly formatted.")
      }
      
    }, error = function(e) {
      stop("Error in featureCounts processing: ", e$message, 
           "\nCheck your BAM files and GRanges object for compatibility.")
    })
  }
  
  # Give warning if all counts are 0
  if(max(counts) == 0){
    warning("All regions have 0 counts, possible causes:\n",
            "- Chromosome names (seqnames) don't match between your GRanges object and BAM files\n",
            "- No reads in your BAM files overlap with the specified regions\n",
            "- Filtering parameters (minMQS, strand, etc.) are too stringent")
  }
  
  # Calculate enrichments
  names(mapped.reads) <- bamNames
  log2enrichments <- matrix(ncol=length(chips), nrow=nrow(counts))
  Avalues <- matrix(ncol=length(chips), nrow=nrow(counts))
  
  if(missing(inputs)){
    message("No input samples provided. Calculating counts per million (CPM) values.")
    # Calculate the cpm for each chip
    for (i in seq_along(chips)){
      log2enrichments[,i] <- (counts[,chips[i]]/mapped.reads[chips[i]]) * 1000000
      Avalues[,i] <- counts[,chips[i]]
    }
  } else {
    message("Input samples provided. Calculating log2(ChIP/Input) enrichments.")
    # Calculate the log2 (Chip/Input) enrichments, normalized by a mapped reads (library size) scaling factor
    for (i in seq_along(chips)){
      chip_counts <- counts[,chips[i]]
      input_counts <- counts[,inputs[i]]
      
      # Check for low counts
      low_count_regions <- which(chip_counts < 3 & input_counts < 3)
      if(length(low_count_regions) > 0 && length(low_count_regions) <= 10) {
        warning("Low counts in both ChIP and Input for ", length(low_count_regions), 
                " regions in pair ", chips[i], "/", inputs[i], 
                ". Results may be unreliable for these regions.")
      } else if(length(low_count_regions) > 10) {
        warning("Low counts in both ChIP and Input for ", length(low_count_regions), 
                " regions in pair ", chips[i], "/", inputs[i], 
                ". Consider increasing sequencing depth or adjusting peak calling parameters.")
      }
      
      # Calculate enrichments
      log2enrichments[,i] <- log2(((chip_counts + pseudocount)/(input_counts + pseudocount))/
                                    ((mapped.reads[chips[i]])/(mapped.reads[inputs[i]])))
      Avalues[,i] <- chip_counts + input_counts
    }
  }
  
  colnames(log2enrichments) <- chips
  rownames(log2enrichments) <- rownames(counts)
  colnames(Avalues) <- chips
  rownames(Avalues) <- rownames(counts)
  
  enrichments <- list(log2enrichments, Avalues)
  message("CountPeakReads completed successfully.")
  return(enrichments)
}