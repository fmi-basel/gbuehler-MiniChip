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
                          useCPM=TRUE,PairedEnd=FALSE, minMQS=255,strand=0, splitOnly=FALSE, nonSplitOnly=FALSE,
                          readExtension3=0,readShiftSize=0,requireBothEndsMapped=FALSE,read2pos=5,mode="F",
                          genome="BSgenome.Mmusculus.UCSC.mm10"){

  if(class(names(peaks)) == "NULL" & class(peaks$name) != "NULL"){
    names(peaks) <- peaks$name
  }
  
  if(class(strand(peaks)) == "NULL" & class(peaks$strand) != "NULL"){
    strand(peaks) <- peaks$strand
  }
  
  if (mode == "Q"){ #use QUASR version 
    #write a table to read in samples for QUASR
    write.table(data.frame(FileName=bamFiles,SampleName=bamNames),file="QUASR.txt",sep="\t",col.names=TRUE,row.names=FALSE,append=FALSE,quote=FALSE)
    
    #translate options
    cl <- makeCluster(20)
    Qpaired <- ifelse(PairedEnd ==TRUE,"fr","no")
    selectReadPosition <- ifelse(read2pos == 3,"end","start")
    orientation <- ifelse(strand == 1, "same", ifelse(strand == 2, "opposite","any"))
    
    #generate project
    proj <- qAlign("QUASR.txt", genome, paired = Qpaired, clObj = cl)
    #generate counts matrices
    prf <- qProfile(proj, resize(peaks, width = 1L, fix = "center"), upstream=span,binSize = step,
                    selectReadPosition= selectReadPosition, orientation = orientation, shift = readShiftSize, 
                    useRead="any", clObj = cl, mapqMin = minMQS)[-1]
    
   if (useCPM == TRUE){ 
     #normalize to million reads (cpm)
    normfactor <- alignmentStats(proj)[,"mapped"]/1000000
    prfcpm <- prf
    for (i in seq_along(prf)){
      prfcpm[[i]] <- prf[[i]]/normfactor[i]
    }
    return(prfcpm)
   } else {
     return(prf)
   }
    unlink("QUASR.txt")
    }  
    
    if (mode == "F"){ #use Feature Counts version 
  
  
  #take the middle of the GRanges region, then define whole heatmap region
  nwindows <- ceiling((span*2)/step)
  if(nwindows %% 2 == 0){
    nwindows <- nwindows + 1
  }
  regionwidth <- step * nwindows
  peaks <- resize(peaks,width=regionwidth,fix="center")

  #remove peaks with negative start values
  peaks <- peaks[start(peaks) >= 0 & width(peaks)== regionwidth]

  #generate window starts and ends across span
  windows <- seq(from=0,to=regionwidth-step,by=step)
  binmids <- windows - regionwidth/2 + step/2

 #generate sliding windows across peaks 
  peaks2plus <- unlist(GenomicRanges::slidingWindows(peaks[strand(peaks)!="-"], width = step, step = step),use.names=FALSE)
  peaks2minus <- rev(unlist(GenomicRanges::slidingWindows(peaks[strand(peaks)=="-"], width = step, step = step),use.names=FALSE))
  peaks2 <- c(peaks2plus,peaks2minus)


  #generate  saf format data frame
    saf <- data.frame(GeneID= names(peaks2), Chr=seqnames(peaks2),
                        Start=start(peaks2), End=end(peaks2),Strand=strand(peaks2))

  #calculate the number of reads in each window for all bam files specified
  f_counts <- featureCounts(bamFiles,annot.ext=saf,useMetaFeatures=FALSE,allowMultiOverlap=TRUE,read2pos=read2pos,
                            minOverlap=minOverlap,readExtension3=readExtension3,countMultiMappingReads=FALSE,fraction=FALSE,
                            minMQS=minMQS,strandSpecific=strand,nthreads=20,verbose=FALSE,isPairedEnd=PairedEnd,
                            splitOnly=splitOnly,nonSplitOnly=nonSplitOnly,readShiftType="downstream",readShiftSize=readShiftSize,requireBothEndsMapped=requireBothEndsMapped)

  #extract gene annotation (for peak names)
#  anno <-  f_counts$annotation[1:(nrow(f_counts$annotation)/length(windows)),]
  GeneID <- matrix(names(peaks2),nrow=length(peaks),ncol=length(windows), byrow=TRUE)[,1]

  #extract number of mapped reads for all bam samples
  if(length(bamNames)>1){
  mapped.reads <- apply(f_counts$stat[c(1,11,12),-1],2,sum)
  } else {
    mapped.reads <-  sum(f_counts$stat[c(1,11,12), -1])
  }

  #extract counts for all samples, name columns and rows, calculate log CPM
  all.counts <- vector("list", length(bamNames))

  for (bam.sample in seq_along(bamNames)){
    counts <- matrix(f_counts$counts[,bam.sample],nrow=length(peaks),ncol=length(windows), byrow=TRUE)

    colnames(counts) <- binmids
   # rownames(counts) <- anno$GeneID
    rownames(counts) <- GeneID

    #prepare counts or cpm for saving as list of matrices
    if (useCPM == TRUE) {
      all.counts[[bam.sample]] <- ((counts)/mapped.reads[bam.sample])*1000000
    } else {
      all.counts[[bam.sample]] <- counts
    }
    #sort the rows by the original peak GRanges object order
    all.counts[[bam.sample]] <- all.counts[[bam.sample]][names(peaks),]
    
  }
  # give warning if all counts are 0
  if(max(all.counts[[1]]) == 0){
    warning("All windows have 0 counts, make sure the chromosome names (seqnames) match between your GRanges object and bam files!")
  }
  
  #return the results
  names(all.counts) <- bamNames
  return(all.counts)
  
    } else {
      warning("mode must be one of Q (QuasR) or F (FeatureCounts)!")
    }
  
}

