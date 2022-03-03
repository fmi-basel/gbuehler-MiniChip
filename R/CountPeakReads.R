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
#' @param requireBothEndsMapped Logical scalar indicating if both ends from the same fragment are required to be successfully aligned before the fragment can be assigned to a feature or meta-feature. 
#' This parameter is only appliable when PairedEnd is TRUE.
#' @param read2pos Specifying whether each read should be reduced to its 5' most base or 3' most base. It has three possible values: NULL, 5 (denoting 5' most base) and 3 (denoting 3' most base).
#' Default value is 5, ie. only the 5' end of the read will be counted. If a read is reduced to a single base, only that base will be considered for the read assignment.
#' Read reduction is performed after read shifting and extension.
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
#'
#' @export
CountPeakReads <- function(peaks,bamFiles,bamNames=bamFiles,chips=bamNames,inputs,width=500,minOverlap = 1,
                           PairedEnd=FALSE, minMQS=255,strand=0, splitOnly=FALSE, nonSplitOnly=FALSE,
                           readExtension3=0,readShiftSize=0,requireBothEndsMapped=FALSE,read2pos=5,pseudocount=8){

#extend  sites to width bp and convert ctcf.sites.gr GRanges object to SAF file 
  if (width==0){
    peaks <- peaks
  } else {
peaks <- resize(peaks,width=width,fix="center")
peaks <- peaks[start(peaks) >= 0 & width(peaks)== width]
  }
  
if(class(names(peaks)) == "NULL" & class(peaks$name) != "NULL"){
  names(peaks) <- peaks$name
}

if(class(strand(peaks)) == "NULL" & class(peaks$strand) != "NULL"){
  strand(peaks) <- peaks$strand
}

saf <- data.frame(GeneID= names(peaks),Chr=seqnames(peaks),
                  Start=start(peaks), End=end(peaks),Strand=strand(peaks))


#calculate the number of reads in each peak for all bam files specified
f_counts <- featureCounts(bamFiles,annot.ext=saf,useMetaFeatures=FALSE,allowMultiOverlap=TRUE,
                          minOverlap=minOverlap,countMultiMappingReads=FALSE,fraction=FALSE,
                          minMQS=minMQS,strandSpecific=strand,nthreads=10,verbose=FALSE,
                          isPairedEnd=PairedEnd,requireBothEndsMapped=requireBothEndsMapped,
                          readShiftType="downstream",readShiftSize=readShiftSize,read2pos = read2pos,
                          splitOnly=splitOnly, nonSplitOnly=nonSplitOnly,readExtension3=readExtension3)


counts <- f_counts$counts
colnames(counts) <- bamNames

# give warning if all counts are 0
if(max(counts) == 0){
  warning("All regions have 0 counts, make sure the chromosome names (seqnames) match between your GRanges object and bam files!")
}

#get the number of mapped reads for each bam sample
if (length(bamFiles) == 1){
  mapped.reads <- sum(f_counts$stat[c(1, 12), -1])
}
else {
  mapped.reads <- apply(f_counts$stat[c(1, 12), -1], 2, sum)   
}
names(mapped.reads) <- bamNames
log2enrichments <- matrix(ncol=length(chips),nrow=nrow(counts))
Avalues <- matrix(ncol=length(chips),nrow=nrow(counts))

if(missing(inputs)){
  #calculate the cpm for each chip
  for (i in seq_along(chips)){
    log2enrichments[,i] <- (counts[,chips[i]]/mapped.reads[chips[i]]) * 1000000
    Avalues[,i] <- counts[,chips[i]]
    
  }
} else {
  #calculate the log2 (Chip/Input) enrichments, normalized by a mapped reads (library size) scaling factor
  for (i in seq_along(chips)){
    log2enrichments[,i] <- log2(((counts[,chips[i]]+pseudocount)/(counts[,inputs[i]]+pseudocount))/((mapped.reads[chips[i]])/(mapped.reads[inputs[i]])))
    Avalues[,i] <- counts[,chips[i]] + counts[,inputs[i]]
  }
  
}

colnames(log2enrichments) <- chips
rownames(log2enrichments) <- rownames(counts)
colnames(Avalues) <- chips
rownames(Avalues) <- rownames(counts)

enrichments <- list(log2enrichments,Avalues)
return(enrichments)
}