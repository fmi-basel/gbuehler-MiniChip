#' @title CountPeakReads
#'
#' @description Generate normalized counts around positions of interest in a genome.
#'
#' @details This function calculates the numebr of reads of different experiments around positions of interest in the genome,
#' for example the summits of ChIP peaks, or transcription start sites. The center of your GRanges object defines the coordinates
#' of these positions of interest (eg peak summits, TSSs). Then the number of reads in a ChIP sample is normalized to the number
#' of reads in the corresponding Input sample, if provided, taking the total library sizes into account. 
#' If no Input sample names are provided, the counts per milliion reads (cpm) for each ChIP samples is returned. 
#'
#' @param peaks A Granges object containing your positions of interest. Must include seqnames (chromosomes), start, end, strand, and name.
#' @param bamFiles The filenames (including full path) of the bam files containing the reads you want to count around the peak summits.
#' @param bamNames A name to describe the bam files you are using (for example: "H3K9me3_reads"). Required argument.
#' @param chips The names of the ChIP samples (corresponding exactly to the ChIP sample names in bamNames).
#' @param inputs The names of the Input samples (corresponding exactly to the Input sample names in bamNames), 
#' provided in an order crresponding to the ChIP samples. If missing, no Input normlaization is performed.
#' @param width The width around the center of the peaks GRanges object provided in which the reads will be counted,default = 500. 
#' If width=0, the provided regions will remain the original size. 
#' @param pseudocount The pseudocount that is added to the counts of Chip and Input before normalization to avoind dividing by 0. Default is 8.
#' @param PairedEnd Are reads in bam files generated with paired-end or single-end sequencing. Default is FALSE (=single-end).
#' @param minMQS The minimum mapping quality. Default is 255, which eliminates multimapping reads when STAR was used to generate bam files.
#' @param strand Strand specific sequencing? If the reads in the bam file are derived from both strands (usually the case for ChIP experiments),
#' use 0 (default). 1 or 2 are the strans specific options (as in feature counts).
#' @param splitOnly Should only split reads (for example spliced reads) be counted? FALSE by default.
#' @param nonSplitOnly Should only non-split reads (for example not spliced reads) be counted? FALSE by default.
#' @param minOverlap The minimum overlap of a read with the window for the read to be counted. Default = 1.
#' @param readExtension3 integer giving the number of bases extended downstream from 3' end of each read. 0 by default. Negative value is not allowed.
#' @param readShiftSize integer specifying the number of bases the reads will be shifted downstream by. 100 by default. Negative value is not allowed.
#' @param requireBothEndsMapped logical indicating if both ends from the same fragment are required to be successfully aligned before the fragment can be assigned to a feature or meta-feature. This parameter is only appliable when isPairedEnd is TRUE.
#' @param read2pos Specifying whether each read should be reduced to its 5' most base or 3' most base. It has three possible values: NULL, 5 (denoting 5' most base) and 3 (denoting 3' most base).
#' Default value is 5, ie. the 5' end of the read will be counted. If a read is reduced to a single base, only that base will be considered for the read assignment.
#' Read reduction is performed after read shifting and extension.
#'
#' @return A list of 2 matrices. Both with a column for each chip sample provided, and the same number of rows as regions provided in the GRanges object.
#' The first matrix contains the (Input) normalized number of reads in a window around each position of interest, which was the center in the supplied GRanges object.
#' The second matrix countains the total number of reads in the the chip sample, or in the chip plus the input sample, if an input sample is provided.
#' 
#'
#' @examples
#' peaks <- GenomicRanges::GRanges(
#' seqnames = Rle(c("chr1", "chr2", "chr1", "chr3"), c(1, 3, 2, 4)),
#' ranges = IRanges(50101:50110, end = 51111:51120),
#' strand = Rle(strand(c("-", "+", "*", "+", "-")), c(1, 2, 2, 3, 2)),
#' summit = 1:10, name = head(letters, 10))
#' bamFiles <- list.files(system.file("extdata", package = "MiniChip"), full.names=TRUE,pattern="*bam$")
#' bamNames <- gsub(paste(system.file("extdata", package = "MiniChip"),"/",sep=""),"",bamFiles)
#' bamNames <- gsub("_chr11.bam","",bamNames)
#' CountPeakReads(peaks=peaks,bamFiles=bamFiles,bamNames=bamNames,chips=bamFiles[1:3],inputs=bamFiles[4:6])
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
CountPeakReads <- function(peaks,bamFiles,bamNames,chips,inputs,width=500,minOverlap = 1,
                           PairedEnd=FALSE, minMQS=255,strand=0, splitOnly=FALSE, nonSplitOnly=FALSE,
                           readExtension3=0,readShiftSize=0,requireBothEndsMapped=FALSE,read2pos=5,pseudocount=8){

#extend  sites to width bp and convert ctcf.sites.gr GRanges object to SAF file 
  if (width==0){
    peaks <- peaks
  } else {
peaks <- resize(peaks,width=width,fix="center")
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

#get the number of mapped reads for each bam sample
mapped.reads <- apply(f_counts$stat[c(1,12),-1],2,sum)
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
  #calculate the log2 (Chip/Input) enrichmnets, normalized by a mapped reads (library size) scaling factor
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