#' @title DifferentialEnrichment
#'
#' @description Calculate differential ChIP enrichment in peaks.
#'
#' @details Calculates read counts for two groups of ChIPs (at least 2 replicates each), normalizes them to total mapped read counts (voom),
#'  and calculates differential enrichment between all groups and a control group using Limma. 
#'
#' @param peaks A GRanges object containing your regions of interest. Must include seqnames (chromosomes), start, end, strand, and name.
#' @param bamFiles Character vector containing the filenames (including the full path) of read alignment files in bam format.
#' @param bamNames Character vector containing the names to describe the \code{bamFiles} you are using (for example: "H3K9me3_reads"). 
#' If no names are supplied, the full \code{bamFiles} names are used.
#' @param group Character vector containing the group names that the ChIP samples belong to, with the order corresponding to the ChIP sample names in \code{bamNames}.
#' @param controlGroup Character scalar giving the name of the control group that all other groups will be compared to. Eg "WT".  If not specified the last element of the group vector will be used. 
#' @param width Integer scalar providing the width around the center of the peaks GRanges object provided in which the reads will be counted,default = 0. 
#' If width=0, the provided regions will remain the original size. 
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
#' @return A data.frame with the logFC (log2 fold change), p-values, and adjusted p-values for each peak and each contrast (each group compared to the control group). 
#'
#'#' @examples
#' peaks <- SimulatePeaks(1000,rep(100,1000),chromosomeSizes=
#' system.file("extdata", "chrNameLength_mm10_chr11.txt", package = "MiniChip"))
#' bamFiles <- list.files(system.file("extdata", package = "MiniChip"), full.names=TRUE,pattern="*bam$")
#' bamNames <- gsub(paste(system.file("extdata", package = "MiniChip"),"/",sep=""),"",bamFiles)
#' bamNames <- gsub("_chr11.bam","",bamNames)
#' DifferentialEnrichment(peaks=peaks,bamFiles=bamFiles,bamNames=bamNames,group=c("A","A","A","IN","IN","IN"))
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
#' @importFrom edgeR DGEList
#' @importFrom edgeR calcNormFactors
#' @importFrom stats model.matrix
#' @importFrom limma voom
#' @importFrom limma lmFit
#' @importFrom limma makeContrasts
#' @importFrom limma contrasts.fit
#' @importFrom limma eBayes
#' @importFrom limma topTable
#' 
#' @export
DifferentialEnrichment <- function(peaks,bamFiles,bamNames=bamFiles,group,controlGroup=group[length(group)], width=0,minOverlap = 1,
                           PairedEnd=FALSE, minMQS=255,strand=0, splitOnly=FALSE, nonSplitOnly=FALSE,
                           readExtension3=0,readShiftSize=0,requireBothEndsMapped=FALSE,read2pos=5){
  
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
  
  #number of mapped reads (within and outside of peaks)
  mapped.reads <- apply(f_counts$stat[c(1,12),-1],2,sum)
  names(mapped.reads) <- bamNames
  
  
#generate DE object and calulte normalisation factors from library size (number of mapped reads)
d <- DGEList(counts,lib.size=mapped.reads)
d <- calcNormFactors(d,method="none")


#generate group info and model (use the control group as the last level)
group <- factor(group,levels=c(unique(as.character(group[group != controlGroup])),controlGroup),
                labels=c(unique(as.character(group[group != controlGroup])),controlGroup))
mm <- model.matrix(~0 + group)

#normalize the counts using voom 
y <- voom(d, mm, plot = F)

#fit linear model
fit <- lmFit(y, mm)
groups2 <- colnames(coef(fit))

#make the list of contrasts
contr <- list()
for (i in 1:(length(groups2)-1)){
  con1 <- paste(groups2[i],groups2[length(groups2)],sep=" - ")
  contr[[i]] <- makeContrasts(con1, levels = colnames(coef(fit)))
  
}
#name the contrasts in the list
contr.names <- character(0)
for (i in 1:(length(groups2)-1)){
contr.names <- c(contr.names,paste(groups2[i],groups2[length(groups2)],sep="_vs_"))
}
names(contr) <- contr.names

#generate a list of DE results 
res <- list()
for (i in seq_along(contr)){
  tmp <- contrasts.fit(fit, contr[[i]])
  tmp <- eBayes(tmp)
  #limma::plotMA(tmp)
  
  top.table <- as.data.frame(topTable(tmp, sort.by = "P", n = Inf))
  top.table$log.adj.P.Val <- -log10(top.table$adj.P.Val)
  top.table$Contrast <- names(contr)[i]
  top.table$ID <- row.names(top.table)
  res[[i]] <- top.table
}

#rbind the list of results across contrasts and return the table
res <- do.call("rbind",res)
return(res)

}
