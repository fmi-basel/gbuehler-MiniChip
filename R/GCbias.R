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
#'
#' @export
GCbias <- function(bamFiles, bamNames=bamFiles, minMQS=255,maxFrag=500,pe="none",restrict="chr11",
                   winWidth=5000, col=inferno, genome, GCprob=TRUE,span=0.1,plot=TRUE,logCPM=TRUE,priorCount=1){
  #define parameters for read extarction
  paramPE <- csaw::readParam(minq=minMQS,max.frag=maxFrag, pe=pe,restrict=restrict)

   #count number of reads in bins
  chip_bins <- csaw::windowCounts(bamFiles,  width=winWidth, param=paramPE, bin = TRUE, filter=0)
  bins_data <- data.frame(SummarizedExperiment::assay(chip_bins))
  colnames(bins_data) <- bamNames

  #calculate cpm
  bins_cpm <- data.frame(csaw::calculateCPM(chip_bins, use.norm.factors=TRUE, use.offsets=FALSE,
                                 log=logCPM, prior.count=priorCount, assay.id="counts"))
  colnames(bins_cpm) <- bamNames

  #extract window coordinates
  bins_ranges <- SummarizedExperiment::rowRanges(chip_bins)

  #extract sequecnes
  bins_seqs <- BSgenome::getSeq(genome,bins_ranges)
  #extract GC content
  bins_cpm$bins_gc <- Biostrings::letterFrequency(bins_seqs,"GC", as.prob = GCprob)

  if (GCprob==FALSE){
    xaxislab <- "number of Gs and Cs"
  } else {
    xaxislab <- "fraction GC"
  }

  #order table by GC content
  bins_cpm <- bins_cpm[order(bins_cpm$bins_gc),]
  
  if(plot==FALSE){
  return(bins_cpm)
  } else {

  #plot GC against cpms
  par(mfrow=c(1,length(bamNames)))
  for (i in seq_along(bamNames)){
  smoothScatter(bins_cpm$bins_gc, bins_cpm[,i], colramp=col,xlab=xaxislab,ylab="cpm",main=bamNames[i])

  #add trendline
  loessGC <- loess(bins_cpm[,i] ~ bins_cpm$bins_gc, span=span)
  lines(bins_cpm$bins_gc,predict(loessGC),lwd=2,col="white")
  }
  }
}
