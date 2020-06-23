#' @title GCbias
#'
#' @description Plot GC content versus read counts.
#'
#' @details This function generates a plot of the number of Gs and Cs against the read count (cpm) in windows across the chosen genome.
#' A seperate plot is generated for each alignment file (bam file) that is supplied. Supports both single and paired-end experiments.
#'
#' @param bamFiles The filenames (including full path) of the bam files containing the reads you want to count around the peak summits.
#' @param bamNames A name to describe the bam files you are using (for example: "H3K9me3_reads"). If no names are supplied, the full bamFiles strings are used.
#' @param minMQS The minimum mapping quality. Default is 255, which eliminates multimapping reads when STAR was used to generate bam files.
#' @param maxFrag An integer scalar, specifying the maximum fragment length corresponding to a read pair. Defaults to 500 base pairs.
#' @param pe A string indicating whether paired-end data is present; set to "none" (the default), "both", "first" or "second".
#' @param restrict A character vector containing the names of allowable chromosomes from which reads will be extracted. Default is "chr1".
#' @param winWidth An integer scalar specifying the width of the window, in which reads are counted and GC content calculated. Default is 5000 base pairs.
#' @param col Color scheme for the plots. Default is "plasma" from the viridis package.
#' @param genome A BSGenome object for the species you are analyzing. Required parameter. For example,use BSgenome.Mmusculus.UCSC.mm10 for mouse.
#' @param GCprob Should the GC content be displayed as absolute counts (GCprob=FALSE) or as fraction of GCs (GCprob=TRUE,default).
#' @param span Span that is used for loess trendline. Default= 0.1
#'
#' @return A smooth Scatter Plot for every suuplied alignment file: The number of Gs and Cs in each window is plotted on the x-axis, the read (fragment if pe="both") counts
#' (in counts/million) is plotted on the y-axis.
#'
#' @examples
#' bamFiles <- c("/work2/gbuehler/deepSeqRepos/bam//HP1a_wt_ChIP_r1_818F1_multi.bam","/work2/gbuehler/deepSeqRepos/bam//HP1a_wt_ChIP_r2_818F3_multi.bam")
#' GCbias(bamFiles=bamFiles)
#'
#' @importFrom csaw readParam
#' @importFrom csaw windowCounts
#' @importFrom csaw calculateCPM
#' @importFrom SummarizedExperiment assay
#' @importFrom csaw calculateCPM
#' @importFrom SummarizedExperiment rowRanges
#' @importFrom BSgenome getSeq
#' @importFrom Biostrings letterFrequency
#' @importFrom viridis inferno
#' @importFrom graphics smoothScatter
#' @importFrom graphics lines
#' @importFrom graphics par
#' @importFrom stats loess
#' @importFrom stats predict
#'
#' @export
GCbias <- function(bamFiles, bamNames=bamFiles, minMQS=255,maxFrag=500,pe="none",restrict="chr1",
                   winWidth=5000, col=inferno, genome, GCprob=TRUE,span=0.1){
  #define parameters for read extarction
  paramPE <- csaw::readParam(minq=minMQS,max.frag=maxFrag, pe=pe,restrict=restrict)

   #count number of reads in bins
  chip_bins <- csaw::windowCounts(bamFiles,  width=winWidth, param=paramPE, bin = TRUE, filter=0)
  bins_data <- data.frame(SummarizedExperiment::assay(chip_bins))
  colnames(bins_data) <- bamNames

  #calculate cpm
  bins_cpm <- data.frame(csaw::calculateCPM(chip_bins, use.norm.factors=TRUE, use.offsets=FALSE,
                                 log=TRUE, prior.count=1, assay.id="counts"))
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

  #plot GC against cpms
  par(mfrow=c(1,length(bamNames)))
  for (i in seq_along(bamNames)){
  smoothScatter(bins_cpm$bins_gc, bins_cpm[,i], colramp=col,xlab=xaxislab,ylab="cpm",main=bamNames[i])

  #add trendline
  loessGC <- loess(bins_cpm[,i] ~ bins_cpm$bins_gc, span=span)
  lines(bins_cpm$bins_gc,predict(loessGC),lwd=2,col="white")
  }
}
