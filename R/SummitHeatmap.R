#' @title SummitHeatmap
#'
#' @description This function generates a heatmap of different experiments (read counts) around the summits of ChIP peaks.
#'
#' @details This function.
#'
#' @param peaks A Granges object containing your ChIP peaks. Must include seqnames (chromosomes), start, end, strand, summit, and name.
#' @param peaknames A character string that describes your peaks. Default is "mypeaks".
#' @param bamFiles The filenames (including full path) of the bam files containing the reads you want to count around the peak summits.
#' @param bamNames A name to describe the bam files you are using (for example: "H3K9me3_reads"). Defaults to "myreads".
#' @param span The distance from the peak summit to the left and right that you want your heatmap to extand to. Default is 2025.
#' @param step The window size in which reads are counted/plotted in the heatmap. Default is 50.
#' @param plotcols The colors for the heatmap. Default is colorRampPalette(brewer.pal(9,"RdPu"))(255).
#' @param plotHM Also plot the heatmap or only count the reads in the windows. Default is TRUE.
#' @param useCPM Normalize the number of reads per window to the total number of reads in the sample (bam file). Default is TRUE.
#' @param PairedEnd Are reads in bam files generated with paired-end or single-end sequencing. Default is FALSE (=single-end).
#'
#'
#' @return A list of matrices of the same length as the number of bam files supllied.
#' Each matrix contains the number of reads (or cpm) in each window (column headers=middle of the window)
#' for each peak (rownames=peaknames). If plotHM is TRUE, a heatmap is plotted for each sample as well (each bam file)
#' and saved in a directory (in your current working directory) that is named after your peaknames object.
#' These heatmaps are clustered individually by euclidian distance.
#'
#' @examples
#' peaks <- GRanges(
#' seqnames = Rle(c("chr1", "chr2", "chr1", "chr3"), c(1, 3, 2, 4)),
#' ranges = IRanges(101:110, end = 111:120,
#' strand = Rle(strand(c("-", "+", "*", "+", "-")), c(1, 2, 2, 3, 2)),
#' summit = 1:10, name = head(letters, 10))
#' bamFiles <- c("/work2/gbuehler/deepSeqRepos/bam//HP1a_wt_ChIP_r1_818F1_multi.bam","/work2/gbuehler/deepSeqRepos/bam//HP1a_wt_ChIP_r2_818F3_multi.bam")
#' SummitHeatmap(peaks,bamFiles,plotHM=FALSE)
#'
#' @importFrom RColorBrewer brewer.pal
#' @importFrom GenomicRanges start
#' @importFrom GenomicRanges end
#' @importFrom GenomicRanges strand
#' @importFrom GenomicRanges elementMetadata
#' @importFrom GenomicRanges seqnames
#' @importFrom Rsubread featureCounts
#' @importFrom gplots heatmap.2
#'
#' @export
SummitHeatmap <- function(peaks,peaknames="mypeaks",bamFiles,bamNames="myreads",span=2025,step=50,
                           plotcols=colorRampPalette(brewer.pal(9,"RdPu"))(255),plotHM=TRUE,useCPM=TRUE,PairedEnd=FALSE){
  #test ifGranges object has all required columns
  if(is.null(peaks$summit)) { stop("Please make sure your peaks GRanges object has a summit metadata culumn")}
  #make the summit of the peak the start and end, then define the whole heatmap region
  start(peaks) <- start(peaks) + peaks$summit - span
  end(peaks) <- start(peaks) + span*2
  #remove peaks with negative start values
  peaks <- peaks[start(peaks) >= 0]

  #generate window starts and ends across span
  windows <- seq(from=0,to=span*2-step,by=step)
  binmids <- windows - span + step/2

  #generate empty saf format data frame
  saf <- data.frame(GeneID= NA,Chr=NA, Start=NA,End=NA,Strand=NA)
  saf <- saf[FALSE, ]

  #generate genomic coordinates for each window across all peaks, concatenate these windows into one big saf format file
  for (w in seq_along(windows)){
    saf.w <- data.frame(GeneID= elementMetadata(peaks)[,"name"], Chr=seqnames(peaks),
                        Start=start(peaks) + windows[[w]], End=start(peaks) + windows[[w]] + step,Strand=strand(peaks))
    saf <- rbind(saf,saf.w)
  }

  #calculate the number of reads in each window for all bam files specified
  f_counts <- featureCounts(bamFiles,annot.ext=saf,useMetaFeatures=FALSE,allowMultiOverlap=TRUE,
                            minOverlap=2,readExtension3=100,countMultiMappingReads=FALSE,fraction=FALSE,
                            minMQS=10,strandSpecific=0,nthreads=4,verbose=FALSE,isPairedEnd=PairedEnd)

  #extract gene annotation (for gene names and strand)
  anno <-  f_counts$annotation[1:(nrow(f_counts$annotation)/length(windows)),]

  #extract number of mapped reads for all bam samples
  mapped.reads <- apply(f_counts$stat[c(1,10),-1],2,sum)

  #extract counts for all samples, name coumns and rows, calculate log CPM
  all.counts <- vector("list", length(bamNames))

  for (bam.sample in seq_along(bamNames)){
    counts <- matrix(f_counts$counts[,bam.sample],nrow=length(peaks),ncol=length(windows), byrow=FALSE)
    colnames(counts) <- binmids
    rownames(counts) <- anno$GeneID
    cpm <-  ((counts+0.1)/mapped.reads[bam.sample])*1000000
    counts.log <- log2(cpm)

    if (plotHM == TRUE){
      #create directory for heatmaps
      ifelse(!dir.exists(sprintf("heatmaps-%s",peaknames)), dir.create(sprintf("heatmaps-%s",peaknames)), FALSE)
      #define distance fucntion for heatmaps
      euc.dist <- function(x){dist(x,method="euclidian")}

      #plot and save heatmaps
      png(sprintf("heatmaps-%s/heatmap-%s-%s.png",peaknames,peaknames,bamNames[bam.sample]),height=1000,width=400)
      heatmap.2(counts.log,col=plotcols,distfun=euc.dist,Colv=FALSE,scale="none",trace="none",dendrogram = "none",
                labRow = "",labCol = binmids,xlab="bins",ylab=sprintf("summits-%s-readcounts-%s",peaknames,bamNames[bam.sample]),
                key=FALSE)
      dev.off()
    } else {
      next
    }

    #prepare counts or cpm for saving as list of matrices
    if (useCPM == TRUE) {
      all.counts[[bam.sample]] <- cpm
    } else {
      all.counts[[bam.sample]] <- counts
    }
  }

  #return the results
  names(all.counts) <- bamNames
  return(all.counts)
}

