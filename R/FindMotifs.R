#' @title FindMotis
#'
#' @description Finds enriched motifs in supplied regions (peaks) using MEME and HOMER.
#'
#' @details Finds enriched motifs in supplied regions (peaks) in the mouse mm10 genome using MEME and HOMER. It will write a fasta file with the peak sequences
#' (for MEME Input) and a bed file with the peak coordinates (for HOMER Input). If runMEME/runHOMER is set to TRUE (the default)
#' it will also run the MEME/HOMER analysis (requires MEME/HOMER to be installed).
#'
#' @param peaks A Granges object containing your positions of interest. Must include seqnames (chromosomes), start, end, unique peak (row) name.
#' @param peaknames A character string that describes your peaks. Default is "mypeaks".
#' @param topn The number of top peaks (counted from the top row of the GRanges object) to use for motif finding. Default is 100.
#' @param width The width of the peak (around the center of each range) to use for motif finding. Default is 100.
#' @param runMEME Should MEME motif finding be run? Default = TRUE.
#' @param runHOMER Should HOMER motif finding be run? Default = TRUE.
#'
#' @return Saves the outputs of MEME and HOMER.
#'
#' @importFrom GenomicRanges start
#' @importFrom GenomicRanges start<-
#' @importFrom GenomicRanges end<-
#' @importFrom GenomicRanges end
#' @importFrom GenomicRanges width
#' @importFrom GenomicRanges strand
#' @importFrom GenomicRanges resize
#' @importFrom BSgenome getSeq
#' @importFrom Biostrings writeXStringSet
#' @importFrom utils write.table
#'
#' @examples
#' peaks <- GenomicRanges::GRanges(
#' seqnames = Rle(c("chr1", "chr2", "chr1", "chr3"), c(1, 3, 2, 4)),
#' ranges = IRanges(50101:50110, end = 51111:51120),
#' strand = Rle(strand(c("-", "+", "*", "+", "-")), c(1, 2, 2, 3, 2)),
#' summit = 1:10, name = head(letters, 10))
#' FindMotifs(peaks)
#'
#' @export
FindMotifs <- function(peaks,peaknames="mypeaks",topn=100,width=100, runMEME=TRUE, runHOMER=TRUE){
peaks <- resize(peaks, width = width, fix = "center")

if(class(names(peaks)) == "NULL" & class(peaks$name) != "NULL"){
  names(peaks) <- peaks$name
}

if(class(strand(peaks)) == "NULL" & class(peaks$strand) != "NULL"){
  strand(peaks) <- peaks$strand
}

#remove peaks with negative start values
peaks <- peaks[start(peaks) >= 0]

#define top n peaks to use
for (i in seq_along(topn)){


  #make fasta file for MEME
  peaks.seq <- getSeq(BSgenome.Mmusculus.UCSC.mm10,peaks[1:topn[i]])
  names(peaks.seq) <- names(peaks)[1:topn[i]]
  writeXStringSet(peaks.seq,sprintf("%s.%sbp_top%s.fasta",peaknames,width,topn[i]),format = "fasta" )

  if(runMEME == TRUE){
  #DEFINE MEME ARGUMENTS
  fastafile <- sprintf("%s.%sbp_top%s.fasta",peaknames,width,topn[i])
  dna <- "-dna"
  mode <- "-mod oops"
  strand <- "-revcomp"
  motiflen <- "-maxw 25"
  nmotifs <- "-nmotifs 3"
  maxsize <- "-maxsize 1000000"
  outputfile <- sprintf("-o %s.%sbp_top%s.motifs.meme",peaknames,width,topn[i])

  #run MEME
  system2(command="/work/gbioinfo/Appz/meme/meme_4.10.0x/bin/meme",args=c(fastafile,dna,mode,strand,motiflen,nmotifs,maxsize,outputfile),wait=FALSE)
  }


  #make bed file for HOMER
  peaks.bed <- data.frame(name=names(peaks)[1:topn[i]],seqnames=seqnames(peaks[1:topn[i]]),
                           start=start(peaks[1:topn[i]]),end=end(peaks[1:topn[i]]), strand=rep("+",length(peaks[1:topn[i]])))
  #peaks.bed <- data.frame(name=names(peaks)[1:topn[i]],seqnames=seqnames(peaks[1:topn[i]]),
  #                        start=start(peaks[1:topn[i]]),end=end(peaks[1:topn[i]]), strand=strand(peaks))
  write.table(peaks.bed,sprintf("%s.%sbp_top%s.bed",peaknames,width,topn[i]),sep="\t",col.names=F,row.names=F,append=FALSE,quote=FALSE)

  if (runHOMER == TRUE){
  #DEFINE HOMER ARGUMENTS
  bedfile <- sprintf("%s.%sbp_top%s.bed",peaknames,width,topn[i])
  genome <- "/work/gbioinfo/DB/genomes/mm10"
  outputfileH <- sprintf("%s.%sbp_top%s.motifs.homer",peaknames,width,topn[i])
  size <- sprintf("-size %d",width)
  cores <- "-p 10"
  noknown <- "-noknown"
  HomerDir <-  "-preparsedDir /tungstenfs/scratch/gbuehler/bioinfo/Annotations/Homer"

  #run HOMER
  #system(command="PATH=$PATH:/work/gbioinfo/Appz/Homer/Homer-current/bin")
  #system(command="PATH=$PATH:/work/gbioinfo/Appz/weblogo/weblogo-current")

  system2(command="/work/gbioinfo/Appz/Homer/Homer-current/bin/findMotifsGenome.pl",args=c(bedfile,genome,outputfileH,size,cores,noknown,HomerDir),wait=FALSE)
 # print(unquote(c("/work/gbioinfo/Appz/Homer/Homer-current/bin/findMotifsGenome.pl",bedfile,genome,outputfileH,size,cores,noknown,HomerDir)))
}
}
}
