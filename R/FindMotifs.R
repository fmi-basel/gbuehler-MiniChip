#' @title FindMotifs
#'
#' @description Finds enriched motifs in supplied regions (peaks) using MEME and HOMER.
#'
#' @details Finds enriched motifs in supplied regions (peaks) in the supplied genome using MEME and HOMER. It will write a fasta file with the peak sequences
#' (for MEME Input) and a bed file with the peak coordinates (for HOMER Input). If \code{runMEME} and/or \code{runHOMER} is set to TRUE
#' it will also run the MEME and/or HOMER analysis (requires MEME and/or HOMER to be installed).
#'
#' @param peaks A Granges object containing your positions of interest. Must include seqnames (chromosomes), start, end, unique peak (row) name.
#' @param peaknames A character scalar that describes your peaks. Default is "mypeaks".
#' @param genome BSGenome object. Required parameter. For example,use BSgenome.Mmusculus.UCSC.mm10 for mouse.
#' @param topn Numeric scalar or vector indicating the number of top peaks (counted from the top row of the GRanges object) to use for motif finding.
#' Default is 100. If a vector is supplied (eg. c(100,1000)), then the function generates a MEME and a HOMER input file (and runs MEME/HOMER) for each element (eg. the 100 top peaks, then the 1000 top peaks).
#' @param width Numeric scalar indicating the width of the peak (around the center of each range) to use for motif finding. Default is 100.
#' @param runMEME Logical scalar, indicating whether the MEME motif finding should be run. Default = FALSE.
#' @param runHOMER Logical scalar, indicating whether the HOMER motif finding should be run. Default = FALSE.
#' @param genomepath Character string giving the full path to the genome that HOMER should use for extracting peak sequences.
#' Only required if \code{runHOMER} = TRUE. Default = "/work/gbioinfo/DB/genomes/mm10".
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
#' library(BSgenome.Mmusculus.UCSC.mm10)
#' peaks <- GenomicRanges::GRanges(
#' seqnames = Rle(c("chr1", "chr2", "chr1", "chr3"), c(1, 3, 2, 4)),
#' ranges = IRanges(50101:50110, end = 51111:51120),
#' strand = Rle(strand(c("-", "+", "*", "+", "-")), c(1, 2, 2, 3, 2)),
#' summit = 1:10, name = head(letters, 10))
#' FindMotifs(peaks,genome=BSgenome.Mmusculus.UCSC.mm10,topn=2)
#' unlink("mypeaks.*")
#'
#' @export
FindMotifs <- function(peaks,peaknames="mypeaks",genome, topn=100,width=100, 
                       runMEME=FALSE, runHOMER=FALSE, genomepath="/work/gbioinfo/DB/genomes/mm10"){
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
  peaks.seq <- getSeq(genome,peaks[1:topn[i]])
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
  outputfileH <- sprintf("%s.%sbp_top%s.motifs.homer",peaknames,width,topn[i])
  size <- sprintf("-size %d",width)
  cores <- "-p 10"
  noknown <- "-noknown"
  HomerDir <-  "-preparsedDir /tungstenfs/scratch/gbuehler/bioinfo/Annotations/Homer"

  #run HOMER
  #system(command="PATH=$PATH:/work/gbioinfo/Appz/Homer/Homer-current/bin")
  #system(command="PATH=$PATH:/work/gbioinfo/Appz/weblogo/weblogo-current")

  system2(command="/work/gbioinfo/Appz/Homer/Homer-current/bin/findMotifsGenome.pl",args=c(bedfile,genomepath,outputfileH,size,cores,noknown,HomerDir),wait=FALSE)
 # print(unquote(c("/work/gbioinfo/Appz/Homer/Homer-current/bin/findMotifsGenome.pl",bedfile,genomepath,outputfileH,size,cores,noknown,HomerDir)))
}
}
}
