#' @title plotTracks
#'
#' @description Plots genome browser like tracks.
#'
#' @details Plots genome browser like tracks for supplied data (stranded or unstranded, eg RNAseq, CHIPseq) in bigwig format, 
#' for a region surrounding a specified gene, or a region specified by a GRanges object. Including a track for all transcripts in that region.
#'
#' @param bwFilesPlus Character vector containing the filenames (including the full path) of read alignment files in bigwig format.
#' Either the only bigwig file per sample if unstranded (eg for ChIPseq) or the bigwig file containing the reads mapped on the plus strand.
#' @param bwFilesMinus Character vector containing the filenames (including the full path) of read alignment files in bigwig format,
#' containing the reads mapped on the minus strand. If missing, only the \code{bwFilesPlus} files will be used. 
#' @param bwNames Character vector containing the names to describe the samples (for example: "H3K9me3"). 
#' @param txdb TxDb object containing all transcripts/genes using Ensembl Gene IDs.
#' @param EnsDb EnsDb object containing an Ensembl annotation database linking GeneIDs to gene symbols. 
#' @param bedFiles Optional character vector containing the filenames (including the full path) of any bed files (eg peaks) to be plotted. 
#' They should contain columns: chrom, start,stop, gene, score, strand (-1 or 1),type (exon).  
#' @param gene Character scalar specifying the gene symbol for the gene around which the tracks should be plotted.
#' If missing, a \code{plotregion} needs to be supplied instead.
#' @param plotregion A GRanges object containing the region of the genome for plotting the tracks (in case \code{gene} is not supplied).
#' @param plotExtension Integer scalar specifying the length (in bp) up and downstream of the gene the plot should be extended. Default = 1000.
#' @param plotranges Numeric vector of the same length as bwNames specifying the y-axis range for each plot. Default=1 for each sample. 
#'
#' @return A plot of genome browser like tracks, including transcripts.
#'
#' @importFrom GenomicRanges start end strand seqnames
#' @importFrom rtracklayer import
#' @importFrom Sushi plotBedgraph SushiColors plotGenes labelgenome
#' @importFrom graphics mtext axis
#' @importFrom stringr str_detect fixed
#' @importFrom GenomicFeatures exons exonsBy genes
#' @importFrom AnnotationDbi select
#' @importFrom dplyr left_join
#'
#' @examples
#' library(TxDb.Mmusculus.UCSC.mm10.ensGene)
#' library(EnsDb.Mmusculus.v79)
#' bwFiles <- list.files(system.file("extdata", 
#' package = "MiniChip"), full.names=TRUE,pattern="*bw$")
#' plotTracks(bwFilesPlus=bwFiles,bwNames=c(rep("Adnp",2)),
#'  txdb=TxDb.Mmusculus.UCSC.mm10.ensGene,EnsDb=EnsDb.Mmusculus.v79,
#' gene="Adnp", plotExtension=500,plotranges=rep(0.8,2))

#' @export
plotTracks <- function(bwFilesPlus,bwFilesMinus,bwNames,txdb,EnsDb,bedFiles,gene,plotregion,
                       plotExtension=1000,plotranges=rep(1,length(bwNames))){

    if (!requireNamespace("Sushi", quietly = TRUE)) {
    stop("Package \"Sushi\" needed for this function to work. Please install it.",
         call. = FALSE)
    }
  
  ######select the gene to plot and generate transcripts bed file#####
  #obtain a GRanges object of genes from the txdb
  genes <- genes(txdb)
  
  #get gene info (gene symbols) from the Ensembl database
  gene_names <- genes(EnsDb, columns=c("gene_name", "gene_biotype"))
 
   #if the txdb gene IDs contain a . (usually . and a number at the end of the ID, for example in txdb derived from GENCODE),
  #but the Ensembl database (gene_names) does not, remove the .# at the end of the Gene IDs
  if(str_detect(genes$gene_id[1], fixed("."))==TRUE & str_detect(gene_names$gene_id[1], fixed("."))==FALSE){
  gene_id <- matrix(unlist(strsplit(genes$gene_id,".",fixed=TRUE)),ncol=2,byrow = TRUE)[,1]
  } else {
  gene_id <- genes$gene_id
  }
  #add the gene symbol to the genes GRanges object
  genes2 <- left_join(data.frame(gene_id=gene_id,GeneID=genes$gene_id),data.frame(gene_names),by="gene_id")
  genes$gene_name <- genes2$gene_name
  
  if(!missing(gene)){
  #generate a GRanges object of only the selected gene
  geneinx <- which(genes$gene_name==gene)
  genes.gr <- genes[geneinx]
  
  } else {
  #if no gene was specified, but instead a genomic region was specified, extract the region into a genes.gr object
  genes.gr <- plotregion
  }
  
  #extend the plot region by the specifed amount of bp
  start(genes.gr) <- start(genes.gr) -plotExtension
  end(genes.gr) <- end(genes.gr) + plotExtension
  
  
  #extract chr, start,end from gene GRanges object
  chrom <- as.character(seqnames(genes.gr))
  chromstart <- start(genes.gr)
  chromend <- end(genes.gr)
  genestrand <- strand(genes.gr)
  
  
  #generate GRanges object of transcripts of selected chromosomal region
  exons <- unlist(exonsBy(txdb,by="tx",use.names=TRUE))
  exons2 <- exons[seqnames(exons) == chrom & start(exons) < chromend & end(exons) > chromstart]
  
  #add the gene symbol to each trasncript name
  exon_gene_IDs <- select(txdb, keys=unique(names(exons2)), columns="GENEID", keytype="TXNAME")
  exon_gene_gr <- unique(data.frame(genes[exon_gene_IDs$GENEID]))
  exon_gene_gr2 <- left_join(exon_gene_IDs,exon_gene_gr,by=c("GENEID"="gene_id"))
  exons2$gene_symbol <- left_join(data.frame(TXNAME=names(exons2)),exon_gene_gr2,by=c("TXNAME"))[,"gene_name"]
  exons2$ID <- paste(exons2$gene_symbol,names(exons2),sep=" ")
  
  #generate bed file of transcripts of selected chromosomal region
  exons.bed <- data.frame(chrom=seqnames(exons2),start=start(exons2),stop=end(exons2),
                          gene=exons2$ID,score=rep(1,length(exons2)), 
                          strand=ifelse(strand(exons2)=="-",-1,1), 
                          type=rep("exon",length(exons2)))
  
  #set the plotting area
  if(!missing(bedFiles)){
  par(mfrow=c(length(bwNames)+1+length(bedFiles),1),mar=c(3,5,2,1))
  } else {
    par(mfrow=c(length(bwNames)+1,1),mar=c(3,5,2,1))
  }
  
  #seq along the bigwig samples provided
  for (i in seq_along(bwNames)){
    #load bigwig data to plot 
    if(!missing(bwFilesMinus)){
   bg.minus <- as.data.frame(import(con=bwFilesMinus[i], format="bw",which=genes.gr))[,c(1:3,6)]
   range=c(-plotranges[i],plotranges[i])
    } else {
      range=c(0,plotranges[i])
    }
   bg.plus <- as.data.frame(import(con=bwFilesPlus[i], format="bw",which=genes.gr))[,c(1:3,6)]
  
   #plot the sample
   plotBedgraph(bg.plus, chrom, chromstart, chromend, range=range,transparency=0.5, color=SushiColors(2)(4)[1])
   if(!missing(bwFilesMinus)){
   plotBedgraph(bg.minus, chrom, chromstart, chromend, range=range,transparency=0.5, color=SushiColors(2)(4)[2], overlay=TRUE, rescaleoverlay=FALSE)
   }
   mtext(bwNames[i],side=2,line=3,cex=1,font=1)
   axis(side=2,las=2,tcl=.2,line=0.5)
  }
  
  #plot the transcripts
  if (nrow(exons.bed) > 0){
  pg = plotGenes(exons.bed,chrom,chromstart,chromend,
                 types = exons.bed$type,
                 colorby=log10(exons.bed$score+0.001),
                 colorbycol= SushiColors(5),colorbyrange=c(0,1.0),
                 labeltext=TRUE,maxrows=50,height=0.4,plotgenetype="box")
  } else {
    message("No transcripts found in plot region!")
  }
  
  #plot additional bed files
  if(!missing(bedFiles)){
  for (i in seq_along(bedFiles)){
  bed_file <- read.table(file=bedFiles[i],sep="\t",header=TRUE)
  bed_file <- bed_file[bed_file$chrom==chrom,]
  pg = plotGenes(bed_file,chrom,chromstart,chromend ,
                 types = bed_file$type,
                 colorby=bed_file$score,
                 colorbycol= SushiColors(5),colorbyrange=c(0,20),
                 labeltext=TRUE,maxrows=50,height=0.4,plotgenetype="box")
  }
  }
  
  #plot the scale bar
  labelgenome( chrom, chromstart,chromend,n=5,scale="Kb")
}
