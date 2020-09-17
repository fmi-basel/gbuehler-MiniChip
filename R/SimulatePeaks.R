#' @title SimulatePeaks
#'
#' @description This function generates a GRanges object with randomly chosen ranges
#'  in the genome.
#'
#' @details This function can be used to shuffle peak regions of your
#'  ChIP data to random locations in the mouse genome. These random peaks
#'  can be used as control regions for further analysis, for example, when
#'  testing the overlap of your peaks with genomic annotations (genes, repeats,...).
#'
#' @param nSites Integer scalar giving the number of desired output ranges.
#' @param peak.widths A numeric vector of the length given in \code{nSites} containing the desired width of the ranges.
#' @param chromosomeSizes A table with two columns: chromosome name, length of the chromosome. 
#' This defines the search space for choosing random regions.
#'
#' @return A GRanges object of randomly chosen genomic ranges of length \code{nSites} with widths \code{peak.widths}.
#'
#' @examples
#' SimulatePeaks(1000,rep(100,1000),chromosomeSizes=system.file("extdata", "chrNameLength_mm10_chr11.txt", package = "MiniChip"))
#'
#' @importFrom data.table data.table
#' @importFrom data.table :=
#' @importFrom data.table .N
#' @importFrom utils read.table
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#'
#' @export
SimulatePeaks <- function(nSites,peak.widths,chromosomeSizes){ 
 
   if (!requireNamespace("data.table", quietly = TRUE)) {
  stop("Package \"data.table\" needed for this function to work. Please install it.",
       call. = FALSE)
   }
  
  chromSizes <- read.table(chromosomeSizes)
  # define chromosomes and regions on them where peaks can fall onto without obtaining negative starting positions in the end
  chromosomes <- data.table(chromosome=chromSizes$V1[chromSizes$V2 > max(peak.widths)], start=max(peak.widths), end=chromSizes$V2[chromSizes$V2 > max(peak.widths)])
  # calculate size of chromosome
  chromosomes[, size := 1 + end-start]
  # Randomly sample chromosome file rows, proportional to the length of each chromosome
  simulated.peaks <- chromosomes[sample(.N, size=nSites, replace=TRUE, prob=chromosomes$size)]
  # Randomly sample uniformly within each chosen chromosome
  simulated.peaks[, position := sample(start:end, size=1), by=1:dim(simulated.peaks)[1]]
  # Remove extra columns, format, and adjust lengths of simulated peaks according to peaks lengths
  simulated.peaks[, end  := position]
  simulated.peaks[, start := position - peak.widths]
  simulated.peaks[, c("size", "position") := NULL]
  # turn it into a Granges object
  random.peaks <- makeGRangesFromDataFrame(data.frame(simulated.peaks),
                                           keep.extra.columns=TRUE,
                                           ignore.strand=TRUE,
                                           seqinfo=NULL,
                                           seqnames.field=c("chromosome"),
                                           start.field=c("start"),
                                           end.field=c("end"),
                                           starts.in.df.are.0based=FALSE)
  names(random.peaks) <- paste(seqnames(random.peaks),start(random.peaks),sep="_")
  return(random.peaks)
}
