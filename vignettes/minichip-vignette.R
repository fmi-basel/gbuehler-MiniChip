## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
suppressPackageStartupMessages({
   library(MiniChip)
  library(ComplexHeatmap)
  library(GenomicFeatures)
  library(ggplot2)
  library(reshape2)
  library(BSgenome.Mmusculus.UCSC.mm10)
  library(TxDb.Mmusculus.UCSC.mm10.knownGene)
})

## -----------------------------------------------------------------------------
peaks1.d <- read.table(system.file("extdata", "Adnp_rep1_chr11_peaks.narrowPeak", package = "MiniChip"),header=FALSE)
peaks2.d <- read.table(system.file("extdata", "Adnp_rep2_chr11_peaks.narrowPeak", package = "MiniChip"),header=FALSE)
peaks3.d <- read.table(system.file("extdata", "Adnp_rep3_chr11_peaks.narrowPeak", package = "MiniChip"),header=FALSE)

names(peaks1.d) <- c("chr","start","end","name","score","empty","foldchange","pvalue","qvalue","summit")
names(peaks2.d) <- c("chr","start","end","name","score","empty","foldchange","pvalue","qvalue","summit")
names(peaks3.d) <- c("chr","start","end","name","score","empty","foldchange","pvalue","qvalue","summit")

peaks1 <- makeGRangesFromDataFrame(peaks1.d,
                                   keep.extra.columns=TRUE,
                                   ignore.strand=TRUE,
                                   seqinfo=NULL,
                                   seqnames.field=c("chr"),
                                   start.field=c("start"),
                                   end.field=c("end"),
                                   starts.in.df.are.0based=FALSE)

peaks2 <- makeGRangesFromDataFrame(peaks2.d,
                                   keep.extra.columns=TRUE,
                                   ignore.strand=TRUE,
                                   seqinfo=NULL,
                                   seqnames.field=c("chr"),
                                   start.field=c("start"),
                                   end.field=c("end"),
                                   starts.in.df.are.0based=FALSE)

peaks3 <- makeGRangesFromDataFrame(peaks3.d,
                                   keep.extra.columns=TRUE,
                                   ignore.strand=TRUE,
                                   seqinfo=NULL,
                                   seqnames.field=c("chr"),
                                   start.field=c("start"),
                                   end.field=c("end"),
                                   starts.in.df.are.0based=FALSE)

## -----------------------------------------------------------------------------
# peaks that occur in at least 2/3 replicates
peaks3.1 <- subsetByOverlaps(peaks3,peaks1)
peaks3.2 <- subsetByOverlaps(peaks3,peaks2)
peaks1.2 <- subsetByOverlaps(peaks1,peaks2)

#keep all 3.1 peaks, find 3.2 and 1.2 peaks that do not overlap 3.1 peaks

# peaks which are found in 3 and 2, but not in 3 and 1, and not in 1 and 2
peaks3.2u <- peaks3.2[overlapsAny(peaks3.2,peaks3.1) ==FALSE & overlapsAny(peaks3.2,peaks1.2) == FALSE]
# peaks which are found in 1 and 2, but not in 3 and 1, and not in 3 and 2
peaks1.2u <- peaks1.2[overlapsAny(peaks1.2,peaks3.1) ==FALSE & overlapsAny(peaks1.2,peaks3.2) == FALSE]
# add together all 3.1 peaks as well as unique 3.2 and 1.2 peaks
peaks <- c(peaks3.1,peaks3.2u,peaks1.2u)

## -----------------------------------------------------------------------------
peaknames <- "Adnp_peaks"
peaks

## -----------------------------------------------------------------------------
#  randomize peak locations
nSites <- length(peaks)
peak.widths <- width(peaks)
random.peaks <-  SimulatePeaks(nSites,peak.widths,chromosomeSizes=system.file("extdata", "chrNameLength_mm10_chr11.txt", package = "MiniChip"))

## -----------------------------------------------------------------------------
#select bam files from Adnp experiment
bamFiles <- list.files(system.file("extdata", package = "MiniChip"), full.names=TRUE,pattern="*bam$")
#all.bamFiles <- list.files("/tungstenfs/scratch/gbuehler/deepSeqRepos/bam/", full.names=TRUE,pattern="*bam$")
#all.bamFiles <- grep("test",all.bamFiles,value=TRUE,invert=TRUE)
#bamFiles <- grep("1950F",all.bamFiles,value=TRUE)[1:8]
bamFiles

## -----------------------------------------------------------------------------
bamNames <- gsub(paste(system.file("extdata", package = "MiniChip"),"/",sep=""),"",bamFiles)
bamNames <- gsub("_chr11.bam","",bamNames)
bamNames

## ----  fig.height=5, fig.width=10---------------------------------------------
GCbias(bamFiles[1:2],bamNames[1:2],pe="none", restrict="chr11",GCprob=TRUE, genome=BSgenome.Mmusculus.UCSC.mm10)

## -----------------------------------------------------------------------------
span <- 3025
step <- 50
summits <- peaks
start(summits) <- start(summits) + summits$summit
end(summits) <- start(summits)
names(summits) <- elementMetadata(summits)[,"name"]

counts <- SummitHeatmap(summits,bamFiles,bamNames,span,step,useCPM=TRUE,readShiftSize=75)

## ----heatmap, fig.height=10, fig.width=8--------------------------------------
ht_list <- DrawSummitHeatmaps(counts[c(1,2,3)], bamNames[c(1,2,3)],orderSample=1, use.log=TRUE,
                              bottomCpm=c(0,0,0), medianCpm = c(2,2,2), topCpm = c(4,4,4), orderWindows=2, TargetHeight=500)
draw(ht_list, padding = unit(c(3, 8, 8, 2), "mm"))

## -----------------------------------------------------------------------------
inx.ChIP <- grep("ChIP",names(counts))
Adnp <- grep("Adnp",names(counts[inx.ChIP]),value=TRUE)
Input <- grep("Input",names(counts),value=TRUE)[1:3]

sampleList <- list(Adnp=Adnp,Input=Input)
meanCounts <- SummarizeHeatmaps(counts,sampleList)

## ---- fig.height=5, fig.width=10----------------------------------------------
#select peaks that overlap a TSS (+/- 1kb)
#txdb=loadDb("/tungstenfs/scratch/gbuehler/michi/Annotations/GENCODE/Mouse/release_M23/gencode.vM23.annotation.txdb.sqlite")
txdb=TxDb.Mmusculus.UCSC.mm10.knownGene
TSSs <- promoters(genes(txdb),upstream=1000,downstream=1000)
summit_TSS_overlap <- overlapsAny(summits,TSSs)
summit_TSS_overlap_names <- names(summits[summit_TSS_overlap == TRUE])

#cumulative plots
CumulativePlots(meanCounts,bamNames = names(meanCounts),
                                   span=3025,step=50,summarizing = "mean",overlapNames = summit_TSS_overlap_names)

## ---- fig.height=10, fig.width=8----------------------------------------------
promoters <- promoters(txdb,upstream=150,downstream=150)
annoname <- "promoters"
promoters.overlap <- AnnotationHeatmap(summits,annotation=promoters,annoname,span,step)
annos <- list(promoters.overlap)
names(annos) <- c("promoters")
allCounts <- c(meanCounts,annos)

cols <- c("darkblue","darkred","darkgreen")
upper.cpm <- c(rep(2,2),rep(1,1))

heatmap_list <- DrawSummitHeatmaps(allCounts, names(allCounts), plotcol= cols, medianCpm = upper.cpm, orderSample = 1, use.log=TRUE, summarizing = "mean", orderWindows=2,MetaScale=c("all","all","individual"), TargetHeight=500)
draw(heatmap_list, padding = unit(c(3, 8, 8, 2), "mm"),show_heatmap_legend=FALSE)

