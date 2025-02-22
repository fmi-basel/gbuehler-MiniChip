% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SimulatePeaks.R
\name{SimulatePeaks}
\alias{SimulatePeaks}
\title{SimulatePeaks}
\usage{
SimulatePeaks(nSites, peak.widths, chromosomeSizes)
}
\arguments{
\item{nSites}{Integer scalar giving the number of desired output ranges.}

\item{peak.widths}{A numeric vector of the length given in \code{nSites} containing the desired width of the ranges.}

\item{chromosomeSizes}{A table with two columns: chromosome name, length of the chromosome.
This defines the search space for choosing random regions.}
}
\value{
A GRanges object of randomly chosen genomic ranges of length \code{nSites} with widths \code{peak.widths}.
}
\description{
This function generates a GRanges object with randomly chosen ranges
in the genome.
}
\details{
This function can be used to shuffle peak regions of your
ChIP data to random locations in the mouse genome. These random peaks
can be used as control regions for further analysis, for example, when
testing the overlap of your peaks with genomic annotations (genes, repeats,...).
}
\examples{
SimulatePeaks(1000,rep(100,1000),chromosomeSizes=
system.file("extdata", "chrNameLength_mm10_chr11.txt", package = "MiniChip"))

}
