\encoding{UTF-8}
\name{geo16S_util}
\alias{geo16S_util}
\alias{taxacomp}
\alias{getrefseq}
\title{Plot chemical metrics of taxa and community proteomes estimated from 16S rRNA gene sequences}
\description{
Functions to plot chemical metrics (\ZC and \nH2O) for estimated microbial community proteomes.
}

\usage{
  getrefseq(filterspecies = TRUE)
  taxacomp(groups = c("Bacteria", "Archaea"), xlim = NULL, ylim = NULL,
    col = seq_along(groups), legend.x = "topleft", identify = FALSE,
    pch = NULL, hline = NULL, filterspecies = TRUE, lcol = NULL)
}

\arguments{
  \item{filterspecies}{logical, remove species with > 20000 sequences?}
  \item{groups}{character, names of taxonomic groups to plot}
  \item{xlim}{numeric, x axis limits}
  \item{ylim}{numeric, y axis limits}
  \item{col}{color for each group}
  \item{legend.x}{position of \code{\link{legend}}}
  \item{identify}{logical, run \code{\link{identify}} for interactive identification of points?}
  \item{pch}{plotting symbol}
  \item{hline}{numeric, add lines to show range of following plot}
  \item{lcol}{color for lines}
}

\details{

\code{getrefseq} retrieves amino acid compositions and taxonomic names from the files \code{extdata/RefSeq/genome_AA.csv.xz} and \code{extdata/RefSeq/taxonomy.csv.xz}.
By default, any species with greater than 20000 reference sequences is omitted; set \code{filterspecies} to FALSE to include all available species.

\code{taxacomp} plots the \ZC and \nH2O of taxonomic groups in RefSeq and their children.
\code{groups} can be the names of one or more taxonomic groups, or one of the following keywords: \samp{majorphyla} (phyla with more than 500 representatives), \samp{majorcellular} (non-viral (cellular) phyla with more than 60 representatives), \samp{proteobacteria}, \samp{acidobacteria}, or \samp{cyanobacteria} (classes within the indicated phylum and their orders).
Note that the keywords are all lower-case; use the upper-case name (e.g. \samp{Proteobacteria}) to plot the phylum and its orders.

}

\examples{
# Default plot for Bacteria and Archaea domains and their children (i.e. phyla)
taxacomp()
# Plot major phyla within Archaea, Bacteria, and viruses
taxacomp("majorphyla")
}
