\encoding{UTF-8}
\name{plot_starburst}
\alias{plot_starburst}
\title{Starburst plot for parent and children taxa}
\description{
Plots chemical metrics of reference proteomes for parent taxa and their children.
Lines are drawn from each parent to all of its children, creating a starburst pattern.
}

\usage{
  plot_starburst(taxa = c("Bacteria", "Archaea"), metrics = c("Zc", "nH2O"),
    refdb = "RefSeq", remove_species_20000 = TRUE, terminal_H2O = 0,
    xlim = NULL, ylim = NULL, pch = NULL, col = seq_along(taxa), lcol = NULL,
    lwd = rep(1, length(taxa)), hline = NULL, legend.x = NA, identify = FALSE)
}

\arguments{
  \item{taxa}{character, names of parent taxa}
  \item{metrics}{character, two chemical metrics (see \code{\link[canprot]{metrics}})}
  \item{refdb}{character, reference database (\samp{RefSeq} or \samp{GTDB})}
  \item{remove_species_20000}{logical, remove species with > 20000 sequences?}
  \item{terminal_H2O}{numeric, how to count terminal groups for calculation of \code{\link[canprot]{nH2O}}}
  \item{xlim}{numeric, x axis limits}
  \item{ylim}{numeric, y axis limits}
  \item{pch}{plotting symbols}
  \item{col}{colors for plotting symbols and lines}
  \item{lcol}{colors for lines; default corresponds to \code{\link{adjustcolor}(col, alpha.f = 0.5)}}
  \item{lwd}{line widths}
  \item{hline}{numeric, coordinates of horizontal lines to add}
  \item{legend.x}{position of \code{\link{legend}}; default is to draw no legend}
  \item{identify}{logical, run \code{\link{identify}} for interactive identification of points?}
}

\details{

\code{plot_starburst} plots the indicated chemical \code{metrics} of the reference proteome for each taxon in \code{taxa}.
Lines are drawn from each taxon to all of their children at the next lowest taxonomic rank.
Use \code{refdb} to change the reference database; the default retrieves amino acid compositions and taxonomic names from the files \code{extdata/RefDB/RefSeq/genome_AA.csv.xz} and \code{extdata/RefDB/RefSeq/taxonomy.csv.xz}.
\code{remove_species_20000} removes all species with more than 20000 reference sequences (applies only for \code{refdb == "RefSeq"}).

\code{terminal_H2O} and \code{hline} are provided to enable reproduction of the published figure for \code{\link{geo16S1}}.
In that figure, the \nH2O values include the contribution of terminal groups (i.e., 1 \H2O per protein, divided by the length of the protein).
Subsequently, in other publications, the contribution of terminal groups to \nH2O has been suppressed.

}

\examples{
# Default plot for Bacteria and Archaea (parents)
# with lines drawn to phyla in each domain (children)
plot_starburst(legend.x = "bottomleft")
}
