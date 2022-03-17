\encoding{UTF-8}
\name{chem16S_plot}
\alias{chem16S_plot}
\alias{taxacomp}
\alias{plotmet}
\alias{addhull}
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
  plotmet(study, cn = FALSE, identify = FALSE, title = TRUE,
    xlim = NULL, ylim = NULL, plot.it = TRUE, points = TRUE, lines = FALSE,
    lineage = NULL, mincount = 200, pch1 = 1, pch2 = 21, dropNA = TRUE, return = "data",
    extracolumn = NULL, add = FALSE, plot.bg = TRUE)
  addhull(x, y, basecol, outline = FALSE, ...)
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
  \item{study}{study key corresponding to file names in \code{extdata/chem16S} (e.g. \samp{HCW+13})}
  \item{cn}{logical, apply copy-number correction to RDP classifier results?}
  \item{title}{character, plot title}
  \item{plot.it}{logical, make a plot?}
  \item{points}{logical, plot points?}
  \item{lines}{logical, plot lines?}
  \item{lineage}{character, restrict plot to this taxonomic lineage}
  \item{mincount}{integer, minimum classification count for each sample}
  \item{pch1}{numeric, symbol for samples in group 1}
  \item{pch2}{numeric, symbol for samples in group 2}
  \item{dropNA}{logical, exclude samples with NA name in metadata?}
  \item{return}{character, indicates whether to return \samp{data} values or group \samp{means}}
  \item{extracolumn}{character, the name of one or more extra columns (from \code{\link{getmdat}}) to include in the output}
  \item{add}{logical, add to existing plot?}
  \item{plot.bg}{logical, plot lines for background correlation from basis species?}
  \item{x}{x values}
  \item{y}{y values}
  \item{basecol}{color used for hull lines}
  \item{outline}{logical, draw hull lines?}
  \item{...}{arguments passed to \code{\link{polygon}} (by \code{addulll}) or \code{getgroup} (by \code{groupmet})}
}

\details{

\code{getrefseq} retrieves amino acid compositions and taxonomic names from the files \code{extdata/refseq/protein_refseq.csv.xz} and \code{extdata/refseq/taxid_names.csv.xz}.
By default, any species with greater than 20000 reference sequences is omitted; set \code{filterspecies} to FALSE to include all available species.

\code{taxacomp} plots the \ZC and \nH2O of taxonomic groups in RefSeq and their children.
\code{groups} can be the names of one or more taxonomic groups, or one of the following keywords: \samp{majorphyla} (phyla with more than 500 representatives), \samp{majorcellular} (non-viral (cellular) phyla with more than 60 representatives), \samp{proteobacteria}, \samp{acidobacteria}, or \samp{cyanobacteria} (classes within the indicated phylum and their orders).
Note that the keywords are all lower-case; use the upper-case name (e.g. \samp{Proteobacteria}) to plot the phylum and its orders.

\code{plotmet} plots the \ZC and \nH2O of community proteomes estimated by combining taxonomic abundances from 16S rRNA gene sequencing with RefSeq proteomes.
If \code{pch1} and \code{pch2} are provided, then two sample groups are identified by their plotting symbols in the metadata (see \code{\link{getmdat}}).
Star-shaped symbols are plotted to indicate the mean values of metrics for each sample group.

\code{addhull} uses \code{\link{polygon}} to draw a convex hull around the points given in \code{x}, \code{y}.
The lines are draw only in \code{outline} is TRUE; otherwise, the area of the hull is filled.
The fill color is a lightened (semi-transparent) shade of the color given in \code{basecol}.

}

\value{
\code{getgroup} returns a list with elements \samp{study}, \samp{param} (study key and chemical metric), \samp{pch1}, \samp{pch2}, \samp{col1}, \samp{col2} (single plotting symbols and color for each sample group), \samp{n1}, \samp{n2} (number of samples in each group), \samp{X1}, \samp{X2} (values of the selected chemical metric for individual taxa within each sample group), \samp{P1}, \samp{P2} (percent abundances of the taxa within each sample group), \samp{Pboth} (percent abundance of each taxon in both sample groups (combined)), \samp{DX}, (contribution of each taxon to the difference between sample groups), \samp{Xall}, (values of the metric for all taxa taken together in sample groups 1 and 2), \samp{taxon} (taxon names).

\samp{DX} is calculated as \code{(X2 - Xall[1]) * P2 / 100 - (X1 - Xall[1]) * P1 / 100}.
Note that \code{sum(DX) / diff(Xall)} equals 1.
}

\seealso{
\code{\link{chem16S}} for the functions used to calculate the chemical metrics.
}

\references{
Herlemann, D. P. R., Lundin, D., Andersson, A. F., Labrenz, M. and Jürgens, K. (2016) Phylogenetic signals of salinity and season in bacterial community composition across the salinity gradient of the Baltic Sea. \emph{Front. Microbiol.} \bold{7}, 1883. \doi{10.3389/fmicb.2016.01883}
}

\examples{
# Default plot for Bacteria and Archaea domains and their children (i.e. phyla)
taxacomp()
# Plot major phyla within Archaea, Bacteria, and viruses
taxacomp("majorphyla")

# Make a plot for the Baltic Sea salinity gradient
# (data from Herlemann et al., 2016)
pm <- plotmet("HLA+16")
# Add a legend
legend <- c("< 6 PSU", "6-20 PSU", "> 20 PSU")
pch <- c(24, 20, 21)
pt.bg <- c(3, NA, 4)
legend("bottomright", legend, pch = pch, col = 1, pt.bg = pt.bg, bg = "white")
# Identify samples with low and high salinity
mdat <- getmdat("HLA+16")
ilo <- mdat$salinity < 6
ihi <- mdat$salinity > 20
# Add convex hulls
addhull(pm$ZC[ilo], pm$nH2O[ilo], "green3")
addhull(pm$ZC[ihi], pm$nH2O[ihi], "blue")

# Some samples (i.e. outliers) have an NA name in the metadata
# to exclude them from further calculations.
# Use dropNA = FALSE to show all samples.
plotmet("UKD+18.sediment", dropNA = FALSE)

# Make a plot for only Proteobacteria
plotmet("HLA+16", lineage = "Proteobacteria")

# Show points for all samples and larger points for
# mean values of high- and low-salinity samples
plotmet("HLA+16", pch1 = 21, pch2 = 24)
}
