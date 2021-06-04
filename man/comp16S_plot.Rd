\encoding{UTF-8}
\name{comp16S_plot}
\alias{comp16S_plot}
\alias{taxacomp}
\alias{plotcomp}
\alias{addhull}
\alias{groupcomp}
\alias{getrefseq}
\title{Plot compositional metrics of taxa and community proteomes inferred from 16S rRNA gene sequences}
\description{
Functions to plot compositional metrics (\ZC and \nH2O) for inferred microbial community proteomes.
}

\usage{
  getrefseq(filterspecies = TRUE)
  taxacomp(groups = c("Bacteria", "Archaea"), xlim = NULL, ylim = NULL,
    col = seq_along(groups), legend.x = "topleft", identify = FALSE,
    pch = NULL, hline = NULL, filterspecies = TRUE)
  plotcomp(study, cn = FALSE, identify = FALSE, title = TRUE,
    xlim = NULL, ylim = NULL, plot.it = TRUE, points = TRUE, lines = FALSE,
    lineage = NULL, pch.up = 21, pch.down = 1, pval = TRUE)
  addhull(x, y, basecol, outline = FALSE, ...)
  groupcomp(study = "XDZ+17", metric = "nH2O", rank = "domain",
    pch.up = 24, pch.down = 21, minpercent = 2, xlim = NULL, ylim = NULL,
    xadj = NULL, yadj = NULL, study2 = NA, mdat = NULL, map = NULL, RDP = NULL)
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
  \item{study}{study key corresponding to file names in \code{extdata/comp16S} (e.g. \samp{HCW+13})}
  \item{cn}{logical, apply copy-number correction to RDP classifier results?}
  \item{title}{character, plot title}
  \item{plot.it}{logical, make a plot?}
  \item{points}{logical, plot points?}
  \item{lines}{logical, plot lines?}
  \item{lineage}{character, restrict plot to this taxonomic lineage}
  \item{pch.up}{numeric, symbol for samples in \dQuote{up-expressed} group}
  \item{pch.down}{numeric, symbol for samples in \dQuote{down-expressed} group}
  \item{pval}{logical, show log10 of p-value in axis labels?}
  \item{x}{x values}
  \item{y}{y values}
  \item{basecol}{color used for hull lines}
  \item{outline}{logical, draw hull lines?}
  \item{...}{other arguments passed to \code{\link{polygon}}}
  \item{metric}{compositional metric, \samp{ZC} or \samp{nH2O}}
  \item{rank}{character, taxonomic rank}
  \item{minpercent}{numeric, include taxa that are at least this abundant}
  \item{xadj}{numeric, x adjustment for text labels}
  \item{yadj}{numeric, y adjustment for text labels}
  \item{study2}{character, calculate differences between two studies instead of groups in one study}
  \item{mdat}{data frame, output from \code{\link{getmdat}}}
  \item{map}{data frame, output from \code{\link{getmap}}}
  \item{RDP}{data frame, output from \code{\link{getRDP}}}
}

\details{

\code{getrefseq} retrieves amino acid compositions and taxonomic names from the files \code{extdata/refseq/protein_refseq.csv.xz} and \code{extdata/refseq/taxid_names.csv.xz}.
By default, any species with greater than 20000 reference sequences is omitted; set \code{filterspecies} to FALSE to include all available species.

\code{taxacomp} plots the \ZC and \nH2O of taxonomic groups in RefSeq and their children.
\code{groups} can be the names of one or more taxonomic groups, or one of the following keywords: \samp{majorphyla} (phyla with more than 500 representatives), \samp{majorcellular} (non-viral (cellular) phyla with more than 60 representatives), \samp{proteobacteria}, \samp{acidobacteria}, or \samp{cyanobacteria} (classes within the indicated phylum and their orders).
Note that the keywords are all lower-case; use the upper-case name (e.g. \samp{Proteobacteria}) to plot the phylum and its orders.

\code{plotcomp} plots the \ZC and \nH2O of community proteomes inferred by combining taxonomic abundances from 16S rRNA gene sequencing with RefSeq proteomes.

\code{addhull} uses \code{\link{polygon}} to draw a convex hull around the points given in \code{x}, \code{y}.
The lines are draw only in \code{outline} is TRUE; otherwise, the area of the hull is filled.
The fill color is a lightened (semi-transparent) shade of the color given in \code{basecol}.

\code{groupcomp} makes composition-abundance plots for sample groups within taxonomic groups.
The sample groups are identified by their plotting symbols (\code{pch.up} and \code{pch.down}) in the metadata (see \code{\link{getmdat}}).

}

\seealso{
\code{\link{comp16S}} for the functions used to calculate the compositional metrics.
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
pc <- plotcomp("HLA+16")
# Add a legend
legend <- c("< 6 PSU", "6-20 PSU", "> 20 PSU")
pch <- c(24, 20, 21)
pt.bg <- c(3, NA, 4)
legend("bottomright", legend, pch = pch, col = 1, pt.bg = pt.bg, bg = "white")
# Add the convex hull - lines then area
addhull(pc$ZC, pc$nH2O, "blue", TRUE)
addhull(pc$ZC, pc$nH2O, "slategray1")

# Make a plot for only Proteobacteria
plotcomp("HLA+16", lineage = "Proteobacteria")

# Plot nH2O differences for individual phyla
groupcomp("HLA+16", rank = "phylum")
}
