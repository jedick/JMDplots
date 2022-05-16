\encoding{UTF-8}
\name{orp16S_plot}
\alias{orp16S_plot}
\alias{plotEZ}
\alias{plotMA}
\title{Plotting functions for orp16S paper}
\description{
Make plots of \ZC vs Eh7 and \ZC vs abundance of Most Abundant Mapped Taxon.
}

\usage{
  plotEZ(study, lineage = NULL, mincount = 100, pch = NULL, col = NULL,
    add = FALSE, type = "p", groupby = NULL, groups = NULL, legend.x = "topleft",
    show = c("lm", "points"), col.line = "gray62", lwd = 1, cex = 1,
    title.line = NA, dxlim = c(0, 0), dylim = c(0, 0), size = NULL)
  plotMA(study, lineage = NULL, mincount = 100, pch = NULL, col = NULL,
    groupby = NULL, groups = NULL, legend.x = "topright")
}

\arguments{
  \item{study}{character, study key}
  \item{lineage}{character, regular expression for filtering lineages, or the keyword \samp{two}}
  \item{mincount}{integer, samples with less than this number of RDP classifications are excluded}
  \item{pch}{point symbol}
  \item{col}{color(s) for points}
  \item{add}{logical, add to existing plot?}
  \item{type}{plot type, as in \code{\link{plot.default}}}
  \item{groupby}{character, column name for grouping samples}
  \item{groups}{character, names of sample groups}
  \item{legend.x}{legend placement (\code{x} argument for \code{\link{legend}})}
  \item{show}{character, show these plot elements}
  \item{col.line}{line color (linear fit)}
  \item{lwd}{line width}
  \item{cex}{character expansion}
  \item{title.line}{numeric, placement of title}
  \item{dxlim}{numeric, adjustment to x-axis limits}
  \item{dylim}{numeric, adjustment to y-axis limits}
  \item{size}{numeric, number of randomly selected samples (see \code{\link{getmdat_orp16S}})}
}

\details{

\code{plotEZ} is used to make \ZC-Eh7 scatterplots for a single study.
By default, \ZC is calculated for estimated proteomes using all bacterial and archaeal sequences.
Set \code{lineage} to \samp{Bacteria} or \samp{Archaea} to make one plot for one of the domains, or to \samp{two} to make two plots for both domains.
The latter option is used in \code{\link{orp16S_S1}} for the Figure S1 in the paper.


\code{plotMA} plots \ZC vs the percent abundance of most abundant mapped taxon (MAMT).
The plot title shows the percentages of MAMT and MAUT (most abundant unmapped taxon).

}

\value{
\code{plotEZ} returns a list with elements \code{study}, \code{envirotype} (environment type), \code{lineage}, \code{metadata} (output of \code{\link{getmdat_orp16S}}), and \code{EZdat}.
\code{EZdat} is a data frame with one row per sample and columns for \code{study}, \code{envirotype}, \code{lineage}, \code{sample} (sample name), \code{Run} (SRA run accession or other database ID), \code{Eh}, \code{Eh7}, and \code{ZC}.
If \code{groupby} and \code{groups} are given, \code{EZdat} additionally contains columns \code{groupby} and \code{group}.
If \code{show} includes \samp{lm}, the output additionally contains \code{EZlm} (the result of \code{\link{lm}(ZC ~ Eh7, EZdat)}), Eh7lim (range of Eh7 values), and ZCpred (fitted \ZC values for this Eh7 range).
}

\seealso{
\code{\link{orp16S}}
}

\examples{
plotEZ("MLL+19", "Bacteria", groupby = "Type", groups = c("Upland", "Paddy", "Sediment"))

plotMA("LMBA21", "Archaea")
}
