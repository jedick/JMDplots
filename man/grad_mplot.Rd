\encoding{UTF-8}
\name{grad_mplot}
\alias{grad_mplot}
\alias{mplot}
\alias{mpage}
\alias{ppage}
\alias{mcomp}
\alias{pcomp}
\title{Plotting functions for chemical composition of metagenomes}
\description{
These are plotting functions for the redox gradients (\code{\link{gradox}}) and salinity gradients (\code{\link{gradH2O}}) papers.
The functions are exported in the package so they are available to user scripts.
}

\usage{
  mplot(study, seqtype, plottype = "bars", ylim = NULL, plot.RNA = TRUE,
    taxid = NULL, dsDNA = TRUE, abbrev = NULL, col = NULL, add.label = TRUE,
    maxdepth = NULL, H2O = FALSE, plot.it = TRUE, add.title = TRUE, yline = 2,
    basis = "QEC", datadir = NULL, mdata = studies)
  mpage(subset = "gradoxSI", H2O = FALSE, plottype = "bars", dsDNA = TRUE,
    set.par = TRUE, add.label = TRUE)
  ppage(subset = "gradoxSI", H2O = FALSE, set.par = TRUE, plot.it = TRUE,
    basis = "QEC", add.label = TRUE)
  mcomp(mout, yvar = "RNA")
  pcomp(mout, pout, seqtype = "MG", type = "ZC", parts = c("plot", "legend"),
    yline = 2, xlim = NULL, ylim = NULL, reorder = TRUE,
    plot.techtype = FALSE, add = FALSE)
}

\arguments{
  \item{study}{character, name of study}
  \item{seqtype}{character, description of sequence data type}
  \item{plottype}{character, type of plot}
  \item{ylim}{numeric, y-axis limits}
  \item{plot.RNA}{logical, plot RNA compositions?}
  \item{taxid}{numeric, select only this taxid}
  \item{dsDNA}{logical, use chemical composition of double-stranded DNA?}
  \item{abbrev}{character, abbreviations for studies}
  \item{col}{character, line color}
  \item{add.label}{logical, label the plot with MG (metagenome) or MT (metatranscriptome)?}
  \item{maxdepth}{numeric, maximum sample depth (meters)}
  \item{H2O}{logical, plot \nH2O instead of \ZC?}
  \item{plot.it}{logical, make a plot (use FALSE to return processed data without plotting)?}
  \item{add.title}{logical, add a title to the plot?}
  \item{yline}{numeric, margin position for y-axis label}
  \item{basis}{character, description of basis species to use}
  \item{subset}{character, name of subset of studies to plot}
  \item{set.par}{logical, set up plot with \code{\link{par}} commands?}
  \item{mout}{list, output from \code{mpage} or \code{ppage}}
  \item{yvar}{character, y-axis variable (\samp{RNA} or \samp{GC})}
  \item{pout}{list, output from \code{ppage}}
  \item{type}{character, type of plot (\samp{ZC} or \samp{both})}
  \item{parts}{character, parts of figure to include}
  \item{xlim}{numeric, x-axis limits}
  \item{reorder}{logical, put points in order of increasing x value?}
  \item{plot.techtype}{logical, identify 454 or Sanger datasets by open symbols?}
  \item{add}{logical, add to existing plot?}
  \item{datadir}{character, location of data directory (default is from package)}
  \item{mdata}{list, metadata for studies (default is from package)}
}

\details{
The \code{study} is composed of two names for the locality separated by an underscore.
The \code{seqtype} contains the name of the database followed by an underscore and MG (metagenome) or MT (metatranscriptome), i.e. \samp{[SRA|IMG|MGRAST|GenBank]_[MG|MT]}.
These can be changed to \samp{_MGP} and \samp{_MTP} to indicate protein sequences.

\code{mplot} plots values of \ZC for a particular combination of \code{study} and \code{seqtype}.
\ZC of both DNA and RNA are plotted, with an offset of -0.28 applied to that of RNA, so that values can be shown on the same scale (see Dick et al., 2019).

\code{mcomp} makes a scatterplot of (\ZC of RNA - \ZC of DNA) vs \ZC of DNA.
Set \code{yvar} to \samp{GC} to change the y-axis to GC ratio of DNA.
(For double-stranded DNA, this should be a linear function of \ZC.)

\code{datadir} and \code{mdata} are used to specify alternate location of compositional data and metadata for metagenomic studies.
This permits the processing of user-provided data sets; see \code{\link{subsurface}} for some examples.
}

\references{
Dick, J. M., Yu, M., Tan, J. and Lu, A. (2019) Changes in carbon oxidation state of metagenomes along geochemical redox gradients. \emph{Front. Microbiol.} \bold{10}, 120. \doi{10.3389/fmicb.2019.00120}
}

\examples{
# plot ZC of DNA and RNA
mplot("Bison_Pool", "IMG_MG")
# plot ZC of protein
mplot("Bison_Pool", "IMG_MGP")

# make a page of plots of ZC of DNA and RNA
mout <- mpage()
# compare ZC of RNA and DNA for these datasets
mcomp(mout)
# compare GC and ZC of DNA for these datasets
mcomp(mout, yvar = "GC")

# make a page of plots of ZC of proteins
pout <- ppage()
# compare ZC of protein and DNA for the metagenomes
pcomp(mout, pout)
# compare ZC of protein and DNA for the metatranscriptomes
pcomp(mout, pout, "MT")

## TODO: compare nH2O of protein and DNA
#mout <- mpage(H2O = TRUE)
#pout <- ppage(H2O = TRUE)
#pcomp(mout, pout, type = "H2O")

# compare ZC and nH2O of protein in a few redox gradients
mout <- ppage("gradoxGS")
pout <- ppage("gradoxGS", H2O = TRUE)
pcomp(mout, pout, type = "both", reorder = FALSE, yline = 2)
# overlay data for protein from the Baltic Sea
mout <- ppage("balticsurface", plot.it = FALSE)
pout <- ppage("balticsurface", H2O = TRUE, plot.it = FALSE)
pcomp(mout, pout, type = "both", reorder = FALSE, add = TRUE)
}
