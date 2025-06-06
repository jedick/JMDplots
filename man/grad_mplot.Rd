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
    datadir = NULL, mdata = studies, add = FALSE,
    pch = 19, var = NULL, srt = 45, ilabel = NULL)
  mpage(subset = "gradoxSI", H2O = FALSE, plottype = "bars", dsDNA = TRUE,
    set.par = TRUE, add.label = TRUE, mfrow = NULL)
  ppage(subset = "gradoxSI", H2O = FALSE, set.par = TRUE, plot.it = TRUE,
    add.label = TRUE, mfrow = NULL)
  mcomp(mout, yvar = "RNA")
  pcomp(mout, pout, seqtype = "MG", vars = NULL, parts = c("plot", "legend"),
    yline = 2, xlim = NULL, ylim = NULL, reorder = TRUE,
    plot.techtype = FALSE, add = FALSE, pch = NULL, lty = 2,
    labels.at = "max", cex.ylab = 1, font = 1, labdx = NULL, labdy = NULL)
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
  \item{H2O}{logical, plot \nH2O instead of \Zc?}
  \item{plot.it}{logical, make a plot (use FALSE to return processed data without plotting)?}
  \item{add.title}{logical, add a title to the plot?}
  \item{yline}{numeric, margin position for y-axis label}
  \item{var}{character, other variable to plot (\samp{GRAVY} or \samp{pI})}
  \item{srt}{numeric, string rotation for non-numeric labels}
  \item{ilabel}{numeric, which labels to include on plot}
  \item{subset}{character, name of subset of studies to plot}
  \item{set.par}{logical, set up plot with \code{\link{par}} commands?}
  \item{mfrow}{numeric, optional setting for \code{\link{par}(mfrow)}}
  \item{mout}{list, output from \code{mpage} or \code{ppage}}
  \item{yvar}{character, y-axis variable (\samp{RNA} or \samp{GC})}
  \item{pout}{list, output from \code{ppage}}
  \item{vars}{character, variables to plot (\samp{ZC}, \samp{H2O-ZC}, \samp{GRAVY}, \samp{pI}, or \samp{pIG})}
  \item{parts}{character, parts of figure to include}
  \item{xlim}{numeric, x-axis limits}
  \item{reorder}{logical, put points in order of increasing x value?}
  \item{plot.techtype}{logical, identify 454 or Sanger datasets by open symbols?}
  \item{add}{logical, add to existing plot?}
  \item{datadir}{character, location of data directory (default is from package)}
  \item{mdata}{list, metadata for studies (default is from package)}
  \item{pch}{numeric, plot symbol (default is determined from sample type)}
  \item{lty}{numeric, line type}
  \item{labels.at}{character, plot study labels at \samp{max} or \samp{min} values (use NA to omit labels)}
  \item{cex.ylab}{numeric, size of y-axis label}
  \item{font}{integer, font to use for study labels (see \code{\link{par}})}
  \item{labdx}{numeric, x-adjustment for label position}
  \item{labdy}{numeric, y-adjustment for label position}
}

\details{
The \code{study} is composed of two names for the locality separated by an underscore.
For \code{mplot}, the \code{seqtype} contains the name of the database followed by an underscore and MG (metagenome) or MT (metatranscriptome), i.e. \samp{[SRA|IMG|MGRAST|GenBank]_[MG|MT]}.
These can be changed to \samp{_MGP} and \samp{_MTP} to indicate protein sequences.

\code{mplot} plots a compositional variable for a particular combination of \code{study} and \code{seqtype}.
The variable plotted is \Zc by default, but this can be changed to \nH2O with \code{H2O = TRUE} or GRAVY or pI (specified in the \code{var} argument).
\Zc of both DNA and RNA are plotted, with an offset of -0.28 applied to that of RNA, so that values can be shown on the same scale (see Dick et al., 2019).
The default value of \code{plottype}, \samp{bars}, specifies a plot with mean values connected by dashed lines and error bars for each sample.
A polygon plot (shaded area instead of error bars) is made if the first character of \code{plottype} is \samp{#}.
In this case, \code{plottype} is interpreted as the fill color for the polygon.

\code{mcomp} makes a scatterplot of (\Zc of RNA - \Zc of DNA) vs \Zc of DNA.
Set \code{yvar} to \samp{GC} to change the y-axis to GC ratio of DNA.
(For double-stranded DNA, this should be a linear function of \Zc.)

For \code{pcomp}, \code{seqtype} is just \samp{MG} or \samp{MT}.
The variables (\code{vars}) are automatically determined based on the values of \code{mout} and \code{pout}.
Use \code{vars = "GRAVY"} or \code{vars = "pI"} to plot the values of \code{\link[canprot]{GRAVY}} or \code{\link[canprot]{pI}} instead of \code{\link[canprot]{Zc}}.
Use \code{vars = "pIG"} to make a \code{\link[canprot]{GRAVY}}-\code{\link[canprot]{pI}} plot.

\code{datadir} and \code{mdata} are used to specify alternate location of compositional data and metadata for metagenomic studies.
This permits the processing of user-provided data sets; see \code{\link{subsurface}} for some examples.
}

\references{
Dick JM, Yu M, Tan J and Lu A (2019) Changes in carbon oxidation state of metagenomes along geochemical redox gradients. \emph{Front. Microbiol.} \bold{10}, 120. \doi{10.3389/fmicb.2019.00120}
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
# plot ZC of DNA vs proteins for the metagenomes
pcomp(mout, pout)
# plot ZC of DNA vs proteins for the metatranscriptomes
pcomp(mout, pout, "MT")

## TODO: plot nH2O of DNA vs proteins
#mout <- mpage(H2O = TRUE)
#pout <- ppage(H2O = TRUE)
#pcomp(mout, pout)

# plot nH2O vs ZC of proteins in a few redox gradients
mout <- ppage("gradoxGS")
pout <- ppage("gradoxGS", H2O = TRUE)
pcomp(mout, pout, reorder = FALSE, yline = 2)
# overlay data for protein from the Baltic Sea
mout <- ppage("balticsurface", plot.it = FALSE)
pout <- ppage("balticsurface", H2O = TRUE, plot.it = FALSE)
pcomp(mout, pout, reorder = FALSE, add = TRUE)

# plot nH2O vs GRAVY of proteins
pcomp(mout, pout, vars = "GRAVY")
}
