\encoding{UTF-8}
\name{grad_mplot}
\alias{grad_mplot}
\alias{mplot}
\alias{mpage}
\alias{ppage}
\title{Plotting functions for chemical composition of metagenomes}
\description{
These are plotting functions for the redox gradients (\code{\link{gradox}}) and salinity gradients (\code{\link{gradH2O}}) papers.
The functions are exported in the package so they are available to user scripts.
}

\usage{
  mplot(study, seqtype, plottype = "bars", ylim = NULL, plot.RNA = TRUE,
    taxid = NULL, dsDNA = TRUE, abbrev = NULL, col = NULL, add.label = TRUE,
    maxdepth = NULL, H2O = FALSE, plot.it = TRUE, add.title = TRUE, yline = 2,
    basis = "QEC")
  mpage(subset = "gradoxSI", H2O = FALSE, plottype = "bars", dsDNA = TRUE, set.par = TRUE)
  ppage(subset = "gradoxSI", H2O = FALSE, set.par = TRUE, plot.it = TRUE, basis = "QEC")
}

\arguments{
  \item{study}{character, name of study}
  \item{seqtype}{character, description of sequence data type}
  \item{plottype}{character, type of plot}
  \item{ylim}{numeric, limits of y-axis}
  \item{plot.RNA}{logical, plot RNA compositions?}
  \item{taxid}{numeric, select only this taxid}
  \item{dsDNA}{logical, use chemical composition of double-stranded DNA?}
  \item{abbrev}{character, abbreviations for studies}
  \item{col}{character, line color}
  \item{add.label}{logical, label the plot?}
  \item{maxdepth}{numeric, maximum sample depth (meters)}
  \item{H2O}{logical, plot \nH2O instead of \ZC?}
  \item{plot.it}{logical, make a plot (use FALSE to return processed data without plotting)?}
  \item{add.title}{logical, add a title to the plot?}
  \item{yline}{numeric, margin position for y-axis label}
  \item{basis}{character, description of basis species to use}
  \item{subset}{character, name of subset of studies to plot}
  \item{set.par}{logical, set up plot with \code{\link{par}} commands?}
}

\examples{
# plot ZC of DNA
mplot("Bison_Pool", "IMG_MG")
# plot ZC of protein
mplot("Bison_Pool", "IMG_MGP")
}
