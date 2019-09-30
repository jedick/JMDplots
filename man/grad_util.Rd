\encoding{UTF-8}
\name{grad_util}
\alias{grad_util}
\alias{mplot}
\title{Utility functions for chemical composition in environmental gradients}
\description{
These are utility functions for the redox gradients (\code{\link{gradox}}) and salinity gradients (\code{\link{gradH2O}}) plots.
The functions are exported in the package so they are available to user scripts.
}

\usage{
  mplot(study, seqtype, plottype="bars", ylim=NULL, plot.RNA=TRUE,
    taxid=NULL, dsDNA=TRUE, abbrev=NULL, col=NULL, add.label=TRUE, maxdepth=NULL)
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
}

\examples{
# plot ZC of DNA
mplot("Bison_Pool", "IMG_MG")
# plot ZC of protein
mplot("Bison_Pool", "IMG_MGP")
}
