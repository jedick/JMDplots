\encoding{UTF-8}
\name{cpcp}
\alias{cpcp}
\alias{Ehplot}
\alias{rankdiff}
\alias{rankplot}
\alias{groupplots}
\alias{mergedplot}
\alias{get_colors}
\title{Chemical portraits of cancer proteomes}
\description{
Plots from the papers by Dick (2016 and 2017).
}

\usage{
  Ehplot(T = 37, pH = 7.4,
    xlim = c(-75, -55), ylim = c(-10, 10), dy = 0.2)
  rankdiff(rank1, rank2, n1 = NULL, n2 = NULL, as.fraction=TRUE)
  rankplot(pdat, T = 37, what = "rankdiff", main = NULL, res = 300,
    plot.it = TRUE, xlim = c(-75, -55), ylim = c(-10, 10), basis = "QEC+")
  groupplots(group = "hypoxia_ZC_down", each100 = FALSE, res = 50,
    plot.it = TRUE)
  mergedplot(gpresult, each100 = FALSE, res = 50)
  get_colors(x, max50 = FALSE)
}

\arguments{
  \item{T}{numeric, temperature in degrees Celsius}
  \item{pH}{numeric, pH}
  \item{xlim}{numeric, range of x axis}
  \item{ylim}{numeric, range of y axis}
  \item{basis}{character, basis species to use (see \code{\link[CHNOSZ]{basis}})}
  \item{dy}{numeric, y spacing between iso-\logfO2 lines}
  \item{rank1}{numeric, ranks in group 1}
  \item{rank2}{numeric, ranks in group 2}
  \item{n1}{numeric, size of group 1}
  \item{n2}{numeric, size of group 2}
  \item{as.fraction}{logical, calculate the fraction of maximum possible difference?}
  \item{pdat}{list, data object generated by a \code{\link{pdat_}} function}
  \item{what}{character, \dQuote{rankdiff} or \dQuote{affinity}}
  \item{main}{character, text to use for title of plot}
  \item{res}{numeric, grid resolution for plot}
  \item{plot.it}{logical, draw a plot?}
  \item{group}{character, description of datasets to include}
  \item{each100}{logical, rescale rank difference of each dataset individually?}
  \item{gpresult}{list, value returned by \code{groupplots}}
  \item{x}{numeric values}
  \item{max50}{logical, use most intense color for all values >= 50?}
}

\details{
  \describe{
    \item{\code{Ehplot}}{
      This function plots selected values of Eh (redox potential) as a function of \logfO2 and \logaH2O.
      The lines are labeled with the Eh value in volts.
      The temperature and pH can be adjusted using the \code{T} and \code{pH} arguments; this affects the lines, but not the positions of labels.
    }
    \item{\code{rankdiff}}{
      Calculate rank-sum difference between two groups, weighted by the sizes of the groups.

      In a combined ranking of two groups, the comparison of sum of ranks has an easy interpretation only for groups of equal size.
      The weighted rank difference is used to compare groups of unequal size.
      The weighting ensures that 1) opposite extreme configurations give weighted rank differences with equal magnitudes, and 2) an evenly distributed (interspersed) ranking of the two groups has a weighted rank difference of zero (Dick, 2016).

      If \code{n1} and \code{n2} are not given, \code{rank1} and \code{rank2} are interpreted as vectors holding the ranks for the two groups.
      If the sizes of the groups are supplied in \code{n1} and \code{n2}, then the single values or higher-dimensional objects in \code{rank1} and \code{rank2} are interpreted as the non-weighted sums of ranks of the two groups.
    }
    \item{\code{rankplot}}{
      Plot ranking of chemical affinities of groups of proteins.

      This function creates a \logaH2O - \logfO2 diagram showing the relative stabilities of the two groups of proteins in the specified \code{dataset}.
      These groups consist of the relatively down- and up-expressed proteins identified by \code{up2} in one of the \code{\link{pdat_}} functions.

      The function generates a colored \code{\link{image}} and \code{\link{contour}} plot showing the weighted difference of sums of ranks of chemical affinities of formation of proteins in the two groups (see \code{\link[CHNOSZ]{affinity}} and \code{\link{rankdiff}}).
      Increasing intensity of blue or red colors (see \code{\link{get_colors}}) represent higher rankings of the down-expressed (\code{up2==FALSE}) or up-expressed (\code{up2==TRUE}) proteins, respectively.
      Alternatively, if \code{what} is \dQuote{affinity}, a maximum affinity diagram is produced (see \code{\link[CHNOSZ]{diagram}} in CHNOSZ), with fields colored red or blue according to the relative expression of the protein.

      If \code{main} is NULL, the title for the plot is taken from the description supplied in \code{pdat}.

      Set \code{plot.it} to FALSE to skip the plotting and instead return a list containing the computed rank differences and x- and y- values and labels.
    }
    \item{\code{groupplots}, \code{mergedplot}}{
      Plot rank difference of chemical affinities for proteins in various datasets and merge the diagrams.

      \code{groupplots} makes weighted rank-difference of affinity (potential) diagrams (see \code{\link{rankplot}}) for each dataset found in the specified \code{group}.
      \code{group} consists of three parts joined by an underscore: the type of experiment (\samp{colorectal}, \samp{pancreatic}, \samp{hypoxia}, or \samp{osmotic}; see \code{\link{pdat_}}), the distinguishing compositional variable (\samp{ZC} or \samp{H2O}), and the direction of change of that variable (\samp{up} or \samp{down}).

      To identify the datasets in any group, compositional summaries for each dataset are read from pre-calculated tables in \code{extdata/summary}.
      Datasets are included for which the absolute mean difference of either \samp{ZC} or \samp{H2O} between up- and down-expressed proteins is greater than 0.01 and the other of \samp{ZC} or \samp{H2O} has \emph{p}-value >= 0.05 and abs(\code{\link[canprot]{CLES}} - 50) < 10.

      \code{groupplots} makes calculations over a large range of \logfO2 and \logaH2O in order to encompass the equipotential lines for most datasets.
      This way, the positions of the median and interquartiles of the equipotential lines can be calculated accurately for the \code{mergedplot}, which covers a smaller range of \logfO2 and \logaH2O.
    }
    \item{\code{get_colors}}{
      \code{get_colors} returns a diverging (blue - light gray - red) color scale.
      Blue and red colors are associated with negative and positive values, respectively.
      The intensity of the color increases with the magnitude of the value.
      For accurate representation, the values should be in a percent scale (i.e. the maximum absolute value is not greater than 100).
      By default, a value of +/- 100 corresponds to greatest intensity.
      Set \code{max50} to TRUE to compress the scale so that greatest intesity is obtained at values of +/- 50 and higher.

      These colors are used in \code{rankplot}.
      The colors were precomputed using \pkg{colorspace}::\code{diverge_hcl}.
    }
  }
}


\section{Files in extdata/Dic16}{
  \describe{
    \item{\code{microbes.csv}, \code{microbial.aa.csv}}{
      Data for microbial occurrence in colorectal cancer and normal tissue, and overall amino acid compositions of proteins in microbial species.
      The files are from the Supporting Information of Dick (2016).
    }
  }
}

\section{Files in extdata/cpcp}{
  \describe{
    \item{\code{colorectal.Rmd}, \code{hyperosmotic.Rmd}, \code{hypoxia.rmd}, \code{pancreatic.Rmd}}{
      Vignettes showing compositional analysis of differential protein expression, based on Figures 1-2 and Tables 1-4 of Dick (2017).
      Updated versions of these plots and the underlying datasets are in Dick (2020).
    }
    \item{\code{summary}}{
      Directory with CSV files summarizing compositional differences, used by \code{groupplots}.
    }
  }
}

\seealso{
Other conversions involving Eh are possible with \code{\link[CHNOSZ]{convert}} (from \pkg{CHNOSZ}).
}

\examples{
Ehplot()

# rankings of H and C in H-H-H-H-C-C-C
rankdiff(1:4, 5:7, as.fraction=FALSE)  # 12
rankdiff(1:4, 5:7)  # 1
# rankings of H and C in C-C-C-H-H-H-H
rankdiff(4:7, 1:3, as.fraction=FALSE)  # -12
rankdiff(4:7, 1:3)  # -1
# rankings of H and C in H-C-H-C-H-C-H
rankdiff(c(1, 3, 5, 7), c(2, 4, 6))  # 0

pdat <- pdat_colorectal("JKMF10")
rankplot(pdat, res=25)
rankplot(pdat, res=25, what="affinity")

\dontrun{
gpresult <- groupplots("osmotic_H2O_down", res = 25)
mergedplot(gpresult, res = 25)

# make plots similar to Figure 3 of Dick, 2017
ZCgroups <- c("colorectal_ZC_up", "pancreatic_ZC_up", "hypoxia_ZC_down")
H2Ogroups <- c("colorectal_H2O_up", "pancreatic_H2O_up", "osmotic_H2O_down")
allgroups <- c(ZCgroups, H2Ogroups)
par(mfrow=c(2, 3))
for(group in allgroups) {
  gpresult <- groupplots(group, plot.it = FALSE)
  mergedplot(gpresult)
  title(main = group)
}}
}

\references{
Dick JM (2016) Proteomic indicators of oxidation and hydration state in colorectal cancer. \emph{PeerJ} \bold{4}, e2238. \doi{10.7717/peerj.2238}

Dick JM (2017) Chemical composition and the potential for proteomic transformation in cancer, hypoxia, and hyperosmotic stress. \emph{PeerJ} \bold{5}, e3421 \doi{10.7717/peerj.3421}

Dick JM (2021) Water as a reactant in the differential expression of proteins in cancer. \emph{Comp. Sys. Onco.} \bold{1}:e1007. \doi{10.1002/cso2.1007}
}
