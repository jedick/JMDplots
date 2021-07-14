\encoding{UTF-8}
\name{MaximAct}
\alias{MaximAct}
\alias{getphyloaa}
\title{Calculate optimal logaH2O and logfO2 for target proteins}
\description{
Calculate the \logaH2O and \logfO2 that maximize the activities of given target proteins on a proteomic background.
}

\usage{
  MaximAct(AA_target, seed = 1:100, nbackground = 2000,
    plot.it = TRUE, xlab = "sample", names = NULL,
    O2 = c(-72, -67), H2O = c(-2, 6), pH = NULL,
    AA_background = NULL)
  getphyloaa(PS_source)
}

\arguments{
  \item{AA_target}{data frame, amino acid composition in the format of \code{\link[CHNOSZ]{thermo}()$protein}}
  \item{seed}{numeric, seeds to use for subsampling background proteins from human proteome}
  \item{nbackground}{numeric, number of background proteins to sample}
  \item{plot.it}{logical, make a plot?}
  \item{xlab}{character or expression, label for \emph{x} axis}
  \item{names}{character, sample names}
  \item{O2}{numeric, range of \logfO2 values}
  \item{H2O}{numeric, range of \logaH2O values}
  \item{pH}{numeric, range of pH values}
  \item{AA_background}{data frame, amino acid compositions of background proteins}
  \item{PS_source}{Source of phylostrata: \samp{TPPG17} or \samp{LMM16}}
}

\details{
\code{MaximAct} is used to compute optimal values of \logaH2O, \logfO2, and optionally pH; that is, those that maximize the equilibrium activity of each of the target proteins against a background of other proteins.
The calculation is described by Dick (2021).

The amino acid compositions of the target and background proteins are given in \code{AA_target} and \code{AA_background}, respectively.
If \code{AA_background} is NULL, the background corresponds to human proteins that have phylostrata assignments in both the Trigos and Liebeskind datasets (see \code{\link[canprot]{PS}}).
The entire human proteome is not used to avoid inclusion of background proteins with unusual amino acid compositions.

The basis species used for the calculation are \samp{QEC} (glutamine - glutamic acid - cysteine - O2 - H2O) (see Dick et al., 2020), or \samp{QEC+} (the same basis species with the addition of \Hplus) if \code{pH} is not NULL.
The NULL default for \code{pH} means that the system consists of uncharged species, and pH is therefore not defined.
The ranges in \code{O2}, \code{H2O}, and \code{pH} can have an optional third value to specify the resolution; e.g. \code{c(-70, -65, 100)} represents 100 points in the range [-70, -65] (see \code{\link[CHNOSZ]{affinity}}).

Due to memory constraints, a single calculation does not include the entire human proteome, but a random subsample of \code{nbackground} human proteins.
The subsampling is performed for each random seed in \code{seed}.
The results are returned in a list with components \code{O2} and \code{H2O}.

By default, a pair of plots is made that are updated during the calculation; at the end, the mean values of \logfO2 and \logaH2O are plotted as red lines.
Change \code{plot.it} to FALSE to prevent making the plots.
The protein names (used to label the x axis of the plots and for the column names in the output) are taken from \code{AA_target$protein} or the \code{names} argument if it is not NULL.

\code{getphyloaa} computes the mean amino acid composition of proteins in each phylostratum using phylostrata from Trigos et al. (2017) (\samp{TPPG17}) or Liebeskind et al. (2016) (\samp{LMM16}).
}

\seealso{
\code{\link{runMaximAct}}, used to set up the calculations reported in Dick (2021).
}

\references{
Dick, J. M., Yu, M. and Tan, J. (2020) Uncovering chemical signatures of salinity gradients through compositional analysis of protein sequences. \emph{Biogeosciences} \bold{17}, 6145--6162. \doi{10.5194/bg-17-6145-2020}

Dick, J. M. (2021) Redox potential linked to water loss from proteins in evolution and developent. bioRxiv. \doi{https://doi.org/10.1101/2021.01.29.428804}

Liebeskind, B. J. and McWhite, C. D. and Marcotte, E. M. (2016) Towards consensus gene ages. \emph{Genome Biol. Evol.} \bold{8}, 1812--1823. \doi{10.1093/gbe/evw113}

Trigos, A. S. and Pearson, R. B. and Papenfuss, A. T. and Goode, D. L. (2017) Altered interactions between unicellular and multicellular genes drive hallmarks of transformation in a diverse range of solid tumors. \emph{Proc. Natl. Acad. Sci.} \bold{114}, 6406--6411. \doi{10.1073/pnas.1617743114}
}

\examples{
# Run MaximAct() for mean amino acid compositions of proteins in phylostrata
AA_target <- getphyloaa("TPPG17")$aa
MA <- MaximAct(AA_target, seed = 1:3, nbackground = 200, O2 = c(-72,-65,15), H2O = c(-2,5,15))
# logfO2 increases and logaH2O decreases between phylostrata 1 and 16
O2 <- colMeans(MA$O2)
stopifnot(O2["16"] - O2["1"] > 0)
H2O <- colMeans(MA$H2O)
stopifnot(H2O["16"] - H2O["1"] < 0)
# The final value of logaH2O is close to 0
stopifnot(round(H2O["16"]) == 0)
}
