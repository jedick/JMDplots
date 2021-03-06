\encoding{UTF-8}
\name{evdevH2O}
\alias{evdevH2O}
\alias{evdevH2O1}
\alias{evdevH2O2}
\alias{evdevH2O3}
\alias{evdevH2O4}
\alias{evdevH2O5}
\alias{evdevH2O6}
\alias{runMaximAct}
\alias{LYSC_example}
\title{Plots from thermodynamics in evolution and development paper}
\description{
Plots from the paper by Dick (2021).
}

\usage{
  evdevH2O1(pdf = FALSE, boot.R = 99)
  evdevH2O2(pdf = FALSE)
  evdevH2O3(pdf = FALSE)
  evdevH2O4(pdf = FALSE)
  evdevH2O5(pdf = FALSE, boot.R = 99)
  evdevH2O6(pdf = FALSE, boot.R = 99)
  runMaximAct(dataset = "TPPG17", seed = 1:100, nbackground = 2000,
    res = 256, bg_organism = "Hsa", writeFiles = TRUE)
  LYSC_example()
}

\arguments{
  \item{pdf}{logical, make a PDF file?}
  \item{boot.R}{integer, number of bootstrap replicates}
  \item{dataset}{character, dataset to use for constructing target proteins}
  \item{seed}{numeric, seeds to use for subsampling background proteins from human proteome}
  \item{nbackground}{numeric, number of background proteins to sample}
  \item{res}{numeric, grid resolution for \logfO2 and \logaH2O}
  \item{writeFiles}{logical, creates CSV and PNG output files?}
  \item{bg_organism}{character, use background proteins from this organism.}
}

\details{
This table gives a brief description of each plotting function.

\tabular{ll}{
  \code{evdevH2O1} \tab Chemical analysis of Trigos and Liebeskind datasets \cr
  \code{evdevH2O2} \tab Thermodynamic analysis of optimal \logaH2O and \logfO2 for target proteins \cr
  \code{evdevH2O3} \tab Optimal \logaH2O and \logfO2 and virtual Eh for target proteins \cr
  \code{evdevH2O4} \tab Ranges of chemical metrics and thermodynamic parameters for different background proteomes \cr
  \code{evdevH2O5} \tab Chemical and thermodynamic analysis of \emph{B. subtilis} biofilm transcriptome and proteome \cr
  \code{evdevH2O6} \tab Organismal water content, proteomic \nH2O, and optimal \logaH2O for fruit fly development \cr
}

\code{runMaximAct} sets up calculations with \code{\link{MaximAct}} to compute the optimal \logaH2O and \logfO2 that maximize the equilibrium activity of target proteins against a proteomic background.
The available \code{dataset}s are described below.
\describe{
  \item{\samp{TPPG17}}{Phylostrata from Trigos et al., 2017; target proteins are mean amino acid compositions for each phylostratum.}
  \item{\samp{LMM16}}{Gene ages from Liebeskind et al., 2016; target proteins are mean amino acid compositions for each gene age.}
  \item{\samp{transcriptome} or \samp{proteome}}{Futo et al., 2021; target proteins are mean amino acid composiitons for each stage of biofilm development, weighted by gene or protein expression values.}
  \item{\samp{fly_embryo} or \samp{fly_adult}}{Differentially expressed proteins from Fabre et al., 2019; target proteins are all proteins (i.e. not mean compositions) with higher expression in embryo or adult flies.}
  \item{\samp{fly_development}}{Fly developmental proteome from Casas-Vila et al., 2017; target proteins are mean amino acid composition for each time point.}
}

\code{bg_organism} can be one of \samp{Hsa} (\emph{Homo sapiens}), \samp{Sce} (\emph{Saccharomyces cerevisiae}), \samp{Eco} (\emph{Escherichia coli}), \samp{Mja} (\emph{Methanocaldococcus jannaschii}), \samp{Dme} (\emph{Drosophila melanogaster}), or \samp{Bsu} (\emph{Bacillus subtilis}).
These are used to retrieve the background proteins used for the \code{AA_background} argument in \code{\link{MaximAct}}.

If \code{writeFiles} is TRUE, a PNG plot file is created and optimal \logfO2 or \logaH2O values for each seed and protein are saved in CSV files named by pasting together the \code{dataset}, \code{bg_organism}, and \samp{O2.csv} or \samp{H2O.csv}, separated by \samp{_}.

\code{LYSC_example} prints the chemical formula, formation reaction, equilibrium constant, and activity product for the example protein (LYSC_CHICK) used in the text.

}

\section{Files in extdata/evdevH2O}{
  \describe{
    \item{\code{runall.sh}}{Script to run \code{runMaximAct} for different datasets and background proteomes.}
    \item{\code{[dataset]_H2O.csv}}{Output file with optimal values of \logaH2O for each random seed.}
    \item{\code{[dataset]_O2.csv}}{Output file with optimal values of \logfO2 for each random seed.}
  }
  The output files are used in \code{evdevH2O3} (thermodynamic analysis of phylostrata and gene ages) and \code{evdevH2O5} (thermodynamic analysis of biofilm developmental stages).
  To reduce package size, the PNG files output by \code{runMaximAct} are not included here.

  Two other files, \code{FOK+21_mean_aa.csv} and \code{CBS+17_mean_aa.csv}, contain mean amino acid compositions of proteins computed using datasets for developmental transcriptomes and proteomes in B. subtilis biofilms (Futo et al., 2021) and proteomes in fruit flies (Casas-Vila et al., 2017).
  These files are used in \code{evdevH2O5} and \code{evdevH2O6}.
}

\references{
Casas-Vilas, N. et al. (2017) The developmental proteome of \emph{Drosophila melanogaster}. \emph{Genome Res.} \bold{27}, 1273--1285. \doi{10.1101/gr.213694.116}

Dick, J. M. (2021) Redox potential linked to water loss from proteins in evolution and developent. bioRxiv. \doi{https://doi.org/10.1101/2021.01.29.428804}

Fabre, B., Korona, D., Lees, J. G., Lazar, I., Livneh, I., Brunet, M., Orengo, C. A., Russell, S. and Lilley, K. S. (2019) Comparison of \emph{Drosophila melanogaster} embryo and adult proteome by SWATH-MS reveals differential regulation of protein synthesis, degradation machinery, and metabolism modules. \emph{J. Proteome Res.} \bold{18}, 2525--2534. \doi{10.1021/acs.jproteome.9b00076}

Futo, M. et al. (2021) Embryo-like features in developing \emph{Bacillus subtilis} biofilms. \emph{Mol. Biol. Evol.} \bold{38}, 31--47. \doi{10.1093/molbev/msaa217}

Liebeskind, B. J. and McWhite, C. D. and Marcotte, E. M. (2016) Towards consensus gene ages. \emph{Genome Biol. Evol.} \bold{8}, 1812--1823. \doi{10.1093/gbe/evw113}

Trigos, A. S. and Pearson, R. B. and Papenfuss, A. T. and Goode, D. L. (2017) Altered interactions between unicellular and multicellular genes drive hallmarks of transformation in a diverse range of solid tumors. \emph{Proc. Natl. Acad. Sci.} \bold{114}, 6406--6411. \doi{10.1073/pnas.1617743114}
}

\examples{
evdevH2O1()

# Start with a fresh plot device to run the next example
# (Something from above messes up split.screen())
dev.off()

MA <- runMaximAct("transcriptome", seed = 1:3, nbackground = 200, res = 15, writeFiles = FALSE)
# The logaH2O decreases at the end of the time course
H2O <- colMeans(MA$H2O)
stopifnot(H2O["14D"] > H2O["1M"])
stopifnot(H2O["1M"] > H2O["2M"])
# The final value of logaH2O is close to 0
stopifnot(round(H2O["2M"]) == 0)
}
