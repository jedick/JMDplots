\encoding{UTF-8}
\name{evdevH2O}
\alias{evdevH2O}
\alias{evdevH2O1}
\alias{evdevH2O2}
\alias{evdevH2O3}
\alias{evdevH2O4}
\alias{evdevH2O5}
\alias{evdevH2O6}
\alias{evdevH2O7}
\alias{evdevH2O8}
\alias{runMaximAct}
\alias{LYSC_example}
\title{Plots from thermodynamics in evolution and development paper}
\description{
Plots from the paper by Dick (2022).
}

\usage{
  evdevH2O1(pdf = FALSE)
  evdevH2O2(pdf = FALSE, boot.R = 99)
  evdevH2O3(pdf = FALSE, H2O = FALSE)
  evdevH2O4(pdf = FALSE)
  evdevH2O5(pdf = FALSE)
  evdevH2O6(pdf = FALSE)
  evdevH2O7(pdf = FALSE, boot.R = 99)
  evdevH2O8(pdf = FALSE, boot.R = 99)
  runMaximAct(dataset = "TPPG17", seed = 1:100, nbackground = 2000,
    res = 256, bg_organism = "Hsa", writeFiles = TRUE)
  LYSC_example()
}

\arguments{
  \item{pdf}{logical, make a PDF file?}
  \item{boot.R}{integer, number of bootstrap replicates}
  \item{H2O}{logical, make plots for \nH2O instead of \ZC?}
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
  \code{evdevH2O1} \tab Comparison of different sets of basis species \cr
  \code{evdevH2O2} \tab Protein length and chemical metrics for phylostratigraphic age groups \cr
  \code{evdevH2O3} \tab Evolution of protein \ZC in eukaryotic lineages \cr
  \code{evdevH2O4} \tab Thermodynamic analysis of optimal \logaH2O and \logfO2 for target proteins \cr
  \code{evdevH2O5} \tab Optimal \logaH2O and \logfO2 and virtual Eh for target proteins \cr
  \code{evdevH2O6} \tab Chemical metrics for and thermodynamic parameters with different background proteomes \cr
  \code{evdevH2O7} \tab Chemical and thermodynamic analysis of \emph{B. subtilis} biofilm transcriptome and proteome \cr
  \code{evdevH2O8} \tab Organismal water content, proteomic \nH2O, and optimal \logaH2O for fruit fly development \cr
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

\section{Files in extdata/evdevH2O/MaximAct}{
  \describe{
    \item{\code{runall.sh}}{Script to run \code{runMaximAct} for different datasets and background proteomes.}
    \item{\code{[dataset]_H2O.csv}}{Output file with optimal values of \logaH2O for each random seed.}
    \item{\code{[dataset]_O2.csv}}{Output file with optimal values of \logfO2 for each random seed.}
  }
  The output files are used in \code{evdevH2O5} (thermodynamic analysis of phylostrata and gene ages) and \code{evdevH2O7} (thermodynamic analysis of biofilm developmental stages).
  To reduce package size, the PNG files output by \code{runMaximAct} are not included here.
}

\section{Files in extdata/evdevH2O/LMM16}{
  \describe{
    \item{\code{mkmetrics.R}}{
      Script to calculate \ZC and \nH2O and get gene ages for proteins in the consensus tables of Liebeskind et al. (2016).
      This script requires files available from UniProt and Zenodo (\doi{10.5281/zenodo.51708}); see comments for details.
    }
    \item{\code{reference_proteomes.csv}}{IDs of UniProt reference proteomes for 31 organisms in Liebeskind et al. (2016) (input file for \code{mkmetrics.R}).}
    \item{\code{modeAges.csv}}{Names of gene age categories for each organism (output by \code{mkmetrics.R}).}
    \item{\code{metrics/*.csv.xz}}{Files with ID, modeAge, \ZC, \nH2O for each protein (output by \code{mkmetrics.R}).}
  }
}

\section{Files in extdata/evdevH2O/JWN+21}{
  \describe{
    \item{\code{mkAA.R}}{
      Script to assemble amino acid compositions from data files of James et al. (2021).
      This script requires files available from Figshare (\doi{10.6084/m9.figshare.12037281.v1}); see comments for details.
    }
    \item{\code{pfam_animal_nontrans_AA.csv.xz}, \code{pfam_animal_trans_AA.csv.xz}, \code{pfam_plant_nontrans_AA.csv.xz}, \code{pfam_plant_trans_AA.csv.xz}}{Output files from the script.}
  }
}

\references{
Casas-Vilas, N. et al. (2017) The developmental proteome of \emph{Drosophila melanogaster}. \emph{Genome Res.} \bold{27}, 1273--1285. \doi{10.1101/gr.213694.116}

Dick, J. M. (2022) A thermodynamic model for water activity and redox potential in evolution and developent. \emph{J. Mol. Evol} \bold{90}, 182--199. \doi{10.1007/s00239-022-10051-7}

Fabre, B., Korona, D., Lees, J. G., Lazar, I., Livneh, I., Brunet, M., Orengo, C. A., Russell, S. and Lilley, K. S. (2019) Comparison of \emph{Drosophila melanogaster} embryo and adult proteome by SWATH-MS reveals differential regulation of protein synthesis, degradation machinery, and metabolism modules. \emph{J. Proteome Res.} \bold{18}, 2525--2534. \doi{10.1021/acs.jproteome.9b00076}

Futo, M. et al. (2021) Embryo-like features in developing \emph{Bacillus subtilis} biofilms. \emph{Mol. Biol. Evol.} \bold{38}, 31--47. \doi{10.1093/molbev/msaa217}

James, J. E., Willis, S. M., Nelson, P. G., Weibel, C., Kosinski, L. J. and Masel, J. (2021) Universal and taxon-specific trends in protein sequences as a function of age. \emph{eLife} \bold{10}, e57347. \doi{10.7554/eLife.57347}

Liebeskind, B. J. and McWhite, C. D. and Marcotte, E. M. (2016) Towards consensus gene ages. \emph{Genome Biol. Evol.} \bold{8}, 1812--1823. \doi{10.1093/gbe/evw113}

Trigos, A. S. and Pearson, R. B. and Papenfuss, A. T. and Goode, D. L. (2017) Altered interactions between unicellular and multicellular genes drive hallmarks of transformation in a diverse range of solid tumors. \emph{Proc. Natl. Acad. Sci.} \bold{114}, 6406--6411. \doi{10.1073/pnas.1617743114}
}

\seealso{
\code{\link{devodata}} for data sources for developmental proteomes.
\code{\link{PS}} for functions used to retrieve phylostrata from UniProt IDs.
}

\examples{
evdevH2O2()

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
