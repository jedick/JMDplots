\encoding{UTF-8}
\name{evdevH2O}
\alias{evdevH2O}
\alias{evdevH2O1}
\alias{evdevH2O2}
\alias{evdevH2O3}
\alias{evdevH2O4}
\alias{evdevH2O5}
\alias{runOptimAct}
\alias{logK_example}
\title{Plots from thermodynamics in evo-devo paper}
\description{
Plots from the paper by Dick (2021).
}

\usage{
  evdevH2O1(pdf = FALSE)
  evdevH2O2(pdf = FALSE)
  evdevH2O3(pdf = FALSE)
  evdevH2O4(pdf = FALSE)
  evdevH2O5(pdf = FALSE)
  runOptimAct(dataset = "TPPG17", seed = 1:100)
  logK_example()
}

\arguments{
  \item{pdf}{logical, make a PDF file?}
  \item{dataset}{character, dataset to use for constructing target proteins}
  \item{seed}{numeric, seeds to use for subsampling background proteins from human proteome}
}

\details{
This table gives a brief description of each plotting function.

\tabular{ll}{
  \code{evdevH2O1} \tab Compositional analysis of Trigos and Liebeskind datasets \cr
  \code{evdevH2O2} \tab Thermodynamic analysis of optimal \logaH2O and \logfO2 for target proteins \cr
  \code{evdevH2O3} \tab Optimal \logaH2O and \logfO2 and effective Eh for target proteins \cr
  \code{evdevH2O4} \tab Compositional and thermodynamic analysis of \emph{B. subtilis} biofilm transcriptome and proteome \cr
  \code{evdevH2O5} \tab Organismal water content, proteomic \nH2O, and optimal \logaH2O for fruit fly development \cr
}

\code{runOptimAct} sets up calculations with \code{\link{OptimAct}} to  compute the optimal \logaH2O and \logfO2 that maximize the equilibrium activity of model proteins for phylostrata.
The possible values of \code{dataset} are described below.
\describe{
  \item{\samp{TPPG17}}{Phylostrata from Trigos et al., 2017; target proteins are mean amino acid compositions for each phylostratum.}
  \item{\samp{LMM16}}{Gene ages from Liebeskind et al., 2016; target proteins are mean amino acid compositions for each gene age.}
  \item{\samp{transcriptome} or \samp{proteome}}{Futo et al., 2021; target proteins are mean amino acid composiitons for each stage of biofilm development, weighted by gene or protein expression values.}
  \item{\samp{fly_embryo} or \samp{fly_adult}}{Differentially expressed proteins from Fabre et al., 2019; target proteins are all proteins (i.e. not mean compositions) with higher expression in embryo or adult flies.}
  \item{\samp{fly_development}}{Fly developmental proteome from Casas-Vila et al., 2017; target proteins are mean amino acid composition for each time point.}
}

\code{logK_example} prints the chemical formula, formation reaction, and equilibrium constant for an example protein (LYSC_CHICK).

}

\section{Files in extdata/evdevH2O}{
  These files were created by running \code{runOptimAct} for each \code{dataset} (see above).
  They are used in \code{evdevH2O3} (phylostrata and gene ages) and \code{evdevH2O4} (biofilm developmental stages).
  To reduce package size, the PNG files output by \code{runOptimAct} are not included here.
  \describe{
    \item{\code{[dataset]_H2O.csv}}{Optimal values of \logaH2O for each random seed.}
    \item{\code{[dataset]_O2.csv}}{Optimal values of \logfO2 for each random seed.}
  }

  Two other files, \code{FOK+21_mean_aa.csv} and \code{CBS+17_mean_aa.csv}, contain mean amino acid compositions of proteins computed using datasets for developmental transcriptomes and proteomes in B. subtilis biofilms (Futo et al., 2021) and proteomes in fruit flies (Casas-Vila et al., 2017).
  These files are used in \code{evdevH2O4} and \code{evdevH2O5}.
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
}
