\encoding{UTF-8}
\name{utegig}
\alias{utegig}
\alias{utegig1}
\alias{utegig2}
\alias{utegig3}
\alias{utegig4}
\alias{utegigS1}
\alias{utegigS2}
\alias{utegigS3}
\alias{calc_logaH2_intermediate}
\title{Plots from geochemistry of genomes paper}
\description{
Plots from the paper by Dick et al. (2022).
}

\usage{
  utegig1(pdf = FALSE)
  utegig2(pdf = FALSE)
  utegig3(pdf = FALSE)
  utegig4(pdf = FALSE)
  utegigS1(pdf = FALSE)
  utegigS2(pdf = FALSE)
  utegigS3(pdf = FALSE)
  calc_logaH2_intermediate(class = NULL)
}

\arguments{
  \item{pdf}{logical, make a PDF file?}
  \item{class}{character, organic compound class}
}

\details{
This table gives a brief description of each plotting function.

\tabular{ll}{
  \code{utegig1} \tab Chemical adaptation to redox conditions in the reference proteomes of methanogens \cr
  \code{utegig2} \tab Relative stabilities of organic compounds depend on redox conditions \cr
  \code{utegig3} \tab Thermodynamic model for methanogen niche partitioning \cr
  \code{utegig4} \tab Chemical and thermodynamic analysis of evolutionary divergence along redox gradients \cr
  \code{utegigS1} \tab \ZC of proteins and lipids in hot springs \cr
  \code{utegigS2} \tab Intermediate \logaH2 calculated for different compound classes \cr
  \code{utegigS3} \tab Comparison of \ZC of proteomes predicted by Glimmer and downloaded from NCBI \cr
}

\code{calc_logaH2_intermediate} is used to find intermediate logaH2, where affinity vs ZC has a slope of 0, for a group of organic compounds.
\code{class} gives the class of compounds to consider (see description of \file{H2_intermediate_*.csv} files below), or NULL to consider compounds in all classes at once.

}

\section{Files in extdata/utegig}{
  \describe{
    \item{\file{Topt.csv}}{Optimal growth temperatures (\degC) of methanogens compiled from the literature.}
    \item{\file{methanogen_AA.csv}}{Amino acid compositions of reference proteomes of methanogens. See \file{R/utegig.R} for the code used to make this file.}
    \item{\file{LG88_Fig1.csv}}{H2 concentrations measured in methanogenic and non-methanogenic sediments, from Lovley and Goodwin (1988).}
    \item{\file{H2_intermediate.csv}}{Intermediate logaH2 (where the slope of affinity vs \ZC is zero) for organic compounds, calculated using \code{calc_logaH2_intermediate}.}
    \item{\file{H2_intermediate_*.csv}}{Intermediate logaH2 (where the slope of affinity vs \ZC is zero) for different classes of organic compounds (\samp{C1andC2}, \samp{Acid}, \samp{Aminoacid}, \samp{Sugar}, \samp{Nucleobase}, \samp{TCA}).}
    \item{\file{methanogen_tree.txt}}{Hierarchical clustering tree for functional genes of Class I and Class I methanogens. Tree was digitized and modified from Fig. 1 of Lyu and Lu (2018) and is provided in Newick format.}
    \item{\file{Thaumarchaeota_predicted_AA.csv}}{Amino acid compositions for predicted proteomes of basal, terrestrial, and shallow-water groups of Thaumarchaeota (habitats are listed in the \samp{protein} column). Proteins were predicted using Glimmer version 3.02b (Delcher et al., 1999) from the genomes downloaded from NCBI (accession numbers in the \samp{organisms} column).}
    \item{\file{Thaumarchaeota_database_AA.csv}}{Amino acid compositions for proteomes of basal, terrestrial, shallow- and deep-water groups of Thaumarchaeota obtained from IMG (numeric TaxonID in the \samp{organisms} column) or NCBI (alphanumeric accession number in the \samp{organisms} column). Taxon IDs and accession numbers were taken from Table S1 of Ren et al. (2019)}.
    \item{\code{Nif_homolog_genomes.R}, \code{Nif_homolog_genomes.csv}, \code{Nif_homolog_AA.csv}}{
      The R script reads Supplemental Table 1A of Poudel et al. (2018) and adds the columns \samp{Refseq.name} and \samp{taxid}, which contain the closest matching genome name in the RefSeq database (see \code{\link[chem16S]{chem16S-package}}) and corresponding taxonomic ID.
      Output file #1: Names of genomes (organisms) containing different nitrogenase (Nif) homologs.
      Output file #2: Amino acid compositions of RefSeq proteomes for each matched genome.
    }
  }
}

\references{
Delcher, A. L., Harmon, D., Kasif, S., White, O. and Salzberg, S. L. (1999) Improved microbial gene identification with Glimmer. \emph{Nucleic Acids Res.} \bold{27}, 4636--4641. \doi{10.1093/nar/27.23.4636}

Dick, J. M., Boyer, G. M., Canovas, P. A. III and Shock, E. L. (2022) Using thermodynamics to extract geochemical information from genomes. Manuscript in preparation.

Lovley, D. R. and Goodwin, S. (1988) Hydrogen concentrations as an indicator of the predominant terminal electron-accepting reactions in aquatic sediments. \emph{Geochim. Cosmochim. Acta} \bold{52}, 2993--3003. \doi{10.1016/0016-7037(88)90163-9}

Lyu, Z. and Lu, Y. (2018) Metabolic shift at the class level sheds light on adaptation of methanogens to oxidative environments. \emph{ISME J.} \bold{12}, 411--423. \doi{10.1038/ismej.2017.173}

Ren, M., Feng, X., ..., and Luo, H. (2019) Phylogenomics suggests oxygen availability as a driving force in Thaumarchaeota evolution. \emph{ISME J.} \bold{13}, 2150--2161. \doi{10.1038/s41396-019-0418-8}
}

\examples{
utegigS1()
}
