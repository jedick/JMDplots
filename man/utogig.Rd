\encoding{UTF-8}
\name{utogig}
\alias{utogig}
\alias{utogig1}
\alias{utogig2}
\alias{utogig3}
\alias{utogig4}
\alias{utogigS1}
\alias{utogigS2}
\alias{utogigS4}
\alias{calc_logaH2_intermediate}
\title{Using thermodynamics to obtain geochemical information from genomes}
\description{
Plots from the paper by Dick et al. (2023).
}

\usage{
  utogig1(pdf = FALSE)
  utogig2(pdf = FALSE, logact = -3)
  utogig3(pdf = FALSE)
  utogig4(pdf = FALSE)
  utogigS1(pdf = FALSE)
  utogigS2(pdf = FALSE)
  utogigS4(pdf = FALSE)
  calc_logaH2_intermediate(class = NULL, logact = -3)
}

\arguments{
  \item{pdf}{logical, make a PDF file?}
  \item{logact}{numeric, logarithm of activity of organic species}
  \item{class}{character, organic compound class}
}

\details{
This table gives a brief description of each plotting function.

\tabular{ll}{
  \code{utogig1} \tab Chemical analysis of reference proteomes of methanogens reveals adaptation to redox conditions \cr
  \code{utogig2} \tab Relative stabilities of organic compounds depend on redox conditions \cr
  \code{utogig3} \tab Thermodynamic model for methanogen niche partitioning \cr
  \code{utogig4} \tab Chemical and thermodynamic analysis of evolutionary divergence along redox gradients \cr
  \code{utogigS1} \tab Comparison of \Zc of proteomes predicted by Glimmer and downloaded from NCBI \cr
  \code{utogigS2} \tab Association between redox gradients and \Zc of proteins and lipids in alkaline Yellowstone hot springs \cr
  \code{utogigS4} \tab \logaH2-\T plots for different organic compound classes \cr
}

\code{utogig3(logact = -6)} is used to make Figure S3.
\code{utogig3() and utogig3(logact -6)} print p-values that are used to make Table S5.
\code{utogig4()} returns an invisible data frame of p-values that is used to make Table S6 and prints values of \logaH2 retrieved from the thermodynamic analysis in Figure 4b.

\code{calc_logaH2_intermediate} is used to find intermediate \logaH2, where affinity vs Zc has a slope of 0, for a group of organic compounds.
\code{class} gives the class of compounds to consider , or NULL to consider compounds in all classes at once.
The output is saved in files named \file{H2_intermediate_*.csv} (see description below), which are used for making some of the plots.

}

\section{Files in extdata/utogig}{
  \describe{
    \item{\file{Topt.csv}}{Optimal growth temperatures (\degC) of methanogens compiled from the literature.}
    \item{\file{methanogen_AA.csv}}{Amino acid compositions of reference proteomes of methanogens. See \file{R/utogig.R} for the code used to make this file.}
    \item{\file{LG88_Fig1.csv}}{H2 concentrations measured in methanogenic and non-methanogenic sediments, from Lovley and Goodwin (1988).}
    \item{\file{H2_intermediate_-3.csv}, \file{H2_intermediate_-6.csv}}{Intermediate \logaH2 (where the slope of affinity vs \Zc is zero) for organic compounds, calculated using \code{calc_logaH2_intermediate} with \code{logact} set to -3 or -6.}
    \item{\file{H2_intermediate_*_-3.csv}}{Intermediate \logaH2 for different classes of organic compounds (\samp{C1andC2}, \samp{Acid}, \samp{Aminoacid}, \samp{Sugar}, \samp{Nucleobase}, \samp{TCA}), all with \code{logact} set to -3.}
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
Delcher AL, Harmon D, Kasif S, White O and Salzberg SL (1999) Improved microbial gene identification with Glimmer. \emph{Nucleic Acids Res.} \bold{27}, 4636--4641. \doi{10.1093/nar/27.23.4636}

Dick JM, Boyer GM, Canovas PA III and Shock EL (2023) Using thermodynamics to obtain geochemical information from genomes. \emph{Geobiology} \bold{21}, 262--273. \doi{10.1111/gbi.12532}

Lovley DR and Goodwin S (1988) Hydrogen concentrations as an indicator of the predominant terminal electron-accepting reactions in aquatic sediments. \emph{Geochim. Cosmochim. Acta} \bold{52}, 2993--3003. \doi{10.1016/0016-7037(88)90163-9}

Lyu Z and Lu Y (2018) Metabolic shift at the class level sheds light on adaptation of methanogens to oxidative environments. \emph{ISME J.} \bold{12}, 411--423. \doi{10.1038/ismej.2017.173}

Ren M, Feng X et al. (2019) Phylogenomics suggests oxygen availability as a driving force in Thaumarchaeota evolution. \emph{ISME J.} \bold{13}, 2150--2161. \doi{10.1038/s41396-019-0418-8}
}

\examples{
utogigS1()
}
