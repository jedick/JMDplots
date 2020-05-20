\encoding{UTF-8}
\name{gradH2O}
\alias{gradH2O}
\alias{gradH2O0}
\alias{gradH2O1}
\alias{gradH2O2}
\alias{gradH2O3}
\alias{gradH2O4}
\alias{gradH2O5}
\alias{gradH2O6}
\alias{gradH2O7}
\alias{gradH2O8}
\alias{NifProteomes}
\title{Plots from salinity gradients paper}
\description{
Plots from the paper by Dick et al. (2020).
}

\usage{
  gradH2O0()
  gradH2O1(pdf = FALSE)
  gradH2O2(pdf = FALSE, vars = "H2O-ZC")
  gradH2O3(pdf = FALSE, var = NULL)
  gradH2O4(pdf = FALSE)
  gradH2O5(pdf = FALSE)
  gradH2O6(pdf = FALSE)
  gradH2O7(pdf = FALSE)
  gradH2O8(pdf = FALSE)
  NifProteomes()
}

\arguments{
  \item{pdf}{logical, make a PDF file?}
  \item{vars}{character, variables to plot (passed to \code{\link{pcomp}})}
  \item{var}{character, variable to plot (passed to \code{\link{mplot}})}
}

\details{
This table gives a brief description of each function.

\tabular{ll}{
  \code{gradH2O0} \tab Number of reactions for each amino acid in \emph{E. coli} metabolic reconstruction (no plot) \cr
  \code{gradH2O1} \tab Derivation of stoichiometric hydration state \cr
  \code{gradH2O2} \tab \nH2O-\ZC scatterplots for redox gradients and the Baltic Sea \cr
  \code{gradH2O3} \tab \nH2O for Baltic Sea metagenome and metatranscriptome in different size fractions \cr
  \code{gradH2O4} \tab \nH2O-\ZC for freshwater, marine, and hypersaline environments \cr
  \code{gradH2O5} \tab Scatterplots of GRAVY vs \ZC and \nH2O \cr
  \code{gradH2O6} \tab \nH2O-\ZC and GRAVY-pI plots for Baltic Sea and Rodriguez-Brito et al. data \cr
  \code{gradH2O7} \tab Differential gene and protein expression, time-course and NaCl vs organic solutes \cr
  \code{gradH2O8} \tab Differentially expressed proteins in halophiles under hypo- and hyperosmotic stress \cr
}

\code{NifProteomes}, which is used in \code{gradH2O2}, calculates the mean and standard deviation of \ZC and \nH2O for proteomes encoding different homologs of nitrogenase (Nif-D, Nif-C, Nif-B, or Nif-A), as listed by Poudel et al. (2018).
The organisms are limited to those that have at least 1000 RefSeq protein sequences, and any duplicated organism names are removed.
}

\section{Files in extdata/gradH2O}{
  \describe{
    \item{\code{Table_S1.xlsx}}{Accession numbers used and sequence processing statistics.}
    \item{\code{MGP}}{Amino acid compositions sampled from metagenomically derived protein sequences, generated as described in \code{\link{gradox}}.}
    \item{\code{mkrds.R}}{Script to make the RDS files from multiple CSV files produced by the \dQuote{ARAST} workflow.}
    \item{\code{Nif_homolog_genomes.csv}, \code{Nif_homolog_genomes.R}}{
      Names of genomes (organisms) containing different nitrogenase (Nif) homologs, extracted from Supplemental Table 1A of Poudel et al. (2018).
      To the source table have been added the columns \samp{Refseq.name} and \samp{taxid}, which contain the closest matching genome name in release 95 of the RefSeq database (see \code{\link{JMDplots-package}}) and corresponding taxonomic ID.
      The R file is the script used to make \code{Nif_homolog_genomes.csv} from the source table.
    }
    \item{\code{reaction_equations.csv}}{
      Equations for metabolic reactions in \emph{E. coli} from the \samp{reactions_GPRs} worksheet from Supplementary Information 1 of Feist et al. (2007).
      This is used in \code{gradH2O0}.
    }
  }
}

\section{Files in extdata/organisms}{
  \describe{
    \item{\code{ecoli.csv.xz}}{
      This data file has amino acid compositions of proteins in the UniProt reference proteome of \emph{E. coli} K12 (\url{https://www.uniprot.org/proteomes/UP000000625}).
      The data frame was created by running \code{ecoli <- CHNOSZ::read.fasta("UP000000625_83333.fasta.gz")} and moving \samp{ECOLI} from the protein name to the \code{organism} column.
      This is used in \code{gradH2O1}.
    }
  }
}

\references{
Dick, J. M., Yu, M. and Tan, J. (2020) Uncovering chemical signatures of salinity gradients through compositional analysis of protein sequences. \emph{Biogeosciences Discussions}. \doi{10.5194/bg-2020-146}

Feist, A. M. et al. (2007) A genome-scale metabolic reconstruction for \emph{Escherichia coli} K-12 MG1655 that accounts for 1260 ORFs and thermodynamic information. \emph{Mol. Syst. Biol.} \bold{3}, 121. \doi{10.1038/msb4100155}

Poudel, S. et al. (2018) Electron transfer to nitrogenase in different genomic and metabolic backgrounds. \emph{J. Bacteriol.} \bold{200}, e00710.1038/ismej.2010.157-17. \doi{10.1128/JB.00757-17}

Rodriguez-Brito, B. et al. (2010) Viral and microbial community dynamics in four aquatic environments. \emph{ISME J.} \bold{4}, 739--751. \doi{10.1038/ismej.2010.1}
}

\examples{
gradH2O2()

# list ZC and nH2O of proteomes for Nif-encoding genomes
NifProteomes()
}
