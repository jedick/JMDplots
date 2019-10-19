\encoding{UTF-8}
\name{gradH2O}
\alias{gradH2O}
\alias{gradH2O1}
\alias{gradH2O2}
\alias{gradH2O3}
\alias{gradH2O4}
\alias{gradH2O5}
\alias{gradH2O6}
\title{Plots from salinity gradients paper}
\description{
Plots from the paper by Dick et al. (2019).
}

\usage{
  gradH2O1(pdf = FALSE)
  gradH2O2(pdf = FALSE)
  gradH2O3(pdf = FALSE)
  gradH2O4(pdf = FALSE)
  gradH2O5(pdf = FALSE)
  gradH2O6(pdf = FALSE)
}

\arguments{
  \item{pdf}{logical, make a PDF file?}
}

\details{
This table gives a brief description of each function.

\tabular{ll}{
  \code{gradH2O1} \tab Basis species comparison \cr
  \code{gradH2O2} \tab \ZC of proteins for selected redox gradients \cr
  \code{gradH2O3} \tab \nH2O-\ZC scatterplots for redox gradients and the Baltic Sea \cr
  \code{gradH2O4} \tab \nH2O for Baltic Sea size fractions \cr
  \code{gradH2O5} \tab \nH2O-\ZC for freshwater, marine, and hypersaline environments \cr
  \code{gradH2O6} \tab Mean differences of \nH2O and \ZC for differentially expressed proteins in hyperosmotic stress \cr
}
}

\section{extdata/gradH2O}{
  \describe{
    \item{\code{MGP}}{Amino acid compositions sampled from metagenomically derived protein sequences, generated as described in \code{\link{gradox}}.}
    \item{\code{Nif_homolog_genomes.csv}, \code{Nif_homolog_genomes.R}}{Names of genomes (organisms) containing different nitrogenase (Nif) homologs, extracted from Supplemental Table 1A of Poudel et al. (2018). To the source table have been added the columns \samp{Refseq.name} and \samp{taxid}, which contain the closest matching genome name in the RefSeq database (see \code{\link{JMDplots-package}}) and corresponding taxonomic ID. The R file is the script used to make \code{Nif_homolog_genomes.csv} from the source table.}
  }
}

\section{extdata/organisms}{
  \describe{
    \item{\code{ecoli.csv.xz}}{
      This data file has amino acid compositions of proteins in the UniProt reference proteome of E. coli K12 (\url{https://www.uniprot.org/proteomes/UP000000625}).
      The data frame was created by running \code{ecoli <- CHNOSZ::read.fasta("UP000000625_83333.fasta.gz")} and moving \samp{ECOLI} from the protein name to the \code{organism} column.
      This is used in \code{gradH2O1}.
    }
  }
}


\references{
Dick, J. M., Tan, J. and Yu, M. (2019) Geochemistry of metagenomes: Distinct compositional signatures of environmental redox and salinity gradients. \emph{(preprint)}

Poudel, S. et al. (2018) Electron transfer to nitrogenase in different genomic and metabolic backgrounds. \emph{J. Bacteriol.} \bold{200}, e00757-17. \doi{10.1128/JB.00757-17}
}

\examples{
gradH2O1()
}
