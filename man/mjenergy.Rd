\encoding{UTF-8}
\name{mjenergy}
\alias{mjenergy}
\alias{mjenergy1}
\alias{mjenergy2}
\alias{mjenergy3}
\alias{calc_affinity}
\title{Plots from energy release in protein synthesis paper}
\description{
Plots from the paper by Dick and Shock (2021).
}

\usage{
  mjenergy1(pdf = FALSE)
  mjenergy2(pdf = FALSE)
  mjenergy3(pdf = FALSE)
  calc_affinity(T = 85, protein = "CSG_METJA")
}

\arguments{
  \item{pdf}{logical, make a PDF file?}
  \item{T}{numeric, temperature in \degC}
  \item{protein}{character, protein name used in CHNOSZ (see \code{\link[CHNOSZ]{pinfo}})}
}

\details{
This table gives a brief description of each plotting function.

\tabular{ll}{
  \code{mjenergy1} \tab Affinities for methanogenesis and amino acid synthesis \cr
  \code{mjenergy2} \tab ZC of amino acids vs frequency in \emph{Mj} proteome \cr
  \code{mjenergy3} \tab Affinities of overall synthesis of proteins in \emph{Mj} proteome \cr
}

\code{calc_affinity} calculates and prints the affinity for amino acid synthesis and polymerization for a single \code{protein} at temperature \code{T}.
The default arguments correspond to the example calculation discussed in the paper (CSG_METJA at 85 \degC).
}

\section{Files in extdata/mjenergy}{
  \code{UP000000805_243232.csv.gz} has amino acid compositions of 1787 proteins in \emph{Methanocaldococcus jannaschii}.
  It was created by processing the file \code{UP000000805_243232.fasta.gz}, which was downloaded from the UniProt reference proteomes FTP site (\url{ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Archaea/}).
  The server timestamp on \code{UP000000805_243232.fasta.gz} was Dec. 2, 2020; the reference proteome was last modified on August 22, 2020 according to \url{https://www.uniprot.org/proteomes/UP000000805}.
  The commands used to create \code{UP000000805_243232.csv.gz} are:
  \verb{
    aa <- CHNOSZ::read.fasta("UP000000805_243232.fasta.gz")
    write.csv(aa, "UP000000805_243232.csv", row.names = FALSE, quote = FALSE)
    system("gzip UP000000805_243232.csv")
  }
}

\references{
Shock, E. L. and Dick, J. M. (2021) The release of energy during protein synthesis at submarine hydrothermal ecosystems. (Preprint)
}

\examples{
mjenergy1()
}
