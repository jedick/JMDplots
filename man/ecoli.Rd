\encoding{UTF-8}
\name{ecoli}
\alias{ecoli}
\docType{data}
\title{Amino acid compositions of proteins in Escherichia coli}
\description{
  This file has the amino acid compositions of proteins in the UniProt reference proteome of E. coli K12 (\url{https://www.uniprot.org/proteomes/UP000000625}).
  The file was created by running \code{ecoli <- CHNOSZ::read.fasta("UP000000625_83333.fasta.gz")} followed by \code{save(ecoli, file = "ecoli.Rdata")}.
}

\usage{
  data(ecoli)
}

\seealso{
  Used in \code{\link{gradH2O1}}.
}
