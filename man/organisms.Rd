\encoding{UTF-8}
\name{organisms}
\alias{organisms}
\alias{ecoli}
\title{Amino acid compositions of proteins in Escherichia coli}
\description{
  This data frame has amino acid compositions of proteins in the UniProt reference proteome of E. coli K12 (\url{https://www.uniprot.org/proteomes/UP000000625}).
  The data file was created by running \code{ecoli <- CHNOSZ::read.fasta("UP000000625_83333.fasta.gz")} followed by \code{saveRDS(ecoli, "ecoli.rds")}.
}

\seealso{
  Used in \code{\link{gradH2O1}}.
}
