\encoding{UTF-8}
\name{aoscp}
\alias{aoscp}
\alias{aoscp1}
\alias{aoscp2}
\alias{aoscp3}
\alias{aoscp4}
\title{Plots from oxidation state of carbon paper}
\description{
Plots from the paper by Dick (2014).
}

\usage{
  aoscp1(pdf = FALSE)
  aoscp2(pdf = FALSE)
  aoscp3(png = FALSE, outline = FALSE)
}

\arguments{
  \item{pdf}{logical, make a PDF file?}
  \item{png}{logical, make a PNG file?}
  \item{outline}{logical, draw only the outlines?}
}

\details{
This table gives a brief description of each function.

\tabular{ll}{
  \code{aoscp1} \tab \ZC of amino acids vs \ZC of RNA codons and hydropathy index of amino acids \cr
  \code{aoscp2} \tab Histograms of \ZC of all human proteins and human membrane proteins \cr
  \code{aoscp3} \tab Drawing of yeast cell color-coded with median ZC of proteins in different locations \cr
  \code{aoscp4} \tab \ZC and Eh ranges in yeast and ER-cytoplasm electron-transfer scheme \cr
}
}

\section{\code{extdata/aoscp}}{
  \describe{
    \item{\code{SGD_associations.csv}}{For yeast genes, this table lists the accessions, SGDID, and the association to cellular components in the Gene Ontology, derived from gene_association.sgd.gz, protein_properties.tab and go_terms.tab downloaded from \url{http://www.yeastgenome.org} on 2013-08-24. All gene associations with the NOT qualifier were removed, as were those without a matching entry in protein_properties.tab (e.g. RNA-coding genes).}
    \item{\code{ZC_HUMAN.csv}, \code{ZC_membrane.csv}}{Compilations of the values of \ZC for human proteins and human membrane proteins. Values in \code{ZC_HUMAN.csv} were calculated from protein sequences in HUMAN.fasta.gz, downloaded from \url{ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/proteomes/HUMAN.fasta.gz} on 2013-08-24 (file dated 2013-07-24). Values in \code{ZC_membrane.csv} were calculated from protein sequences in all *.fa files in Additional File 2 of Almén et al. (2009)}
    \item{\code{codons.csv}}{In the first column, the three-letter abbreviations for each of the RNA codons; in the second column, the names of the corresponding amino acids.}
    \item{\code{rubisco.csv}}{UniProt IDs for Rubisco and optimal growth temperatures of organisms.}
    \item{\code{cell/*.png}}{PNG images for each of the cellular components used in \code{aoscp3}.}
    \item{\code{protein_refseq61.csv.xz}}{Amino acid compositions of proteins in different organisms, based on RefSeq release 61 (see \code{extdata/refseq} for an updated file). A trimmed version of this file was previously present in \pkg{CHNOSZ} (to version 1.3.3). Here, the columns are \code{protein} (\dQuote{refseq}), \code{organism} (taxonomic id), \code{ref} (names of RefSeq source files (followed by number of sequences and number of amino acids) [organism name in square brackets]), \code{abbrv} (total number of amino acids), \code{chains} (1). The remaining columns have the average frequency of each amino acid in the proteins.}
    \item{\code{prep.R}}{Code use to prepare the data files from various sources.}
  }
}

\references{
Dick, J. M. (2014) Average oxidation state of carbon in proteins. \emph{J. R. Soc. Interface} \bold{11}, 20131095. \doi{10.1098/rsif.2013.1095}

Almén, M. S., Nordström, K. J. V., Fredriksson, R. and Schiöth, H. B. (2009) Mapping the human membrane proteome: a majority of the human membrane proteins can be classified according to function and evolutionary origin. \emph{BMC Biol.} \bold{7}, 50. \doi{10.1186/1741-7007-7-50}
}

\examples{
aoscp1()
}