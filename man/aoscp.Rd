\encoding{UTF-8}
\name{aoscp}
\alias{aoscp}
\alias{aoscp1}
\alias{aoscp2}
\alias{aoscp3}
\alias{aoscp4}
\alias{aoscp5}
\alias{aoscp6}
\alias{aoscp99}
\title{Average oxidation state of carbon in proteins}
\description{
Plots from the paper by Dick (2014).
}

\usage{
  aoscp1(pdf = FALSE)
  aoscp2(pdf = FALSE)
  aoscp3(png = FALSE, pdf = FALSE, outline = FALSE)
  aoscp4(pdf = FALSE)
  aoscp5(pdf = FALSE, file = NULL)
  aoscp6(pdf = FALSE)
  aoscp99(pdf = FALSE)
}

\arguments{
  \item{pdf}{logical, make a PDF file?}
  \item{png}{logical, make the base PNG file (no labels)?}
  \item{outline}{logical, draw only the outlines?}
  \item{file}{character, file with amino acid compositions of organisms}
}

\details{
This table gives a brief description of each function.

\tabular{ll}{
  \code{aoscp1} \tab \ZC of amino acids vs \ZC of RNA codons and hydropathy index of amino acids \cr
  \code{aoscp2} \tab Histograms of \ZC of all human proteins and human membrane proteins \cr
  \code{aoscp3} \tab Drawing of yeast cell color-coded with median ZC of proteins in different locations \cr
  \code{aoscp4} \tab \ZC and Eh ranges in yeast and ER-cytoplasm electron-transfer scheme \cr
  \code{aoscp5} \tab Average oxidation state of carbon in proteins from different organisms \cr
  \code{aoscp6} \tab \ZC and \Topt of different rubiscos and thermodynamic comparison \cr
  \code{aoscp99} \tab \ZC of ferredoxin, thioredoxin, and glutaredoxin vs midpoint reduction potential (unpublished) \cr
}

NOTE: for the default \code{file} = NULL in \code{aoscp5}, values of \ZC are taken from the packaged file \code{ZC_refseq.csv.xz} (see below).
If \code{file} is not NULL, then value of \ZC are calculated from the amino acid compositions and the file \code{ZC_refseq.csv} is written to the working directory.
}

\section{Files in extdata/aoscp}{
  \describe{
    \item{\code{SGD_associations.csv.xz}}{
      For yeast genes, this table lists the accessions, SGDID, and the association to cellular components in the Gene Ontology, derived from gene_association.sgd.gz, protein_properties.tab and go_terms.tab downloaded from \url{https://www.yeastgenome.org/} on 2013-08-24.
      All gene associations with the NOT qualifier were removed, as were those without a matching entry in protein_properties.tab (e.g. RNA-coding genes).
    }
    \item{\code{ZC_HUMAN.csv.xz}, \code{ZC_membrane.csv.xz}}{
      Compilations of the values of \ZC for human proteins and human membrane proteins.
      Values in \code{ZC_HUMAN.csv} were calculated from protein sequences in HUMAN.fasta.gz, downloaded from \url{https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/} on 2013-08-24 (file dated 2013-07-24).
      Values in \code{ZC_membrane.csv} were calculated from protein sequences in all *.fa files in Additional File 2 of Almén et al. (2009).
      NOTE: To save space in the package, the accession and length columns have been removed from these files; the complete files are present in the Supplementary data of Dick (2014).
    }
    \item{\code{codons.csv}}{In the first column, the three-letter abbreviations for each of the RNA codons; in the second column, the names of the corresponding amino acids.}
    \item{\code{rubisco.csv}}{UniProt IDs for Rubisco and optimal growth temperatures of organisms.}
    \item{\code{rubisco/*.fasta}}{FASTA sequence files for each Rubisco identified in table 3.}
    \item{\code{cell/*.png}}{PNG images for each of the cellular components used in \code{aoscp3}.}
    \item{\code{midpoint.csv}}{
      List of protein names, UniProt IDs and standard midpoint reduction potentials of ferredoxin and thioredoxin in spinach and glutaredoxin in \emph{E. coli}.
      Start and stop positions, taken from UniProt, identify the protein chain excluding initiator methionines or signal peptides.
    }
    \item{\code{midpoint/*.fasta}}{FASTA sequence files for the redoxins.}
    \item{\code{ZC_refseq.csv.xz}}{
      This file has the taxid, organism name, number of amino acids, and calculated \ZC in 6758 microbial organisms.
      This file was created using \code{aoscp5(file = "protein_refseq.csv")}, where \code{"protein_refseq.csv"} has the average amino acid compositions of proteins in different organisms.
      That file, which is based on RefSeq release 61, is in the Supplementary data of Dick (2014) (\url{https://royalsocietypublishing.org/doi/suppl/10.1098/rsif.2013.1095}); a trimmed version of the file was previously present in \pkg{CHNOSZ} (to version 1.3.3).
      In that file, the columns are \code{protein} (\dQuote{refseq}), \code{organism} (taxonomic id), \code{ref} (names of RefSeq source files (followed by number of sequences and number of amino acids) [organism name in square brackets]), \code{abbrv} (total number of amino acids), \code{chains} (1).
      The remaining columns have the average frequency of each amino acid in the proteins.
    }
    \item{\code{prep.R}}{Code use to prepare the data files from various sources.}
  }
}

\references{
Dick JM (2014) Average oxidation state of carbon in proteins. \emph{J. R. Soc. Interface} \bold{11}, 20131095. \doi{10.1098/rsif.2013.1095}

Almén MS, Nordström KJV, Fredriksson R and Schiöth HB (2009) Mapping the human membrane proteome: a majority of the human membrane proteins can be classified according to function and evolutionary origin. \emph{BMC Biol.} \bold{7}, 50. \doi{10.1186/1741-7007-7-50}
}

\examples{
aoscp1()
}
