\encoding{UTF-8}
\name{PS}
\alias{PS}
\title{Retrieve Phylostrata for Given UniProt IDs}
\description{
Retrieves the phylostrata for protein-coding genes according to Liebeskind et al. (2016) or Trigos et al. (2017).
}

\usage{
  PS(uniprot, source = "TPPG17")
}

\arguments{
  \item{uniprot}{character, UniProt accession numbers}
  \item{source}{character, \samp{TPPG17} or \samp{LMM16}}
}

\details{
The phylostratum for each protein is found by matching the UniProt ID in one of these data files:

\describe{
  \item{\code{extdata/evdevH2O/phylostrata/TPPG17.csv.xz}}{
    This file has columns \samp{GeneID} (gene name), \samp{Entrez}, \samp{Entry}, and \samp{Phylostrata}.
    Except for \samp{Entry}, the values are from Dataset S1 of Trigos et al. (2017).
    UniProt acession numbers in \samp{Entry} were generated using the UniProt mapping tool first for \samp{Entrez}, followed by \samp{GeneID} for the unmatched genes.
    \samp{Entry} is NA for genes that remain unmatched to any proteins after both mapping steps.
  }
  \item{\code{extdata/evdevH2O/phylostrata/LMM16.csv.xz}}{
    This file has columns \samp{UniProt}, \samp{modeAge}, and \samp{PS}.
    The data are from file \code{main_HUMAN.csv} in Gene-Ages v1.0 (\url{https://zenodo.org/record/51708}; Liebeskind et al. (2016)).
    The modeAges were converted to phylostrata values 1-8 (\samp{PS} column) in this order: Cellular_organisms, Euk_Archaea, Euk+Bac, Eukaryota, Opisthokonta, Eumetazoa, Vertebrata, Mammalia.
  }
}
}

\references{
Liebeskind, B. J. and McWhite, C. D. and Marcotte, E. M. (2016) Towards consensus gene ages. \emph{Genome Biol. Evol.} \bold{8}, 1812--1823. \doi{10.1093/gbe/evw113}

Trigos, A. S. and Pearson, R. B. and Papenfuss, A. T. and Goode, D. L. (2017) Altered interactions between unicellular and multicellular genes drive hallmarks of transformation in a diverse range of solid tumors. \emph{Proc. Natl. Acad. Sci.} \bold{114}, 6406--6411. \doi{10.1073/pnas.1617743114}
}

\seealso{
These data, but not this function, are used in the \code{\link{evdevH2O}} plots.
}

\examples{
# Get protein expression data for one dataset
pd <- canprot::pdat_colorectal("JKMF10")
# Get phylostrata
PS(pd$pcomp$uniprot)

# Compare the two sources
PSdir <- system.file("extdata/evdevH2O/phylostrata", package = "JMDplots")
TPPG17 <- read.csv(file.path(PSdir, "TPPG17.csv.xz"))
LMM16 <- read.csv(file.path(PSdir, "LMM16.csv.xz"))
IDs <- intersect(TPPG17$Entry, LMM16$UniProt)
PS_TPPG17 <- TPPG17$Phylostrata[match(IDs, TPPG17$Entry)]
PS_LMM16 <- LMM16$PS[match(IDs, LMM16$UniProt)]
plot(jitter(PS_TPPG17), jitter(PS_LMM16), pch = ".")
}
