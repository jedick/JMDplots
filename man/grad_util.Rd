\encoding{UTF-8}
\name{grad_util}
\alias{grad_util}
\alias{H2OAA}
\alias{ZCAA}
\alias{NifProteomes}
\alias{GRAVY}
\title{Utility functions for chemical composition of proteins}
\description{
These functions calculate compositional metrics of proteins given a data frame of amino acid compositions.
}

\usage{
  ZCAA(AAcomp, nothing = NULL)
  H2OAA(AAcomp, basis = "rQEC")
  GRAVY(AAcomp)
  NifProteomes()
}

\arguments{
  \item{AAcomp}{data frame, amino acid compositions}
  \item{nothing}{dummy argument}
  \item{basis}{character, basis species}
}

\details{
These functions are written to enable rapid calculations of carbon oxidation state (\ZC) and per-residue stoichiometric water content (\nH2O) of proteins.
For calculations of \ZC from chemical formulas of organic molecules (not only proteins), use the \code{\link[CHNOSZ]{ZC}} function in \pkg{CHNOSZ}.

The default for \code{basis} stands for relative water content, calculated as residuals of a linear model for \nH2O vs \ZC of the 20 common amino acids using the basis species glutamine - glutamic acid - cysteine - \H2O - \O2.
This can be changed to \samp{QEC} to use the stoichiometric water content calculated directly from the basis species (not the residuals of the model).
See Dick et al. (2019) for more information.

Note that \ZC is a per-carbon average, but \nH2O is a per-residue average.
The contribution of \H2O from the terminal groups of proteins is counted, so shorter proteins have slightly greater \nH2O.

\code{nothing} is a dummy argument so that these functions can be called with the same numbers of arguments (used internally by some plotting functions).

\code{GRAVY} calculates the grand average of hydropathicity using values of the hydropathy index for individual amino acids given by Kyte and Doolittle (1982).
The resulting averages are equal to those calculated with the ProtParam tool (\url{https://web.expasy.org/protparam/}).

\code{NifProteomes} calculates the mean and standard deviation of \ZC and \nH2O for proteomes encoding different homologs of nitrogenase (Nif-D, Nif-C, Nif-B, or Nif-A), as listed by Poudel et al. (2018).
The organisms are limited to those that have at least 1000 RefSeq protein sequences, and any duplicated organism names are removed.
}

\section{Files in extdata/gradH2O}{
  \describe{
    \item{\code{Nif_homolog_genomes.csv}, \code{Nif_homolog_genomes.R}}{
      Names of genomes (organisms) containing different nitrogenase (Nif) homologs, extracted from Supplemental Table 1A of Poudel et al. (2018).
      To the source table have been added the columns \samp{Refseq.name} and \samp{taxid}, which contain the closest matching genome name in release 95 of the RefSeq database (see \code{\link{JMDplots-package}}) and corresponding taxonomic ID.
      The R file is the script used to make \code{Nif_homolog_genomes.csv} from the source table.
    }
  }
}

\references{
Dick, J. M., Tan, J. and Yu, M. (2019) Geochemistry of metagenomes: Distinct compositional signatures of environmental redox and salinity gradients. \emph{(preprint)}

Kyte, J. and Doolittle, R. F. (1982) A simple method for displaying the hydropathic character of a protein. \emph{J. Mol. Biol.} \bold{157}, 105--132. \doi{10.1016/0022-2836(82)90515-0}

Poudel, S. et al. (2018) Electron transfer to nitrogenase in different genomic and metabolic backgrounds. \emph{J. Bacteriol.} \bold{200}, e00757-17. \doi{10.1128/JB.00757-17}
}

\examples{
# for reference, compute ZC of alanine and glycine "by hand"
ZC.Gly <- ZC("C2H5NO2")
ZC.Ala <- ZC("C3H7NO2")
# define the composition of a Gly-Ala-Gly tripeptide
AAcomp <- data.frame(Gly = 2, Ala = 1)
# calculate the ZC of the tripeptide (value: 0.571)
ZC.GAG <- ZCAA(AAcomp)
# this is equal to the carbon-number-weighted average of the amino acids
nC.Gly <- 2 * 2
nC.Ala <- 1 * 3
ZC.average <- (nC.Gly * ZC.Gly + nC.Ala * ZC.Ala) / (nC.Ala + nC.Gly)
stopifnot(all.equal(ZC.GAG, ZC.average))

# compute the per-residue nH2O of Gly-Ala-Gly
basis("QEC")
nH2O.GAG <- species("Gly-Ala-Gly")$H2O
# divide by the length to get residue average (we keep the terminal H-OH)
nH2O.residue <- nH2O.GAG / 3
# compare with the value calculated by H2OAA() (-0.2)
nH2O.H2OAA <- H2OAA(AAcomp, "QEC")
stopifnot(all.equal(nH2O.residue, nH2O.H2OAA))

# list ZC and nH2O of proteomes for Nif-encoding genomes
NifProteomes()

# calculate GRAVY for lysozyme and ribonuclease
# first get the protein index in CHNOSZ's list of proteins
iprotein <- pinfo(c("LYSC_CHICK", "RNAS1_BOVIN"))
# then get the amino acid compositions
AAcomp <- pinfo(iprotein)
# then calculate GRAVY
Gcalc <- as.numeric(GRAVY(AAcomp))
# these are equal to values obtained with ProtParam
# (uniprot.org --> LYSC_CHICK and RNAS1_BOVIN)
Gref <- c(-0.472, -0.663)
stopifnot(all.equal(round(Gcalc, 3), Gref))
}
