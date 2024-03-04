\encoding{UTF-8}
\name{mjenergy}
\alias{mjenergy}
\alias{mjenergy1}
\alias{mjenergy2}
\alias{mjenergy3}
\alias{calc_affinity}
\alias{mjenergy_Table_S3}
\alias{mjenergy_Dataset_S1}
\title{Energy release during protein synthesis at ultramafic-hosted hydrothermal systems}
\description{
Plots from the paper by Dick and Shock (2021).
Affinities are calculated for each protein in the \emph{Methanocaldococcus jannaschii} (\emph{Mj}) proteome.
Amino acid compositions are taken from the file \code{UP000000805_243232.csv.xz} (see \code{\link{JMDplots-package}}).
}

\usage{
  mjenergy1(pdf = FALSE)
  mjenergy2(pdf = FALSE)
  mjenergy3(pdf = FALSE, write.csv = FALSE)
  calc_affinity(T = 85, protein = "CSG_METJA")
  mjenergy_Table_S3()
  mjenergy_Dataset_S1(write.csv = FALSE)
}

\arguments{
  \item{pdf}{logical, make a PDF file?}
  \item{T}{numeric, temperature in \degC}
  \item{protein}{character, protein name used in CHNOSZ (see \code{\link[CHNOSZ]{pinfo}})}
  \item{write.csv}{logical, save a CSV file with the calculation results?}
}

\details{
This table gives a brief description of each plotting function.

\tabular{ll}{
  \code{mjenergy1} \tab Affinities for methanogenesis and amino acid synthesis \cr
  \code{mjenergy2} \tab ZC of amino acids vs frequency in the \emph{Mj} proteome \cr
  \code{mjenergy3} \tab Affinities of overall synthesis of proteins in the \emph{Mj} proteome \cr
}

\code{calc_affinity} calculates and prints the affinity for amino acid synthesis and polymerization for a single \code{protein} at temperature \code{T}.
The default arguments correspond to the example calculation discussed in the paper (CSG_METJA at 85 \degC).

The files described below correspond to Tables S1 and S2 in the supporting information for the paper.
\code{mjenergy_Table_S3} outputs the amino acid synthesis reactions that are listed in Table S3.
\code{mjenergy_Dataset_S1} outputs the UniProt accession number, elemental composition, sequence length, and \Zc for each protein in the \emph{Mj} proteome.
Use \code{mjenergy_Dataset_S1(write.csv = TRUE)} to instead save the results in a file named \samp{Dataset_S1.csv}, which corresponds to Dataset S1 in the paper.

Use \code{mjenergy3(write.csv = TRUE)} to save the affinities (MJ / mol protein) in two files named \samp{Dataset_S2.csv} and \samp{Dataset_S3.csv} for Rainbow and Endeavor vent fields.
The column names are fluid temperatures (\degC) and the rownames are UniProt accession numbers of the proteins.
Regardless of the value of \code{write.csv}, \code{mjenergy} returns the affinities of the proteins invisibly in a list with two elements named \samp{Rainbow} and \samp{Endeavor}

}

\section{Files in extdata/mjenergy}{
  \describe{
    \item{\code{SC10_Rainbow.csv}}{Values of temperature (\degC), pH and logarithms of activity of \CO2, \H2, \NH4plus, \H2S and \CH4 for mixing of seawater and hydrothermal fluid at Rainbow field (Mid-Atlantic Ridge), taken from Shock and Canovas (2010).}
    \item{\code{SC10_Endeavor.csv}}{Values of temperature (\degC), pH and logarithms of activity of \CO2, \H2, \NH4plus, \H2S and \CH4 for mixing of seawater and hydrothermal fluid at Endeavor field (Juan de Fuca Ridge), taken from Shock and Canovas (2010).}
  }
}

\references{
Shock E and Canovas P (2010) The potential for abiotic organic synthesis and biosynthesis at seafloor hydrothermal systems. \emph{Geofluids} \bold{10}, 161--192. \doi{10.1111/j.1468-8123.2010.00277.x}

Dick JM and Shock EL (2021) The release of energy during protein synthesis at ultramafic-hosted submarine hydrothermal ecosystems. \emph{J. Geophys. Res.: Biogeosciences} \bold{126}, e2021JG006436. \doi{10.1029/2021JG006436}
}

\examples{
mjenergy1()
}
