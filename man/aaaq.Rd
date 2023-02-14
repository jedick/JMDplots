\encoding{UTF-8}
\name{aaaq}
\alias{aaaq}
\alias{aaaq4}
\alias{aaaq10}
\alias{aaaq13}
\title{Amino acid group additivity}
\description{
Plots from the paper by Dick et al. (2006).
}

\usage{
  aaaq4(pdf = FALSE)
  aaaq10(pdf = FALSE)
  aaaq13(pdf = FALSE)
}

\arguments{
  \item{pdf}{logical, make a PDF file?}
}

\details{
This table gives a brief description of each function.

\tabular{ll}{
  \code{aaaq4} \tab \pKa of sidechain and amino acid backbone groups \cr
  \code{aaaq10} \tab Net charge of proteins as a function of pH and \T \cr
  \code{aaaq13} \tab Eh-pH diagram for extracellular alpha-amylases \cr
}
}

\section{Files in extdata/aaaq}{
  \describe{
    \item{\code{BKM60_Fig7.csv}}{Eh-pH values for normal, wet and waterlogged soils from Fig. 7 of Baas Becking et al. (1960).}
    \item{\code{NT67.csv}}{Hydrogen ion titration curve for unfolded ribonuclease (\samp{RNAS1_BOVN}) from Nozaki and Tanford (1967).}
    \item{\code{RT71.csv}}{Hydrogen ion titration curve for unfolded lysozyme (\samp{LYSC_CHICK}) from Roxby and Tanford (1971).}
    \item{\code{WG00.csv}}{Hydrogen ion titration curve for unfolded staphylococcal nuclease (\samp{NUC_STAAU}) from Whitten and García-Morena E. (2000).}
  }
}

\references{
Baas Becking LGM, Kaplan IR and Moore D (1960) Limits of the natural environment in terms of pH and oxidation-reduction potentials. \emph{Journal of Geology} \bold{68}, 243--284. \url{https://www.jstor.org/stable/30059218}

Dick JM, LaRowe DE and Helgeson HC (2006) Temperature, pressure, and electrochemical constraints on protein speciation: Group additivity calculation of the standard molal thermodynamic properties of ionized unfolded proteins. \emph{Biogeosciences} \bold{3}, 311--336. \doi{10.5194/bg-3-311-2006}

Nozaki Y and Tanford C (1967) Proteins as random coils. II. Hydrogen ion titration curve of ribonuclease in 6 M guanidine hydrochloride \emph{J. Amer. Chem. Soc.} \bold{89}, 742--749. \doi{10.1021/ja00980a003}

Roxby R and Tanford C (1971) Hydrogen ion titration curve of lysozyme in 6 M guanidine hydrochloride. \emph{Biochemistry} \bold{10}, 3348--3352. \doi{10.1021/bi00794a005}

Whitten ST and García-Morena E. B (2000) pH dependence of stability of staphylococcal nuclease: Evidence of substantial electrostatic interactions in the denatured state. \emph{Biochemistry} \bold{39}, 14292--14304. \doi{10.1021/bi001015c}
}

\examples{
aaaq4()
}
