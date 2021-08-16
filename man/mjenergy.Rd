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
  \describe{
    \item{\code{SC10_Rainbow.csv}}{Values of temperature (\degC), pH and logarithms of activity of \CO2, \H2, \NH4plus, \H2S and \CH4 for mixing of seawater and hydrothermal fluid at Rainbow field (Mid-Atlantic Ridge), taken from Shock and Canovas (2010).}
    \item{\code{SC10_Endeavor.csv}}{Values of temperature (\degC), pH and logarithms of activity of \CO2, \H2, \NH4plus, \H2S and \CH4 for mixing of seawater and hydrothermal fluid at Endeavor field (Juan de Fuca Ridge), taken from Shock and Canovas (2010).}
  }
}

\references{
Shock, E. and Canovas, P. (2010) The potential for abiotic organic synthesis and biosynthesis at seafloor hydrothermal systems. \emph{Geofluids} \bold{10}, 161--192. \doi{10.1111/j.1468-8123.2010.00277.x}

Shock, E. L. and Dick, J. M. (2021) The release of energy during protein synthesis at submarine hydrothermal ecosystems. (Preprint)
}

\examples{
mjenergy1()
}
