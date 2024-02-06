\encoding{UTF-8}
\name{carrollite}
\alias{carrollite}
\alias{add_Co_aqueous}
\alias{calc_carrollite}
\alias{carrollite_5}
\alias{carrollite_8}
\alias{carrollite_S3}
\alias{carrollite_S4}
\alias{carrollite_S5}
\alias{carrollite_S6}
\title{Carrollite stability}
\description{
Plots from the paper by von der Heyden et al. (2024).
}

\usage{
  calc_carrollite()
  add_Co_aqueous()
  carrollite_5(res, pdf = FALSE)
  carrollite_8(res, pdf = FALSE)
  carrollite_S3(pdf = FALSE)
  carrollite_S4(complexes, pdf = FALSE)
  carrollite_S5(pdf = FALSE)
  carrollite_S6(res, pdf = FALSE)
}

\arguments{
  \item{res}{numeric, plot resolution}
  \item{pdf}{logical, make a PDF file?}
  \item{complexes}{list, output from \code{add_Co_aqueous}}
}

\details{
This table gives a brief description of each function.

\tabular{ll}{
  \code{calc_carrollite} \tab Calculate thermodynamic parameters for carrollite \cr
  \code{add_Co_aqueous} \tab Add thermodynamic parameters for aqueous Co complexes \cr
  \code{carrollite_5} \tab Mineral stability diagram with carrollite \cr
  \code{carrollite_8} \tab Comparison of Cu-Co and Fe-Cu diagrams \cr
  \code{carrollite_S3} \tab Experimental and fitted \Cp for carrollite \cr
  \code{carrollite_S4} \tab Experimental and fitted formation constants for Cl\S{-} complexes \cr
  \code{carrollite_S5} \tab \logK of reactions showing temperatures of linnaeite in and carrollite out \cr
  \code{carrollite_S6} \tab Compare mosaic stack (with carrollite) to \code{\link[CHNOSZ]{solubility}} calculation for Co or Cu without carrollite \cr
}
}

\references{
von der Heyden BP, Dick J, Rosenfels RC, Carlton L, Lilova K, Navrotsky A, Subramani T, Woodfield BF and Gibson A (2024) Growth and stability of stratiform carrollite (CuCo\s{2}S\s{4}) in the Tenke-Fungurume ore district, Central African Copperbelt. \emph{The Canadian Journal of Mineralogy and Petrology} \bold{62}, 77--93. \doi{10.3389/feart.2019.00180}
}

\examples{
# This adds a new species named carrollite_test to OBIGT.
# The vignette checks that the parameters of the new species are
# consistent with those of carrollite in the default OBIGT database.
calc_carrollite()
}
