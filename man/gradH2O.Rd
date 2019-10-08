\encoding{UTF-8}
\name{gradH2O}
\alias{gradH2O}
\alias{gradH2O1}
\alias{gradH2O2}
\alias{gradH2O3}
\alias{gradH2O4}
\alias{gradH2O5}
\alias{gradH2O6}
\title{Plots from salinity gradients paper}
\description{
Plots from the paper by Dick et al. (2019).
}

\usage{
  gradH2O1(pdf = FALSE)
  gradH2O2(pdf = FALSE)
  gradH2O3(pdf = FALSE)
  gradH2O4(pdf = FALSE)
  gradH2O5(pdf = FALSE)
  gradH2O6(pdf = FALSE)
}

\arguments{
  \item{pdf}{logical, make a PDF file?}
}

\details{
This table gives a brief description of each function.

\tabular{ll}{
  \code{gradH2O1} \tab Basis species comparison \cr
  \code{gradH2O2} \tab \ZC of proteins for selected redox gradients \cr
  \code{gradH2O3} \tab \nH2O-\ZC scatterplots for redox gradients and the Baltic Sea \cr
  \code{gradH2O4} \tab \nH2O for Baltic Sea size fractions \cr
  \code{gradH2O5} \tab \nH2O-\ZC for freshwater, marine, and hypersaline environments \cr
  \code{gradH2O6} \tab Mean differences of \nH2O and \ZC for differentially expressed proteins in hyperosmotic stress \cr
}
}

\references{
Dick, J. M., Tan, J. and Yu, M. (2019) Geochemistry of metagenomes: Distinct compositional signatures of environmental redox and salinity gradients. \emph{(preprint)}
}

\examples{
gradH2O1()
}
