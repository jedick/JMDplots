\encoding{UTF-8}
\name{bison}
\alias{bison}
\alias{bison1}
\alias{bison2}
\alias{bison3}
\alias{bison4}
\alias{bison5}
\alias{bison6}
\alias{bison7}
\alias{bison8}
\title{Plots from Bison Pool hot spring papers}
\description{
Plots from the papers by Dick and Shock (2011 and 2013).
}

\usage{
  bison1()
  bison2()
  bison3()
  bison4()
  bison5()
  bison6()
  bison7(equil.results)
  bison8(equil.results)
}

\arguments{
  \item{equil.results}{list, value returned by \code{bison6}}
}

\details{
This table gives a brief description of each function.

\tabular{ll}{
  \code{bison1} \tab Measured \T and pH \cr
  \code{bison2} \tab Carbon oxidation state of proteins \cr
  \code{bison3} \tab Chemical affinities (table) \cr
  \code{bison4} \tab Relative stabilities along a temperature and chemical gradient \cr
  \code{bison5} \tab Comparing old and new methionine sidechain parameters \cr
  \code{bison6} \tab Metastable equilibrium model for relative abundances \cr
  \code{bison7} \tab Activity of hydrogen comparison \cr
  \code{bison8} \tab Comparison of model and observed abundances \cr
}
}

\section{Files in extdata/bison}{
  \describe{
    \item{\code{DS11.csv}, \code{DS13.csv}}{
      These files contain average amino acid compositions of groups of proteins inferred from metagenomic sequences.
      The groups correspond to proteins with a common functional annotation (DS11) or taxonomic affiliation (DS13).
      The data are from Dick and Shock, 2011 and 2013.
    }
  }
}

\references{
Dick, J. M. and Shock, E. L. (2011) Calculation of the relative chemical stabilities of proteins as a function of temperature and redox chemistry in a hot spring. \emph{PLOS One} \bold{6}, e22782. \doi{10.1371/journal.pone.0022782}

Dick, J. M. and Shock, E. L. (2013) A metastable equilibrium model for the relative abundance of microbial phyla in a hot spring. \emph{PLOS One} \bold{8}, e72395. \doi{10.1371/journal.pone.0072395}
}

\examples{
bison1()
}
