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
\alias{protein.equil}
\title{Relative stabilities of proteins inferred from a hot-spring metagenome}
\description{
Plots from the papers by Dick and Shock (2011 and 2013).
}

\usage{
  bison1()
  bison2(plots = 1:2, add.titles = TRUE)
  bison3()
  bison4()
  bison5()
  bison6(plot.it = TRUE)
  bison7(equil.results)
  bison8(equil.results)
  protein.equil(protein, T=25, loga.protein = 0, digits = 4)
}

\arguments{
  \item{plots}{numbers of plots to make}
  \item{add.titles}{logical, add titles to plots?}
  \item{plot.it}{logical, make a plot?}
  \item{equil.results}{list, value returned by \code{bison6}}
  \item{protein}{character, names of proteins; numeric, species index of proteins; data frame; amino acid composition of proteins}
  \item{T}{numeric, temperature in \degC}
  \item{loga.protein}{numeric, decimal logarithms of reference activities of proteins}
  \item{digits}{integer, number of significant digits (see \code{\link{signif}})}
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

\code{protein.equil} produces a series of messages showing step-by-step a calculation of the chemical activities of proteins in metastable equilibrium.
For the first protein, it shows the standard Gibbs energies of the reaction to form the nonionized protein from the basis species and of
the ionization reaction of the protein (if \samp{H+} is in the basis),
then the standard Gibbs energy/RT of the reaction to form the (possibly ionized) protein per residue.
The per-residue values of \samp{logQstar} and \samp{Astar/RT} are also shown for the first protein.
Equilibrium calculations are then performed, only if more than one protein is specified.
This calculation applies the Boltzmann distribution to the calculation of the equilibrium degrees of formation of the residue equivalents of the proteins,
then converts them to activities of proteins taking account of \code{loga.protein} and protein length.
If the \code{protein} argument is numeric (indicating rownumbers in \code{thermo()$protein}),
the values of \samp{Astar/RT} are compared with the output of \code{\link[CHNOSZ]{affinity}},
and those of the equilibrium degrees of formation of the residues and the chemical activities of the proteins with the output of \code{\link[CHNOSZ]{diagram}}.
If the values in any of these tests are are not \code{\link{all.equal}} an error is produced indicating a bug. 
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
Dick JM and Shock EL (2011) Calculation of the relative chemical stabilities of proteins as a function of temperature and redox chemistry in a hot spring. \emph{PLOS One} \bold{6}, e22782. \doi{10.1371/journal.pone.0022782}

Dick JM and Shock EL (2013) A metastable equilibrium model for the relative abundance of microbial phyla in a hot spring. \emph{PLOS One} \bold{8}, e72395. \doi{10.1371/journal.pone.0072395}
}

\examples{
bison1()
}
