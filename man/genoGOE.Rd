\encoding{UTF-8}
\name{genoGOE}
\alias{genoGOE}
\alias{genoGOE_1}
\alias{genoGOE_2}
\alias{genoGOE_3}
\alias{genoGOE_4}
\alias{genoGOE_5}
\title{Earth's oxygenation recorded in genomes}
\description{
Plots from the manuscript by Dick (2025).
}

\usage{
  genoGOE_1(pdf = FALSE, panel = NULL)
  genoGOE_2(pdf = FALSE, metric = "Zc")
  genoGOE_3(pdf = FALSE)
  genoGOE_4(pdf = FALSE, panel = NULL)
  genoGOE_5(pdf = FALSE)
}

\arguments{
  \item{pdf}{logical, make a PDF file?}
  \item{panel}{character, create only this panel (disables figure layout, PDF output, and panel labeling)}
  \item{metric}{character, which metric to plot (see \code{\link[canprot]{metrics}})}
}

\details{

\tabular{ll}{
  \code{genoGOE_1} \tab Genome-wide differences of oxidation state between two lineages of methanogens \cr
  \code{genoGOE_2} \tab Carbon oxidation state of proteins in eukaryotic gene age groups \cr
  \code{genoGOE_3} \tab Evolutionary oxidation of ancestral Rubiscos and relative stabilities as a function of redox and pH \cr
  \code{genoGOE_4} \tab Evolutionary oxidation and relative stabilities for genomes with S-cycling genes \cr
  \code{genoGOE_5} \tab Carbon oxidation state of ancestral and extant nitrogenases \cr
}
}

\references{
Dick JM (2025) Manuscript in preparation.
}

\examples{
genoGOE_2()
}
