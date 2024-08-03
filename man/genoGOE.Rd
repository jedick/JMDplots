\encoding{UTF-8}
\name{genoGOE}
\alias{genoGOE}
\alias{genoGOE_1}
\alias{genoGOE_2}
\alias{genoGOE_3}
\title{Genomes record the Great Oxidation Event}
\description{
Plots from the manuscript by Dick (2024).
}

\usage{
  genoGOE_1(pdf = FALSE)
  genoGOE_2(pdf = FALSE, metric = "Zc")
  genoGOE_3(pdf = FALSE)
}

\arguments{
  \item{pdf}{logical, make a PDF file?}
  \item{metric}{character, which metric to plot (see \code{\link[canprot]{metrics}})}
}

\details{

\tabular{ll}{
  \code{genoGOE_1} \tab Genome-wide differences of oxidation state between two lineages of methanogens \cr
  \code{genoGOE_2} \tab Carbon oxidation state of proteins in eukaryotic gene age groups \cr
  \code{genoGOE_3} \tab Evolutionary oxidation of ancestral Rubiscos and thermodynamic prediction of redox boundaries around the GOE \cr
}
}

\references{
Dick J (2024) Manuscript in preparation.
}

\examples{
genoGOE_2()
}
