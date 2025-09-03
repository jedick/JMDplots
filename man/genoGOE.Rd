\encoding{UTF-8}
\name{genoGOE}
\alias{genoGOE}
\alias{genoGOE_1}
\alias{genoGOE_2}
\alias{genoGOE_3}
\alias{genoGOE_4}
\alias{genoGOE_5}
\title{Evolutionary oxidation of proteins in Earth's history}
\description{
Plots from the manuscript by Dick et al. (2025).
}

\usage{
  genoGOE_1(pdf = FALSE)
  genoGOE_2(pdf = FALSE, panel = NULL)
  genoGOE_3(pdf = FALSE, metric = "Zc")
  genoGOE_4(pdf = FALSE)
  genoGOE_5(pdf = FALSE, panel = NULL)
}

\arguments{
  \item{pdf}{logical, make a PDF file?}
  \item{panel}{character, create only this panel (disables figure layout, PDF output, and panel labeling)}
  \item{metric}{character, which metric to plot (see \code{\link[canprot]{metrics}})}
}

\details{

\tabular{ll}{
  \code{genoGOE_1} \tab Ranges of carbon oxidation state for organic compounds, amino acids, and proteins \cr
  \code{genoGOE_2} \tab Genome-wide differences of oxidation state between two lineages of methanogens \cr
  \code{genoGOE_3} \tab Carbon oxidation state of proteins in eukaryotic gene age groups \cr
  \code{genoGOE_4} \tab Carbon oxidation state of reconstructed ancestral sequences and extant proteins \cr
  \code{genoGOE_5} \tab From carbon oxidation state to relative stability diagrams \cr
}
}

\references{
Dick JM et al. (2025) Manuscript in preparation.
}

\examples{
genoGOE_3()
}
