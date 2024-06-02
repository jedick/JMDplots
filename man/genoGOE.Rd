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
  genoGOE_3(pdf = FALSE, x = "pH", y = "Eh", basis = "QEC")
}

\arguments{
  \item{pdf}{logical, make a PDF file?}
  \item{metric}{character, which metric to plot (see \code{\link[canprot]{metrics}})}
  \item{x}{character, x variable (\samp{pH} or \samp{O2})}
  \item{y}{character, y variable (\samp{Eh}, \samp{O2}, or \samp{H2O})}
  \item{basis}{character, basis species (\samp{QEC} or \samp{CHNOS})}
}

\details{

\tabular{ll}{
  \code{genoGOE_1} \tab Protein oxidation is pervasive in genomes of methanogens that diverged after the GOE \cr
  \code{genoGOE_2} \tab Gene age groups in eukaryotic genomes show protein oxidation around the GOE and NOE \cr
  \code{genoGOE_3} \tab Oxidation of reconstructed ancestral Rubisco sequences around the GOE corresponds to high-Eh conditions \cr
}
}

\references{
Dick J (2024) Manuscript in preparation.
}

\examples{
genoGOE_2()
}
