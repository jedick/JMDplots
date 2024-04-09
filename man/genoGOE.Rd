\encoding{UTF-8}
\name{genoGOE}
\alias{genoGOE}
\alias{genoGOE_1}
\alias{genoGOE_2}
\title{Genomic record of Earth's surface oxygenation}
\description{
Plots from the manuscript by Dick (2024).
}

\usage{
  genoGOE_1(pdf = FALSE, metric = "Zc")
  genoGOE_2(pdf = FALSE, x = "pH", y = "Eh", basis = "QEC")
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
  \code{genoGOE_1} \tab Carbon oxidation state of proteins as a function of gene age in two lineages \cr
  \code{genoGOE_2} \tab Chemical analysis and thermodynamic calculations for ancenstral Rubiscos \cr
}
}

\references{
Dick J (2024) Manuscript in preparation.
}

\examples{
genoGOE_1()
}
