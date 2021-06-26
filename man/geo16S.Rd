\encoding{UTF-8}
\name{geo16S}
\alias{geo16S}
\alias{geo16S1}
\alias{geo16S2}
\alias{geo16S3}
\alias{geo16S4}
\alias{geo16S5}
\alias{geo16S6}
\alias{geo16S_S1}
\title{Plots from geobiochemistry of 16S sequences paper}
\description{
Plots from the paper by Dick and Tan (2021).
}

\usage{
  geo16S1(pdf = FALSE)
  geo16S2(pdf = FALSE)
  geo16S3(pdf = FALSE, plot.it = TRUE)
  geo16S4(pdf = FALSE)
  geo16S5(pdf = FALSE)
  geo16S6(pdf = FALSE)
  geo16S_S1(pdf = FALSE)
}

\arguments{
  \item{pdf}{logical, make a PDF file?}
  \item{plot.it}{logical, make a plot?}
}

\details{
This table gives a brief description of each plotting function.

\tabular{ll}{
  \code{geo16S1} \tab Distinct chemical parameters of predicted proteomes for major taxonomic groups \cr
  \code{geo16S2} \tab Inferred community proteomes from different environments have distinct chemical signatures \cr
  \code{geo16S3} \tab Lower carbon oxidation state is tied to oxygen depletion in water columns \cr
  \code{geo16S4} \tab Changes of abundance and chemical composition for individual taxonomic groups \cr
  \code{geo16S5} \tab Decreased carbon oxidation state of inferred proteomes for communities affected by unconventional oil and gas extraction \cr
  \code{geo16S6} \tab Contributions of classes to overall \ZC difference between oxidizing and reducing conditions \cr
  \code{geo16S_S1} \tab Abundance and \ZC of classes in oxidizing and reducing conditions \cr
}

For \code{geo16S3}, the \code{plot.it} argument can be used to prevent plotting, and instead return the \ZC and \nH2O values for all the points in the plot.
This is used by \code{code16S2} to make the convex hull for the datasets plotted by \code{geo16S3}.

}

\section{Files in extdata/geo16S}{
  \describe{
    \item{\code{geo16S_S1.csv}}{File created by running \code{geo16S_S1}, used by \code{geo16S6}.}
  }
}

\references{
Dick, J. M. and Tan, J. (2021) Genomic evidence for a chemical link between redox conditions and microbial community composition. bioRxiv. \doi{https://doi.org/10.1101/2021.05.31.446500}
}

\examples{
geo16S1()
}
