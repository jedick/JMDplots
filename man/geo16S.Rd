\encoding{UTF-8}
\name{geo16S}
\alias{geo16S}
\alias{geo16S1}
\alias{geo16S2}
\alias{geo16S3}
\alias{geo16S4}
\alias{geo16S5}
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
}

\arguments{
  \item{pdf}{logical, make a PDF file?}
  \item{plot.it}{logical, make a plot?}
}

\details{
This table gives a brief description of each plotting function.

\tabular{ll}{
  \code{geo16S1} \tab Chemical compositions of taxonomic groups \cr
  \code{geo16S2} \tab Natural environment datasets \cr
  \code{geo16S3} \tab Stratified lakes and seawater \cr
  \code{geo16S4} \tab Compositional differences at different taxonomic levels \cr
  \code{geo16S5} \tab Shale gas datasets \cr
}

For \code{geo16S3}, the \code{plot.it} argument can be used to prevent plotting, and instead return the \ZC and \nH2O values for all the points in the plot.
This is used by \code{code16S2} to make the convex hull for the datasets plotted by \code{geo16S3}.

}

\examples{
geo16S1()
}
