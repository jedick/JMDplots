\encoding{UTF-8}
\name{aoscp}
\alias{aoscp}
\alias{aoscp1}
\alias{aoscp2}
\alias{aoscp3}
\alias{yeast.ZC}
\title{Plots from oxidation state of carbon paper}
\description{
Plots from the paper by Dick (2014).
}

\usage{
  aoscp1(pdf = FALSE)
  aoscp2(pdf = FALSE)
  aoscp3(png = FALSE, outline = FALSE)
  yeast.ZC(location)
}

\arguments{
  \item{pdf}{logical, make a PDF file?}
  \item{png}{logical, make a PNG file?}
  \item{outline}{logical, draw only the outlines?}
  \item{location}{character, one-letter code for subcellular location from the \emph{Saccharomyces} Genome Database}
}

\details{
This table gives a brief description of each function.

\tabular{ll}{
  \code{aoscp1} \tab \ZC of amino acids vs \ZC of RNA codons and hydropathy index of amino acids \cr
  \code{aoscp2} \tab Histograms of \ZC of all human proteins and human membrane proteins \cr
  \code{aoscp3} \tab Drawing of yeast cell color-coded with median ZC of proteins in different locations \cr
}
}

\references{
Dick, J. M. (2014) Average oxidation state of carbon in proteins. \emph{J. R. Soc. Interface} \bold{11}, 20131095. \doi{10.1098/rsif.2013.1095}
}

\examples{
aoscp1()
}
