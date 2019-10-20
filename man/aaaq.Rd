\encoding{UTF-8}
\name{aaaq}
\alias{aaaq}
\alias{aaaq4}
\alias{aaaq13}
\title{Plots from amino acid group additivity paper}
\description{
Plots from the paper by Dick et al. (2006).
}

\usage{
  aaaq4(pdf = FALSE)
  aaaq13(pdf = FALSE)
}

\arguments{
  \item{pdf}{logical, make a PDF file?}
}

\details{
This table gives a brief description of each function.

\tabular{ll}{
  \code{aaaq4} \tab \pKa of sidechain and amino acid backbone groups \cr
  \code{aaaq13} \tab Eh-pH diagram for extracellular alpha-amylases \cr
}
}

\section{extdata/aaaq}{
  \describe{
    \item{\code{BKM60_Fig7.csv}}{Eh-pH values for normal, wet and waterlogged soils from Fig. 7 of Baas Becking et al. (1960).}
  }
}

\references{
Baas Becking, L. G. M., Kaplan, I. R. and Moore, D. (1960) Limits of the natural environment in terms of pH and oxidation-reduction potentials. \emph{Journal of Geology} \bold{68}, 243--284. \url{http://www.jstor.org/stable/30059218}

Dick, J. M., LaRowe, D. E. and Helgeson, H. C. (2006) Temperature, pressure, and electrochemical constraints on protein speciation: Group additivity calculation of the standard molal thermodynamic properties of ionized unfolded proteins. \emph{Biogeosciences} \bold{3}, 311--336. \doi{10.5194/bg-3-311-2006}
}

\examples{
aaaq4()
}
