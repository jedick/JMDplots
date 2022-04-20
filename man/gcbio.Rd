\encoding{UTF-8}
\name{gcbio}
\alias{gcbio}
\alias{gcbio1}
\alias{gcbio2}
\alias{gcbio3}
\alias{gcbio4}
\alias{gcbio5}
\alias{gcbioS1}
\alias{calc_logaH2_intermediate}
\title{Plots from frontier of geochemical biology paper}
\description{
Plots from the paper by Dick et al. (2022).
}

\usage{
  gcbio1(pdf = FALSE)
  gcbio2(pdf = FALSE)
  gcbio3(pdf = FALSE)
  gcbio4(pdf = FALSE)
  gcbio5(pdf = FALSE)
  gcbioS1(pdf = FALSE)
  calc_logaH2_intermediate(class = NULL)
}

\arguments{
  \item{pdf}{logical, make a PDF file?}
  \item{class}{character, organic compound class}
}

\details{
This table gives a brief description of each plotting function.

\tabular{ll}{
  \code{gcbio1} \tab \ZC of proteins and lipids in hot springs \cr
  \code{gcbio2} \tab Overlay \ZC on phylogenetic tree of methanogens \cr
  \code{gcbio3} \tab Affinities of organic synthesis reactions \cr
  \code{gcbio4} \tab Affinities for methanogen proteomes and \logaH2 of habitable niches \cr
  \code{gcbio5} \tab Interdisciplinary fields in a biological hierarchy \cr
  \code{gcbioS1} \tab Intermediate \logaH2 calculated for different compound classes \cr
}

\code{calc_logaH2_intermediate} is used to find intermediate logaH2, where affinity vs ZC has a slope of 0, for a group of organic compounds.
\code{class} gives the class of compounds to consider (see description of \file{H2_intermediate_*.csv} files below), or NULL to consider compounds in all classes at once.

}

\section{Files in extdata/gcbio}{
  \describe{
    \item{\file{Topt.csv}}{Optimal growth temperatures (\degC) of methanogens compiled from the literature.}
    \item{\file{methanogen_AA.csv}}{Amino acid compositions of reference proteomes of methanogens. See \file{R/gcbio.R} for the code used to make this file.}
    \item{\file{LG88_Fig1.csv}}{H2 concentrations measured in methanogenic and non-methanogenic sediments, from Lovley and Goodwin (1988).}
    \item{\file{H2_intermediate.csv}}{Intermediate logaH2 (where the slope of affinity vs \ZC is zero) for organic compounds, calculated using \code{calc_logaH2_intermediate}.}
    \item{\file{H2_intermediate_*.csv}}{Intermediate logaH2 (where the slope of affinity vs \ZC is zero) for different classes of organic compounds (\samp{C1andC2}, \samp{Acid}, \samp{Aminoacid}, \samp{Sugar}, \samp{Nucleobase}, \samp{TCA}).}
  }
}

\references{
Dick, J. M., Boyer, G. M., Canovas, P. A. III and Shock, E. L. (2022) The frontier of geochemical biology: Ecological applications of carbon oxidation state and the next steps for thermodynamic models. Manuscript in preparation.

Lovley, D. R. and Goodwin, S. (1988) Hydrogen concentrations as an indicator of the predominant terminal electron-accepting reactions in aquatic sediments. \emph{Geochim. Cosmochim. Acta} \bold{52}, 2993--3003. \doi{10.1016/0016-7037(88)90163-9}
}

\examples{
gcbio5()
}
