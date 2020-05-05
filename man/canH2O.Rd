\encoding{UTF-8}
\name{canH2O}
\alias{canH2O}
\alias{canH2O1}
\alias{canH2O2}
\alias{canH2O3}
\alias{canH2O4}
\alias{canH2OT1}
\alias{canH2OT2}
\alias{canH2OT3}
\alias{canH2OS1}
\alias{canH2OS2}
\alias{canH2OS3}
\alias{canH2OS4}
\alias{canH2OS5}
\title{Plots from water in cancer paper}
\description{
Plots from the paper by Dick (2020).
Function names ending with \samp{T*} and \samp{S*} refer to supplementary tables and figures.
}

\usage{
  canH2O1(pdf = FALSE)
  canH2O2(pdf = FALSE)
  canH2O3(pdf = FALSE)
  canH2O4(pdf = FALSE)
  canH2OT1()
  canH2OT2()
  canH2OS1(pdf = FALSE)
  canH2OS2(pdf = FALSE)
  canH2OS3(pdf = FALSE)
  canH2OS4(pdf = FALSE)
  canH2OS5(pdf = FALSE)
}

\arguments{
  \item{pdf}{logical, make a PDF file?}
}

\details{
This table gives a brief description of each function.

\tabular{ll}{
  \code{canH2O1} \tab Study overview \cr
  \code{canH2O2} \tab Median differences of protein length, \nH2O and \ZC for cell culture and cancer tissue \cr
  \code{canH2O3} \tab \nH2O-\ZC and phylostrata plots for TCGA and HPA datasets \cr
  \code{canH2O4} \tab Compositional analysis of phylostrata and phylostrata-nAA plots for cancer tissue \cr
  \code{canH2OT1} \tab Stoichiometric matrix for amino acids with QEC basis species \cr
  \code{canH2OT2} \tab Mean differences and \emph{p}-values across all datasets \cr
  \code{canH2OS1} \tab rQEC derivation \cr
  \code{canH2OS2} \tab HPA-TCGA scatterplots for \ZC and \nH2O \cr
  \code{canH2OS3} \tab Trigos and Liebeskind phylostrata against \ZC and \nH2O for cancer tissue \cr
  \code{canH2OS4} \tab Scatterplots of \ZC and \nH2O vs phylostrata for HPA datasets  \cr
  \code{canH2OS5} \tab Scatterplots of hypoxia scores and \ZC or \nH2O for TCGA or HPA datasets \cr
}
}

\section{Files in extdata/canH2O}{
  \describe{
    \item{\code{BHL+19_Fig1a.csv}}{Median hypoxia scores for different cancer types, digitized from Figure 1a of Bhandari et al. (2019).}
  }
}


\references{
Bhandari, V. et al. (2019) Molecular landmarks of tumor hypoxia across cancer types. \emph{Nat. Genet.} \bold{51}, 308--318. \doi{10.1038/s41588-018-0318-2}

Dick, J. M. (2020) Water as a reactant in the differential expression of proteins in cancer. \emph{(preprint)}
}

\examples{
canH2O1()
}
