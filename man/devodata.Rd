\encoding{UTF-8}
\name{devodata}
\alias{devodata}
\alias{getCBS17}
\title{Data from developmental proteomics studies}
\description{
Data and chemical parameters for fruit fly and biofilm developmental proteomes.
}

\usage{
  getCBS17(X, boot.R)
}

\arguments{
  \item{X}{character, name of parameter (\samp{nH2O})}
  \item{boot.R}{integer, number of bootstrap replicates}
}

\details{
\code{getCBS17} calculates the weighted mean and bootstrap confidence interval for chemical parameters in the \emph{Drosophila melanogaster} developmental proteome (Casas-Vila et al., 2017).
}

\section{Files in extdata/devodata}{
  \describe{
    \item{\code{CBS+17_abundance.csv.xz}}{Processed protein abundance data from Supplemental Table S1 of Casas-Vila et al. (2017). Values in the columns for "imputed.log2.LFQ.intensity" were exponentiated, and means of all 4 replicates at each time point were taken to give the values in this file. Protein IDs were mapped to UniProt IDs in the "Entry" column.}
    \item{\code{CBS+17_aa.csv.xz}}{Amino acid compositions obtained from the canonical UniProt sequences for all proteins.}
    \item{\code{CBS+17.R}, \code{CBS+17_mean_aa.csv}}{R script and output for abundance-weighted mean amino acid compositions for proteins at each developmental time point.}
  }

}

\seealso{
\code{\link{evdevH2O4}} and \code{\link{evdevH2O5}} are functions that use these files.
}

\references{
Casas-Vilas, N. et al. (2017) The developmental proteome of \emph{Drosophila melanogaster}. \emph{Genome Res.} \bold{27}, 1273--1285. \doi{10.1101/gr.213694.116}

Futo, M. et al. (2021) Embryo-like features in developing \emph{Bacillus subtilis} biofilms. \emph{Mol. Biol. Evol.} \bold{38}, 31--47. \doi{10.1093/molbev/msaa217}
}

\examples{
getCBS17()
}
