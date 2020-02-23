\encoding{UTF-8}
\name{pdat_HPA}
\alias{pdat_HPA}
\title{Get Protein Expression Data}
\description{
  Calculate chemical composition for protein expression data from the Human Protein Atlas (Uhlen et al., 2015).
  This function is used in the \code{\link{canH2O}} paper.
}

\usage{
  pdat_HPA(dataset = 2020, basis = "rQEC")
}

\arguments{
  \item{dataset}{character, dataset name}
  \item{basis}{character, keyword for basis species to use}
}

\details{
This function is used to calculate compositional metrics (including carbon oxidation state and stoichiometric hydration state) for differentially expressed proteins identified in the Human Protein Atlas.
See \code{\link[canprot]{pdat_}} for more information about the calculations and results returned by the \code{pdat_} functions.
}

\section{Files in extdata/expression/pancan}{
  \describe{
    \item{\code{HPA.csv.xz}}{
      This file has gene names and Ensembl Gene IDs, UniProt IDs, and expression level values for cancer and normal tissue derived from the \code{pathology.tsv.zip} and \code{normal_tissue.tsv.zip} data files in the Human Protein Atlas, version 19 (\url{https://proteinatlas.org}).
      To calculate the expression level values, antibody staining intensities were converted to a numeric scale (not detected: 0, low: 1, medium: 3, high: 5), and values for all available samples were averaged (including all cell types for normal tissues).
      Ensembl Gene IDs were converted to UniProt IDs using the UniProt mapping tool (\url{https://www.uniprot.org/mapping/}).
    }
  }
}

\references{
Uhlen, M. et al. (2015) A pathology atlas of the human cancer transcriptome. \emph{Science} \bold{357}, eaan2507. \doi{10.1126/science.aan2507}
}

\seealso{
\code{\link{pdat_TCGA}}, \code{\link{canH2O}}
}

\examples{
# list names of datasets for cancer types
pdat_HPA(2020)
# process one dataset (breast cancer)
pdat_HPA("HPA19_1")
}
