\encoding{UTF-8}
\name{pdat_TCGA}
\alias{pdat_TCGA}
\title{Get Protein (from Gene) Expression Data}
\description{
  Calculate chemical composition for proteins corresponding to differentially expressed genes listed by GEPIA2 (Tang et al., 2019).
  Expression levels for normal tissue and cancer are from the Genotype-Tissue Expression project (GTEx Consortium, 2017) and The Cancer Genome Atlas (TCGA Research Network et al., 2013; Hutter and Zenklusen, 2018).
  This function is used in the \code{\link{canH2O}} paper.
}

\usage{
  pdat_TCGA(dataset = 2020, basis = "rQEC")
}

\arguments{
  \item{dataset}{character, dataset name}
  \item{basis}{character, keyword for basis species to use}
}

\details{
This function is used to calculate compositional metrics (including carbon oxidation state and stoichiometric hydration state) for proteins corresponding to differentially expressed genes.
See \code{\link[canprot]{pdat_}} for more information about the calculations and results returned by the \code{pdat_} functions.
}

\section{Files in extdata/expression/pancan}{
  \describe{
    \item{\code{GEPIA2.csv.xz}}{
      This file has gene names, Ensembl Gene IDs, UniProt IDs, and log2 fold change values for differentially expressed genes in each cancer type available on the GEPIA2 server (\url{http://gepia2.cancer-pku.cn/#degenes}) run using default settings (ANOVA, log2FC cutoff = 1, q-value cutoff = 0.01).
      Ensembl Gene IDs were converted to UniProt IDs using the UniProt mapping tool (\url{https://www.uniprot.org/mapping/}).
    }
  }
}

\references{
GTEx Consortium (2017) Genetic effects on gene expression across human tissues. \emph{Nature} \bold{550}, 204--213. \doi{10.1038/nature24277}

Hutter, C. and Zenklusen, J. C. (2018) The Cancer Genome Atlas: Creating lasting value beyond its data. \emph{Cell} \bold{173}, 283--285. \doi{10.1016/j.cell.2018.03.042}

Tang, Z., Kang, B., Li, C., Chen, T. and Zhang, Z. (2019) GEPIA2: an enhanced web server for large-scale expression profiling and interactive analysis. \emph{Nucleic Acids Res.} \bold{47}, W556--W560. \doi{10.1093/nar/gkz430}

The Cancer Genome Atlas Research Network et al. (2013) The Cancer Genome Atlas Pan-Cancer analysis project. \emph{Nat. Genetics} \bold{45}, 1113--1120. \doi{10.1038/ng.2764}
}

\seealso{
\code{\link{pdat_HPA}}, \code{\link{canH2O}}
}

\examples{
# list names of datasets for each cancer type
pdat_TCGA(2020)
# process one dataset (adrenocortical carcinoma)
pdat_TCGA("GEPIA2_ACC")
}
