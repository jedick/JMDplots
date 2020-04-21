\encoding{UTF-8}
\name{pdat_}
\alias{pdat_TCGA}
\alias{pdat_HPA}
\alias{pdat_saltygenes}
\title{Get Protein Expression Data}
\description{
  These functions are used to calculate compositional metrics, including carbon oxidation state and stoichiometric hydration state, for differentially expressed proteins or proteins corresponding to differentially expressed genes.
  See the help page for \code{\link[canprot]{pdat_}} in the \pkg{canprot} package for more information about the calculations and results returned by the \code{pdat_} functions.
  These functions are used in the \code{\link{canH2O}} paper.
}

\usage{
  pdat_TCGA(dataset = 2020, basis = "rQEC")
  pdat_HPA(dataset = 2020, basis = "rQEC")
  pdat_saltygenes(dataset = 2020, basis = "rQEC")
}

\arguments{
  \item{dataset}{character, dataset name}
  \item{basis}{character, keyword for basis species to use}
}

\details{
\code{pdat_TCGA} calculates chemical composition for proteins corresponding to differentially expressed genes listed by GEPIA2 (Tang et al., 2019).
Expression levels for normal tissue and cancer are from the Genotype-Tissue Expression project (GTEx Consortium, 2017) and The Cancer Genome Atlas (TCGA Research Network et al., 2013; Hutter and Zenklusen, 2018).

\code{pdat_HPA} calculates chemical composition for protein expression data from the Human Protein Atlas (Uhlen et al., 2015).

\code{pdat_saltygenes} gets data for differentially expressed genes in transcriptomic experiments for bacteria.
}

\section{Files in extdata/expression/pancan}{
  \describe{
    \item{\code{GEPIA2.csv.xz}}{
      This file has gene names, Ensembl Gene IDs, UniProt IDs, and log2 fold change values for differentially expressed genes in each cancer type available on the GEPIA2 server (\url{http://gepia2.cancer-pku.cn/#degenes}) run using default settings (ANOVA, log2FC cutoff = 1, q-value cutoff = 0.01).
      Ensembl Gene IDs were converted to UniProt IDs using the UniProt mapping tool (\url{https://www.uniprot.org/mapping/}).
    }

    \item{\code{HPA.csv.xz}}{
      This file has gene names and Ensembl Gene IDs, UniProt IDs, and expression level values for cancer and normal tissue derived from the \code{pathology.tsv.zip} and \code{normal_tissue.tsv.zip} data files in the Human Protein Atlas, version 19 (\url{https://proteinatlas.org}).
      To calculate the expression level values, antibody staining intensities were converted to a numeric scale (not detected: 0, low: 1, medium: 3, high: 5), and values for all available samples were averaged (including all cell types for normal tissues).
      Ensembl Gene IDs were converted to UniProt IDs using the UniProt mapping tool (\url{https://www.uniprot.org/mapping/}).
    }

  }
}

\references{
GTEx Consortium (2017) Genetic effects on gene expression across human tissues. \emph{Nature} \bold{550}, 204--213. \doi{10.1038/nature24277}

Hutter, C. and Zenklusen, J. C. (2018) The Cancer Genome Atlas: Creating lasting value beyond its data. \emph{Cell} \bold{173}, 283--285. \doi{10.1016/j.cell.2018.03.042}

Tang, Z., Kang, B., Li, C., Chen, T. and Zhang, Z. (2019) GEPIA2: an enhanced web server for large-scale expression profiling and interactive analysis. \emph{Nucleic Acids Res.} \bold{47}, W556--W560. \doi{10.1093/nar/gkz430}

The Cancer Genome Atlas Research Network et al. (2013) The Cancer Genome Atlas Pan-Cancer analysis project. \emph{Nat. Genetics} \bold{45}, 1113--1120. \doi{10.1038/ng.2764}

Uhlen, M. et al. (2015) A pathology atlas of the human cancer transcriptome. \emph{Science} \bold{357}, eaan2507. \doi{10.1126/science.aan2507}
}

\seealso{
\code{\link{canH2O}}
}

\examples{
# list datasets for cancer types in TCGA
pdat_TCGA(2020)
# process one dataset (adrenocortical carcinoma)
pdat_TCGA("GEPIA2_ACC")
# show the median values and differences of
# ZC and nH2O between groups of proteins
# corresponding to up- and down-regulated genes
(get_comptab(pdat_TCGA("GEPIA2_ACC")))

# list datasets for cancer types in HPA
pdat_HPA(2020)
# process one dataset (breast cancer)
pdat_HPA("HPA19_1")
}
