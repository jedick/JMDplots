\encoding{UTF-8}
\name{pdat_CH16}
\alias{pdat_CH16}
\title{Get Protein (from Gene) Expression Data}
\description{
  Calculate chemical composition for proteins corresponding to 302 genes with common expression changes across cancer types, defined as genes with unidirectional expression changes in at least 13 of 32 cancer data sets (Chen and He, 2016).
}

\usage{
  pdat_CH16(dataset = 2020, basis = "rQEC")
}

\arguments{
  \item{dataset}{character, dataset name}
  \item{basis}{character, keyword for basis species to use}
}

\details{
See \code{\link[canprot]{pdat_}} for more information about the calculations and results returned by the \code{pdat_} functions.
}

\section{Files in extdata/expression/pancan}{
  \describe{
    \item{\code{CH16.csv}}{
      Ensembl Gene IDs, UniProt IDs, and regulation in cancer, from Supplementary Table S3 of Chen and He (2016).
      Ensembl Gene IDs were converted to UniProt IDs using the UniProt mapping tool (\url{https://www.uniprot.org/mapping/}).
    }
  }
}

\references{
Chen, H. and He, X. (2016) The convergent cancer evolution toward a single cellular destination. \emph{Mol. Biol. Evol.} \bold{33}, 4--12. \doi{10.1093/molbev/msv212}
}

\seealso{
This function is used in \code{\link{canH2O2}}.
}

\examples{
# show the median values and differences of
# ZC and nH2O between groups of proteins
# corresponding to up- and down-regulated genes
(get_comptab(pdat_CH16("CH16")))
}
