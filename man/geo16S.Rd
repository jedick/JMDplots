\encoding{UTF-8}
\name{geo16S}
\alias{geo16S}
\alias{geo16S1}
\alias{geo16S2}
\alias{geo16S3}
\alias{geo16S4}
\alias{geo16S5}
\alias{geo16S6}
\alias{geo16S_S1}
\alias{geo16S_S2}
\alias{geo16S_S3}
\alias{geo16S_S4}
\alias{geo16S_S5}
\title{Plots from geobiochemistry of 16S sequences paper}
\description{
Plots from the paper by Dick and Tan (2021).
}

\usage{
  geo16S1(pdf = FALSE)
  geo16S2(pdf = FALSE)
  geo16S3(pdf = FALSE)
  geo16S4(pdf = FALSE)
  geo16S5(pdf = FALSE)
  geo16S6(pdf = FALSE)
  geo16S_S1(pdf = FALSE)
  geo16S_S2(pdf = FALSE)
  geo16S_S3(pdf = FALSE)
  geo16S_S4(pdf = FALSE)
  geo16S_S5(pdf = FALSE)
}

\arguments{
  \item{pdf}{logical, make a PDF file?}
}

\details{
This table gives a brief description of each plotting function.

\tabular{ll}{
  \code{geo16S1} \tab Distinct chemical parameters of reference proteomes for major taxonomic groups \cr
  \code{geo16S2} \tab Estimated community proteomes from different environments have distinct chemical signatures \cr
  \code{geo16S3} \tab Lower carbon oxidation state is tied to oxygen depletion in water columns \cr
  \code{geo16S4} \tab Decreased carbon oxidation state of estimated proteomes for communities affected by shale gas extraction \cr
  \code{geo16S5} \tab Contributions of classes to overall \ZC difference between oxidizing and reducing conditions \cr
  \code{geo16S6} \tab Comparison of protein \ZC from metagenomic or metatranscriptomic data with estimates from 16S and reference sequences \cr
  \code{geo16S_S1} \tab RefSeq and 16S rRNA data processing outline \cr
  \code{geo16S_S2} \tab Scatterplots of \ZC and \nH2O for genera vs higher taxonomic levels \cr
  \code{geo16S_S3} \tab Changes of abundance and chemical composition at different taxonomic levels \cr
  \code{geo16S_S4} \tab Abundance and \ZC of classes in oxidizing and reducing conditions \cr
  \code{geo16S_S5} \tab Correlations between \ZC estimated from metagenomes and 16S rRNA sequences \cr
}

}

\section{Files in extdata/geo16S}{
  \describe{
    \item{\code{geo16S_S4.csv}}{File created by running \code{write.csv(geo16S_S4(), "geo16S_S4.csv", row.names = FALSE, quote = FALSE)}, used by \code{geo16S5}.}
    \item{\code{metadata/*.csv}}{Sample metadata for each study.}
    \item{\code{RDP/*.csv.xz}}{RDP Classifier results combined into a single CSV file for each study, created with the \code{mkRDP} function in the \code{process.R} script (\code{\link{chem16S}}).}
    \item{\code{AWDM15.csv}}{Accession numbers and sample descriptions for the Human Microbiome Project, Soils, and Mammalian Guts datasets taken from Supplementary Material of Aßhauer et al. (2015) and used in \code{geo16S6}.}
  }
}

\section{Files in extdata/geo16S/ARAST}{
  \describe{
    \item{\code{*_AA.csv}}{Mean amino acid compositions of protein fragments inferred from metagenomes, produced by \code{runARAST.R} and used in \code{geo16S6}.}
    \item{\code{scripts/}}{Directory with script files \code{ARAST.R} for the analysis pipeline, \code{runARAST.R} to run the pipeline for each dataset, and supporting Bash and Perl script files. The \code{ARAST.R} file is a modified version of the pipeline used in a previous study (see \code{\link{gradox}}).}
  }
}

\references{
Aßhauer, K. P., Wemheuer, B., Daniel, R. and Meinicke, P. (2015) Tax4Fun: predicting functional profiles from metagenomic 16S rRNA data. \emph{Bioinformatics} \bold{31}, 2882--2884. \doi{10.1093/bioinformatics/btv287}

Dick, J. M. and Tan, J. (2021) Genomic evidence for a chemical link between redox conditions and microbial community composition. bioRxiv. \doi{https://doi.org/10.1101/2021.05.31.446500}
}

\examples{
geo16S1()
}
