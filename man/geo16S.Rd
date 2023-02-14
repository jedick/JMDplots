\encoding{UTF-8}
\name{geo16S}
\alias{geo16S}
\alias{geo16S1}
\alias{geo16S2}
\alias{geo16S3}
\alias{geo16S4}
\alias{geo16S5}
\alias{geo16S_S1}
\alias{geo16S_S2}
\alias{geo16S_S3}
\alias{geo16S_S4}
\alias{geo16S_S5}
\alias{geo16S_S6}
\alias{getmdat_geo16S}
\alias{getmetrics_geo16S}
\alias{plotmet_geo16S}
\title{Chemical links between redox conditions and community reference proteomes}
\description{
Plots from the paper by Dick and Tan (2022).
}

\usage{
  geo16S1(pdf = FALSE)
  geo16S2(pdf = FALSE)
  geo16S3(pdf = FALSE)
  geo16S4(pdf = FALSE)
  geo16S5(pdf = FALSE)
  geo16S_S1(pdf = FALSE)
  geo16S_S2(pdf = FALSE)
  geo16S_S3(pdf = FALSE)
  geo16S_S4(pdf = FALSE)
  geo16S_S5(pdf = FALSE, H2O = FALSE)
  geo16S_S6(pdf = FALSE)
  getmdat_geo16S(study, metrics = NULL, dropNA = TRUE)
  getmetrics_geo16S(study, quiet = TRUE, ...)
  plotmet_geo16S(study, quiet = TRUE, ...)
}

\arguments{
  \item{pdf}{logical, make a PDF file?}
  \item{H2O}{logical, make plots for \nH2O instead of \ZC?}
  \item{study}{character, study name}
  \item{metrics}{data frame, output of \code{\link[chem16S]{getmetrics}}}
  \item{dropNA}{logical, exclude samples with NA name in metadata?}
  \item{quiet}{logical, change to FALSE to print details about data processing}
  \item{...}{additional arguments passed to \code{\link[chem16S]{readRDP}} (for \code{getmetrics_geo16S}) or \code{\link[chem16S]{plotmet}} (for \code{plotmet_geo16S})}
}

\details{
This table gives a brief description of each plotting function.

\tabular{ll}{
  \code{geo16S1} \tab Distinct chemical parameters of reference proteomes for major taxonomic groups \cr
  \code{geo16S2} \tab Estimated community proteomes from different environments have distinct chemical signatures \cr
  \code{geo16S3} \tab Lower carbon oxidation state is tied to oxygen depletion in water columns \cr
  \code{geo16S4} \tab Common trends of carbon oxidation state of community reference proteomes for shale gas wells and hydrothermal systems \cr
  \code{geo16S5} \tab Comparison of protein \ZC from metagenomic or metatranscriptomic data with estimates from 16S and reference sequences \cr
  \code{geo16S_S1} \tab RefSeq and 16S rRNA data processing outline \cr
  \code{geo16S_S2} \tab Scatterplots of \ZC and \nH2O for genera vs higher taxonomic levels \cr
  \code{geo16S_S3} \tab \nH2O-\ZC plots for major phyla and their genera \cr
  \code{geo16S_S4} \tab Venn diagrams for phylum and genus names in the RefSeq (NCBI), RDP, and SILVA taxonomies \cr
  \code{geo16S_S5} \tab Correlations between \ZC estimated from metagenomes and 16S rRNA sequences \cr
  \code{geo16S_S6} \tab Correlation of \ZC with GC content of metagenomic and 16S amplicon reads \cr
}

\code{getmdat_geo16S} gets metadata for the indicated \code{study} and adds columns for plot parameters (\code{pch}, \code{col}).
For some datasets, sample subsets are indicated by appending a suffix to the study name separated by an underscore.
The default for \code{dropNA} means to exclude samples with NA name in the metadata file, which is used to exclude outliers from the analysis.
If \code{metrics} is supplied, the samples are sorted in the same order as the metadata file, and the function returns a list with both \samp{metadata} and the sorted \samp{metrics}.

\code{getmetrics_geo16S} calculates chemical metrics (\ZC and \nH2O) for the indicated study.
\code{...} is used to supply arguments to \code{\link[chem16S]{readRDP}}.

}

\section{Files in extdata/geo16S}{
  \describe{
    \item{\file{pipeline.R}}{Pipeline for sequence data processing (uses external programs fastq-dump, vsearch, seqtk, RDP Classifier).}
    \item{\file{metadata/*.csv}}{Sample metadata for each study.}
    \item{\file{RDP/*.csv.xz}}{RDP Classifier results combined into a single CSV file for each study, created with the \code{mkRDP} function in \file{pipeline.R}.}
    \item{\file{AWDM15.csv}}{Accession numbers and sample descriptions for the Human Microbiome Project, Soils, and Mammalian Guts datasets taken from Supplementary Material of Aßhauer et al. (2015) and used in \code{geo16S5}. GC content for metagenome and 16S amplicon reads was taken from \code{https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR*******}.}
  }
}

\section{Files in extdata/geo16S/ARAST}{
  \describe{
    \item{\code{*_AA.csv}}{Sum of amino acid compositions of protein fragments inferred from metagenomes, produced by \code{runARAST.R} and used in \code{geo16S5}. The \code{chains} column has the number of protein fragments.}
    \item{\code{scripts/}}{Directory with script files \code{ARAST.R} for the analysis pipeline, \code{runARAST.R} to run the pipeline for each dataset, and supporting Bash and Perl script files. The \code{ARAST.R} file is a modified version of the pipeline used in a previous study (see \code{\link{gradox}}).}
  }
}

\section{Files in extdata/geo16S/taxonomy}{
  \describe{
    \item{\code{process.R}}{Script used for extracting taxonomic names from RDP and SILVA sequence files. See comments for details about data sources.}
    \item{\code{[RDP|SILVA][phyla|genera].txt}}{Files with archaeal and bacterial phylum and genus names in RDP and SILVA, used in \code{geo16S_S4}.}
  }
}

\references{
Aßhauer KP, Wemheuer B, Daniel R and Meinicke P (2015) Tax4Fun: predicting functional profiles from metagenomic 16S rRNA data. \emph{Bioinformatics} \bold{31}, 2882--2884. \doi{10.1093/bioinformatics/btv287}

Dick JM and Tan J (2022) Chemical links between redox conditions and estimated community proteomes from 16S rRNA and reference protein sequences. \emph{Microb. Ecol.} \doi{10.1007/s00248-022-01988-9}
}

\examples{
# Make Figure 1
geo16S1()

# Get metrics of community reference proteomes and paired metadata for one study
metrics <- getmetrics_geo16S("BCA+21")
mdat <- getmdat_geo16S("BCA+21", metrics = metrics)

# Get *all* available metadata
metadata <- getmdat_geo16S("BCA+21", dropNA = FALSE)
stopifnot(nrow(metadata) > nrow(mdat$metadata))

# Make a nH2O-ZC plot with lots of messages printed
# Symbols are coded in getmetrics_geo16S
# (blue: oxic, black: suboxic, red: euxinic)
plotmet_geo16S("SVH+19", quiet = FALSE)

# List datasets used in geo16S paper
mdatdir <- system.file("extdata/geo16S/metadata", package = "JMDplots")
gsub(".csv", "", dir(mdatdir))
# Get metadata for one study
mdat <- getmdat_geo16S("BGPF13")
}
