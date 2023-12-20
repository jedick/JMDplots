\encoding{UTF-8}
\name{microhum}
\alias{microhum}
\alias{microhum_1}
\alias{microhum_2}
\alias{microhum_3}
\alias{microhum_4}
\alias{microhum_S1}
\alias{microhum_S2}
\alias{dataset_metrics}
\title{Water and oxygen in the genomic adaptation of human microbiomes}
\description{
Plots from the paper by Dick (2023).
}

\usage{
  microhum_1(pdf = FALSE)
  microhum_2(pdf = FALSE)
  microhum_3(pdf = FALSE)
  microhum_4(pdf = FALSE)
  microhum_S1(pdf = FALSE)
  microhum_S2(pdf = FALSE)
  dataset_metrics()
}

\arguments{
  \item{pdf}{logical, make a PDF file?}
}

\details{
This table briefly describes each plotting function.

\tabular{ll}{
  \code{microhum_1} \tab Chemical variation of microbial proteins across body sites and multi-omics comparison \cr
  \code{microhum_2} \tab Differences of chemical metrics between controls and COVID-19/IBD patients \cr
  \code{microhum_3} \tab Overview of chemical variation of human microbiomes inferred from multi-omics datasets \cr
  \code{microhum_4} \tab Oxygen tolerance of genera in body sites, COVID-19, and IBD \cr
  \code{microhum_S1} \tab Changes of chemical metrics for microbiomes associated with viral inactivation \cr
  \code{microhum_S2} \tab Differences of \nO2 and \nH2O between obligate anaerobic and aerotolerant genera \cr
}

\code{dataset_metrics} is used to precompute mean values of chemical metrics for controls and COVID-19/IBD patients in each study (chemical metrics are for community reference proteomes; source data are 16S rRNA sequences).

The data files listed below are stored in \file{extdata/microhum}:

}

\section{Data and scripts for 16S rRNA datasets}{
  \describe{
    \item{\file{pipeline.R}}{Pipeline for sequence data processing (uses external programs fastq-dump, vsearch, seqtk, RDP Classifier, and GNU Parallel (Tange, 2023))}
    \item{\file{RDP-GTDB/*.tab.xz}}{
      RDP Classifier results combined into a single CSV file for each study, created with the \code{classify} and \code{mkRDP} functions in \file{pipeline.R}.
      GTDB 16S SSU rRNA training files for RDP Classifier are available at \url{https://doi.org/10.5281/zenodo.7633100}.
    }
    \item{\file{metadata/*.csv}}{Sample metadata for each study}
  }
}

\section{Data and scripts for metagenomes}{
  \describe{
    \item{\file{ARAST}}{Directory with scripts and output for metagenomes}
    \item{\file{ARAST/ARAST.R}}{Processing pipeline (modified from Dick et al. (2019) to implement removal of human sequences)}
    \item{\file{ARAST/arast_sortme_rna.pl}, \file{ARAST/sort_helper.sh}}{Helper scripts for ARAST.R}
    \item{\file{ARAST/runARAST.R}}{Script to run pipeline with specific settings for each dataset}
    \item{\file{ARAST/*_aa.csv}}{Output of runARAST.R with amino acid compositions of protein sequences inferred from metagenomes}
  }
}

\section{Data and scripts for metaproteomes}{
  \describe{
    \item{\file{metaproteome/*/mkaa.R}}{
      Scripts to process metaproteomic data (5 directories).
      See the comments in each \file{mkaa.R} for the required files that contain source data.
      Required files are available from databases or SI tables and are not included here.
    }
    \item{\file{metaproteome/*/*_aa.csv}}{Output of mkaa.R with amino acid composition of proteins}
  }
}

\section{Data and scripts for metagenome-assembled genomes (MAGs)}{
  \describe{
    \item{\file{KWL22}}{Directory with scripts and data for analyzing MAGs from Ke et al. (2022)}
    \item{\file{KWL22/COVID19_metadata.txt}}{Downloaded from \url{https://github.com/Owenke247/COVID-19/blob/main/Pre-processed_Files/COVID19_metadata.txt}}
    \item{\file{KWL22/mkaa.R}}{Script to get amino acid compositions for proteins predicted by Prodigal}
    \item{\file{KWL22/KWL22_MAGs_prodigal_aa.csv.xz}}{Output of mkaa.R with amino acid compositions of proteins for 5403 MAGs}
  }
}

\section{Other files}{
  \describe{
    \item{\file{MR18_Table_S1_modified.csv}}{List of strictly anaerobic and aerotolerant genera modified from Table S1 of Million and Raoult (2018)}
    \item{\file{dataset_metrics.csv}}{Mean values of chemical metrics for samples in 16S studies, created with \code{dataset_metrics}}
    \item{\file{images/body.png}, \file{images/testtube.png}}{Images used for \code{microhum_3}}
  }
}

\references{
Dick JM, Yu M, Tan J and Lu A (2019) Changes in carbon oxidation state of metagenomes along geochemical redox gradients. \emph{Front. Microbiol.} \bold{10}, 120. \doi{10.3389/fmicb.2019.00120}

Dick JM (2023) Two dimensions of chemical variation of the human microbiome across body sites and in COVID-19 patients. bioRxiv. \doi{10.1101/2023.02.12.528246}

Ke S, Weiss ST and Liu Y-Y (2022) Dissecting the role of the human microbiome in COVID-19 via metagenome-assembled genomes. \emph{Nat. Commun.} \bold{13}, 5253. \doi{10.1038/s41467-022-32991-w}

Million M and Raoult D (2018) Linking gut redox to human microbiome. \emph{Human Microbiome Journal} \bold{10}, 27--32. \doi{10.1016/j.humic.2018.07.002}

Tange O (2023) GNU Parallel 20230222 ('Gaziantep'). \doi{10.5281/zenodo.7668338}
}

\examples{
# Overview figure
microhum_3()
}
