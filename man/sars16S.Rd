\encoding{UTF-8}
\name{sars16S}
\alias{sars16S}
\alias{sars16S_1}
\alias{sars16S_2}
\alias{sars16S_3}
\alias{sars16S_4}
\alias{COVID_means}
\title{Plots from human microbiome paper}
\description{
Plots from the paper by Dick (2023).
}

\usage{
  sars16S_1(pdf = FALSE)
  sars16S_2(pdf = FALSE)
  sars16S_3(pdf = FALSE)
  sars16S_4(pdf = FALSE)
  COVID_means()
}

\arguments{
  \item{pdf}{logical, make a PDF file?}
}

\details{
This table briefly describes each plotting function.

\tabular{ll}{
  \code{sars16S_1} \tab Chemical metrics of reference proteomes as a function of oxygen tolerance and comparison with metaproteomes \cr
  \code{sars16S_2} \tab Chemical variation of microbial proteins across body sites, after viral inactivation treatment, and multi-omics comparison \cr
  \code{sars16S_3} \tab Differences of chemical metrics of nasopharyngeal, oropharyngeal, and gut microbiomes between controls and COVID-19 patients \cr
  \code{sars16S_4} \tab Overview of chemical variation of the human microbiome inferred from multi-omics datasets \cr
}

\code{COVID_means} is used to precompute mean values of chemical metrics for controls and COVID-19 patients in each study (chemical metrics are for community reference proteomes; source data are 16S rRNA sequences).

}

\section{Data and scripts for 16S rRNA datasets}{
  \describe{
    \item{\file{pipeline.R}}{Pipeline for sequence data processing (uses external programs fastq-dump, vsearch, seqtk, RDP Classifier).}
    \item{\file{RDP-GTDB/*.csv.xz}}{
      RDP Classifier results combined into a single CSV file for each study, created with the \code{classify} and \code{mkRDP} functions in \file{pipeline.R}.
      GTDB 16S SSU rRNA training files for RDP Classifier are available at \url{https://doi.org/10.5281/zenodo.7633100}.
    }
    \item{\file{metadata/*.csv}}{Sample metadata for each study.}
  }
}

\section{Data and scripts for metagenomes}{
  These files are in \file{extdata/sars16S/ARAST}.
  \describe{
    \item{\file{ARAST}}{Directory with scripts and output for metagenomes}
    \item{\file{ARAST/ARAST.R}}{Processing pipeline (modified from Dick et al. (2019) to implement removal of human sequences)}
    \item{\file{ARAST/arast_sortme_rna.pl}, \file{ARAST/sort_helper.sh}}{Helper scripts for ARAST.R}
    \item{\file{ARAST/runARAST.R}}{Script to run pipeline with specific settings for each dataset}
    \item{\file{ARAST/*_AA.csv}}{Output of runARAST.R with amino acid composition of proteins (3 files)}
  }
}

\section{Data and scripts for metaproteomes}{
  These files are in \file{extdata/sars16S/metaproteome}.
  See the comments in each \file{mkaa.R} for the required files that contain source data.
  Required files are available from databases or SI tables and are not included here.
  \describe{
    \item{\file{metaproteome/*/mkaa.R}}{Scripts to process metaproteomic data (5 directories)}
    \item{\file{metaproteome/*/*_aa.csv}}{Output of mkaa.R with amino acid composition of proteins}
  }
}

\section{Data and scripts for metagenome-assembled genomes (MAGs)}{
  These files are in \file{extdata/sars16S/KWL22}.
  \describe{
    \item{\file{KWL22}}{Directory with scripts and data for analyzing MAGs from Ke et al. (2022)}
    \item{\file{KWL22/COVID19_metadata.txt}}{Downloaded from \url{https://github.com/Owenke247/COVID-19/blob/main/Pre-processed_Files/COVID19_metadata.txt}}
    \item{\file{KWL22/mkaa.R}}{Script to get amino acid compositions for proteins predicted by Prodigal}
    \item{\file{KWL22/KWL22_MAGs_prodigal_aa.csv.xz}}{Output of mkaa.R with amino acid compositions of proteins for 5403 MAGs}
  }
}

\section{Other files}{
  \describe{
    \item{\file{extdata/orp16S/MR18_Table_S1.csv}}{List of strictly anaerobic and aerotolerant genera from Table S1 of Million and Raoult (2018).}
    \item{\file{extdata/sars16S/COVID_means_GTDB.csv}}{Mean values of chemical metrics for samples in 16S studies, created with \code{COVID_means}.}
  }
}

\references{
Dick JM, Yu M, Tan J and Lu A (2019) Changes in carbon oxidation state of metagenomes along geochemical redox gradients. \emph{Front. Microbiol.} \bold{10}, 120. \doi{10.3389/fmicb.2019.00120}

Dick JM (2023) Two dimensions of chemical variation of the human microbiome across body sites and in COVID-19 patients. Manuscript in preparation.

Ke S, Weiss ST and Liu Y-Y (2022) Dissecting the role of the human microbiome in COVID-19 via metagenome-assembled genomes. \emph{Nat. Commun.} \bold{13}, 5253. \doi{10.1038/s41467-022-32991-w}

Million M and Raoult D (2018) Linking gut redox to human microbiome. \emph{Human Microbiome Journal} \bold{10}, 27--32. \doi{10.1016/j.humic.2018.07.002}
}

\examples{
# Summary figure
sars16S_4()
}
