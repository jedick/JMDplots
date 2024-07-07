\encoding{UTF-8}
\name{microhum}
\alias{microhum}
\alias{microhum_1}
\alias{microhum_2}
\alias{microhum_3}
\alias{microhum_4}
\alias{microhum_5}
\alias{microhum_6}
\alias{microhum_1_1}
\alias{microhum_3_1}
\alias{microhum_5_1}
\alias{microhum_6_1}
\alias{dataset_metrics}
\title{Adaptations of microbial genomes to human body chemistry}
\description{
Plots from the manuscript by Dick (2024).
}

\usage{
  microhum_1(pdf = FALSE)
  microhum_2(pdf = FALSE)
  microhum_3(pdf = FALSE)
  microhum_4(pdf = FALSE)
  microhum_5(pdf = FALSE)
  microhum_6(pdf = FALSE)
  microhum_1_1(pdf = FALSE)
  microhum_3_1(pdf = FALSE)
  microhum_5_1(pdf = FALSE)
  microhum_6_1(pdf = FALSE)
  dataset_metrics()
}

\arguments{
  \item{pdf}{logical, make a PDF file?}
}

\details{
This table briefly describes each plotting function.

\tabular{ll}{
  \code{microhum_1} \tab Consistency between shotgun metagenomes and community reference proteomes \cr
  \code{microhum_2} \tab Chemical metrics are broadly different among genera and are similar between GTDB and low-contamination genomes from UHGG \cr
  \code{microhum_3} \tab Chemical variation of microbial proteins across body sites (multi-omics comparison) \cr
  \code{microhum_4} \tab Differences of chemical metrics between controls and COVID-19 or IBD patients \cr
  \code{microhum_5} \tab Differences of relative abundances of genera between controls and patients \cr
  \code{microhum_6} \tab Oxygen tolerance of genera in body sites, COVID-19, and IBD \cr
  \code{microhum_1_1} \tab Amount of putative human DNA removed from HMP metagenomes in screening step \cr
  \code{microhum_3_1} \tab Differences of \nO2 and \nH2O between untreated and viral-inactivated samples \cr
  \code{microhum_5_1} \tab Chemical metrics of reference proteomes for genera with known oxygen tolerance \cr
  \code{microhum_5_1} \tab Differences of \nO2 and \nH2O between subcommunities of obligate anaerobes and aerotolerant genera in controls and patients \cr
}

\code{dataset_metrics} is used to precompute mean values of chemical metrics for controls and COVID-19/IBD patients in each study (chemical metrics are for community reference proteomes; source data are 16S rRNA sequences).

The data files listed below are stored in \file{extdata/microhum}:

}

\section{Data and scripts for 16S rRNA datasets}{
  \describe{
    \item{\file{16S/pipeline.R}}{Pipeline for sequence data processing (uses external programs fastq-dump, vsearch, seqtk, RDP Classifier, and GNU Parallel (Tange, 2023))}
    \item{\file{16S/RDP-GTDB/*.tab.xz}}{
      RDP Classifier results combined into a single CSV file for each study, created with the \code{classify} and \code{mkRDP} functions in \file{pipeline.R}.
      RDP Classifier training files for 16S rRNA sequences from GTDB are available at \url{https://doi.org/10.5281/zenodo.7633099}.
    }
    \item{\file{16S/metadata/*.csv}}{Sample metadata for each study}
    \item{\file{16S/dataset_metrics_*.csv}}{Mean values of chemical metrics for samples in 16S studies, created with \code{dataset_metrics}}
  }
}

\section{Data and scripts for metagenomes}{
  \describe{
    \item{\file{ARAST}}{Directory with scripts and output for metagenomes}
    \item{\file{ARAST/ARAST.R}}{Processing pipeline (modified from Dick et al. (2019) to implement screening for human sequences)}
    \item{\file{ARAST/arast_sortme_rna.pl}, \file{ARAST/sort_helper.sh}}{Helper scripts for ARAST.R}
    \item{\file{ARAST/runARAST.R}}{Script to run pipeline with specific settings for each dataset}
    \item{\file{ARAST/*_aa.csv}}{Summed amino acid compositions of protein sequences inferred from each metagenomic sequencing run (produced by runARAST.R)}
    \item{\file{ARAST/*_stats.csv}}{Metagenome sequence processing statistics (produced by runARAST.R)}
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
    \item{\file{KWL22/BioSample_metadata.txt}}{BioSample metadata obtained from NCBI BioProjects PRJNA624223 and PRJNA650244}
    \item{\file{KWL22/mkaa.R}}{Script to get amino acid compositions for proteins predicted by Prodigal}
    \item{\file{KWL22/KWL22_MAGs_prodigal_aa.csv.xz}}{Output of mkaa.R with amino acid compositions of proteins for 5403 MAGs}
  }
}

\section{Other files}{
  \describe{
    \item{\file{MR18_Table_S1_modified.csv}}{List of strictly anaerobic and aerotolerant genera modified from Table S1 of Million and Raoult (2018)}
    \item{\file{Figure_5_genera.txt}}{List of genera in Figure 5, created from the value invisibly returned by \code{microhum_5}}
  }
}

\references{
Dick JM, Yu M, Tan J and Lu A (2019) Changes in carbon oxidation state of metagenomes along geochemical redox gradients. \emph{Front. Microbiol.} \bold{10}, 120. \doi{10.3389/fmicb.2019.00120}

Dick JM (2024) Adaptations of microbial genomes to human body chemistry. bioRxiv. \doi{10.1101/2023.02.12.528246}

Ke S, Weiss ST and Liu Y-Y (2022) Dissecting the role of the human microbiome in COVID-19 via metagenome-assembled genomes. \emph{Nat. Commun.} \bold{13}, 5253. \doi{10.1038/s41467-022-32991-w}

Million M and Raoult D (2018) Linking gut redox to human microbiome. \emph{Human Microbiome Journal} \bold{10}, 27--32. \doi{10.1016/j.humic.2018.07.002}

Tange O (2023) GNU Parallel 20230222 ('Gaziantep'). \doi{10.5281/zenodo.7668338}
}

\examples{
# Figure 1
microhum_1()
}
