\encoding{UTF-8}
\name{JMDplots-package}
\alias{JMDplots-package}
\alias{JMDplots}
\docType{package}
\title{Plots from papers by Jeffrey M. Dick}

\description{
This package contains data and code used to make the plots in different papers by the author.
Each paper is identified with a project name, as listed in the table below.
The plots are available in the corresponding vignettes in the package.

Some papers are not included in this package.
Many plots from the papers by Dick and Shock (2011 and 2013) are shown in the vignette \dQuote{Hot-spring proteins in CHNOSZ} (\url{../../CHNOSZ/doc/hotspring.pdf}).
Most plots from the papers by Dick (2016 and 2017) are in the vignettes for \pkg{canprot} (\url{../../canprot/doc/index.html}).

Note that the code in this package is not identical to that used in the papers, but has been minimally modified to be included in this package.
Please refer to the supporting information of the papers for the original scripts.
}

\details{
  \itemize{
    \item \code{\link{gradH2O}} - Hydration state of metagenomes in salinity gradients (Dick et al., 2019b)
    \item \code{\link{chnszten}} - CHNOSZ, ten years after first CRAN submission (Dick, 2019)
    \item \code{\link{gradox}} - Carbon oxidation state of metagenomes in redox gradients (Dick et al., 2019a)
    \item \code{\link{aoscp}} - Average oxidation state of carbon in proteins (Dick, 2014)
    \item \code{aaaq} - Amino acid group additivity for ionized proteins (Dick et al., 2006)
  }
  These names are used for the vignettes and functions; the function names have figure numbers appended, as in \code{\link{gradox1}}.
  Note that \code{aaaq} has a vignette but no functions.
  Data for each of the papers are stored in the corresponding directories under \code{extdata}; see the documentation page for each paper for more details.
  There are some additional directories, described below, with other data.

  \describe{
    \item{\code{aaaq/BKM60_Fig7.csv}}{Eh-pH values for normal, wet and waterlogged soils from Fig. 7 of Baas Becking et al. (1960).}
    \item{\code{subsurface}}{
      Metagenome-derived amino acid compositions for two subsurface environments.
      These files are provided to illustrate the usage of user-provided data sets; see \code{\link{subsurface}}.
    }

    \item{\code{organisms/Sce.csv.xz}}{
      Data frame of amino acid composition of 6716 proteins from the \emph{Saccharomyces} Genome Database (SGD).
      Values in the first three columns are the \code{ORF} names of proteins, \code{SGDID}, and \code{GENE} names.
      The remaining twenty columns (\code{ALA}..\code{VAL}) contain the numbers of the respective amino acids in each protein.
      The sources of data for \samp{Sce.csv} are the files \samp{protein_properties.tab} and \samp{SGD_features.tab} (for the gene names), downloaded from \url{http://www.yeastgenome.org} on 2013-08-24.
      A shorter version of this file was previously present in \pkg{CHNOSZ} (to version 1.3.3).
      See \code{\link{yeast.aa}}
    }

    \item{\code{organisms/yeastgfp.csv.xz}}{
      Has 28 columns; the names of the first five are \code{yORF}, \code{gene name}, \code{GFP tagged?}, \code{GFP visualized?}, and \code{abundance}.
      The remaining columns correspond to the 23 subcellular localizations considered in the YeastGFP project (Huh et al., 2003 and Ghaemmaghami et al., 2003) and hold values of either \code{T} or \code{F} for each protein.
      \samp{yeastgfp.csv} was downloaded on 2007-02-01 from http://yeastgfp.ucsf.edu using the Advanced Search, setting options to download the entire dataset and to include localization table and abundance, sorted by orf number.
      See \code{\link{yeastgfp}}.
    }

    \item{\code{organisms/ecoli.rds}}{
      Amino acid compositions of proteins from \emph{E. coli}, used in \code{\link{gradH2O1}}.
      See \code{\link{organisms}}.
    }

    \item{\code{refseq}}{
      Total amino acid compositions of all proteins for each taxon (bacteria, archaea, and viral) in the NCBI Reference Sequence database, based on RefSeq release 95 of 2019-07-12 (O'Leary et al., 2016).
      This directory contains the scripts used to process the RefSeq files and resulting CSV file with amino acid compositions.
      \describe{
        \item{\code{protein_refseq.csv.xz}}{The columns are arranged as follows: \code{protein} (\dQuote{refseq}), \code{organism} (taxonomic id), \code{ref} (organism name), \code{abbrv} (empty), \code{chains} (number of protein sequences for this organism). The remaining columns have the total counts of amino acids.}
      }
    }
  }

}



\section{Requirements}{
Requires \pkg{CHNOSZ} 1.3.2 for functionality used in \code{\link{chnszten}}.
Requires \pkg{canprot} > 0.1.2 for updated plotting style in \code{\link[canprot]{diffplot}}.
}

\references{
Baas Becking, L. G. M., Kaplan, I. R. and Moore, D. (1960) Limits of the natural environment in terms of pH and oxidation-reduction potentials. \emph{Journal of Geology} \bold{68}, 243--284. \url{http://www.jstor.org/stable/30059218}

Dick, J. M., LaRowe, D. E. and Helgeson, H. C. (2006) Temperature, pressure, and electrochemical constraints on protein speciation: Group additivity calculation of the standard molal thermodynamic properties of ionized unfolded proteins. \emph{Biogeosciences} \bold{3}, 311--336. \doi{10.5194/bg-3-311-2006}

Dick, J. M. (2014) Average oxidation state of carbon in proteins. \emph{J. R. Soc. Interface} \bold{11}, 20131095. \doi{10.1098/rsif.2013.1095}

Dick, J. M. (2019) CHNOSZ: Thermodynamic calculations and diagrams for geochemistry. \emph{Front. Earth Sci.} \bold{7}:180. \doi{10.3389/feart.2019.00180}

Dick, J. M., Yu, M., Tan, J. and Lu, A. (2019a) Changes in carbon oxidation state of metagenomes along geochemical redox gradients. \emph{Front. Microbiol.} \bold{10}, 120. \doi{10.3389/fmicb.2019.00120}

Dick, J. M., Tan, J. and Yu, M. (2019b) Geochemistry of metagenomes: Distinct compositional signatures of environmental redox and salinity gradients. \emph{(preprint)}

O'Leary, N. A. et al. (2016) Reference sequence (RefSeq) database at NCBI: current status, taxonomic expansion, and functional annotation. \emph{Nucleic Acids Res.} \bold{44}, D733-D745. \doi{10.1093/nar/gkv1189}
}
