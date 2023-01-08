\encoding{UTF-8}
\name{JMDplots-package}
\alias{JMDplots-package}
\alias{JMDplots}
\docType{package}
\title{Plots from papers by Jeffrey M. Dick}

\description{
This package contains data and code used to make the plots in various papers.
Each paper is identified with a project name, as listed in the table below.
The plots are available in the corresponding vignettes in the package.

Note that the code in this package is not identical to that used in the papers, but has been minimally modified to be included in this package.
Please refer to the supporting information of the papers for the original scripts.
}

\details{
  \itemize{
    \item \code{\link{utogig}} - Using thermodynamics to obtain geochemical information from genomes (Dick et al., 2022)
    \item \code{\link{orp16S}} - Geochemical influence on protein evolution (Dick and Meng, 2021)
    \item \code{\link{geo16S}} - Chemical links between redox conditions and estimated community proteomes (Dick and Tan, 2022)
    \item \code{\link{evdevH2O}} - Thermodynamic model for water activity and redox potential in evolution and development (Dick, 2022)
    \item \code{\link{mjenergy}} - Energy release in protein synthesis (Dick and Shock, 2021)
    \item \code{\link{canH2O}} - Water as a reactant in the differential expression of proteins in cancer (Dick, 2021)
    \item \code{\link{gradH2O}} - Stoichiometric hydration state of metagenomes in salinity gradients (Dick et al., 2020)
    \item \code{\link{chnosz10}} - CHNOSZ, ten years after first CRAN submission (Dick, 2019)
    \item \code{\link{gradox}} - Carbon oxidation state of metagenomes in redox gradients (Dick et al., 2019)
    \item \code{\link{cpcp}} - Potential diagrams for cancer proteomes (Dick, 2016, 2017)
    \item \code{\link{aoscp}} - Average oxidation state of carbon in proteins (Dick, 2014)
    \item \code{\link{bison}} - Bison Pool hot spring (Dick and Shock, 2011, 2014)
    \item \code{\link{scsc}} - Subcellular locations of \emph{Saccharomyces cerevisiae} (Dick, 2009)
    \item \code{\link{aaaq}} - Amino acid group additivity for ionized proteins (Dick et al., 2006)
  }
  These names are used for the vignettes and functions; the function names have figure numbers appended, as in \code{\link{gradox1}}.
  Data for each of the papers are stored in the corresponding directories under \code{extdata}; see the documentation page for each paper for more details.
  There are some other directories, described below.

  \describe{

    \item{\code{subsurface}}{
      Metagenome-derived amino acid compositions for two subsurface environments.
      These files are provided to illustrate the usage of user-provided data sets; see \code{\link{subsurface}}.
    }

    \item{\code{organisms}}{

      \describe{

        \item{\code{Sce.csv.xz}}{
          Data frame of amino acid composition of 6716 proteins from the \emph{Saccharomyces} Genome Database (SGD).
          Values in the first three columns are the \code{ORF} names of proteins, \code{SGDID}, and \code{GENE} names.
          The remaining twenty columns (\code{ALA}..\code{VAL}) contain the numbers of the respective amino acids in each protein.
          The sources of data for \samp{Sce.csv} are the files \samp{protein_properties.tab} and \samp{SGD_features.tab} (for the gene names), downloaded from \url{https://www.yeastgenome.org/} on 2013-08-24.
          A shorter version of this file was previously present in \pkg{CHNOSZ} (to version 1.3.3).
          Used in \code{\link{yeast.aa}}.
        }

        \item{\code{UP000000805_243232.csv.xz}}{
          Amino acid compositions of 1787 proteins in \emph{Methanocaldococcus jannaschii}.
          It was created by processing the file \code{UP000000805_243232.fasta.gz}, which was downloaded from the UniProt reference proteomes FTP site (\url{ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Archaea/}).
          The server timestamp on \code{UP000000805_243232.fasta.gz} was Dec. 2, 2020; the reference proteome was last modified on August 22, 2020 according to \url{https://www.uniprot.org/proteomes/UP000000805}.
          The commands used to create \code{UP000000805_243232.csv.gz} are:
          \verb{
            aa <- CHNOSZ::read.fasta("UP000000805_243232.fasta.gz")
            write.csv(aa, "UP000000805_243232.csv", row.names = FALSE, quote = FALSE)
            system("xz UP000000805_243232.csv")
          }
          Used in \code{\link{mjenergy}}.
        }

        \item{\code{UP000000625_83333.csv.xz}}{
          This data file has amino acid compositions of 4392 proteins in the UniProt reference proteome of \emph{Escherichia coli} K12 (\url{https://www.uniprot.org/proteomes/UP000000625}; last modified 2021-03-07, file timestamp 2021-06-16, accessed on 2021-07-12).
          The data frame was created by running \code{CHNOSZ::read.fasta("UP000000625_83333.fasta.gz")}.
        }

        \item{\code{UP000000803_7227.csv.xz}}{
          This data file has amino acid compositions of 4392 proteins in the UniProt reference proteome of \emph{Drosophila melanogaster} (\url{https://www.uniprot.org/proteomes/UP000000803}; last modified 2021-03-07, file timestamp 2021-06-16, accessed on 2021-07-12).
          Used in \code{\link{evdevH2O}}.
        }

        \item{\code{UP000001570_224308.csv.xz}}{
          This data file has amino acid compositions of 4392 proteins in the UniProt reference proteome of \emph{Bacillus subtilis} strain 168 (\url{https://www.uniprot.org/proteomes/UP000001570}; last modified 2021-03-09, file timestamp 2021-06-16, accessed on 2021-07-12).
          Used in \code{\link{evdevH2O}}.
        }


      }

      This directory also has subcellular location data for yeast; see \code{\link{yeast}} for more information.

    }

    \item{\code{vignettes}}{
      This directory has vignettes for differential expression data: TCGA.Rmd, HPA.Rmd, and osmotic_gene.Rmd.
      The CSV files generated by these vignettes are also kept here; they are used for plots in \code{\link{canH2O}}.
    }

    \item{\code{OBIGT}}{
      The \code{OldAA.csv} file has thermodynamic data for glycine and methionine, [Gly] and [Met] sidechain groups, and the protein backbone group from Dick et al. (2006).
      These data have been superseded in the default OBIGT database in CHNOSZ and are kept here in order to reproduce calculations from some papers (\code{\link{aaaq}}, \code{\link{scsc}}, and \code{\link{bison}}).
    }

    \item{\code{RefSeq}}{
      This directory has two files, \file{protein_refseq.csv.xz} (amino acid compositions of species-level archaeal, bacterial, and viral taxa in the RefSeq database) and \file{taxid_names.csv.xz} (taxonomic names for each of those species).
      The scripts to produce these files, and further processed data files, are in the \code{extdata/RefSeq} directory of \pkg{chem16S} (see \code{\link[chem16S]{chem16S-package}}).
    }

  }

}

\references{
Dick JM, LaRowe DE and Helgeson HC (2006) Temperature, pressure, and electrochemical constraints on protein speciation: Group additivity calculation of the standard molal thermodynamic properties of ionized unfolded proteins. \emph{Biogeosciences} \bold{3}, 311--336. \doi{10.5194/bg-3-311-2006}

Dick JM (2009) Calculation of the relative metastabilities of proteins in subcellular compartments of \emph{Saccharomyces cerevisiae}. \emph{BMC Syst. Biol.} \bold{3}, 75. \doi{10.1186/1752-0509-3-75}

Dick JM and Meng D (2021) Local and global chemical shaping of bacterial communities by redox potential. bioRxiv. \doi{10.1101/2021.10.12.464155}

Dick JM and Shock EL (2011) Calculation of the relative chemical stabilities of proteins as a function of temperature and redox chemistry in a hot spring. \emph{PLOS One} \bold{6}, e22782. \doi{10.1371/journal.pone.0022782}

Dick JM and Shock EL (2013) A metastable equilibrium model for the relative abundance of microbial phyla in a hot spring. \emph{PLOS One} \bold{8}, e72395. \doi{10.1371/journal.pone.0072395}

Dick JM (2014) Average oxidation state of carbon in proteins. \emph{J. R. Soc. Interface} \bold{11}, 20131095. \doi{10.1098/rsif.2013.1095}

Dick JM (2016) Proteomic indicators of oxidation and hydration state in colorectal cancer. \emph{PeerJ} \bold{4}, e2238. \doi{10.7717/peerj.2238}

Dick JM (2017) Chemical composition and the potential for proteomic transformation in cancer, hypoxia, and hyperosmotic stress. \emph{PeerJ} \bold{5}, e3421 \doi{10.7717/peerj.3421}

Dick JM, Yu M, Tan J and Lu A (2019) Changes in carbon oxidation state of metagenomes along geochemical redox gradients. \emph{Front. Microbiol.} \bold{10}, 120. \doi{10.3389/fmicb.2019.00120}

Dick JM (2019) CHNOSZ: Thermodynamic calculations and diagrams for geochemistry. \emph{Front. Earth Sci.} \bold{7}:180. \doi{10.3389/feart.2019.00180}

Dick JM, Yu M and Tan J (2020) Uncovering chemical signatures of salinity gradients through compositional analysis of protein sequences. \emph{Biogeosciences} \bold{17}, 6145--6162. \doi{10.5194/bg-17-6145-2020}

Dick JM (2021) Water as a reactant in the differential expression of proteins in cancer. \emph{Comp. Sys. Onco.} \bold{1}:e1007. \doi{10.1002/cso2.1007}

Dick JM and Shock EL (2021) The release of energy during protein synthesis at ultramafic-hosted submarine hydrothermal ecosystems. \emph{J. Geophys. Res.: Biogeosciences} \bold{126}, e2021JG006436. \doi{10.1029/2021JG006436}

Dick JM and Tan J (2022) Chemical links between redox conditions and estimated community proteomes from 16S rRNA and reference protein sequences. \emph{Microb. Ecol.} \doi{10.1007/s00248-022-01988-9}

Dick JM (2022) A thermodynamic model for water activity and redox potential in evolution and developent. \emph{J. Mol. Evol} \bold{90}, 182--199. \doi{10.1007/s00239-022-10051-7}

Dick JM, Boyer GM, Canovas PA III and Shock EL (2022) Using thermodynamics to obtain geochemical information from genomes. Manuscript in preparation.
}
