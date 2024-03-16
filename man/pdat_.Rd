\encoding{UTF-8}
\name{pdat_}
\alias{pdat_}
\alias{pdat_breast}
\alias{pdat_colorectal}
\alias{pdat_lung}
\alias{pdat_liver}
\alias{pdat_pancreatic}
\alias{pdat_prostate}
\alias{pdat_hypoxia}
\alias{pdat_secreted}
\alias{pdat_3D}
\alias{pdat_glucose}
\alias{pdat_osmotic_bact}
\alias{pdat_osmotic_euk}
\alias{pdat_osmotic_halo}
\alias{.pdat_multi}
\alias{.pdat_osmotic}
\alias{pdat_TCGA}
\alias{pdat_HPA}
\alias{pdat_osmotic_gene}
\alias{pdat_aneuploidy}
\alias{pdat_yeast_stress}
\alias{pdat_fly}
\title{Get protein expression data}
\description{
  These functions are used to retrieve the amino acid compositions of differentially expressed proteins or proteins corresponding to differentially expressed genes.
}

\usage{
  pdat_breast(dataset = 2020)
  pdat_colorectal(dataset = 2020)
  pdat_liver(dataset = 2020)
  pdat_lung(dataset = 2020)
  pdat_pancreatic(dataset = 2020)
  pdat_prostate(dataset = 2020)
  pdat_hypoxia(dataset = 2020)
  pdat_secreted(dataset = 2020)
  pdat_3D(dataset = 2020)
  pdat_glucose(dataset = 2020)
  pdat_osmotic_bact(dataset = 2020)
  pdat_osmotic_euk(dataset = 2020)
  pdat_osmotic_halo(dataset = 2020)
  .pdat_multi(dataset = 2020)
  .pdat_osmotic(dataset = 2017)
  pdat_osmotic_gene(dataset = 2020)
  pdat_TCGA(dataset = 2020)
  pdat_HPA(dataset = 2020)
  pdat_aneuploidy(dataset = 2020)
  pdat_yeast_stress(dataset = 2020)
  pdat_fly(dataset = NULL)
}

\arguments{
  \item{dataset}{character, dataset name}
}

\details{

The \code{pdat_} functions assemble lists of up- and down-regulated proteins and retrieve their amino acid compositions using \code{\link[canprot]{human_aa}}.
The result can be used with \code{\link{get_comptab}} to make a table of chemical metrics that can then be plotted with \code{\link{diffplot}}.

If \code{dataset} is \samp{2020} (the default) or \samp{2017}, the function returns the names of all datasets in the compilation for the respective year.

Each dataset name starts with a reference key indicating the study (i.e. paper or other publication) where the data were reported.
The reference keys are made by combining the first characters of the authors' family names with the 2-digit year of publication.

If a study has more than one dataset, the reference key is followed by an underscore and an identifier for the particular dataset.
This identifier is saved in the variable named \code{stage} in the functions, but can be any descriptive text.

To retrieve the data, provide a single dataset name in the \code{dataset} argument.
Protein expression data is read from the CSV files stored in \code{extdata/expression/} or \code{extdata/canH2O}, under the subdirectory corresponding to the name of the \code{pdat_} function.
Some of the functions also read amino acid compositions (e.g. for non-human proteins) from the files in \code{extdata/aa/}.

The following functions are used in the \code{\link{canH2O}} paper:

\itemize{
\item \code{pdat_colorectal}, \code{pdat_pancreatic}, \code{pdat_breast}, \code{pdat_lung}, \code{pdat_prostate}, and \code{pdat_liver} retrieve data for protein expression in different cancer types.

\item \code{pdat_hypoxia} gets data for cellular extracts in hypoxia and \code{pdat_secreted} gets data for secreted proteins (e.g. exosomes) in hypoxia.

\item \code{pdat_3D} retrieves data for 3D (e.g. tumor spheroids and aggregates) compared to 2D (monolayer) cell culture.

\item \code{.pdat_osmotic} retrieves data for hyperosmotic stress, for the 2017 compilation only.
In 2020, this compilation was expanded and split into \code{pdat_osmotic_bact} (bacteria), \code{pdat_osmotic_euk} (eukaryotic cells) and \code{pdat_osmotic_halo} (halophilic bacteria and archaea).

\item \code{pdat_glucose} gets data for high-glucose experiments in eukaryotic cells.

\item \code{.pdat_multi} retrieves data for studies that have multiple types of datasets (e.g. both cellular and secreted proteins in hypoxia), and is used internally by the specific functions (e.g. \code{pdat_hypoxia} and \code{pdat_secreted}).
}

\code{pdat_osmotic_gene}, which is used in the \code{\link{gradH2O}} paper, gets data for differentially expressed genes in transcriptomic experiments for bacteria.

\code{pdat_fly}, which is used in the \code{\link{evdevH2O}} paper, gets data for differentially expressed genes or proteins in adult flies compared to embryos (Fabre et al., 2019).

The following functions are used in the \code{\link{canH2O}} paper:

\code{pdat_TCGA} calculates chemical composition for proteins corresponding to differentially expressed genes listed by GEPIA2 (Tang et al., 2019).
Expression levels for normal tissue and cancer are from the Genotype-Tissue Expression project (GTEx Consortium, 2017) and The Cancer Genome Atlas (TCGA Research Network et al., 2013; Hutter and Zenklusen, 2018).

\code{pdat_HPA} calculates chemical composition for protein expression data from the Human Protein Atlas (Uhlen et al., 2015).

\code{pdat_aneuploidy} gets data for differentially expressed genes in aneuploid vs haploid yeast cells (Tsai et al., 2019).

\code{pdat_yeast_stress} gets data for differentially expressed genes in hyper- and hypo-osmotic stress in yeast (Gasch et al., 2000).

Except for \code{pdat_aneuploidy} and \code{pdat_fly}, each of the functions has a corresponding vignette.
These vignettes are located in \code{inst/vignettes}, so they are not built with the package, but can be run using \code{\link{mkvig}}.
The output files (CSV) from the vignettes are included in the package, and they are used in some figures: \code{\link{gradH2O7}} (osmotic_gene), \code{\link{canH2O3}} (TCGA, HPA), \code{\link{canH2O5}} (yeast_stress).
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

Fabre B, Korona D, Lees JG, Lazar I, Livneh I, Brunet M, Orengo CA, Russell S and Lilley KS (2019) Comparison of \emph{Drosophila melanogaster} embryo and adult proteome by SWATH-MS reveals differential regulation of protein synthesis, degradation machinery, and metabolism modules. \emph{J. Proteome Res.} \bold{18}, 2525--2534. \doi{10.1021/acs.jproteome.9b00076}

Gasch AP, Spellman PT, Kao CM, Carmel-Harel O, Eisen MB, Storz G, Botstein D and Brown PO (2000) Genomic expression programs in the response of yeast cells to environmental changes. \emph{Mol. Biol. Cell} \bold{11}, 4241--4257. \doi{10.1091/mbc.11.12.4241}

Hutter C and Zenklusen JC (2018) The Cancer Genome Atlas: Creating lasting value beyond its data. \emph{Cell} \bold{173}, 283--285. \doi{10.1016/j.cell.2018.03.042}

Tang Z, Kang B, Li C, Chen T and Zhang Z (2019) GEPIA2: an enhanced web server for large-scale expression profiling and interactive analysis. \emph{Nucleic Acids Res.} \bold{47}, W556--W560. \doi{10.1093/nar/gkz430}

The Cancer Genome Atlas Research Network et al. (2013) The Cancer Genome Atlas Pan-Cancer analysis project. \emph{Nat. Genetics} \bold{45}, 1113--1120. \doi{10.1038/ng.2764}

Tsai H-J et al. (2019) Hypo-osmotic-like stress underlies general cellular defects of aneuploidy. \emph{Nature} \bold{570}, 117--121. \doi{10.1038/s41586-019-1187-2}

Uhlen M et al. (2015) A pathology atlas of the human cancer transcriptome. \emph{Science} \bold{357}, eaan2507. \doi{10.1126/science.aan2507}
}

\value{
A list consisting of:
\describe{
  \item{\code{dataset}}{Name of the dataset}
  \item{\code{description}}{Descriptive text for the dataset, used for making the tables in the vignettes}
  \item{\code{pcomp}}{UniProt IDs together with amino acid compositions obtained using \code{\link[canprot]{human_aa}}}
  \item{\code{up2}}{Logical vector with length equal to the number of proteins; TRUE for up-regulated proteins and FALSE for down-regulated proteins}
}
}

\seealso{
\code{\link{canprot}}
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

# List datasets in the 2017 complilation for colorectal cancer
pdat_colorectal(2017)
# Get proteins and amino acid compositions for one dataset
pdat_colorectal("JKMF10")

}
