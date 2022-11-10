\encoding{UTF-8}
\name{orp16S}
\alias{orp16S}
\alias{orp16S_1}
\alias{orp16S_2}
\alias{orp16S_3}
\alias{orp16S_4}
\alias{orp16S_5}
\alias{orp16S_6}
\alias{orp16S_S1}
\alias{orp16S_S2}
\alias{orp16S_info}
\alias{orp16S_T1}
\alias{orp16S_D3}
\alias{getmdat_orp16S}
\alias{getmetrics_orp16S}
\title{Plots from chemical shaping of bacterial communities paper}
\description{
Plots from the paper by Dick and Meng (2021).
}

\usage{
  orp16S_1(pdf = FALSE)
  orp16S_3(pdf = FALSE)
  orp16S_4(pdf = FALSE, EMP_primers = FALSE)
  orp16S_5(pdf = FALSE)
  orp16S_6(pdf = FALSE)

  orp16S_S1(pdf = FALSE)
  orp16S_S2(pdf = FALSE)

  orp16S_info(study)
  orp16S_T1(samesign = FALSE)
  orp16S_D3(mincount = 100)

  getmdat_orp16S(study, metrics = NULL, dropNA = TRUE, size = NULL, quiet = TRUE)
  getmetrics_orp16S(study, mincount = 100, quiet = TRUE, ...)
}

\arguments{
  \item{pdf}{logical, make a PDF file?}
  \item{EMP_primers}{logical, include only datasets using Earth Microbime Project primers (515F/806R)?}
  \item{study}{character, study name}
  \item{samesign}{logical, only count slopes with same sign of minimum and maximum values in the 95\% CI?}
  \item{mincount}{integer, samples with less than this number of RDP classifications are excluded}
  \item{metrics}{data frame, output of \code{\link[chem16S]{getmetrics}}}
  \item{dropNA}{logical, exclude samples with NA name in metadata?}
  \item{size}{numeric, number of samples to randomly select}
  \item{quiet}{logical, change to FALSE to print details about data processing}
  \item{...}{additional arguments passed to \code{\link[chem16S]{readRDP}}}
}

\details{
This table gives a brief description of each plotting function.

\tabular{ll}{
  \code{orp16S_1} \tab Methods overview and chemical depth profiles in Winogradsky columns \cr
  \code{orp16S_2} \tab Sample locations on world map and Eh-pH diagram (outline is from Baas Becking et al., 1960) \cr
  \code{orp16S_3} \tab Associations between Eh7 and \ZC at local scales \cr
  \code{orp16S_4} \tab Associations between Eh7 and \ZC at a global scale. Set \code{EMP_primers} to TRUE to make Figure S3. \cr
  \code{orp16S_5} \tab Comparison of 16S_-based community reference proteomes with metaproteomes \cr
  \code{orp16S_6} \tab Most abundant genera at low and high Eh7 in geothermal and hyperalkaline areas \cr

  \code{orp16S_S1} \tab \ZC-Eh scatterplots for each dataset \cr
  \code{orp16S_S2} \tab Comparison of Eh7, Eh, and \O2 as predictors of carbon oxidation state \cr
}

\code{orp16S_S1} creates the files \file{EZdat.csv} and \file{EZlm.csv} that are used by the other plotting functions and correspond to Datasets S1 and D2 in the paper.
\file{EZdat.csv} has sample data and metadata: study key, environment type, lineage (Bacteria or Archaea), sample name, SRA or MG-RAST accession number, grouping description, sample group, values of T (°C), pH, \O2 concentration (\umol/L), Eh (mV), Eh7 (mV), and \ZC.
\file{EZlm.csv} has regression results for each dataset and domain: study key, environment type, lineage, number of samples, Eh7 range, slope (V\S{-1}) and intercept (\ZC) of of linear regression, and Pearson correlation coefficient.

\code{orp16S_T1} returns a table summarizing regression slopes for all datasets (based on \file{EZlm.csv}), corresponding to Table 1a in the paper.
\code{orp16S_T1(samesign = TRUE)} generates the values listed in Table 1b in the paper.
\code{orp16S_D3} creates a file with percent of RDP classifications to genus level mapped to the NCBI taxonomy, and genus names and \code{ZC} for the groups of samples corresponding to the first and fourth quartile of Eh7 values (i.e., reducing and oxidizing conditions), corresponding to Dataset S3 in the paper.

\code{getmdat_orp16S} gets metadata for the indicated \code{study}, calculates Eh7 from Eh and pH and T (if available), classifies samples into 3 clusters according to their Eh values, and adds columns for plot parameters (\samp{pch}, \samp{col}).
The function also adds a column named \samp{O2_umol_L}, with \O2 concentration in (\umol/L), which is converted from source units of \O2 concentration as needed.
The default for \code{dropNA} means to exclude samples with NA name in the metadata file.
If \code{metrics} is supplied, metadata are kept only for samples with available metrics, the sample metrics are placed in the same order as the remaining metadata, and the function returns a list with both \samp{metadata} and the sorted \samp{metrics}.
\code{size} can be used to specify a number of samples (from among those with both metadata and \code{metrics}, if given) to be randomly selected.

\code{getmetrics_orp16S} calculates chemical metrics (\ZC and \nH2O) for the indicated study.
\code{...} is used to supply argument to \code{\link[chem16S]{readRDP}}; note that the default of \code{mincount = 100} was used for data processing in the paper (see \code{\link{plotEZ}}).

\code{orp16S_info} prints general information about a study (more specifically, a dataset): name, number of samples with Bacteria and Archaea (with \code{mincount = 100}), bibligraphic key, range of T, pH, Eh, and Eh7, and whether the slope of the \ZC-Eh7 correlations for Bacteria and Archaea (if available) are negative or positive.
The function also returns information in a data frame, which is used to fill in Table S1 of the paper.

}

\section{Files in extdata/orp16S}{
  \describe{
    \item{\file{pipeline.R}}{Pipeline for sequence data processing (uses external programs fastq-dump, vsearch, seqtk, RDP Classifier).}
    \item{\file{metadata/*.csv}}{Sample metadata for each study.}
    \item{\file{RDP/*.csv.xz}}{
      RDP Classifier results combined into a single CSV file for each study, created with the \code{classify} and \code{mkRDP} function2 in \file{pipeline.R}.
    }
    \item{\file{hydro_p}}{Shapefiles for the North American Great Lakes, downloaded from USGS (2010).}
    \item{\file{EZdat.csv}}{Data for all samples, created by \code{orp16S_S1}.}
    \item{\file{EZlm.csv}}{Linear fits between Eh7 and \ZC for each dataset created by \code{orp16S_S1}.}
    \item{\file{BKM60.csv}}{Outline of Eh-pH range of natural environments, digitized from Fig. 32 of Baas Becking et al. (1960).}
    \item{\file{HWLH22_Jun_2018_aa.csv}, \file{HWLH22_Sep_2018_aa.csv}}{Sum of amino acid composition of all peptides in each sample from the metaproteomic study of Han et al. (2022).}
    \item{\file{mkaa.R}}{Script used to process source metaproteomic data to make the \file{*aa.csv} files.}
    \item{\file{Dataset_S3.csv}}{File created by \code{orp16S_D3}.}
  }
}

\section{Other files}{
  \describe{
    \item{\code{doc/orp16S.bib}}{Reference keys for all datasets (all of them are listed here; only some of them are cited in the vignette orp16S.Rmd).}
  }
}

\seealso{
\code{\link{plotEZ}} is used to make the \ZC-Eh7 scatterplots and linear fits in \code{orp16S_3} and \code{orp16S_S1}.
}

\references{
Baas Becking LGM, Kaplan IR and Moore D (1960) Limits of the natural environment in terms of pH and oxidation-reduction potentials. \emph{Journal of Geology} \bold{68}, 243--284. \url{https://www.jstor.org/stable/30059218}

Dick JM and Meng D (2021) Local and global chemical shaping of bacterial communities by redox potential. bioRxiv. \doi{10.1101/2021.10.12.464155}

Han Y, Wang Q, Li Q and Hu C (2022) Active metabolism and biomass dynamics of biocrusts are shaped by variation in their successional state and seasonal energy sources. \emph{Science of the Total Environment} \bold{831}, 154756. \doi{10.1016/j.scitotenv.2022.154756}

USGS (2010) Great Lakes and Watersheds Shapefiles. ScienceBase Catalog, U.S. Geological Survey. \doi{https://www.sciencebase.gov/catalog/item/530f8a0ee4b0e7e46bd300dd}
}

\examples{
# Methods outline and Winogradsky columns
orp16S_1()

# Sample locations and Eh-pH diagram
orp16S_2()

# Summarize results for individual datasets
orp16S_T1()

# "Raw" metrics with lots of messages printed
# NOTE: the messages include values for the last two columns of Table S2 of the paper
# (% Classification to genus level and % Mapping to NCBI taxonomy)
getmetrics_orp16S("MLL+19", quiet = FALSE)
# Metrics with associated metadata
metrics <- getmetrics_orp16S("SPA+21")
mdat <- getmdat_orp16S("SPA+21", metrics = metrics)
# This has 'metadata' and 'metrics' data frames with same number of rows
str(mdat)
stopifnot(length(unique(sapply(mdat, nrow))) == 1)

# Printing info about one dataset
orp16S_info("SPA+21")
# Some keys have suffixes to retrieve subsets of data (see code of getmdat_orp16S)
orp16S_info("PCL+18_Acidic")
orp16S_info("PCL+18_Alkaline")

# There's no plotmet_orp16S, but this is what it would look like in one line
plotmet(getmdat_orp16S("MLL+19", getmetrics_orp16S("MLL+19")))

\dontrun{
# Get unexported object with datasets in each environment type
envirotype <- JMDplots:::envirotype
# Function to get info for all datasets in one environment type
envinfo <- function(ienv) do.call(rbind, lapply(envirotype[[ienv]], orp16S_info))
# Get info for all datasets in all environment types
allenv <- lapply(1:7, envinfo)
# Put in names of environments
cbindenv <- function(a, b) cbind(envirotype = a, b)
allenv <- mapply(cbindenv, names(envirotype)[1:7], allenv, SIMPLIFY = FALSE)
allenv <- do.call(rbind, c(allenv, make.row.names = FALSE))
# Now write to CSV file for copying to Table S1
write.csv(allenv, "Table_S1_values.csv", row.names = FALSE)
}

# Check that references for all studies are available
mdatdir <- system.file("extdata/orp16S/metadata", package = "JMDplots")
studies <- gsub(".csv", "", dir(mdatdir))
bibfile <- system.file("doc/orp16S.bib", package = "JMDplots")
bibentry <- bibtex::read.bib(bibfile)
stopifnot(all(studies \%in\% names(bibentry)))
}
