\encoding{UTF-8}
\name{orp16S}
\alias{orp16S}
\alias{orp16S1}
\alias{orp16S2}
\alias{orp16S3}
\alias{orp16S4}
\alias{orp16S5}
\alias{orp16S6}
\alias{orp16S_S1}
\alias{getmdat_orp16S}
\alias{getmetrics_orp16S}
\alias{orp16S_info}
\title{Plots from chemical shaping of bacterial communities paper}
\description{
Plots from the paper by Dick and Meng (2021).
}

\usage{
  orp16S1(pdf = FALSE)
  orp16S2(pdf = FALSE)
  orp16S3(pdf = FALSE)
  orp16S4(pdf = FALSE)
  orp16S5(pdf = FALSE)
  orp16S6(pdf = FALSE)
  orp16S_S1(pdf = FALSE)
  getmdat_orp16S(study, metrics = NULL, dropNA = TRUE, size = NULL, quiet = TRUE)
  getmetrics_orp16S(study, mincount = 100, quiet = TRUE, ...)
  orp16S_info(study)
}

\arguments{
  \item{pdf}{logical, make a PDF file?}
  \item{study}{character, study name}
  \item{metrics}{data frame, output of \code{\link[chem16S]{getmetrics}}}
  \item{dropNA}{logical, exclude samples with NA name in metadata?}
  \item{size}{numeric, number of samples to randomly select}
  \item{quiet}{logical, change to FALSE to print details about data processing}
  \item{mincount}{integer, samples with less than this number of RDP classifications are excluded}
  \item{...}{additional arguments passed to \code{\link[chem16S]{readRDP}}}
}

\details{
This table gives a brief description of each plotting function.

\tabular{ll}{
  \code{orp16S1} \tab Geobiochemical predictive framework \cr
  \code{orp16S2} \tab Chemical and geobiochemical depth profiles in Winogradsky columns \cr
  \code{orp16S3} \tab Sample locations on world map \cr
  \code{orp16S4} \tab Analysis of selected datasets for each environment type \cr
  \code{orp16S5} \tab Local and global analysis of \ZC and Eh7 data \cr
  \code{orp16S6} \tab Distinctions in carbon oxidation state estimated for different hot springs, and global fits for all environments. \cr
  \code{orp16S_S1} \tab \ZC-Eh scatterplots for all studies \cr
}

\code{orp16S_S1} creates the files EZdat.csv and EZlm.csv that are used by the other functions.

\code{getmdat_orp16S} gets metadata for the indicated \code{study}, calculates Eh7 from Eh and pH and T (if available), classifies samples into 3 clusters according to their Eh values, and adds columns for plot parameters (\code{pch}, \code{col}).
The default for \code{dropNA} means to exclude samples with NA name in the metadata file.
If \code{metrics} is supplied, metadata are kept only for samples with available metrics, the sample metrics are placed in the same order as the remaining metadata, and the function returns a list with both \samp{metadata} and the sorted \samp{metrics}.
\code{size} can be used to specify a number of samples (from among those with both metadata and \code{metrics}, if given) to be randomly selected.

\code{getmetrics_orp16S} calculates chemical metrics (\ZC and \nH2O) for the indicated study.
\code{...} is used to supply argument to \code{\link[chem16S]{readRDP}}; note that \code{mincount = 50} is needed to reproduce the data processing used in the paper (see \code{\link{plotEZ}}).

\code{orp16S_info} prints general information about a study (more specifically, a dataset): name, number of samples with Bacteria and Archaea (with \code{mincount = 50}), bibligraphic key, range of T, pH, Eh, and Eh7, and whether the slope of the ZC-Eh7 correlations for Bacteria and Archaea (if available) are negative (< -0.01 V-1), positive (> 0.01 V-1), or close to zero.
This information is used to fill in Table S1 of the paper.

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
  }
}

\section{Other files}{
  \describe{
    \item{\code{doc/orp16S.bib}}{Reference keys for all datasets (all of them are listed here; only some of them are cited in the vignette orp16S.Rmd).}
  }
}

\seealso{
\code{\link{plotEZ}} is used to make the \ZC-Eh7 scatterplots and linear fits in \code{orp16S4} and \code{orp16S_S1}.
}

\references{
Dick, J. M. and Meng, D. (2021) Local and global chemical shaping of bacterial communities by redox potential. bioRxiv. \doi{10.1101/2021.10.12.464155}

USGS (2010) Great Lakes and Watersheds Shapefiles. ScienceBase Catalog, U.S. Geological Survey. \doi{https://www.sciencebase.gov/catalog/item/530f8a0ee4b0e7e46bd300dd}
}

\examples{
orp16S1()

orp16S2()

orp16S3()

# "Raw" metrics with lots of messages printed
getmetrics_orp16S("SVH+19", mincount = 50, quiet = FALSE)
# Metrics with associated metadata
metrics <- getmetrics_orp16S("SPA+21", mincount = 50)
mdat <- getmdat_orp16S("SPA+21", metrics = metrics)
# This has 'metadata' and 'metrics' data frames with same number of rows
str(mdat)
stopifnot(length(unique(sapply(mdat, nrow))) == 1)

# Printing info about one dataset
orp16S_info("SPA+21")
# Some keys have suffixes to retrieve subsets of data (see code of getmdat_orp16S)
orp16S_info("PCL+18_Acidic")
orp16S_info("PCL+18_Alkaline")

# Check that references for all studies are available
mdatdir <- system.file("extdata/orp16S/metadata", package = "JMDplots")
studies <- gsub(".csv", "", dir(mdatdir))
bibfile <- system.file("doc/orp16S.bib", package = "JMDplots")
bibentry <- bibtex::read.bib(bibfile)
stopifnot(all(studies \%in\% names(bibentry)))
}
