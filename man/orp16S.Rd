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
}

\arguments{
  \item{pdf}{logical, make a PDF file?}
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

To set up the \code{\link{chem16S}} environment, change \code{options(chem16Sdir)} before running these functions (see example).

\code{orp16S_S1} creates the files EZdat.csv and EZlm.csv that are used by the other functions.

}

\section{Files in extdata/orp16S}{
  \describe{
    \item{\code{hydro_p}}{Shapefiles for the North American Great Lakes, downloaded from USGS (2010).}
    \item{\code{EZdat.csv}}{Data for all samples, created by \code{orp16S_S1}.}
    \item{\code{EZlm.csv}}{Linear fits between Eh7 and \ZC for each dataset created by \code{orp16S_S1}.}
    \item{\code{metadata/*.csv}}{Sample metadata for each study.}
    \item{\code{RDP/*.csv.xz}}{RDP Classifier results combined into a single CSV file for each study, created with the \code{mkRDP} function in the \code{process.R} script (\code{\link{chem16S}}).}
    \item{\code{OTU2RDP.R}}{Script to construct \code{RDP/OFY+19.csv} from Table S2 of Omae et al. (2019).}
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
Dick, J. M. and Meng, D. (2021) Local and global chemical shaping of bacterial communities by redox potential.

USGS (2010) Great Lakes and Watersheds Shapefiles. ScienceBase Catalog, U.S. Geological Survey. \doi{https://www.sciencebase.gov/catalog/item/530f8a0ee4b0e7e46bd300dd}

Omae, K., Fukuyama, Y., Yasuda, H., Mise, K., Yoshida, T. and Sako, Y. (2019) Diversity and distribution of thermophilic hydrogenogenic carboxydotrophs revealed by microbial community analysis in sediments from multiple hydrothermal environments in Japan. \emph{Archives of Microbiology} \bold{201}, 969-982. \doi{10.1007/s00203-019-01661-9}
}

\examples{
orp16S1()

options(chem16Sdir = system.file("extdata/orp16S", package = "JMDplots"))

orp16S2()

orp16S3()

# Check that references for all studies are available
mdatdir <- system.file("extdata/orp16S/metadata", package = "JMDplots")
studies <- gsub(".csv", "", dir(mdatdir))
bibfile <- system.file("doc/orp16S.bib", package = "JMDplots")
bibentry <- bibtex::read.bib(bibfile)
stopifnot(all(studies \%in\% names(bibentry)))
}
