\encoding{UTF-8}
\name{subsurface}
\alias{subsurface}
\title{User-provided compositional data for metagenomes}
\description{
Metadata for two subsurface metagenomes, used as an example of user-provided data.
}

\usage{
  data(subsurface)
}

\details{
To use user-provided data, do one of the following:
Call \code{\link{ppage}} with an argument for \code{subset} that is a list with the metadata, and have the CSV files of DNA, RNA, or protein composition under the \code{data} directory of the working directory.
Or, call \code{\link{mplot}} with the path of the data directory given in \code{datadir}.

For the example here, the \code{extdata/subsurface/data} directory has additional files of amino acid compositions for two subsurface environments: North Pond metagenome (Tully et al., 2018) and Peru Margin metatranscriptome (Orsi et al., 2013).
}

\references{
Orsi, W. D, Edgcomb, V. P., Christman, G. D. and Biddle, J. F. (2013) Gene expression in the deep biosphere. \emph{Nature} \bold{499}, 205-208. \doi{10.1038/nature12230}

Tully, B. J., Wheat, C. G., Glazer, B. T. and Huber, J. A. (2018) A dynamic microbial community with high functional redundancy inhabits the cold, oxic subseafloor aquifer \emph{ISME J.} \bold{12}, 1-16. \doi{10.1038/ismej.2017.187}
}

\examples{
# load and summarize the metadata
data(subsurface)
str(subsurface)

# make plot for one data set
datadir <- system.file("extdata/subsurface/data", package = "JMDplots")
mplot("North_Pond", "SRA_MGP", datadir = datadir, mdata = subsurface$studies)

# make plots for multiple data sets
WD <- setwd(system.file("extdata/subsurface", package = "JMDplots"))
ppage(subsurface)
setwd(WD)
}
