\encoding{UTF-8}
\name{subsurface}
\alias{subsurface}
\title{User-provided compositional data for metagenomes}
\description{
Metadata for two subsurface metagenomes, used as an example of user-provided data.
}

\details{
To use user-provided data, do one of the following:
Call \code{\link{ppage}} with an argument for \code{subset} that is a list with the metadata, and have the CSV files of DNA, RNA, or protein composition under the \code{data} directory of the working directory.
Or, call \code{\link{mplot}} with the path of the data directory given in \code{datadir}.

For the example here, the \code{extdata/subsurface/data} directory has additional files of amino acid compositions for two subsurface environments: North Pond metagenome (Tully et al., 2018) and Peru Margin metatranscriptome (Orsi et al., 2013).
}

\references{
Orsi WD, Edgcomb VP, Christman GD and Biddle JF (2013) Gene expression in the deep biosphere. \emph{Nature} \bold{499}, 205-208. \doi{10.1038/nature12230}

Tully BJ, Wheat CG, Glazer BT and Huber JA (2018) A dynamic microbial community with high functional redundancy inhabits the cold, oxic subseafloor aquifer \emph{ISME J.} \bold{12}, 1-16. \doi{10.1038/ismej.2017.187}
}

\examples{
# create the metadata
subsurface <- list(
  Peru_Margin = list( # 20180313 [OECB13]
    depth = c("5m", "30m", "50m", "70m", "91m", "159m"),
    SRA_MT = c("SRX187180","SRX330851","SRX331883","SRX331884","SRX331956","SRX331958"),
    MT_range = c(1.64, 1.72),
    MTP_range = c(-0.18, -0.13),
    abbrev = "PM",
    group = "rock",
    plot_real_x = FALSE
  ),
  North_Pond = list( # 20180320 [TWGH18]
    # time points of samples from U1832A 
    sample = c("TP1", "TP2", "TP3", "TP4", "TP5", "TP6", "TP7", "TP8", "TP9", "CTD"),
    SRA_MG = c("SRX3143893", "SRX3143894", "SRX3143895", "SRX3143896", "SRX3143897",
               "SRX3143898", "SRX3143899", "SRX3143900", "SRX3143886", "SRX3143890"),
    MG_range = c(1.63, 1.71),
    MGP_range = c(-0.18, -0.14),
    abbrev = "NP",
    group = c("oxygenated", "oxygenated", "oxygenated", "oxygenated", "oxygenated",
              "oxygenated", "oxygenated", "oxygenated", "oxygenated", "oxygenated0"),
    plot_real_x = FALSE
  )
)

# make plot for one data set
datadir <- system.file("extdata/subsurface/data", package = "JMDplots")
mplot("North_Pond", "SRA_MGP", datadir = datadir, mdata = subsurface)

# make plots for multiple data sets
WD <- setwd(system.file("extdata/subsurface", package = "JMDplots"))
ppage(subsurface)
setwd(WD)
}
