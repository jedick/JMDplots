\encoding{UTF-8}
\name{yeast}
\alias{yeast}
\alias{yeast.aa}
\alias{yeastgfp}
\alias{yeast.ZC}
\title{Composition, Localization, and Abundances of Proteins in Yeast}
\description{
Retrieve the amino acid compositions of one or more proteins from \emph{Saccharomyces cerevisiae} and get localizations and abundances reported by the YeastGFP project or locations from the SGD gene associations.
}

\usage{
  yeast.aa(protein = NULL)
  yeastgfp(location, exclusive = TRUE)
  yeast.ZC(location)
}

\arguments{
  \item{protein}{character, name of protein}
  \item{location}{character, name (\code{yeastgpf}) or one-letter code (\code{yeast.ZC}) for subcellular location}
  \item{exclusive}{logical, report only proteins exclusively localized to a compartment?}
}

\details{
\code{yeast.aa} retrieves the amino acid composition(s) of the indicated proteins in \emph{Saccharomyces cerevisiae}.
The calculation depends on the data file \code{\link{extdata}/organisms/Sce.csv.xz}, which contains the amino acid compositions of the proteins.
The \code{protein} argument should be a vector or a list of vectors of one or more SGD IDs, Open Reading Frame (ORF) or gene names that are found in these files.
The output data frame contains rows with NA compositions for names that are not matched.
If \code{protein} is NULL, all available proteins are returned except those with NA amino acid composition.

\code{yeastgfp} returns the identities and abundances of proteins with the requested subcellular localization(s) (specified in \code{location}) using data from the YeastGFP project that is stored in \code{\link{extdata}/abundance/yeastgfp.csv.xz}.
If \code{exclusive} is \code{FALSE}, the function grabs all proteins that are localized to a compartment even if they are also localized to other compartments.
If \code{exclusive} is \code{TRUE} (the default), only those proteins that are localized exclusively to the requested compartments are identified, unless there are no such proteins, then the non-exclusive localizations are used (applies to the \samp{bud} localization).

\code{yeast.ZC} returns the \code{\link[CHNOSZ:util.formula]{ZC}} of proteins in the subcellular location identified by the one-letter abbreviation.
This abbreviation identifies a \dQuote{Cellular Component} in the Gene Ontology and the associatied proteins (see \code{\link{aoscp}}).
}

\value{
For \code{yeast.aa}, a data frame, or list of data frames, containing the amino acid composition(s) of the specified protein(s) in the format of \code{\link{thermo}$protein}.

For \code{yeastgfp}, a list with elements named \code{protein} (names of proteins) and \code{abundance} (counts or concentrations without any conversion from the units in the data file).
If \code{location} is NULL, \code{yeastgfp} returns the names of all known locations, and if the length of \code{location} is >1, the \code{protein} and \code{abundance} values are lists of the results for each location.
}

\section{Files in extdata/organisms}{
  \describe{
    \item{\code{yeastgfp.csv.xz}}{
      Has 28 columns; the names of the first five are \code{yORF}, \code{gene name}, \code{GFP tagged?}, \code{GFP visualized?}, and \code{abundance}.
      The remaining columns correspond to the 23 subcellular localizations considered in the YeastGFP project (Huh et al., 2003 and Ghaemmaghami et al., 2003) and hold values of either \code{T} or \code{F} for each protein.
      \samp{yeastgfp.csv} was downloaded on 2007-02-01 from http://yeastgfp.ucsf.edu using the Advanced Search, setting options to download the entire dataset and to include localization table and abundance, sorted by orf number.
      Used in \code{yeastgfp}.
    }

    \item{\code{TBD+05.csv}}{
      lists genes with transcriptomic expression changes in carbon limitation stress response experiments in yeast (Tai et al., 2005). Used in an example in this help page.
    }
  }
}


\seealso{
\code{\link{scsc}}
}

\examples{
\dontshow{reset()}
# the first few names in UniProt for "aminotransferase yeast"
genes <- c("AATC", "ARO8", "BCA1", "AMPL", "BCA2", "ARO9")
# the corresponding ORF names
ORF <- c("YLR027C", "YGL202W", "YHR208W", "YKL103C", "YJR148W", "YHR137W")
# we only match two of them by gene name, but all by ORF name
aa <- yeast.aa(genes)
aa <- yeast.aa(ORF)
# what are their formulas and average oxidation states of carbon
protein.formula(aa)
ZC(protein.formula(aa))

## potential fields for overall protein compositions 
## transcriptionally induced and repressed in aerobic
## and anaerobic carbon limitation
## (experiments of Tai et al., 2005)
# the activities of ammonium and sulfate used here
# are similar to the non-growth-limiting concentrations
# used by Boer et al., 2003
basis(c("glucose", "H2O", "NH4+", "hydrogen", "SO4-2", "H+"),
  c(-1, 0, -1.3, 999, -1.4, -7))
# the names of the experiments in TBD+05.csv
expt <- c("Clim.aerobic.down", "Clim.aerobic.up",
  "Clim.anaerobic.down", "Clim.anaerobic.up")
file <- system.file("extdata/organisms/TBD+05.csv", package="JMDplots")
dat <- read.csv(file, as.is=TRUE)
# yeast.aa: get the amino acid compositions
# aasum: average them together
for(thisexpt in expt) {
  p <- dat$protein[dat[, thisexpt]]
  aa <- yeast.aa(p)
  aa <- aasum(aa, average=TRUE, protein=thisexpt)
  add.protein(aa)
}
species(expt, "Sce")
a <- affinity(C6H12O6=c(-30, 0), H2=c(-20, 0))
d <- diagram(a, normalize=TRUE, fill=NULL)
title(main=paste("Formation potential of proteins associated with\n",
  "transcriptional response to carbon limitation in yeast"))
# the affinity of formation favors the proteins upregulated 
# by carbon limitation at low chemical potentials of C6H12O6 ...
res <- length(a$vals[[1]])
stopifnot(c(d$predominant[1, 1], d$predominant[1, res])==grep("up", expt))
# ... and favors proteins downregulated by aerobic conditions
# at high hydrogen fugacities
stopifnot(c(d$predominant[res, res], d$predominant[res, 1])==grep("down", expt))

## overall oxidation state of proteins exclusively localized 
## to cytoplasm of S. cerevisiae with/without abundance weighting
y <- yeastgfp("cytoplasm")
aa <- yeast.aa(y$protein)
aaavg <- aasum(aa, average=TRUE)
ZC(protein.formula(aaavg))
# the average composition weighted by abundance
waaavg <- aasum(aa, abundance=y$abundance, average=TRUE)
ZC(protein.formula(waaavg))
}

\references{
Boer, V. M., de Winde, J. H., Pronk, J. T. and Piper, M. D. W. (2003) The genome-wide transcriptional responses of \emph{Saccharomyces cerevisiae} grown on glucose in aerobic chemostat cultures limited for carbon, nitrogen, phosphorus, or sulfur. \emph{J. Biol. Chem.} \bold{278}, 3265--3274. \doi{10.1074/jbc.M209759200}

Ghaemmaghami, S., Huh, W., Bower, K., Howson, R. W., Belle, A., Dephoure, N., O'Shea, E. K. and Weissman, J. S. (2003) Global analysis of protein expression in yeast. \emph{Nature} \bold{425}(6959), 737--741. \doi{10.1038/nature02046}

Huh, W. K., Falvo, J. V., Gerke, L. C., Carroll, A. S., Howson, R. W., Weissman, J. S. and O'Shea, E. K. (2003) Global analysis of protein localization in budding yeast. \emph{Nature} \bold{425}(6959), 686--691. \doi{10.1038/nature02026}

Tai, S. L., Boer, V. M., Daran-Lapujade, P., Walsh, M. C., de Winde, J. H., Daran, J.-M. and Pronk, J. T. (2005) Two-dimensional transcriptome analysis in chemostat cultures: Combinatorial effects of oxygen availability and macronutrient limitation in \emph{Saccharomyces cerevisiae}. \emph{J. Biol. Chem.} \bold{280}, 437--447. \doi{10.1074/jbc.M410573200}

YeastGFP project. Yeast GFP Fusion Localization Database, http://yeastgfp.ucsf.edu; Current location: \url{https://yeastgfp.yeastgenome.org}
}
