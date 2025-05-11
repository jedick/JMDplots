\encoding{UTF-8}
\name{yeast}
\alias{yeast}
\alias{yeast.aa}
\alias{yeastgfp}
\alias{yeast.ZC}
\title{Composition, localization, and abundances of proteins in yeast}
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
The calculation depends on the data file \code{\link[CHNOSZ]{extdata}/organisms/Sce.csv.xz}, which contains the amino acid compositions of the proteins.
The \code{protein} argument should be a vector or a list of vectors of one or more SGD IDs, Open Reading Frame (ORF) or gene names that are found in these files.
The output data frame contains rows with NA compositions for names that are not matched.
If \code{protein} is NULL, all available proteins are returned except those with NA amino acid composition.

\code{yeastgfp} returns the identities and abundances of proteins with the requested subcellular localization(s) (specified in \code{location}) using data from the YeastGFP project that is stored in \code{\link[CHNOSZ]{extdata}/abundance/yeastgfp.csv.xz}.
If \code{exclusive} is \code{FALSE}, the function grabs all proteins that are localized to a compartment even if they are also localized to other compartments.
If \code{exclusive} is \code{TRUE} (the default), only those proteins that are localized exclusively to the requested compartments are identified, unless there are no such proteins, then the non-exclusive localizations are used (applies to the \samp{bud} localization).

\code{yeast.ZC} returns the \code{\link[CHNOSZ:util.formula]{ZC}} of proteins in the subcellular location identified by the one-letter abbreviation.
This abbreviation identifies a \dQuote{Cellular Component} in the Gene Ontology and the associatied proteins (see \code{\link{aoscp}}).
}

\value{
For \code{yeast.aa}, a data frame, or list of data frames, containing the amino acid composition(s) of the specified protein(s) in the format of \code{\link[CHNOSZ]{thermo}$protein}.

For \code{yeastgfp}, a list with elements named \code{protein} (names of proteins) and \code{abundance} (counts or concentrations without any conversion from the units in the data file).
If \code{location} is NULL, \code{yeastgfp} returns the names of all known locations, and if the length of \code{location} is >1, the \code{protein} and \code{abundance} values are lists of the results for each location.
}

\seealso{
\code{\link{scsc}}
}

\examples{
\dontshow{reset()}
# The first few names in UniProt for "aminotransferase yeast"
genes <- c("AATC", "ARO8", "BCA1", "AMPL", "BCA2", "ARO9")
# the corresponding ORF names
ORF <- c("YLR027C", "YGL202W", "YHR208W", "YKL103C", "YJR148W", "YHR137W")
# We only match two of them by gene name, but all by ORF name
aa <- yeast.aa(genes)
aa <- yeast.aa(ORF)
# What are their formulas and average oxidation states of carbon
protein.formula(aa)
ZC(protein.formula(aa))
}

\references{
Ghaemmaghami S, Huh W, Bower K, Howson RW, Belle A, Dephoure N, O'Shea EK and Weissman JS (2003) Global analysis of protein expression in yeast. \emph{Nature} \bold{425}(6959), 737--741. \doi{10.1038/nature02046}

Huh WK, Falvo JV, Gerke LC, Carroll AS, Howson RW, Weissman JS and O'Shea EK (2003) Global analysis of protein localization in budding yeast. \emph{Nature} \bold{425}(6959), 686--691. \doi{10.1038/nature02026}

YeastGFP project. Yeast GFP Fusion Localization Database, http://yeastgfp.ucsf.edu; Current location: \url{https://yeastgfp.yeastgenome.org}
}
