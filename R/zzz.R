# JMDplots/R/zzz.R
# Read data from RDS files 20200828

# The 'JMDplots' environment is made here in open code
# https://stackoverflow.com/questions/41954302/where-to-create-package-environment-variables
JMDplots <- new.env()

.onAttach <- function(libname, pkgname) {
  with(JMDplots, {
    gradox_MGP <- readRDS(system.file("extdata/gradox/MGP.rds", package = "JMDplots"))
    gradox_MGR <- readRDS(system.file("extdata/gradox/MGR.rds", package = "JMDplots"))
    gradox_MGD <- readRDS(system.file("extdata/gradox/MGD.rds", package = "JMDplots"))
    gradH2O_MGP <- readRDS(system.file("extdata/gradH2O/MGP.rds", package = "JMDplots"))
  })
}

# For geo16S, orp16S, and microhum: read taxon_AA.csv once, to make things faster 20221017
taxon_AA <- list(
  RefSeq = read.csv(system.file("extdata/RefSeq/taxon_AA.csv.xz", package = "chem16S")),
  GTDB = read.csv(system.file("extdata/GTDB/taxon_AA.csv.xz", package = "chem16S"))
)

# For pdat_* functions: re-implement protcomp() here (no longer provided by canprot) 20240303
protcomp <- function(uniprot = NULL, aa_file = NULL) {
  aa <- human.aa(uniprot, aa_file, stop.if.missing = TRUE, warn.if.duplicated = TRUE)
  # Return list with UniProt IDs and amino acid composition
  list(uniprot = uniprot, aa = aa)
}
