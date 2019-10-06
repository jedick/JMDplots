# JMDplots/R/zzz.R
# load amino acid compositions of proteins
# 20191006

# read amino acid compositions of E. coli proteins and print some information
ecoli <- readRDS(system.file("/extdata/organisms/ecoli.rds", package = "JMDplots"))

.onAttach <- function(libname,pkgname) {
  packageStartupMessage(paste("JMDplots::ecoli has amino acid compositions of", nrow(ecoli), "proteins"))
}
