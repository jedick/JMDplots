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

# Set 'chem16Sdir' option 20210607
# Adapted from R/src/library/grDevices/zzz.R
.onLoad <- function(libname, pkgname) {
  op.JMDplots <- list(chem16Sdir = system.file("extdata/chem16S", package = "JMDplots"))
  toset <- !(names(op.JMDplots) %in% names(.Options))
  if(any(toset)) options(op.JMDplots[toset])
}
