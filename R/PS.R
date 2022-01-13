# JMDplots/PS.R
# get phylostrata for given UniProt IDs
# 20191127 jmd
# 20220113 moved to JMDplots

PS <- function(uniprot, source = "TPPG17") {
  file <- system.file(paste0("extdata/evdevH2O/phylostrata/", source, ".csv.xz"), package = "JMDplots")
  dat <- read.csv(file, as.is = TRUE)
  if(source == "TPPG17") {
    iPS <- match(uniprot, dat$Entry)
    out <- dat$Phylostrata[iPS]
  } else {
    iPS <- match(uniprot, dat$UniProt)
    out <- dat$PS[iPS]
  }
  out
}
