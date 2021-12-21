# JMDplots/devodata
# Get abundance-weighted means and bootstrap confidence intervals
# for chemical metrics (protein length, ZC, nH2O) of developmental proteomes
# 20210708 jmd first version

getCBS17 <- function(metric = "H2O", boot.R = 99) {

  # Read abundance data
  datadir <- system.file("extdata/evdevH2O/devodata", package = "JMDplots")
  dat <- read.csv(file.path(datadir, "CBS+17_abundance.csv.xz"))
  # Read amino acid compositions of proteins
  aa <- read.csv(file.path(datadir, "CBS+17_aa.csv.xz"))

  # Names of time points
  tp <- c("e02", "e06", "e12", "e20", "L1", "L2", "L3", "L3c", "p1", "p2", "p3", "p4", "p5", "Ayf", "Aym", "Af", "Am")
  # Columns with the proteome data
  itp <- match(tp, colnames(dat))

  # Match UniProt IDs
  aa$protein <- sapply(strsplit(aa$protein, "\\|"), tail, 1)
  iaa <- match(dat$Entry, aa$protein)
  # Reorder aa data frame
  aa <- aa[iaa, ]
  stopifnot(all(aa$protein == dat$Entry))
  # Get protein length
  pl <- protein.length(aa)

  # Calculate the metric for each protein
  if(metric == "ZC") X <- ZCAA(aa)
  if(metric == "H2O") X <- H2OAA(aa)
  if(metric == "length") X <- protein.length(aa)

  # https://stackoverflow.com/questions/46231261/bootstrap-weighted-mean-in-r
  samplewmean <- function(d, i, j) {
      d <- d[i]
      w <- j[i]
      weighted.mean(d, w)
  }

  # Run boot() for each timepoint
  boot.list <- lapply(tp, function(TP) {
    weights <- dat[, TP]
    boot::boot(X, samplewmean, R = boot.R, j = weights)
  })

  # Get the means and confidence intervals
  mean.X <- sapply(boot.list, "[[", "t0")
  ci.list <- lapply(boot.list, boot::boot.ci, conf = 0.95, type = "perc")
  low.X <- sapply(lapply(ci.list, "[[", "percent"), "[", 4)
  high.X <- sapply(lapply(ci.list, "[[", "percent"), "[", 5)

  list(mean = mean.X, low = low.X, high = high.X)

}

getFOK21 <- function(metric = "H2O", boot.R = 99) {

  # Read abundance data
  datadir <- system.file("extdata/evdevH2O/devodata", package = "JMDplots")
  dat <- read.csv(file.path(datadir, "FOK+21_abundance.csv.xz"))
  # Read amino acid compositions of proteins
  aa <- read.csv(file.path(datadir, "FOK+21_aa.csv.xz"))

  # Columns with the transcriptome and proteome data
  iT <- grep("T_", colnames(dat))
  iP <- grep("P_", colnames(dat))
  iTP <- c(iT, iP)

  # Replace NA values with 0 (not-measured proteins)
  dat[, 5:20][is.na(dat[, 5:20])] <- 0
  # Match UniProt IDs
  aa$protein <- sapply(strsplit(aa$protein, "\\|"), tail, 1)
  iaa <- match(dat$UniProt, aa$protein)
  # Reorder aa data frame
  aa <- aa[iaa, ]
  stopifnot(all(aa$protein == dat$UniProt))
  # Get protein length
  pl <- protein.length(aa)

  # Calculate the metric for each protein
  if(metric == "ZC") X <- ZCAA(aa)
  if(metric == "H2O") X <- H2OAA(aa)
  if(metric == "length") X <- protein.length(aa)

  # https://stackoverflow.com/questions/46231261/bootstrap-weighted-mean-in-r
  samplewmean <- function(d, i, j) {
      d <- d[i]
      w <- j[i]
      weighted.mean(d, w)
  }

  # Run boot() for each timepoint
  boot.list <- lapply(iTP, function(i) {
    weights <- dat[, i]
    boot::boot(X, samplewmean, R = boot.R, j = weights)
  })

  # Get the means and confidence intervals
  mean.X <- sapply(boot.list, "[[", "t0")
  ci.list <- lapply(boot.list, boot::boot.ci, conf = 0.95, type = "perc")
  low.X <- sapply(lapply(ci.list, "[[", "percent"), "[", 4)
  high.X <- sapply(lapply(ci.list, "[[", "percent"), "[", 5)

  list(mean = mean.X, low = low.X, high = high.X)

}
