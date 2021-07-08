# JMDplots/devodata
# Get abundance-weighted means and bootstrap confidence intervals
# for chemical parameters (protein length, ZC, nH2O) of developmental proteomes
# 20210708 jmd first version

getCBS17 <- function(X = "H2O", boot.R = 99) {

  # Read abundance data
  datadir <- system.file("extdata/devodata", package = "JMDplots")
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

  # Calculate nH2O for each protein
  if(X == "H2O") {
    X <- H2OAA(aa)
    # Remove the +1 H2O from terminal groups (it shouldn't be counted in the weights)
    X <- ( X * pl - 1 ) / pl
  }

  # Names of time points
  tp <- c("e02", "e06", "e12", "e20", "L1", "L2", "L3", "L3c", "p1", "p2", "p3", "p4", "p5", "Ayf", "Aym", "Af", "Am")

  # https://stackoverflow.com/questions/46231261/bootstrap-weighted-mean-in-r
  samplewmean <- function(d, i, j) {
      d <- d[i]
      w <- j[i]
      weighted.mean(d, w)
  }

  # Run boot() for each timepoint
  boot.list <- lapply(tp, function(TP) {
#    # Weight by abundance *and* length
#    abundances <- dat[, TP]
#    weights <- abundances * pl

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
