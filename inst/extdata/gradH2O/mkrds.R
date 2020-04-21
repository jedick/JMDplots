# mkrds.R
# make RDS file from multiple CSV files
# (CSV files have nucleotide or amino acid compositions and are made
#  using mprep() in the ARAST pipeline: https://doi.org/10.5281/zenodo.2314933)
# 20191022

mkrds <- function() {
  # loop over directories
  for(section in c("MGP")) {
    files <- dir(section)
    out <- vector("list", length(files))
    names(out) <- files
    for(file in files) {
      path <- file.path(section, file)
      dat <- read.csv(path)
      out[[file]] <- dat
    }
    saveRDS(out, paste0(section, ".rds"), version = 2, compress = "bzip2")
  }
}

# the reverse of mkrds: export files
rds2file <- function() {
  alldat <- readRDS("MGP.rds")
  files <- names(alldat)
  for(file in files) {
    dat <- alldat[[file]]
    write.csv(dat, file, row.names = FALSE, quote = FALSE)
  }
}
