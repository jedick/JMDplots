# orp16S/metaproteome/KTS+17/mkaa.R
# Sum amino acid composition of protein sequences for each sample
# 20221029 jmd version 1
# 20221227 revised to get protein sequences

# REQUIRED FILES:
# SodaLakes_AllCombined_Cluster95ID_V2.fasta.xz
#   - Downloaded from https://ftp.pride.ebi.ac.uk/pride/data/archive/2017/05/PXD006343/SodaLakes_AllCombined_Cluster95ID_V2.fasta
# GEM-(01).pep.xml
#   - Downloaded from https://ftp.pride.ebi.ac.uk/pride/data/archive/2017/05/PXD006343/GEM-(01).pep.xml
# LCM-(01).pep.xml
#   - Downloaded from https://ftp.pride.ebi.ac.uk/pride/data/archive/2017/05/PXD006343/LCM-(01).pep.xml

# Get amino acid composition of proteins in reference database
refaa <- canprot::read_fasta("SodaLakes_AllCombined_Cluster95ID_V2.fasta.xz")
# Remove control character
refaa$protein <- gsub("\r", "", refaa$protein)
# Just use first identifier (e.g., OTU_X_325115443|emb|CBZ50998.1| -> OTU_X_325115443)
refaa$protein <- sapply(strsplit(refaa$protein, "\\|"), "[", 1)

samples <- c("GEM", "LCM")

out <- lapply(samples, function(sample) {

  ## Read source data
  file <- paste0(sample, "-(01).pep.xml.xz")
  print(file)
  dat <- readLines(file)

  ## Extract data from XML file
  # Get line numbers of spectrum query results
  iquery <- grep("<spectrum_query", dat, fixed = TRUE)
  # Get line numbers of first search hits for each query
  isearch <- iquery + 2
  # Get protein IDs
  proteins <- gsub('\\".*', "", gsub('.*protein\\=\\"', "", dat[isearch]))
  # Just use first identifier
  proteins <- sapply(strsplit(proteins, "\\|"), "[", 1)
  # Exclude contaminants
  proteins <- proteins[-grep("CRAP_", proteins)]

  ## Sum amino acid composition of proteins
  iref <- match(proteins, refaa$protein)
  stopifnot(all(!is.na(iref)))
  canprot::sum_aa(refaa[iref, ])

})

## Assemble one data frame
out <- do.call(rbind, out)
out$protein <- "metaproteome-proteins"
out$organism <- samples
out$ref <- "KTS+17"
write.csv(out, "KTS+17_aa.csv", row.names = FALSE, quote = FALSE)
