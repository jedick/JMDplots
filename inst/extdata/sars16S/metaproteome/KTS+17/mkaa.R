# chem16S/metaproteome/KTS+17/protein/mkaa.R
# Sum amino acid composition of protein sequences for each sample
# 20221029 jmd version 1
# 20221227 revised to get protein sequences

# Get amino acid composition of proteins in reference database
refaa <- CHNOSZ::read.fasta("SodaLakes_AllCombined_Cluster95ID_V2.fasta.xz")
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
  CHNOSZ::aasum(refaa[iref, ])

})

## Assemble one data frame
out <- do.call(rbind, out)
out$protein <- "metaproteome-proteins"
out$organism <- samples
out$ref <- "KTS+17"
write.csv(out, "KTS+17_aa.csv", row.names = FALSE, quote = FALSE)
