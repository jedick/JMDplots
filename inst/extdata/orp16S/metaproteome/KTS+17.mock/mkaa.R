# chem16S/metaproteome/KTS+17.mock/protein/mkaa.R
# Sum amino acid composition of protein sequences for each sample
# 20221029 jmd version 1
# 20221222 revised to get protein sequences

# Get amino acid composition of proteins in reference database
refaa <- CHNOSZ::read.fasta("Mock_Comm_RefDB_V3.fasta.xz")
# Remove functional descriptions after ID (as in "AK199_peg.3220\tSEED:fig|6666666.200881.peg.3220\telements")
refaa$protein <- sapply(strsplit(refaa$protein, "\t"), "[", 1)
# Make sure IDs are unique
stopifnot(all(!duplicated(refaa$protein)))

samples <- c("C1", "C2", "C3", "C4", "P1", "P2", "P3", "P4", "U1", "U2", "U3", "U4")

out <- lapply(samples, function(sample) {

  ## Read source data
  # Use files for 260 min 1D-LC-MS/MS run
  # https://ftp.pride.ebi.ac.uk/pride/data/archive/2017/05/PXD006118/
  file <- paste0("Run4and5_", sample, ".pep.xml.xz")
  print(file)
  dat <- readLines(file)

  ## Extract data from XML file
  # Get line numbers of spectrum query results
  iquery <- grep("<spectrum_query", dat, fixed = TRUE)
  # Get line numbers of first search hits for each query
  isearch <- iquery + 2
  # Get protein IDs
  proteins <- gsub('\\".*', "", gsub('.*protein\\=\\"', "", dat[isearch]))
  # Remove functional descriptions after ID (as in "CV_peg.930&#x9;SEED:fig")
  proteins <- sapply(strsplit(proteins, "&#x9;"), "[", 1)

  ## Sum amino acid composition of proteins
  iref <- match(proteins, refaa$protein)
  stopifnot(all(!is.na(iref)))
  aa <- CHNOSZ::aasum(refaa[iref, ])
  aa$protein <- "metaproteome-proteins"
  aa$organism <- sample
  aa$ref <- "KTS+17.mock"
  aa

})

## Put data from both files in one data frame
out <- do.call(rbind, out)
write.csv(out, "KTS+17.mock_aa.csv", row.names = FALSE, quote = FALSE)
