# chem16S/PMM+18/mkaa.R
# Add amino acids of all proteins
# 20220829 jmd v1
# 20221227 Get protein IDs and spectral counts

# Read protein database
refaa <- CHNOSZ::read.fasta("140826_StM_Chimney39surface53surface_fw_rev_cont.fasta.xz")
refaa$protein <- gsub("\r", "", refaa$protein)

# Files with PRT block extracted from the mztab.gz file
files <- c("Mudpit_121026_V2_P6_AO_SM_Chimney_1a_10.pride.mztab_PRT.tab", "Mudpit_121028_V2_P6_AO_SM_Chimney_2a_1.pride.mztab_PRT.tab")
samples <- c("DffDRAFT", "FcsdDRAFT")

# Loop over files
out <- lapply(files, function(file) {

  dat <- read.csv(file, sep = "\t")
  # Exclude contaminants and decoy sequences
  dat <- dat[!grepl("^REVERSE", dat$accession), ]
  ## This database doesn't have sequences labeled with "contamination"
  #dat <- dat[!grepl("contamination", dat$accession), ]
  dat <- dat[!grepl("HUMAN", dat$accession), ]
  dat <- dat[!grepl("PIG", dat$accession), ]

  # Match protein IDs
  irefaa <- match(dat$accession, refaa$protein)
  aa <- refaa[irefaa, ]
  # Multiply amino acid compositions by PSMs
  aa[, 5:25] <- aa[, 5:25] * dat$num_psms_ms_run.1.
  # Sum amino acid compositions
  CHNOSZ::aasum(aa)

})

out <- do.call(rbind, out)
out$protein <- "metaproteome"
out$organism <- samples
out$ref <- "PMM+18"
write.csv(out, "PMM+18_aa.csv", row.names = FALSE, quote = FALSE)
