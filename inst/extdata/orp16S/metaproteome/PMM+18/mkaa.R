# orp16S/metaproteome/PMM+18/mkaa.R
# Add amino acids of all proteins
# 20220829 jmd v1
# 20221227 Get protein IDs and spectral counts

# REQUIRED FILES:
# 140826_StM_Chimney39surface53surface_fw_rev_cont.fasta.xz
#   - Downloaded from https://ftp.pride.ebi.ac.uk/pride/data/archive/2018/04/PXD009105/140826_StM_Chimney39surface53surface_fw_rev_cont.fasta
# Mudpit_121026_V2_P6_AO_SM_Chimney_1a_10.pride.mztab_PRT.tab
#   - Downloaded from https://ftp.pride.ebi.ac.uk/pride/data/archive/2018/04/PXD009105/generated/Mudpit_121026_V2_P6_AO_SM_Chimney_1a_10.pride.mztab.gz
#   - Save the PRH header line and PRT block extracted from the mztab.gz file
# Mudpit_121028_V2_P6_AO_SM_Chimney_2a_1.pride.mztab_PRT.tab
#   - Downloaded from https://ftp.pride.ebi.ac.uk/pride/data/archive/2018/04/PXD009105/generated/Mudpit_121028_V2_P6_AO_SM_Chimney_2a_1.pride.mztab.gz
#   - Save the PRH header line and PRT block extracted from the mztab.gz file

# 20221028 NOTE:
# Sample descriptions from Bach et al. (2011) http://nbn-resolving.de/urn:nbn:de:gbv:46-00102250-11
# Dff (diffuse): RMR-D = 053-ROV-03 --> "shimmering orange chimney for Eoghan"
# Fcsd (focused): RMR5 = 039-ROV-01 --> "chimney tip, venting gray smoke"

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
