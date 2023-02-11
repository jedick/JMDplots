# chem16S/MPB+19/mkaa.R
# Add amino acids of all proteins
# 20220829 jmd v1
# 20221227 Get protein IDs and spectral counts

# Read protein database
refaa <- CHNOSZ::read.fasta("DEAD_Chimneys_20170720_nr_fw_rev_cont.fasta.xz")
refaa$protein <- gsub("\r", "", refaa$protein)

# The PRT block extracted from the mztab.gz file
dat <- read.csv("Mudpit_121204_V2_P6_CH_SM_Chimney_3a_1.pride.mztab_PRT.tab", sep = "\t")
# Exclude contaminants and decoy sequences
dat <- dat[!grepl("^REVERSE", dat$accession), ]
dat <- dat[!grepl("contamination", dat$accession), ]

# Match protein IDs
irefaa <- match(dat$accession, refaa$protein)
aa <- refaa[irefaa, ]
# Multiply amino acid compositions by PSMs
aa[, 5:25] <- aa[, 5:25] * dat$num_psms_ms_run.1.
# Sum amino acid compositions
aa <- CHNOSZ::aasum(aa)
aa$protein <- "metaproteome"
aa$organism <- "StM-R1"
aa$ref <- "MPB+19"
write.csv(aa, "MPB+19_aa.csv", row.names = FALSE, quote = FALSE)
