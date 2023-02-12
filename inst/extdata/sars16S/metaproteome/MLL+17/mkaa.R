# sars16S/metaproteome/MLL+17/mkaa.R
# Calculate amino acid compositions from metaproteome data and UniProt sequences
# 20220829

# REQUIRED FILES:
# Maier_Proteomics_Analysis.csv
#   - Downloaded from https://zenodo.org/record/838741/files/Maier_Proteomics_Analysis.xlsx?download=1
#   - Merged two rows of header text
#   - Moved UniProt ID to first column
#   - Export selected columns (UniProtID and 24 runs) as CSV file
# uniparc.tsv and uniparc.fasta
#   - Created from https://www.uniprot.org/id-mapping on 2022-08-29
#   - Mapped unique_IDs.txt with (From database: UniProtKB AC/ID) and (To database: UniParc)
#   - Choose Download -> Format: TSV with default options -> save file as uniparc_in.tsv
#   - Choose Download -> Format: FASTA with Compressed: No -> save file as uniparc.fasta

# Additional processing for uniparc_in.tsv to make uniparc.tsv
#   - Save selected columns and first UniProt ID:
# dat <- read.csv("uniparc_in.tsv", sep = "\t")
# dat$UniProtKB <- sapply(strsplit(dat$UniProtKB, ";"), "[", 1)
# dat <- dat[, c("From", "Entry", "UniProtKB")]
# write.table(dat, "uniparc.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# INTERMEDIATE FILE:
# unique_IDs.txt
#   - List of unique UniProt IDs without organism strings
#   - Created with these commands:
# dat <- read.csv("Maier_Proteomics_Analysis.csv", na.strings="-")
# IDs <- na.omit(unique(sapply(strsplit(dat$UniProt.ID, "_"), "[", 1)))
# writeLines(IDs, "unique_IDs.txt")

# Read data
dat <- read.csv("Maier_Proteomics_Analysis.csv")
uniparc <- CHNOSZ::read.fasta("uniparc.fasta")
# Normalize protein IDs
dat$UniProt.ID <- sapply(strsplit(dat$UniProt.ID, "_"), "[", 1)
#uniparc$protein <- sapply(strsplit(uniparc$protein, "\\|"), "[", 2)
# Get ID mapping results
# For mapping, use UniParc instead of UniProtKB to get obsolete sequences
map <- read.csv("uniparc.tsv", sep = "\t")
imap <- match(dat$UniProt.ID, map$From)
Entry <- map$Entry[imap]
# Match protein IDs
iuni <- match(Entry, uniparc$protein)
uniparc <- uniparc[iuni, ]
# Put in UniProt IDs for checking
UniProtKB <- sapply(strsplit(map$UniProtKB, ";"), "[", 1)
UniProtKB <- UniProtKB[imap]
uniparc$abbrv <- UniProtKB
# Exclude NA IDs (gene names and "-")
ina <- is.na(iuni)
dat <- dat[!ina, ]
uniparc <- uniparc[!ina, ]
# Make sure most of the IDs are identical
# (Some are not because of different first UniProt IDs in UniParc mapping)
stopifnot(sum(dat$UniProt.ID == uniparc$abbrv) / nrow(dat) > 0.6)

# Get amino acid composition for each sample
aa <- lapply(2:ncol(dat), function(icol) {
  aacomp <- colSums(uniparc[5:25] * dat[, icol])
  # Build data frame
  cbind(protein = "metaproteome", organism = colnames(dat)[icol], ref = "MLL+17", abbrv = NA, as.data.frame(t(aacomp)))
})

aa <- do.call(rbind, aa)
# Check that we did matrix multiplication correctly
stopifnot(all(colSums(dat[, 2:ncol(dat)]) == aa$chains))
aa$organism <- gsub("^X", "", aa$organism)
write.csv(aa, "MLL+17_aa.csv", row.names = FALSE, quote = FALSE)
