# sars16S/metaproteome/JZW+22/mkaa.R
# Calculate amino acid composition from proteins
# 20230207 jmd

# REQUIRED FILE:
# proteinGroups.txt
#   - Downloaded from https://ftp.pride.ebi.ac.uk/pride/data/archive/2022/11/PXD026727/results_Total.zip
#   - Extracted proteinGroups.txt from ZIP file

# List faa files dowloaded from HOMD
faafiles <- dir("/home/sequence/HOMD/JZW+22_PROKKA")
# Get organism IDs from file names
orgids <- sapply(strsplit(faafiles, "\\."), "[", 1)

# Read protein groups table of PXD026727
dat <- read.csv("proteinGroups.txt", sep = "\t")
# Keep proteins with HOMD IDs
dat <- dat[grep("^SEQF", dat$Majority.protein.IDs), ]
IDs <- dat$Majority.protein.IDs

# Get the first listed protein
prots <- sapply(strsplit(IDs, ";"), "[", 1)
# List the organism IDs
orgs <- substr(prots, 1, 8)
# There are some organism IDs that arent't available in HOMD
nothave <- which(!orgs %in% orgids)
if(length(nothave) > 0) {
  # Use the second Majority protein ID for these
  prots[nothave] <- strsplit(IDs[nothave], ";")[[1]][2]
  orgs <- substr(prots, 1, 8)
}
# Check that all organisms are available
stopifnot(all(orgs %in% orgids))

# Loop over proteins
aaprot <- lapply(prots, function(prot) {
  org <- substr(prot, 1, 8)
  # Read the faa file
  faafile <- paste0("/home/sequence/HOMD/JZW+22_PROKKA/", org, ".faa.xz")
  aa <- CHNOSZ::read.fasta(faafile)
  iaa <- match(prot, aa$protein)
  # Return the amino acid composition for this protein
  aa[iaa, ]
})
# Make data frame of amino acid composition for each protein
aaprot <- do.call(rbind, aaprot)

# Identify columns with LFQ intensity for each sample
icol <- grep("LFQ.intensity", colnames(dat))
# Loop over samples
aa <- lapply(icol, function(i) {
  aa <- aasum(aaprot, abundance = dat[, i])
  # Set 'chains' to number of proteins with non-zero LFQ intensity
  aa$chains <- sum(dat[, i] > 0)
  aa
})
# Make data frame and add sample names
aa <- do.call(rbind, aa)
aa$protein <- "HOMD_PROKKA"
aa$organism <- gsub("LFQ.intensity.", "", colnames(dat[icol]))
aa$ref <- "JZW+22"
# Write result
write.csv(aa, "JZW+22_aa.csv", row.names = FALSE, quote = FALSE)
