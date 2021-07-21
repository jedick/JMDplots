# JMDplots/devodata/FOK+20.R
# Make proteins with abundance-weighted amino acid composition from transcriptomes and proteomes 20201221

# Read abundance data
dat <- read.csv("FOK+21_abundance.csv.xz")
# Read amino acid compositions of proteins
aa <- read.csv("FOK+21_aa.csv.xz")

# Columns with the transcriptome and proteome data
iT <- grep("T_", colnames(dat))
iP <- grep("P_", colnames(dat))
iTP <- c(iT, iP)

# Initialize data frame with output amino acid compositions
aaout <- thermo()$protein[rep(1, length(iTP)), ]
aaout[, 6:25] <- 0
aaout$chains <- 1
aaout$abbrv <- NA
aaout$ref <- "FOK+21"
aaout$organism[grepl("T_", colnames(dat)[iTP])] <- "transcriptome"
aaout$organism[grepl("P_", colnames(dat)[iTP])] <- "proteome"
aaout$protein <- sapply(strsplit(colnames(dat)[iTP], "_"), "[", 2)

# Replace NA values with 0 (not-measured proteins)
dat[, 5:20][is.na(dat[, 5:20])] <- 0
# Match UniProt IDs
aa$protein <- sapply(strsplit(aa$protein, "\\|"), tail, 1)
iaa <- match(dat$UniProt, aa$protein)
# Reorder aa data frame
aa <- aa[iaa, ]
stopifnot(all(aa$protein == dat$UniProt))

# Loop over transcriptome and proteome measurements
for(i in seq_along(iTP)) {
  values <- dat[, iTP[i]]
  aarow <- colSums(aa[, 6:25] * values)
  aaout[i, 6:25] <- aarow
}
aaout[, 6:25] <- round(aaout[, 6:25], 3)
write.csv(aaout, "FOK+21_mean_aa.csv", row.names = FALSE, quote = FALSE)
