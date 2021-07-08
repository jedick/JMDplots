# jedick/devodata/CBS+17.R
# Make proteins with abundance-weighted amino acid composition from proteomes
# 20210114

# Read abundance data
dat <- read.csv("CBS+17_abundance.csv")
# Read amino acid compositions of proteins
aa <- read.csv("CBS+17_aa.csv")

# Names of time points
tp <- c("e02", "e06", "e12", "e20", "L1", "L2", "L3", "L3c", "p1", "p2", "p3", "p4", "p5", "Ayf", "Aym", "Af", "Am")
# Columns with the proteome data
itp <- match(tp, colnames(dat))

# Initialize data frame with output amino acid compositions
aaout <- thermo()$protein[rep(1, length(itp)), ]
aaout[, 6:25] <- 0
aaout$chains <- 1
aaout$abbrv <- NA
aaout$ref <- "CBS+17"
aaout$organism <- "DROME"
aaout$protein <- tp

# Match UniProt IDs
aa$protein <- sapply(strsplit(aa$protein, "\\|"), tail, 1)
iaa <- match(dat$Entry, aa$protein)
# Reorder aa data frame
aa <- aa[iaa, ]
stopifnot(all(aa$protein == dat$Entry))

# Loop over time points
for(i in seq_along(itp)) {
  values <- dat[, itp[i]]
  aarow <- colSums(aa[, 6:25] * values)
  aaout[i, 6:25] <- aarow
}
write.csv(aaout, "CBS+17_mean_aa.csv", row.names = FALSE, quote = FALSE)
