# chem16S/OFY+19/OTU2RDP.R
# 20210917

# Script to read Table S2, take most abundant OTUs, assign taxonomic ranks, and save in RDP tab format
dat <- read.csv("SI/Table_S2.csv", check.names = FALSE)
# Keep most abundant OTUs
idat <- dat$"Total read" > 500
dat <- dat[idat, ]
# Get taxonomic names in RefSeq
names <- read.csv(system.file("extdata/refseq/taxid_names.csv.xz", package = "JMDplots"))
names <- reshape2::melt(names, id = "taxid", na.rm = TRUE)

# Find taxonomic rank for lowest-level names in OTU table
OTUname <- sapply(strsplit(dat$"Taxonomy path", ";"), "tail", 1)
inames <- match(OTUname, names$value)
OTUrank <- names$variable[inames]

# Find taxonomic ranks for 2nd-lowest-level names in OTU table
iNA <- is.na(OTUrank)
# For each OTU, reverse the list of taxonomic names and take the second one
OTUname2 <- sapply(sapply(strsplit(dat$"Taxonomy path"[iNA], ";"), "rev"), "[", 2)
inames2 <- match(OTUname2, names$value)
OTUrank2 <- names$variable[inames2]
OTUrank[iNA] <- OTUrank2
OTUname[iNA] <- OTUname2

# Assemble data frame
counts <- dat[, -(1:6)]
out <- data.frame(taxid = dat$OTU, lineage = dat$"Taxonomy path", name = OTUname, rank = OTUrank, counts, check.names = FALSE)
write.table(out, "OFY+19.tab", row.names = FALSE, quote = FALSE, sep = "\t")
