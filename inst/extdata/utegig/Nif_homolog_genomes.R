# Nif_homolog_genomes.R
# 20191014 jmd

# This is the R code used to make Nif_homolog_genomes.csv.
# Input file PCF+18_SuppTable1A.csv is prepared by exporting Supplemental Table 1A
# of Poudel et al. (2018) (doi:10.1128/JB.00757-17) as a CSV file
# (excluding the first line, which has the name and title of the table).

# Modifications for 2022 utegig paper (made after publication of 2020 gradH2O paper):
# 20210527 Updated Nif_homolog_genomes.csv for RefSeq release 206
# 20220531 Updated this script to read protein_refseq.csv.xz from the chem16S package
# 20220629 Moved this script to extdata/utegig, but we still read ../gradH2O/PCF+18_SuppTable1A.csv

# Read the organism names
org <- read.csv("../gradH2O/PCF+18_SuppTable1A.csv", as.is = TRUE)

# Read RefSeq amino acid compositions and taxon names
refseq <- read.csv(system.file("extdata/refseq/protein_refseq.csv.xz", package = "chem16S"), as.is = TRUE)
taxa <- read.csv(system.file("extdata/refseq/taxid_names.csv.xz", package = "chem16S"), as.is = TRUE)
# Make sure the data tables have consistent taxids
stopifnot(all(refseq$organism == taxa$taxid))
# Keep taxids classified at species level 20220104
isspecies <- !is.na(taxa$species)
refseq <- refseq[isspecies, ]
taxa <- taxa[isspecies, ]
# Exclude Bacteria and Archaea species with less than 500 sequences
ivirus <- taxa$superkingdom == "Viruses"
ivirus[is.na(ivirus)] <- FALSE
ilow <- refseq$chains < 500 & !ivirus
refseq <- refseq[!ilow, ]
taxa <- taxa[!ilow, ]

# Find matching organism name in refseq data
findorg <- function(orgname, refseq) {
  i <- NA
  while(is.na(i)) {
    i <- match(orgname, refseq$ref)
    # If we don't get a match, remove characters from end
    orgname <- substr(orgname, 1, nchar(orgname) - 1)
    if(orgname == "") break
  }
  # Return the matching row number
  i
}
# Get all the matches
irefseq <- sapply(org$Genome.name, findorg, refseq = refseq)
# Get the matching organism names and taxids
Refseq.name <- refseq$ref[irefseq]
taxid <- refseq$organism[irefseq]
# Construct output data frame
Nif <- cbind(org[, c(1, 6:9)], Refseq.name = Refseq.name, taxid = taxid)
write.csv(Nif, "Nif_homolog_genomes.csv", row.names = FALSE, quote = c(2, 6))

# Moved from NifProteomes() 20220629
# Drop NA taxids
Nif <- Nif[!is.na(Nif$taxid), ]
# The Nif types, arranged from anaerobic to aerobic
types <- c("Nif-D", "Nif-C", "Nif-B", "Nif-A")
aalist <- list()
for(i in 1:length(types)) {
  type <- types[i]
  # Get the taxids for genomes with this type of Nif
  iNif <- Nif$Type == type
  taxid <- Nif$taxid[iNif]
  # Remove duplicated taxids 20191018
  taxid <- taxid[!duplicated(taxid)]
  # Get the row number in the refseq data frame
  irefseq <- match(taxid, refseq$organism)
  # Get the amino acid composition from refseq
  AAcomp <- refseq[irefseq, ]
  # Store amino acid compositions 20220603
  AAcomp$protein <- type
  aalist[[i]] <- AAcomp
  AA <- do.call(rbind, aalist)
}
write.csv(AA, "Nif_homolog_AA.csv", row.names = FALSE, quote = FALSE)
