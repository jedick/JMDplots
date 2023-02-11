# chem16S/TWC+22/mkaa.R
# Extract UniProt IDs from protein database accessions
# and sum amino acid compositions for each sample
# 20220830

samples <- c("H1", "H2", "H3", "H4", "H5", "H6", "H7", "H8",
"UC1", "UC2", "UC4", "UC5", "UC6", "UC8", "UC9", "UC11", "UC15", "UC23F")

# Run the following chunk to extract UniProt IDs to files
if(FALSE) {
# Get lists of UniProt IDs
allIDs <- lapply(samples, function(sample) {
  dat <- readLines(paste0(sample, "_proteinID.txt"))
  # Exclude accessions with "." and "_"
  idot <- grepl("\\.", dat)
  isub <- grepl("_", dat)
  dat <- dat[!isub & !idot]
  writelines(dat, paste0(sample, "_uniprotID.txt"))
  dat
})

# Save unique UniProt IDs for mapping on uniprot.org
allIDs <- do.call(c, allIDs)
uniqueIDs <- unique(allIDs)
writeLines(uniqueIDs, "uniqueIDs.txt")
}

# uniprot.fasta and uniprot.tsv are FASTA (canonical sequence) and TSV files obtained
# by mapping uniqueIDs.txt on https://www.uniprot.org/id-mapping (UniProtKB_AC-ID -> UniProtKB)
# uniprot_bacteria.fasta is filtered by taxonomy (in ID mapping tool) to include only sequences from bacteria
uniprot <- CHNOSZ::read.fasta("uniprot_bacteria.fasta")
# Remove "sp|" prefix
uniprot$protein <- sapply(strsplit(uniprot$protein, "\\|"), "[", 2)
# Read mapping table
map <- read.csv("uniprot.tsv", sep = "\t")

# Sum amino acid composition of all proteins for each sample
aa <- lapply(samples, function(sample) {
  ID <- readLines(paste0(sample, "_uniprotID.txt"))
  imap <- match(ID, map$From)
  Entry <- map$Entry[imap]
  iuni <- match(Entry, uniprot$protein)
  # Sum amino acid composition
  CHNOSZ::aasum(uniprot[iuni, ])
})

aa <- do.call(rbind, aa)
aa$protein <- "metaproteome"
aa$organism <- samples
aa$ref <- "TWC+22"
aa$abbrv <- NA
write.csv(aa, "TWC+22_aa.csv", row.names = FALSE, quote = FALSE)
