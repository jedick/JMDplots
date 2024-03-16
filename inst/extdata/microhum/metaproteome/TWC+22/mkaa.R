# microhum/metaproteome/TWC+22/mkaa.R
# Extract UniProt IDs from protein database accessions
# and sum amino acid compositions for each sample
# 20220830

# REQUIRED FILES:
# *_proteinID.txt (18 files)
#   - UniProt IDs extracted from each mzid file with the following commands:
# for i in H1 H2 H3 H4 H5 H6 H7 H8 UC1 UC2 UC4 UC5 UC6 UC8 UC9 UC11 UC15 UC23F; do
#  wget https://ftp.pride.ebi.ac.uk/pride/data/archive/2021/06/PXD022433/$i.mzid.gz
#  zcat $i.mzid.gz | grep \<PeptideEvidence\ | grep isDecoy\=\"false\"  | sed -e s/.*\"DBSeq_//g | sed -e s/\".*//g > ${i}_proteinID.txt
# done

# *_uniprotID.txt (18 files)
#   - UniProt IDs extracted from *_proteinID.txt with the following code

samples <- c("H1", "H2", "H3", "H4", "H5", "H6", "H7", "H8",
"UC1", "UC2", "UC4", "UC5", "UC6", "UC8", "UC9", "UC11", "UC15", "UC23F")

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
}

# uniqueIDs.txt
#   - Save unique UniProt IDs for mapping on uniprot.org
if(FALSE) {
allIDs <- do.call(c, allIDs)
uniqueIDs <- unique(allIDs)
writeLines(uniqueIDs, "uniqueIDs.txt")
}

# uniprot_bacteria.fasta, uniprot.tsv
#   - Mapped uniqueIDs.txt on https://www.uniprot.org/id-mapping (UniProtKB_AC-ID -> UniProtKB)
#   - Download TSV -> save as uniprot.tsv
#   - Use ID mapping tool to filter by taxonomy: bacteria
#   - Download FASTA (canonical sequences) -> save as uniprot_bacteria.fasta

uniprot <- canprot::read_fasta("uniprot_bacteria.fasta")
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
  canprot::sum_aa(uniprot[iuni, ])
})

aa <- do.call(rbind, aa)
aa$protein <- "metaproteome"
aa$organism <- samples
aa$ref <- "TWC+22"
aa$abbrv <- NA
write.csv(aa, "TWC+22_aa.csv", row.names = FALSE, quote = FALSE)
