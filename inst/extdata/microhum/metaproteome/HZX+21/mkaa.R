# microhum/metaproteome/HZX+21/mkaa.R
# Sum amino acid composition of bacterial UniProt sequences for each sample
# 20220830

# REQUIRED FILES:
# - Table_S2.1.csv
#   - From preprint on ResearchSquare https://doi.org/10.21203/rs.3.rs-208797/v1
#   - Download https://assets.researchsquare.com/files/rs-208797/v1/1bc449eb8d5174ad6b7f414c.rar
#   - Extract Table S2. Global metaproteome.xlsx
#   - Export selected columns from Sheet S2.1 to CSV file
# The exported columns are:
# Protein.Group Accession Ctrl16 Ctrl17 Ctrl18 Ctrl19 Ctrl22 Ctrl31 Ctrl33 Ctrl34 Ctrl36 Ctrl37 Ctrl38
# Ctrl42 Ctrl60 Ctrl61 Ctrl69 Ctrl90 Ctrl95 Ctrl96 Ctrl113 Ctrl114 P1_0322 P1_0327 P2_0324 P3_0324
# P3_0401 P4_0324 P4_0329 P4_0401 P5_0422 P5_0402 P5_0328 P5_0413 P5_0404 P5_0419 P5_0325 P6_0324
# P6_0329 P6_0404 P6_0419 P7_0324 P7_0329 P7_0401 P8_0324 P8_0329 P9_0324 P9_0329 P9_0401 P10_0322
# P11_0322 P11_0404 P11_0408 P11_0419 P11_0422 P11_0501 P11_0505 P12_0327 P12_0404 P12_0408 P12_0423
# P12_0426 P12_0428 P12_0501 P12_0505 P13_0330 P13_0404 P13_0408 P13_0419 P13_0426 P13_0428 P13_0501
# P13_0503 P13_0505 QC1 QC2 QC3 QC4 QC5

# 20230212 NOTE:
# The metaproteomics data are also available on iProX (this file was not used in this analysis):
# https://download.iprox.cn/IPX0002453000/IPX0002453001/Metaproteomics.csv

# uniprotIDs.txt
#   - File created with this code
if(FALSE) {
  dat <- read.csv("Table_S2.1.csv")
  # Get accessions starting with tr| or sp|
  prefix <- sapply(strsplit(dat$Accession, "\\|"), "[", 1)
  dat <- dat[prefix %in% c("tr", "sp"), ]
  # Get UniProt IDs
  uniprotIDs <- sapply(strsplit(dat$Accession, "\\|"), "[", 2)
  # Write UniProt IDs for mapping
  writeLines(uniprotIDs, "uniprotIDs.txt")
}

# uniprot.tsv, uniprot.fasta, uniprot_bacteria.fasta
#   - Mapped uniqueIDs.txt on https://www.uniprot.org/id-mapping (UniProtKB_AC-ID -> UniProtKB)
#   - Download TSV -> save as uniprot.tsv
#   - Download FASTA (canonical sequences) -> save as uniprot.fasta
#   - Use ID mapping tool to filter by taxonomy: bacteria
#   - Download FASTA (canonical sequences) -> save as uniprot_bacteria.fasta

for(taxonomy in c("", "_bacteria")) {

  dat <- read.csv("Table_S2.1.csv")
  uniprotIDs <- sapply(strsplit(dat$Accession, "\\|"), "[", 2)

  # Get amino acid composition from FASTA file
  uniprot <- canprot::read.fasta(paste0("uniprot", taxonomy, ".fasta"))
  # Remove "tr|" and "sp|" prefixes
  uniprot$protein <- sapply(strsplit(uniprot$protein, "\\|"), "[", 2)
  # Read mapping table
  map <- read.csv("uniprot.tsv", sep = "\t")

  # Identify mapped IDs
  imap <- match(uniprotIDs, map$From)
  Entry <- map$Entry[imap]
  iuni <- match(Entry, uniprot$protein)
  uniprot <- uniprot[iuni, ]
  # Take out NA mappings
  ina <- is.na(iuni)
  dat <- dat[!ina, ]
  uniprot <- uniprot[!ina, ]
  # Check that accessions are correct
  accession <- sapply(strsplit(dat$Accession, "\\|"), "[", 2)
  # Mapping has no effect on 99% of accessions
  stopifnot(sum(accession == uniprot$protein) / length(accession) > 0.99)

  # Get amino acid composition for each sample
  samples <- colnames(dat)[3:74]
  aa <- lapply(samples, function(sample) {
    # Multiply amino acid composition by spectral counts
    aa <- uniprot[, 5:25] * dat[, sample]
    aa <- cbind(uniprot[, 1:4], aa)
    canprot::aasum(aa)
  })

  aa <- do.call(rbind, aa)
  aa$protein <- "metaproteome"
  aa$organism <- samples
  aa$ref <- "HZX+21"
  aa$abbrv <- NA
  write.csv(aa, paste0("HZX+21", taxonomy, "_aa.csv"), row.names = FALSE, quote = FALSE)
}
