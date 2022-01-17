# geo16S/taxonomy/process.R
# Get names of phyla and genera from RDP and SILVA
# 20220116

## Where to get source files for this script:

## For RDP (see https://sourceforge.net/p/rdp-classifier/news/2020/07/rdp-classifier-213-july-2020-release-note/)
# wget https://sourceforge.net/projects/rdp-classifier/files/RDP_Classifier_TrainingData/RDPClassifier_16S_trainsetNo18_rawtrainingdata.zip
# unzip RDPClassifier_16S_trainsetNo18_rawtrainingdata.zip
# grep "^>" RDPClassifier_16S_trainsetNo18_rawtrainingdata/trainset18_062020_speciesrank.fa | grep -v "Eukaryota" > RDP_headers.txt

## For SILVA
# wget https://www.arb-silva.de/fileadmin/silva_databases/release_138_1/Exports/SILVA_138.1_SSURef_NR99_tax_silva.fasta.gz
# zcat SILVA_138.1_SSURef_NR99_tax_silva.fasta.gz | grep "^>" | grep -v "Eukaryota" | grep -v "uncultured" > SILVA_headers.txt

## Read RDP headers
RDP <- readLines("RDP_headers.txt")
# Make sure all headers have a phylum
stopifnot(all(grepl("phylum__", RDP)))
# Get the phylum
phylumtxt <- sapply(strsplit(sapply(strsplit(RDP, "phylum__"), "[", 2), ";"), "[", 1)
# Get the unique phyla
RDPphyla <- unique(phylumtxt)
# Make sure all headers have a genus
stopifnot(all(grepl("genus__", RDP)))
# Get the genus
genustxt <- sapply(strsplit(sapply(strsplit(RDP, "genus__"), "[", 2), ";"), "[", 1)
# Get the unique genera
RDPgenera <- unique(genustxt)

## Read SILVA headers
SILVA <- readLines("SILVA_headers.txt")
phylatxt <- sapply(strsplit(SILVA, ";"), "[", 2)
SILVAphyla <- unique(phylatxt)
# Genus try 1: from position in taxonomy
genustxt1 <- sapply(strsplit(SILVA, ";"), "[", 6)
SILVAgenera1 <- unique(genustxt1)
# Genus try 2: from first name of species binomial
speciestxt <- sapply(strsplit(SILVA, ";"), "tail", n = 1)
genustxt2 <- sapply(strsplit(speciestxt, " "), "[", 1)
SILVAgenera2 <- unique(genustxt2)
# True genera should be present in both places
SILVAgenera <- intersect(SILVAgenera1, SILVAgenera2)

# Save lists of phyla and genera
writeLines(RDPphyla, "RDPphyla.txt")
writeLines(RDPgenera, "RDPgenera.txt")
writeLines(SILVAphyla, "SILVAphyla.txt")
writeLines(SILVAgenera, "SILVAgenera.txt")
