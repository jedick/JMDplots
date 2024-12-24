# genoGOE/process_GTDB.R
# - Get DNA and protein sequences of marker genes for methanogens
# - Get amino acid compositions of all proteins for each methanogens
# 20240528

# Read methanogen genomes information
mg <- read.csv("methanogen_genomes.csv")

### MARKER GENES: get DNA and protein sequences

# Set sequence type
seqtype <- "faa"
dir.create(seqtype)

# List files in GTDB directory
# (Source: https://data.gtdb.ecogenomic.org/releases/release220/220.0/genomic_files_reps/ar53_marker_genes_reps_r220.tar.gz)
sourcedir <- file.path("~/tmp/ar53_marker_genes_reps_r220", seqtype)
files <- dir(sourcedir)

# Loop over GTDB files
for(file in files) {
  # Read headers
  headers <- read_fasta(file.path(sourcedir, file), type = "headers")
  # Remove ">" from headers
  headers <- gsub(">", "", headers)
  # Find matching headers for methanogen genomes
  iseq <- na.omit(match(mg$Genome, headers))
  # Read lines (headers + sequences)
  lines <- read_fasta(file.path(sourcedir, file), iseq = iseq, type = "lines")
  # Save lines
  writeLines(lines, file.path("methanogen", "marker", seqtype, file))
}

### ALL PROTEINS: get amino acid composition

# Loop over genomes
lapply(mg$Genome, function(genome) {
  # Source of protein.faa.gz files:
  # https://data.gtdb.ecogenomic.org/releases/release220/220.0/genomic_files_reps/gtdb_proteins_aa_reps_r220.tar.gz)
  file <- file.path("tmp/genome/faa", paste0(genome, "_protein.faa.gz"))
  aa <- read_fasta(file)
  aa$organism <- gsub("_protein.faa", "", aa$organism)
  # Get GC ratio from headers and put it into 'abbrv' column
  headers <- read_fasta(file, type = "headers")
  aa$abbrv <- as.numeric(sapply(strsplit(headers, "gc_cont="), "[", 2))
  # Save file
  outfile <- file.path("methanogen", "aa", paste0(genome, "_aa.csv"))
  write.csv(aa, outfile, row.names = FALSE, quote = FALSE)
})
