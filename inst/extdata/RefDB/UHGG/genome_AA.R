# Read and sum amino acid composition of all proteins for each genome 20221022
# Script modified for MGnify, based on GTDB/genome_AA.R in chem16S 20231230
genome_AA <- function() {

  # Read taxonomy file
  # ./taxonomy.csv.xz has 2350 genomes with Contamination < 2 & Completeness > 95
  # fullset/taxonomy.csv.xz has all 4744 genomes in UHGG v2.0.1
  taxonomy <- read.csv("taxonomy.csv.xz")
  # Get genome names (accessions)
  genome <- taxonomy$genome
  # List xz-ipped files previously downloaded using commands in getMGnify.R
  files <- file.path("genomes", paste0(genome, ".faa.xz"))

  # Loop over FASTA files (one for each genome)
  ifile <- seq_along(files)
  aa <- lapply(ifile, function(i) {
    # Print progress message
    if(i %% 100 == 0) print(i)
    # Read amino acid composition 
    aa <- suppressMessages(read.fasta(files[i]))
    # Sum amino acid composition
    aasum(aa)
  })
  aa <- do.call(rbind, aa)

  # Put in full genome names (with version suffix .1, .2, etc.)
  aa$organism <- genome
  # Use generic protein name
  aa$protein <- "MGnify"

  # Save result
  write.csv(aa, "genome_AA.csv", row.names = FALSE, quote = FALSE)

}
