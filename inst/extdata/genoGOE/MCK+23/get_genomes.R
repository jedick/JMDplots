# MCK+23/get_genomes.R
# Get genome IDs from bootstrap files
# 20240913 jmd

# For read.tree()
library(ape)
# Place to store genome and organism name lists
genomes <- list()
names <- list()

# The "bootstraps" directory contains *.bootstraps files from
# https://doi.org/10.6084/m9.figshare.20369313.v1
files <- dir("bootstraps", pattern = "bootstraps", full.names = TRUE)
for(file in files) {
  gene <- substr(basename(file), 1, 4)
  tree <- read.tree(file)
  label <- tree[[1]]$tip.label
  accessions <- gsub("..", "_", substr(label, 1, 16), fixed = TRUE)
  organisms <- gsub("_.$", "", substr(label, 19, 10000))
  print(paste(gene, "has", length(accessions), "genomes and", length(unique(accessions)), "unique genomes"))
  idup <- duplicated(accessions)
  accessions <- accessions[!idup]
  organisms <- organisms[!idup]
  genomes[[gene]] <- accessions
  names[[gene]] <- organisms
}

# List genes in reverse evolutionary order
genes <- c("dmdA", "mddA", "dmsA", "aprA", "aprB", "soxA", "soxB", "soxC", "soxX", "soxY", "soxZ", "dsrA", "dsrB")
# Place to store genome table
outdf <- NULL
for(gene in genes) {
  df <- data.frame(genome = genomes[[gene]], organism = names[[gene]], gene = TRUE)
  colnames(df)[3] <- gene
  if(is.null(outdf)) outdf <- df else
    outdf <- merge(outdf, df, all = TRUE, sort = FALSE)
}

# Loop over genes to group younger ones at top
for(gene in rev(genes)) {
  has <- outdf$genome %in% genomes[[gene]]
  outdf <- rbind(outdf[has, ], outdf[!has, ])
}
write.csv(outdf, "genomes.csv", row.names = FALSE, quote = FALSE, na = "")


# Download genomes
dat <- read.csv("genomes.csv")
for(accession in dat$genome) {
  outfile <- file.path(".", paste0(accession, ".zip"))
  if(!file.exists(outfile)) {
    print(outfile)
    # 'datasets' command from NCBI's command line tools
    cmd <- paste("datasets download genome accession", accession, "--include protein")
    system(cmd)
    file.rename("ncbi_dataset.zip", outfile)
  }
}

# Get protein FASTA
aalist <- lapply(dat$genome, function(accession) {
  infile <- file.path(".", paste0(accession, ".zip"))
  print(infile)
  cmd <- paste("unzip", infile)
  system(cmd)
  faafile <- file.path("ncbi_dataset/data/", accession, "protein.faa")
  outfile <- paste0(accession, ".faa")
  myaa <- NULL
  if(file.exists(faafile)) {
    myaa <- canprot::read_fasta(faafile)
    myaa <- canprot::sum_aa(myaa)
    myaa$protein <- "genome"
    myaa$organism <- accession
  }
  system("rm -rf ncbi_dataset/ README.md")
  myaa
})
aa <- do.call(rbind, aalist)
write.csv(aa, "genome_aa.csv", row.names = FALSE, quote = FALSE)
