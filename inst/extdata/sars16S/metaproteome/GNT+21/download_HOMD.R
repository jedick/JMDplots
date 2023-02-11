# chem16S/metaproteome/GNT+21/download_HOMD.R
# Download selected faa files from Human Oral Microbome Database
# 20221105 jmd

# Get bacterial protein IDs from SI Tables of GNT+21
dat <- readLines("TableS9-S12_Majority_Protein_IDs.txt")
## All IDs
#IDs <- unlist(strsplit(dat, ";"))
# First IDs
IDs <- sapply(strsplit(dat, ";"), "[", 1)

# Get unique organism IDs
orgs <- unique(sapply(strsplit(IDs, "_"), "[", 1))
# Skip already downloaded files
orgs <- orgs[!orgs %in% gsub(".faa", "", dir())]
# Download faa files
for(org in orgs) {
  URL <- paste0("https://homd.org/ftp/genomes/PROKKA/current/faa/", org, ".faa")
  cmd <- paste("wget", URL)
  print(cmd)
  system(cmd)
}
