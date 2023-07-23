# microhum/metaproteome/JZW+22/download_HOMD.R
# Download selected faa files from Human Oral Microbome Database
# 20230207 jmd

# REQUIRED FILE:
# proteinGroups.txt
#   - Downloaded from https://ftp.pride.ebi.ac.uk/pride/data/archive/2022/11/PXD026727/results_Total.zip
#   - Extracted proteinGroups.txt from ZIP file

# Read IDs from protein groups table
dat <- read.csv("proteinGroups.txt", sep = "\t")
protein.IDs <- dat$Majority.protein.IDs
# Get first IDs
IDs <- sapply(strsplit(protein.IDs, ";"), "[", 1)
# Keep HOMD IDs
IDs <- IDs[grep("^SEQF", IDs)]

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

# The following organisms are not available in HOMD 20230207
# SEQF1058 SEQF3075 SEQF1068 SEQF1063 SEQF2480 SEQF2762 SEQF3069 SEQF2791
