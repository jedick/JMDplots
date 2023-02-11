# chem16S/metaproteome/HTZ+17/download.R
# Download pepXML.gz files for selected samples
# 20221222 jmd v1

# Read sample information
dat <- read.csv("HTZ+17_TableS2_selected.csv", check.names = FALSE)
# Read PRIDE file names
readme <- read.csv("README.txt", sep = "\t")
# Use only R1 pepXML.gz files
readme <- readme[grep("pepXML.gz", readme$NAME), ]
readme <- readme[grep("_R1", readme$NAME), ]
# Get URIs
ireadme <- sapply(sapply(dat$"MetaP Pride File Prefix", grep, x = readme$NAME), "[", 1)
URIs <- readme$URI[ireadme]
URIs <- gsub("ftp:", "https:", URIs)
for(URI in URIs) {
  cmd <- paste("wget", URI)
  print(cmd)
  system(cmd)
}

# Append file names to sample information
File <- readme$NAME[ireadme]
dat <- cbind(dat, File = File)
write.csv(dat, "HTZ+17_TableS2_selected_filenames.csv", row.names = FALSE)
