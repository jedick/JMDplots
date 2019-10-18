# Nif_homolog_genomes.R
# 20191014 jmd

# This is the R code used to make Nif_homolog_genomes.csv.
# Input file PCF+18_SuppTable1A.csv is prepared by exporting Supplemental Table 1A
# of Poudel et al. (2018) (doi:10.1128/JB.00757-17) as a CSV file
# (excluding the first line, which has the name and title of the table).

# read the organism names and refseq data
org <- read.csv("PCF+18_SuppTable1A.csv", as.is = TRUE)
RSfile <- system.file("extdata/refseq/protein_refseq.csv.xz", package = "JMDplots")
refseq <- read.csv(RSfile, as.is = TRUE)
# find matching organism name in refseq data
findorg <- function(orgname, refseq) {
  i <- NA
  while(is.na(i)) {
    i <- match(orgname, refseq$ref)
    # if we don't get a match, remove characters from end
    orgname <- substr(orgname, 1, nchar(orgname) - 1)
    if(orgname == "") break
  }
  # return the matching row number
  i
}
# get all the matches
irefseq <- sapply(org$Genome.name, findorg, refseq = refseq)
# get the matching organism names and taxids
Refseq.name <- refseq$ref[irefseq]
taxid <- refseq$organism[irefseq]
# construct output data frame
out <- cbind(org[, c(1, 6:9)], Refseq.name = Refseq.name, taxid = taxid)
write.csv(out, "Nif_homolog_genomes.csv", row.names = FALSE, quote = c(2, 6))
