# getMGnify.R
# Script to download faa files and taxonomy for MGnify genomes
# 20231229

# Download faa files for MGnify genomes
# First change to destination directory
#setwd("genomes")
# Get list of genomes
dat <- read.csv("../MGnify_genomes.csv")
# Loop over accessions
for(i in seq_along(dat$Accession)) {
  # Skip if file already exists
  if(file.exists(paste0(dat$Accession[i], ".faa"))) next
  # Make wget command
  cmd <- paste0("wget https://www.ebi.ac.uk/metagenomics/api/v1/genomes/", dat$Accession[i], "/downloads/", dat$Accession[i], ".faa")
  # Run command
  system(cmd)
}

# Download taxonomy
setwd("taxonomy")
dat <- read.csv("../MGnify_genomes.csv")
# https://community.rstudio.com/t/web-scraping-with-r-and-rvest-when-the-web-page-uses-javascript/128110/3
library(RSelenium)
# Running rsDriver() the first time downloads a lot of stuff
# To continue, start selenium server manually
# java -jar .local/share/binman_seleniumserver/generic/4.0.0-alpha-2/selenium-server-standalone-4.0.0-alpha-2.jar
remDr <- rsDriver(browser='chrome', port=4444L)
browser <- remDr$client
browser$open()

for(accession in dat$Accession) {

  taxonomy_file <- paste0(accession, ".csv")
  # Skip if taxonomy file already exists
  if(file.exists(taxonomy_file)) {
    #print(paste("Skipping", accession, "(already downloaded)"))
    next
  }

  browser$navigate(file.path("https://www.ebi.ac.uk/metagenomics/genomes", accession))
  # Wait for the page to load
  Sys.sleep(1)
  pagesource <- browser$getPageSource()
  taxonomy_text <- strsplit(strsplit(pagesource[[1]], "Taxonomic lineage:</b> ", fixed = TRUE)[[1]][2], "</p>", fixed = TRUE)[[1]][1]
  # Skip if the taxonomy wasn't loaded (might need to re-run the script or increase sleep)
  if(is.na(taxonomy_text)) {
    print(paste("Skipping", accession, "(page load incomplete?)"))
    next
  }
  taxonomy_names <- strsplit(taxonomy_text, " &gt; ", fixed = TRUE)[[1]]
  taxonomy <- data.frame(
    genome = accession,
    kingdom = taxonomy_names[1],
    phylum = taxonomy_names[2],
    class = taxonomy_names[3],
    order = taxonomy_names[4],
    family = taxonomy_names[5],
    genus = taxonomy_names[6],
    species = taxonomy_names[7]
  )
  print(paste("Downloaded", accession))
  write.csv(taxonomy, taxonomy_file, row.names = FALSE)

}
