# chem16S/metaproteome/GPM+22/download.R
# 20220906 modified from:

# https://genomicislands.wordpress.com/2013/03/11/downloading-multiple-sequences-from-genbank-quickly-and-easily-using-ape-in-r/
library(ape)

if(FALSE) {
  # Loop over mzid.gz files and get sequence IDs
  gzfiles <- dir(pattern = "*.mzid.gz")
  ID <- lapply(gzfiles, function(file) {
    cmd <- paste("zcat", file, '| grep \\<DBSequence | sed -e s/.*accession\\=\\"//g | sed -e s/\\".*//g')
    system(cmd, intern = TRUE)
  })
  # Write unique IDs to file
  ID <- unique(unlist(ID))
  writeLines(ID, "all_IDs.txt")
}

# Use EDirect to retrieve the sequences from NCBI
# https://www.ncbi.nlm.nih.gov/books/NBK179288/

if(FALSE) {
  ## Retrieve all sequences from NCBI
  # cat all_IDs.txt | efetch -db protein -format fasta > all_proteins.fasta

  # Get amino acid composition of all proteins
  aa <- CHNOSZ::read.fasta("all_proteins.fasta")
  write.csv(aa, "all_proteins_AA.csv", row.names = FALSE, quote = FALSE)
  writeLines(aa$protein, "all_proteins_ID.txt")

  ## Get XML file for sequences retrieved above (includes lineage)
  IDs <- readLines("all_proteins_ID.txt")
  # Split IDs into groups of 100
  # https://statisticsglobe.com/split-vector-into-chunks-in-r
  n <- 100
  ID_list <- split(IDs, ceiling(seq_along(IDs) / n))
  for(i in 1:length(ID_list)) {
    # Write this group of IDs to a file
    writeLines(ID_list[[i]], "cur_ID.txt")
    # Name of output GenBank file padded with zeros
    file <- paste0(sprintf("%04d", i), ".gb")
    # Get GenBank records for this group
    cmd <- paste("cat cur_ID.txt | efetch -db protein -format gb >", file)
    system(cmd)
    # Print group number and number of ORGANISM lines
    cmd <- paste("grep ORGANISM", file, "| wc -l")
    nORG <- system(cmd, intern = TRUE)
    print(paste("Group", i, "has", nORG, "organisms"))
  }

  # Combine GenBank files
  # cat *.db > all_proteins_gb.txt
}

if(FALSE) {
  ## Get only bacterial sequences
  # Read accession list
  IDs <- readLines("all_IDs.txt")
  # Get the bacterial sequences for each ID
  for(i in 1:length(IDs)) {
    cmd <- paste0('esearch -db protein -query "', IDs[i], '" -organism bacteria | efetch -format fasta >> bacterial_proteins.fasta')
    if(i==1) cmd <- paste0('esearch -db protein -query "', IDs[i], '" -organism bacteria | efetch -format fasta > bacterial_proteins.fasta')
    system(cmd)
    print(paste(i, IDs[i]))
  }
}

if(FALSE) {
  # Figure out restart location for interrupted download
  fa <- readLines("bacterial_proteins.fasta")
  iend <- tail(grep("^>", fa), 1)
  id <- gsub(">", "", strsplit(fa[iend], " ")[[1]][1])
  IDs <- readLines("all_IDs.txt")
  iid <- which(IDs == id)
  # Start at the next ID
  iid + 1
}
