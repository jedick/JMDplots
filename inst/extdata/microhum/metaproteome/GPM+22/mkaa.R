# chem16S/metaproteome/GPM+22/mkaa.R
# Sum amino acid composition of evidence proteins for each sample
# 20220906 jmd

for(type in c("all", "bacterial")) {

  # List *.dat files
  datfiles <- dir(pattern = "*.dat")
  # Get run names from files
  runs <- gsub(".dat", "", datfiles)
  # Get amino acid compositions of all referenced database proteins
  # (created by download.R)
  all_aa <- read.csv("all_proteins_AA.csv")

  # Get compositions of only bacterial proteins 20220907
  if(type == "bacterial") {
    gb <- readLines("all_proteins_gb.txt.xz")
    iVERSION <- grep("^VERSION", gb)
    iORGANISM <- grep("^\ \ ORGANISM", gb)
    # Make sure every accession (VERSION) has an organism line
    stopifnot(all(iORGANISM - iVERSION > 0))
    # Get accessions and first lineage line
    accession <- sapply(strsplit(gb[iVERSION], " "), "[", 6)
    lineage <- gb[iORGANISM + 1]
    kingdom <- gsub("\\.", "", gsub(" ", "", sapply(strsplit(lineage, ";"), "[", 1)))
    # > table(kingdom)
    # kingdom
    # Archaea     Bacteria    Eukaryota Unclassified
    #    1268        94134        13322          272
    isbacteria <- kingdom == "Bacteria"
    all_aa <- all_aa[isbacteria, ]
  }

  # Get amino acid composition for each run
  aa <- lapply(runs, function(run) {
    # Get peptide sequences from mzid.gz file
    cmd <- paste0('zcat ', run, '_p0.05.mzid.gz | grep \\<PeptideEvidence\\ | grep isDecoy\\=\\"false\\" | sed -e s/.*dBSequence_ref\\=\\"//g | sed -e s/\\".*//g')
    DBref <- system(cmd, intern = TRUE)
    IDs <- gsub("DBSeq_1_", "", DBref)
    # Match IDs to data frame of protein amino acid composition
    iall <- match(IDs, all_aa$protein)
    # Print run and percent matches / percent unique
    perc_match <- sum(!is.na(iall)) / length(iall) * 100
    perc_unique <- length(unique(iall)) / length(iall) * 100
    print(paste0(run, " ", round(perc_match, 1), "% matched ", round(perc_unique, 1), "% unique"))
    aa <- canprot::aasum(all_aa[iall, ])
    # Get the patient ID from the dat file
    dat <- readLines(paste0(run, ".dat"))
    pid <- gsub(".mgf", "", strsplit(dat[10], "Feces_")[[1]][2])
    # Change "PO" to "P0"
    pid <- gsub("PO", "P0", pid)
    # Make sample ID
    sample <- paste(pid, run, sep = "-")
    aa$organism <- sample
    aa
  })

  aa <- do.call(rbind, aa)
  aa$protein <- "metaproteome"
  aa$ref <- "GPM+22"

  # Assign PCR negative/positive to samples
  pid <- sapply(strsplit(aa$organism, "_"), "[", 1)
  aa$abbrv <- "positive"
  # Data from Figure S2A of Grenga et al. (2022)
  aa$abbrv[pid %in% c("P23", "P27", "P02", "P24", "P10", "P25", "P31", "P28", "P15", "P32", "P30", "P17")] <- "negative"

  if(type == "all") write.csv(aa, "GPM+22_aa.csv", row.names = FALSE, quote = FALSE)
  if(type == "bacterial") write.csv(aa, "GPM+22_bacteria_aa.csv", row.names = FALSE, quote = FALSE)

}
