# Calculate amino acid composition from proteins
# 20221112 jmd

# List faa files dowloaded from HOMD
faafiles <- dir("/home/sequence/HOMD/GNT+21_PROKKA")
# Get organism IDs from file names
orgids <- sapply(strsplit(faafiles, "\\."), "[", 1)

# Use two SI Tables of Granato et al. (2021)
for(table in c("S9", "S11")) {

  # Read list of Majority protein IDs
  file <- paste0("Table_", table, ".csv")
  dat <- read.csv(file)
  # Keep proteins with HOMD IDs
  dat <- dat[grep("^SEQF", dat$Majority.protein.IDs), ]
  # Get Majority protein IDs
  IDs <- dat$Majority.protein.IDs
  # Get the first listed protein
  prots <- sapply(strsplit(IDs, ";"), "[", 1)
  # List the organism IDs
  orgs <- substr(prots, 1, 8)
  # There is one organism ID that isn't available in HOMD (SEQF2791)
  nothave <- which(!orgs %in% orgids)
  if(length(nothave) > 0) {
    # Use the second Majority protein ID for this one
    prots[nothave] <- strsplit(IDs[nothave], ";")[[1]][2]
    orgs <- substr(prots, 1, 8)
  }
  # Check that all organisms are available
  stopifnot(all(orgs %in% orgids))

  # Loop over proteins
  aaprot <- lapply(prots, function(prot) {
    org <- substr(prot, 1, 8)
    # Read the faa file
    faafile <- paste0("/home/sequence/HOMD/GNT+21_PROKKA/", org, ".faa.xz")
    aa <- CHNOSZ::read.fasta(faafile)
    iaa <- match(prot, aa$protein)
    # Return the amino acid composition for this protein
    aa[iaa, ]
  })
  # Make data frame and sum amino acid composition
  aaprot <- do.call(rbind, aaprot)

  # Identify columns with LFQ intensity for each sample
  icol <- grep("LFQ.intensity", colnames(dat))
  # Loop over samples
  aa <- lapply(icol, function(i) {
    # Get LFQ intensity and set NA to 0
    LFQ <- dat[, i]
    LFQ[is.na(LFQ)] <- 0
    # Multiply amino acid composition by LFQ and take the sum
    aa <- aasum(aaprot, abundance = LFQ)
    # Set 'chains' to number of proteins with non-zero LFQ intensity
    aa$chains <- sum(dat[, i] > 0)
    aa
  })

  # Make data frame and add sample names
  aa <- do.call(rbind, aa)
  aa$protein <- "HOMD_PROKKA"
  aa$organism <- gsub("LFQ.intensity.", "", colnames(dat[icol]))
  # Write result
  if(table == "S9") {
    aa$ref <- "GNT+21_cells"
    write.csv(aa, "GNT+21_cells_aa.csv", row.names = FALSE, quote = FALSE)
  }
  if(table == "S11") {
    aa$ref <- "GNT+21_supernatant"
    write.csv(aa, "GNT+21_supernatant_aa.csv", row.names = FALSE, quote = FALSE)
  }
}
