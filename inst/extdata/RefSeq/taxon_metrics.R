# Compute chemical metrics for each RefSeq group 20200927
taxon_metrics <- function() {
  # Read amino acid compositions of all groups
  AA <- system.file("extdata/RefSeq/taxon_AA.csv.xz", package = "chem16S")
  # Build output data frame; rename columns for better readability
  out <- data.frame(rank = AA$protein, group = AA$organism, ntaxa = AA$ref, parent = AA$abbrv, nH2O = NA, Zc = NA, nC = NA)
  # Calculate metrics
  out$nH2O <- round(canprot::nH2O(AA), 6)
  out$Zc <- round(canprot::Zc(AA), 6)
  out$nC <- round(CAA(AA), 6)
  write.csv(out, "taxon_metrics.csv", row.names = FALSE, quote = FALSE)
}

# Function used in taxon_metrics() to calculate number of carbon atoms in amino acid compositions 20200927
CAA <- function(AAcomp) {
  # the number of carbons of the amino acids
  nC_AA <- c(Ala = 3, Cys = 3, Asp = 4, Glu = 5, Phe = 9, Gly = 2, His = 6, 
    Ile = 6, Lys = 6, Leu = 6, Met = 5, Asn = 4, Pro = 5, Gln = 5, 
    Arg = 6, Ser = 3, Thr = 4, Val = 5, Trp = 11, Tyr = 9)
  # find columns with names for the amino acids
  isAA <- colnames(AAcomp) %in% c("Ala", "Cys", "Asp", "Glu", "Phe", "Gly", "His", "Ile", "Lys", 
    "Leu", "Met", "Asn", "Pro", "Gln", "Arg", "Ser", "Thr", "Val", "Trp", "Tyr")
  iAA <- match(colnames(AAcomp)[isAA], names(nC_AA))
  # calculate the nC for all occurrences of each amino acid
  multC <- t(t(AAcomp[, isAA]) * nC_AA[iAA])
  # calculate the total nC, then the per-residue nC
  nCtot <- rowSums(multC)
  nCtot / rowSums(AAcomp[, isAA])
}
al
