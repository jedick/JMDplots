# evdevH2O/JWN+21/mkAA.R
# Assemble amino acid compositions from data files of James et al. (2021)
# 20211221

# Data files are in the AACompFiles directory extracted from HomologyDictionaryFiles.zip available at
# https://figshare.com/articles/dataset/Data_from_Universal_and_taxon-specific_trends_in_protein_sequences_as_a_function_of_age/12037281

set <- "pfam_animal_nontrans"

for(set in c("pfam_plant_nontrans", "pfam_plant_trans", "pfam_animal_nontrans", "pfam_animal_trans")) {

  # Read one file to get ages
  dat <- read.csv(paste0("AACompFiles/AAcomp_", set, "_notransform_o.csv"))

  # Create blank amino acid data frame
  AAtmp <- structure(list(
    protein = NA, organism = NA, ref = NA, abbrv = NA, chains = 1,
    Ala = NA, Cys = NA, Asp = NA, Glu = NA, Phe = NA,
    Gly = NA, His = NA, Ile = NA, Lys = NA, Leu = NA,
    Met = NA, Asn = NA, Pro = NA, Gln = NA, Arg = NA,
    Ser = NA, Thr = NA, Val = NA, Trp = NA, Tyr = NA), row.names = 1L, class = "data.frame")
  AA <- AAtmp[rep(1, nrow(dat)), ]
  AA$protein <- dat$Age
  AA$organism <- "plant_trans"
  rownames(AA) <- 1:nrow(dat)

  # Corresponding 3-letter and 1-letter abbreviations
  AA3 <- CHNOSZ::aminoacids(3)
  AA1 <- CHNOSZ::aminoacids(1)
  aa1 <- tolower(AA1)

  # Loop over amino acids
  for(i in 1:20) {
    file <- paste0("AACompFiles/AAcomp_", set, "_notransform_", aa1[i], ".csv")
    dat <- read.csv(file)
    # Check that ages are the same
    stopifnot(all(dat$Age == AA$protein))
    # Put amino acid composition into data frame
    AA[, AA3[i]] <- dat$Metric
  }

  # Round value to reduce size
  AA[, 6:25] <- round(AA[, 6:25], 6)
  # Save the file
  outfile <- paste0(set, "_AA.csv")
  write.csv(AA, outfile, row.names = FALSE, quote = FALSE)

}
