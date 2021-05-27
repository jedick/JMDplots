# comp16S/process.R

# File preparation functions
# mkRDP("KGP+12")     # Combine output of RDP Classifier into one file --> RDP/KGP+12.tab
# mkAA()              # Use RefSeq sequences to make amino acid composition of taxonomic groups at genus and higher ranks --> groupAA.csv
# mkmetrics()         # Calculate compositional metrics (nH2O, ZC) for each RefSeq group --> RefSeq_metrics.csv

# Combine output files of RDP Classifier into a single CSV file for each study 20200903
# Use RDP Classifier merge-count command 20200910
# Read study metadata to combine samples from different sources 20200921
mkRDP <- function(study) {
  # Get metadata for this study
  mdat <- getmdat(study, dropNA = FALSE)
  # Change to temporary directory
  olddir <- setwd(tempdir())
  on.exit(setwd(olddir))
  # Remove any existing .txt files
  file.remove(Sys.glob("*.txt"))

  # Loop over runs
  for(i in 1:nrow(mdat)) {
    # The RDP result file
    RDPfile <- file.path("/home/comp16S", mdat$study[i], "RDP", paste0(mdat$Run[i], ".txt"))
    print(RDPfile)
    # Copy the RDP result file to here (temporary directory)
    file.copy(RDPfile, ".")
  }

  # Now list all the RDP result files
  files <- dir(pattern = "txt")
  # Name of output file
  outfile <- file.path(olddir, "RDP", paste0(study, ".tab"))
  # Remove output file in case it exists
  if(file.exists(outfile)) file.remove(outfile)
  # RDP jar file
  RDPjar <- "/opt/rdp_classifier_2.13/dist/classifier.jar"
  # Run merge-count command
  print(cmd <- paste("java -jar", RDPjar, "merge-count", outfile, paste(files, collapse = " ")))
  system(cmd)
}

# Calculate average amino acid composition of genus and higher ranks in RefSeq 20200911
# Note: column names are chosen to be compatible with 'protein' data in CHNOSZ
mkAA <- function(ranks = c("genus", "family", "order", "class", "phylum", "superkingdom")) {

  # Read RefSeq amino acid compositions and taxon names
  refseq <- read.csv(system.file("extdata/refseq/protein_refseq.csv.xz", package = "JMDplots"), as.is = TRUE)
  taxa <- read.csv(system.file("extdata/refseq/taxid_names.csv.xz", package = "JMDplots"), as.is = TRUE)
  # Make sure the data tables have consistent taxids
  stopifnot(all(refseq$organism == taxa$taxid))

  # Make a list to hold the output
  out <- vector("list", length(ranks))
  names(out) <- ranks

  # Loop over ranks
  for(rank in ranks) {

    # Find the column corresponding to this rank
    icol <- grep(rank, colnames(taxa))

    # Get names of all taxa at this rank
    names <- na.omit(unique(taxa[, icol]))
    print(paste(rank, length(names)))

    # Create blank amino acid data frame
    AAtmp <- structure(list(
      protein = NA, organism = NA, ref = NA, abbrv = NA, chains = NA,
      Ala = NA, Cys = NA, Asp = NA, Glu = NA, Phe = NA,
      Gly = NA, His = NA, Ile = NA, Lys = NA, Leu = NA,
      Met = NA, Asn = NA, Pro = NA, Gln = NA, Arg = NA,
      Ser = NA, Thr = NA, Val = NA, Trp = NA, Tyr = NA), row.names = 1L, class = "data.frame")
    AA <- AAtmp[rep(1, length(names)), ]
    AA$protein <- rank
    AA$organism <- names
    rownames(AA) <- 1:length(names)

    # Loop over names
    for(i in 1:length(names)) {
      # Find all RefSeq taxa that have this taxon name
      taxon <- names[i]
      istax <- taxa[, icol] == taxon
      istax[is.na(istax)] <- FALSE
      if(any(istax)) {
        # Sum the number of sequences ("chains" column) and amino acid composition for this taxon
        sumAA <- colSums(refseq[istax, 5:25])
        AA[i, 5:25] <- sumAA
        # Put the number of taxa into the "ref" column
        AA$ref[i] <- sum(istax)
        # Put the parent taxon into the "abbrv" column
        parent <- "Root"
        if(icol < ncol(taxa)) {
          parent <- unique(taxa[istax, icol+1])
          if(length(parent) > 1) parent <- paste(parent, collapse = ";")
        }
        AA$abbrv[i] <- parent
      }
    }

    if(rank == "genus") {
      # Add individual taxids that are used for RDP-NCBI mappings 20200922
      addspecies <- refseq$ref %in% c("Candidatus Marinimicrobia bacterium")
      adds <- refseq[addspecies, ]
      adds$organism <- adds$ref
      adds$ref <- 1
      adds$protein <- "species"
      AA <- rbind(adds, AA)
    }

    # Make per-protein average
    AA[, 5:25] <- round(AA[, 5:25] / AA$chains, 1)
    out[[rank]] <- AA

  }

  # Combine the data frames for all ranks
  out <- do.call(rbind, out)
  # Replace NA parent with ""
  out$abbrv[is.na(out$abbrv)] <- ""
  write.csv(out, "groupAA.csv", row.names = FALSE, quote = FALSE)
}

# Compute compositional metrics for each RefSeq group 20200927
mkmetrics <- function() {
  # Read amino acid compositions of all groups
  AA <- read.csv("groupAA.csv", as.is = TRUE)
  # Build output data frame; rename columns for better readability
  out <- data.frame(rank = AA$protein, group = AA$organism, ntaxa = AA$ref, parent = AA$abbrv, nH2O = NA, ZC = NA, nC = NA)
  # Calculate metrics
  out$nH2O <- canprot::H2OAA(AA)
  out$ZC <- canprot::ZCAA(AA)
  out$nC <- CAA(AA)
  write.csv(out, "RefSeq_metrics.csv", row.names = FALSE, quote = FALSE)
}

# Function used in mkmetrics() to calculate number of carbon atoms in amino acid compositions 20200927
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

