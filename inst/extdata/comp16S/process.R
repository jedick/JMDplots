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
  #  # Include a Cyanobacteria "class" (following RDP Classifer; NCBI has NA class for all Cyanobacteria)
  #  names <- c(names, "Cyanobacteria")
  #  # Also include an Acidobacteria "class"
  #  names <- c(names, "Acidobacteria")

    # Create blank amino acid data frame
    AA <- thermo()$protein[rep(1, length(names)), ]
    AA[] <- NA
    AA$protein <- rank
    AA$organism <- names
    rownames(AA) <- 1:length(names)

    # Loop over names
    for(i in 1:length(names)) {
      # Find all RefSeq taxa that have this taxon name
      taxon <- names[i]
      istax <- taxa[, icol] == taxon
  #    if(taxon %in% c("Cyanobacteria", "Acidobacteria")) istax <- taxa$phylum == taxon
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

#    if(rank == "genus") {
#      # Add individual taxids that are used for RDP-NCBI mappings 20200922
#      addspecies <- refseq$ref %in% c("Planktothrix agardhii", "Candidatus Marinimicrobia bacterium")
#      addspecies <- refseq$ref %in% c("Spartobacteria bacterium LR76")
#      adds <- refseq[addspecies, ]
#      adds$organism <- adds$ref
#      adds$ref <- 1
#      adds$protein <- "species"
#      AA <- rbind(adds, AA)
#    }

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
  out$nH2O <- H2OAA(AA)
  out$ZC <- ZCAA(AA)
  out$nC <- CAA(AA)
  write.csv(out, "RefSeq_metrics.csv", row.names = FALSE, quote = FALSE)
}

