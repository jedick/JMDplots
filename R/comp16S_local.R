# JMDplots/comp16S_local.R
# Separated from comp16S.R to make it easier to define local versions 20210605

# Utility functions
# getmdat("KGP+12")     # Metadata for this study (study, name, Run, BioSample, sample, type, cohort)
# getRDP("KGP+12")      # RDP results for all samples in this study

# Get metadata for a study, appending columns for pch, col, minuend and subtrahend (pairs for difference calculation) 20200914
getmdat <- function(study, dropNA = TRUE) {
  # Read metadata file
  # Remove suffix after underscore 20200929
  studyfile <- gsub("_.*", "", study)
  datadir <- system.file("extdata/comp16S", package = "JMDplots")
  file <- file.path(datadir, "metadata", paste0(studyfile, ".csv"))
  mdat <- read.csv(file, as.is = TRUE, check.names = FALSE)

  if(dropNA) {
    # Exclude samples with NA name (e.g. very low ZC in JHM+16 - outlier?) 20200916
    noname <- is.na(mdat$name)
    if(any(noname)) {
      print(paste0("getmdat [", study, "]: dropping ", sum(noname), " samples with NA name"))
      mdat <- mdat[!is.na(mdat$name), ]
    }
  }
  # Default point symbol: red filled circle for tumor, black open circle for normal, NA otherwise
  pch <- sapply(mdat$type, switch, tumor = 21, normal = 1, NA)
  col <- sapply(mdat$type, switch, tumor = 2, normal = 1, NA)
  # Default differences: tumor minus normal
  isminuend <- mdat$type == "tumor"
  issubtrahend <- mdat$type == "normal"
  subject <- NULL

  ## Identify samples for computing differences in each study

  # Natural environments
  if(study == "BGPF13") { # Heart Lake Geyser Basin, Yellowstone
    pch <- sapply(mdat$cohort, switch, Bacteria = 22, Archaea = 23)
    col <- sapply(mdat$cohort, switch, Bacteria = 5, Archaea = 6)
  }
  if(study == "MPB+17") { # Manus Basin
    type <- mdat$type
    iwater <- type == "water/fluid"
    type[iwater][mdat$T[iwater] > 50] <- "highT"
    type[iwater][mdat$T[iwater] < 50] <- "lowT"
    pch <- sapply(type, switch, lowT = 21, highT = 23, "fauna surface" = 8, "rock/chimney" = 20, NA)
    col <- sapply(type, switch, lowT = 4, highT = 2, "fauna surface" = "yellow4", "rock/chimney" = 1, NA)
  }
  if(study == "SVH+19") { # Black Sea
    pch <- sapply(mdat$type, switch, oxic = 24, suboxic = 20, euxinic = 25, NA)
    col <- sapply(mdat$type, switch, oxic = 4, suboxic = 1, euxinic = 2, NA)
  }
  if(study == "XDZ+17") { # Qarhan Salt Lake
    pch <- sapply(mdat$type, switch, normal = 24, saline = 21)
    col <- sapply(mdat$type, switch, normal = 3, saline = 4)
  }
  if(study == "JHM+16") { # Lake Fryxell microbial mats
    pch <- sapply(mdat$type, switch, oxic = 24, transition = 20, anoxic = 25)
    col <- sapply(mdat$type, switch, oxic = 4, transition = 1, anoxic = 2)
  }
  if(study == "HLA+16") { # Baltic Sea
    #pch <- sapply(mdat$type, switch, Oligohaline = 24, Mesohaline = 20, Marine = 21)
    #col <- sapply(mdat$type, switch, Oligohaline = 3, Mesohaline = 1, Marine = 4)
    # Use near-surface samples 20210601
    mdat <- mdat[mdat$depth <= 20, ]
    type <- rep("moderate", nrow(mdat))
    type[mdat$salinity < 6] <- "low"
    type[mdat$salinity > 20] <- "high"
    pch <- sapply(type, switch, low = 24, moderate = 20, high = 21)
    col <- sapply(type, switch, low = 3, moderate = 1, high = 4)
  }
  if(study == "ZLM+16") { # Tibetan Plateau Lakes
    type <- rep("moderate", nrow(mdat))
    type[mdat$lake %in% c("Keluke", "Qing")] <- "low"
    type[mdat$lake == "Gasikule"] <- "high"
    pch <- sapply(type, switch, low = 24, moderate = 20, high = 21)
    col <- sapply(type, switch, low = 3, moderate = 1, high = 4)
  }
  if(study == "HCW+13") { # Guerrero Negro microbial mat
    pch <- sapply(mdat$zone, switch, A = 24, B = 20, C = 25)
    col <- sapply(mdat$zone, switch, A = 4, B = 1, C = 2)
  }

  # Unconvential oil and gas environments
  if(study == "UKD+18.water") {
    pch <- sapply(mdat$type, switch, "MSA+" = 21, "MSA-" = 1)
    col <- sapply(mdat$type, switch, "MSA+" = 2, "MSA-" = 1)
  }
  if(study == "UKD+18.sediment") {
    pch <- sapply(mdat$type, switch, "MSA+" = 21, "MSA-" = 1)
    col <- sapply(mdat$type, switch, "MSA+" = 2, "MSA-" = 1)
  }
  if(grepl("UKD\\+18.*2012", study)) {
    mdat <- mdat[mdat$year == 2012, ]
    pch <- sapply(mdat$type, switch, "MSA+" = 21, "MSA-" = 1)
    col <- sapply(mdat$type, switch, "MSA+" = 2, "MSA-" = 1)
  }
  if(grepl("UKD\\+18.*2013", study)) {
    mdat <- mdat[mdat$year == 2013, ]
    pch <- sapply(mdat$type, switch, "MSA+" = 21, "MSA-" = 1)
    col <- sapply(mdat$type, switch, "MSA+" = 2, "MSA-" = 1)
  }
  if(grepl("UKD\\+18.*2014", study)) {
    mdat <- mdat[mdat$year == 2014, ]
    pch <- sapply(mdat$type, switch, "MSA+" = 21, "MSA-" = 1)
    col <- sapply(mdat$type, switch, "MSA+" = 2, "MSA-" = 1)
  }
  if(grepl("UKD\\+18.*2015", study)) {
    mdat <- mdat[mdat$year == 2015, ]
    pch <- sapply(mdat$type, switch, "MSA+" = 21, "MSA-" = 1)
    col <- sapply(mdat$type, switch, "MSA+" = 2, "MSA-" = 1)
  }
  if(grepl("UKD\\+18.*2016", study)) {
    mdat <- mdat[mdat$year == 2016, ]
    pch <- sapply(mdat$type, switch, "MSA+" = 21, "MSA-" = 1)
    col <- sapply(mdat$type, switch, "MSA+" = 2, "MSA-" = 1)
  }
  if(study == "CUN+18") {
    pch <- sapply(mdat$type, switch, "UOG+" = 21, "UOG-" = 1, 0)
    col <- sapply(mdat$type, switch, "UOG+" = 2, "UOG-" = 1, 1)
  }
  if(study == "MMA+20") {
    # Exclude AMD streams
    mdat <- mdat[!(grepl("Bark_Camp_Sed", mdat$sample) | grepl("Boone_Sed", mdat$sample) | grepl("Boone_Dup_Sed", mdat$sample)), ]
#    # Include only streams categorized as "low/lowest" or "high/highest" disturbance intensity
#    mdat <- mdat[mdat$sDII < 15 | mdat$sDII > 30, ]
    # Include only streams categorized as "lowest" or "highest" disturbance intensity
    mdat <- mdat[mdat$sDII < 5 | mdat$sDII > 40, ]
    pch <- ifelse(mdat$sDII >= 20, 21, 1)
    col <- ifelse(mdat$sDII >= 20, 2, 1)
  }
  if(study == "MMA+20_spring") {
    # Exclude AMD streams
    mdat <- mdat[!(grepl("Bark_Camp_Sed", mdat$sample) | grepl("Boone_Sed", mdat$sample) | grepl("Boone_Dup_Sed", mdat$sample)), ]
    # Include only streams categorized as "lowest" or "highest" disturbance intensity
    mdat <- mdat[mdat$sDII < 5 | mdat$sDII > 40, ]
    # Include only spring samples
    ispring <- grep("^04", sapply(strsplit(mdat$sample, "_"), tail, 1))
    mdat <- mdat[ispring, ]
    pch <- ifelse(mdat$sDII >= 20, 21, 1)
    col <- ifelse(mdat$sDII >= 20, 2, 1)
  }
  if(study == "MMA+20_fall") {
    # Exclude AMD streams
    mdat <- mdat[!(grepl("Bark_Camp_Sed", mdat$sample) | grepl("Boone_Sed", mdat$sample) | grepl("Boone_Dup_Sed", mdat$sample)), ]
    # Include only streams categorized as "lowest" or "highest" disturbance intensity
    mdat <- mdat[mdat$sDII < 5 | mdat$sDII > 40, ]
    # Include only fall samples
    ifall <- grep("^09", sapply(strsplit(mdat$sample, "_"), tail, 1))
    mdat <- mdat[ifall, ]
    pch <- ifelse(mdat$sDII >= 20, 21, 1)
    col <- ifelse(mdat$sDII >= 20, 2, 1)
  }
  if(study == "CHM+14") {
    # Injected fluids and later produced water
    mdat <- mdat[mdat$day == 0 | mdat$day >= 49, ]
    pch <- ifelse(mdat$day >= 49, 21, 1)
    col <- ifelse(mdat$day >= 49, 2, 1)
  }
  if(study == "HRR+18") {
    mdat <- mdat[mdat$day %in% c(0, 130, 220), ]
    pch <- ifelse(mdat$day > 10, 21, 1)
    col <- ifelse(mdat$day > 10, 2, 1)
  }
  if(study == "ZLF+19") {
    # Source water and flowback day 18
    mdat <- mdat[mdat$day %in% c(-1, 18), ]
    pch <- ifelse(mdat$day >= 1, 21, 1)
    col <- ifelse(mdat$day >= 1, 2, 1)
  }

  # Stratified water datasets
  if(study == "MZG+20") {
    # Identify shallowest and deepest samples from each lake
    newdat <- lapply(unique(mdat$lake), function(lake) {
      ilake <- mdat$lake == lake
      range <- range(mdat$depth[ilake])
      iext <- mdat$depth[ilake] %in% range
      extdat <- mdat[ilake, ][iext, ]
      extdat$type <- ifelse(extdat$depth == range[1], "shallowest", "deepest")
      extdat
    })
    newdat <- do.call(rbind, newdat)
    # Keep all samples, but plot only shallowest and deepest samples
    mdat$type <- newdat$type[match(mdat$Run, newdat$Run)]
    pch <- ifelse(mdat$type == "shallowest", 24, 25)
    col <- ifelse(mdat$type == "shallowest", 4, 2)
  }
  if(study == "HXZ+20") {
    pch <- ifelse(mdat$depth < 70, 24, 25)
    col <- ifelse(mdat$depth < 70, 4, 2)
  }
  if(study == "GBL+15") {
    pch <- ifelse(mdat$depth < 85, 24, 25)
    col <- ifelse(mdat$depth < 85, 4, 2)
  }
  if(study == "BCA+20") {
    type <- rep("transition", nrow(mdat))
    type[mdat$depth < 3] <- "oxic"
    type[mdat$depth > 4] <- "anoxic"
    pch <- sapply(type, switch, oxic = 24, transition = 20, anoxic = 25)
    col <- sapply(type, switch, oxic = 4, transition = 1, anoxic = 2)
  }

  minuend <- subtrahend <- rep(NA, nrow(mdat))
  if(!is.null(subject)) {
    # Enumerate minuend and subtrahend 20200915
    # We can't deal with subjects that appear more than twice
    if(max(table(subject[isminuend | issubtrahend])) > 2) stop("one or more subjects have ambiguous samples for difference")
    # List subjects in each group
    sminuend <- subject[isminuend]
    ssubtrahend <- subject[issubtrahend]
    # Find subjects with paired data
    paired <- intersect(sminuend, ssubtrahend)
    ispaired <- subject %in% paired
    # Identify the minuend and subtrahend in all pairs
    iminuend <- which(ispaired & isminuend)
    isubtrahend <- which(ispaired & issubtrahend)
    # Number the pairs sequentially
    pairs <- seq_along(isubtrahend)
    minuend[iminuend] <- pairs
    # Find the paired sample for each subject
    isubject <- match(subject[isubtrahend], subject[iminuend])
    subtrahend[isubtrahend] <- pairs[isubject]
  } else subject <- rep(NA, nrow(mdat))
  mdat <- cbind(mdat, pch, col, subject, minuend, subtrahend)
  mdat
}

# Get RDP results for all samples in a study 20200912
getRDP <- function(study, cn = FALSE, mdat = NULL, lineage = NULL) {
  # Handle missing arguments
  if(is.null(mdat)) mdat <- getmdat(study)
  # Get the sample Runs (SRR numbers)
  Run <- mdat$Run
  # Remove suffix after underscore 20200929
  studyfile <- gsub("_.*", "", study)
  # Read output of RDP classifer
  datadir <- system.file("extdata/comp16S", package = "JMDplots")
  file <- file.path(datadir, "RDP", paste0(studyfile, ".tab.xz"))
  dat <- read.table(file, sep = "\t", header = TRUE)
  # Get counts for each sample
  icol <- match(Run, colnames(dat))
  # Keep the "lineage", "rank", "name", and counts columns
  dat <- dat[, c(2, 4, 3, icol)]

  if(grepl("NLF+20", study, fixed = TRUE)) {
    # For NLF+20, just use all counts as-is
    out <- dat
    totalcounts <- colSums(out[, -(1:3), drop = FALSE])
  } else {
    # Find the ranks for the "unclassified" (i.e. classified at higher rank than genus)
    iunclass <- grepl("unclassified_", dat$name)
    if(any(iunclass)) {
      unclassname <- gsub("unclassified_", "", dat$name[iunclass])
      # Split the lineage text and find the length of each one
      slineage <- strsplit(dat$lineage[iunclass], ";")
      sllength <- vapply(slineage, length, 1)
      # The rank is a certain number of positions before the end
      irank <- sllength - 2
      unclassrank <- mapply("[", slineage, irank)
      dat$rank[iunclass] <- unclassrank
      dat$name[iunclass] <- unclassname
    }
    # Keep the rows for counts at genus and higher ranks ("unclassified")
    igenus <- dat$rank == "genus"
    out <- dat[igenus | iunclass, ]
    # Get the total counts
    totalcounts <- colSums(out[, -(1:3), drop = FALSE])
    # Get the "rootrank" counts
    rootcounts <- dat[1, -(1:3)]
    # Stop if total counts doesn't equal "rootrank" counts
    stopifnot(all(abs(totalcounts - rootcounts) < 0.1))
  }


  # Keep specified lineage 20200924
  if(!is.null(lineage)) {
    precount <- sum(out[, -(1:3)])
    irow <- grepl(lineage, out$lineage)
    if(!any(irow)) stop(paste("nothing available for lineage =", lineage))
    out <- out[irow, ]
    postcount <- sum(out[, -(1:3)])
    lpercent <- formatC(postcount / precount * 100, 1, format = "f")
    print(paste0("getRDP [", study, "]: keeping ", lineage, " lineage (", lpercent, "%)"))
  }

  # Keep the rows with any counts > 0
  groupcounts <- rowSums(out[, -(1:3), drop = FALSE])
  out <- out[groupcounts > 0, ]

  # Get the number of counts classified at genus level
  igenus <- out$rank == "genus"
  genuscounts <- colSums(out[igenus, -(1:3), drop = FALSE])
  # Print median percent genus counts
  genuspercent <- round(100 * median(genuscounts / totalcounts))
  print(paste0("getRDP [", study, "]: ", genuspercent, "% of classifications at genus level"))

  # Remove classifications at root and domain level (Bacteria and Archaea),
  # and Chlorophyta, Chloroplast and Bacillariophyta 20200922
  RDPgroups <- paste(out$rank, out$name, sep = "_")
  rmgroups <- c("rootrank_Root", "domain_Bacteria", "domain_Archaea", "class_Chloroplast", "family_Chloroplast", "genus_Chlorophyta", "genus_Bacillariophyta")
  isrm <- RDPgroups %in% rmgroups
  if(any(isrm)) {
    irm <- which(isrm)
    rmpercent <- round(rowSums(out[irm, -(1:3)]) / sum(totalcounts) * 100, 1)
    for(i in seq_along(irm)) {
      # Only print message if removed group is >= 0.1% 20200927
      if(rmpercent[i] >= 0.1) print(paste0("getRDP [", study, "]: removing ", RDPgroups[irm[i]], " (", rmpercent[i], "%)"))
    }
    out <- out[!isrm, ] 
  }
  # Recalculate total counts
  totalcounts <- colSums(out[, -(1:3), drop = FALSE])

  if(is.null(lineage)) {
    # Discard samples with < 200 total counts 20201001
    ismall <- totalcounts < 200
    if(any(ismall)) {
      print(paste0("getRDP [", study, "]: discarding ", sum(ismall), " samples with < 200 total counts"))
      out <- out[, c(TRUE, TRUE, TRUE, !ismall)]
    }
  }
  # Recalculate total counts
  totalcounts <- colSums(out[, -(1:3), drop = FALSE])
  # Report the median number of counts 20200917
  # Change this to range 20200924
  print(paste0("getRDP [", study, "]: count range is [", paste(round(range(totalcounts)), collapse = " "), "]"))

  # Adjust counts for 16S rRNA gene copy number 20200927
  if(cn) {
    # Data from rdp_classifier_2.13/src/data/classifier/16srrna/bergeyTrainingTree.xml (20200720)
    bergey <- read.csv("bergeyTrainingTree.csv")
    # Paste together rank and name
    RDPgroups <- paste(out$rank, out$name, sep = "_")
    bergeygroups <- paste(bergey$rank, bergey$name, sep = "_")
    # Get the copy number for these groups
    ibergey <- match(RDPgroups, bergeygroups)
    cpNumber <- bergey$cpNumber[ibergey]
    # Divide RDP Classifier counts by copy number
    # Round to 3 decimal places following output of RDP Classifier with copy-number adjustment
    out[, -(1:3)] <- round(out[, -(1:3)] / cpNumber, 3)
  }

  out
}

