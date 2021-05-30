# JMDplots/comp16S.R
# Calculate compositional metrics based on 16S data and RefSeq proteins 20200902
# Revised to include "unclassified" groups in RDP (i.e. classified above genera) 20200911
# Moved to JMDplots 20210416

# Utility functions
# getmdat("KGP+12")     # Metadata for this study (study, name, Run, BioSample, sample, type, cohort)
# getRDP("KGP+12")      # RDP results for all samples in this study
# getmap("KGP+12")      # Map RDP to RefSeq taxonomy (match to rows of groupAA.csv)
# getmetrics("KGP+12")  # Calculate compositional metrics (nH2O, ZC) for each sample

# Plotting functions
# taxacomp()                # nH2O-ZC plot for taxa (default: Bacteria, Archaea) and their children
# plotcomp("KGP+12")        # nH2O-ZC plot for specified study
# diffcomp("KGP+12")        # DnH2O-DZC plot for specified study

#####################
# Utility functions #
#####################

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
    pch <- sapply(mdat$type, switch, marine = 21, hydrothermal = 23, fauna = 8, rock = 20, NA)
    col <- sapply(mdat$type, switch, marine = 4, hydrothermal = 2, fauna = "yellow4", rock = 1, NA)
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
  rmgroups <- c("rootrank_Root", "domain_Bacteria", "domain_Archaea", "genus_Chlorophyta", "class_Chloroplast", "genus_Bacillariophyta")
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

# Map RDP to RefSeq taxonomy 20200912
getmap <- function(study, RDP = NULL, lineage = NULL) {
  # Handle missing arguments
  if(is.null(RDP)) {
    mdat <- getmdat(study)
    RDP <- getRDP(study, mdat = mdat, lineage = lineage)
  }
  # Make group names by combining rank and name
  RDPgroups <- paste(RDP$rank, RDP$name, sep = "_")
  # Calculate group abundances for displaying messages
  groupcounts <- rowSums(RDP[, -(1:3), drop = FALSE])

  NCBIgroups <- vapply(RDPgroups, switch, "",
    # 20200920 Lots of Escherichia in urine [WZZ+18]
    "genus_Escherichia/Shigella" = "genus_Escherichia",
    "phylum_Cyanobacteria/Chloroplast" = "phylum_Cyanobacteria",
    # 20200924 Manus Basin [MPB+17]
    "genus_Marinimicrobia_genera_incertae_sedis" = "species_Candidatus Marinimicrobia bacterium",
    # 20200929 Unclassified Cyanobacteria are just Cyanobacteria
    "class_Cyanobacteria" = "phylum_Cyanobacteria",
    "genus_Spartobacteria_genera_incertae_sedis" = "class_Spartobacteria",
    # 20210502 Processing Guerrero Negro
    "class_Planctomycetacia" = "class_Planctomycetia",
    # 20210526 NCBI taxonomy no longer has an Actinobacteria "class"
    "class_Actinobacteria" = "phylum_Actinobacteria",
    "order_Rhizobiales" = "order_Hyphomicrobiales",
    # 20210530 Acidobacteria
    "genus_Gp1" = "genus_Acidobacterium",
    "genus_Gp6" = "genus_Luteitalea",
    # 20210530 Cyanobacteria
    "genus_GpI" = "genus_Nostoc",
    "genus_GpIIa" = "genus_Synechococcus",
    "genus_GpVI" = "genus_Pseudanabaena",
    "family_Family II" = "family_Synechococcaceae",

    ## NOT USED

    # 20200929 Yellowstone [BGPF13]
    # Not used because there's not much information about this group
    #"genus_Armatimonadetes_gp7" = "phylum_Armatimonadetes",

    # 20210530 Marcellus Shale [CHM+14]
    # https://lpsn.dsmz.de/family/arcobacteraceae
    # Not used because this causes a large low-ZC deviation in Blue Hole 100m sample
    #"family_Arcobacteraceae" = "family_Campylobacteraceae",

  NA_character_)
  iswitch <- !is.na(NCBIgroups)
  if(any(iswitch)) {
    # Print message(s) about switched names and abundance
    from <- names(NCBIgroups)[iswitch]
    to <- NCBIgroups[iswitch]
    switchcounts <- groupcounts[iswitch]
    switchpercent <- round(switchcounts / sum(groupcounts) * 100, 1)
    # Only print message for mappings of groups at least 0.1% abundant 20200927
    if(any(switchpercent >= 0.1)) {
      print(paste0("getmap [", study, "]: using the following RDP --> NCBI mapping(s):"))
      for(i in seq_along(from)) {
        if(switchpercent[i] >= 0.1) message(paste0(from[i], " --> ", to[i], " (", switchpercent[i], "%)"))
      }
    }
    # Actually make the switch!
    RDPgroups[iswitch] <- NCBIgroups[iswitch]
  }

  # Read amino acid composition of all taxonomic groups in RefSeq
  datadir <- system.file("extdata/comp16S", package = "JMDplots")
  AA <- read.csv(file.path(datadir, "groupAA.csv"), as.is = TRUE)
  AAgroups <- paste(AA$protein, AA$organism, sep = "_")
  iAA <- match(RDPgroups, AAgroups)
  # Print summary of missing groups
  naAA <- is.na(iAA)
  nacounts <- groupcounts[naAA]
  napercent <- nacounts / sum(groupcounts) * 100
  naorder <- order(napercent, decreasing = TRUE)
  napercent <- round(napercent[naorder], 2)
  nagroups <- RDPgroups[naAA][naorder]
  if(sum(naAA) > 0) namsg <- paste0("getmap [", study, "]: can't map RDP group ", nagroups[1], " (", napercent[1], "%)")
  if(sum(naAA) > 1) namsg <- paste0("getmap [", study, "]: can't map RDP groups ", nagroups[1], " (", napercent[1], "%), ",
                                    nagroups[2], " (", napercent[2], "%)")
  if(sum(naAA) > 2) namsg <- paste0("getmap [", study, "]: can't map RDP groups ", nagroups[1], " (", napercent[1], "%), ",
                                    nagroups[2], " (", napercent[2], "%), ", nagroups[3], " (", napercent[3], "%)")
  if(sum(naAA) > 3) namsg <- paste0("getmap [", study, "]: can't map RDP groups ", nagroups[1], " (", napercent[1], "%), ",
                                    nagroups[2], " (", napercent[2], "%), ", sum(naAA) - 2, " others (", sum(napercent[-(1:2)]), "%)")
  if(sum(naAA) > 0) print(namsg)
  # Print message about total mapped percent 20200927
  mappedpercent <- formatC(100 - sum(napercent), 1, format = "f")
  print(paste0("getmap [", study, "]: mapped ", mappedpercent, "% of RDP classifications to NCBI taxonomy"))
  # Return result
  iAA
}

# Get compositional metrics for all samples in a study 20200927
getmetrics <- function(study, cn = FALSE, mdat = NULL, RDP = NULL, map = NULL, lineage = NULL, metrics = NULL) {
  # Handle missing arguments
  if(is.null(mdat)) mdat <- getmdat(study)
  if(is.null(RDP)) RDP <- getRDP(study, cn = cn, mdat = mdat, lineage = lineage)
  if(is.null(map)) map <- getmap(study, RDP = RDP, lineage = lineage)
  # Keep metadata only for samples with >= 100 counts 20201001
  mdat <- mdat[mdat$Run %in% colnames(RDP), ]
  # Exclude NA mappings
  RDP <- RDP[!is.na(map), ]
  map <- na.omit(map)
  if(length(map) == 0) stop("no mappings to available RefSeq taxa!")

  # Get compositional metrics of RefSeq groups
  datadir <- system.file("extdata/comp16S", package = "JMDplots")
  if(is.null(metrics)) metrics <- read.csv(file.path(datadir, "RefSeq_metrics.csv"), as.is = TRUE)
  metrics <- metrics[map, ]
  # Make sure the mapping is correct
  equalrank <- RDP$rank == metrics$rank
  # Don't test particular RDP-NCBI mappings that cross ranks
  iclassCyano <- RDP$rank == "class" & RDP$name == "Cyanobacteria"
  igenusSparto <- RDP$rank == "genus" & RDP$name == "Spartobacteria_genera_incertae_sedis"
  iclassActino <- RDP$rank == "class" & RDP$name == "Actinobacteria"
  igenusMarini <- RDP$rank == "genus" & RDP$name == "Marinimicrobia_genera_incertae_sedis"
  igenusGpI <- RDP$rank == "genus" & RDP$name == "GpI"
  igenusGpIIa <- RDP$rank == "genus" & RDP$name == "GpIIa"
  equalrank <- equalrank[!(iclassCyano | igenusSparto | iclassActino | igenusMarini | igenusGpI | igenusGpIIa)]
  stopifnot(all(equalrank))

  # Get classification matrix (rows = taxa, columns = samples)
  RDPmat <- RDP[, -(1:3)]
  # Calculate average nH2O for all samples
  nH2O <- colSums(RDPmat * metrics$nH2O) / colSums(RDPmat)
  # To calculate ZC, we need to compute the sum of charge (ZC * nC) and the sum of carbon atoms
  sumZ <- colSums(RDPmat * metrics$ZC * metrics$nC)
  sumC <- colSums(RDPmat * metrics$nC)
  ZC <- sumZ / sumC

  # Create output data frame
  out <- data.frame(Run = colnames(RDPmat), sample = mdat$sample, nH2O = nH2O, ZC = ZC)
  rownames(out) <- 1:nrow(out)
  out
}

######################
# Plotting functions #
######################

# Make a nH2O-ZC plot for selected taxa and all their children 20200911
taxacomp <- function(which = c("Bacteria", "Archaea"), xlim = NULL, ylim = NULL,
  col = seq_along(taxa), legend.x = "topleft", identify = FALSE, pch = NULL, hline = NULL) {

  # Read compositional metrics of all taxa
  datadir <- system.file("extdata/comp16S", package = "JMDplots")
  metrics <- read.csv(file.path(datadir, "RefSeq_metrics.csv"), as.is = TRUE)
  # Default point symbols
  taxa <- which
  if(is.null(pch)) pch <- rep(21, length(taxa))
  lty <- 1

  # For "majorphyla", get names of phyla with more than 500 representatives
  if(identical(which, "majorphyla")) {
    phyla <- metrics[metrics$rank == "phylum", ]
    phyla <- phyla[phyla$ntaxa > 500, ]
    phyla <- phyla[order(phyla$ntaxa, decreasing = TRUE), ]
    # Move viruses to end 20200926
    ivirus <- phyla$parent == "Viruses"
    phyla <- rbind(phyla[!ivirus, ], phyla[ivirus, ])
    taxa <- phyla$group
    col <- seq_along(taxa)
    if(is.null(ylim)) ylim <- c(-0.9, -0.65)
    pch <- ifelse(phyla$parent == "Bacteria", 21, ifelse(phyla$parent == "Archaea", 24, 23))
  }

  # "majorcellular" is like "majorphyla" but excludes Viruses
  if(identical(which, "majorcellular")) {
    phyla <- metrics[metrics$rank == "phylum" & metrics$parent != "Viruses", ]
    phyla <- phyla[phyla$ntaxa > 60, ]
    phyla <- phyla[order(phyla$ntaxa, decreasing = TRUE), ]
    # Swap Chloroflexi and Crenarchaeota so latter doesn't have same color as Euryarchaeota 20210527
    phyla[14:15, ] <- phyla[15:14, ]
    taxa <- phyla$group
    col <- seq_along(taxa)
    if(is.null(ylim)) ylim <- c(-0.81, -0.68)
    pch <- rep(22, length(taxa))
    pch[1:8] <- 21
    pch[phyla$parent == "Archaea"] <- 24
  }

  # Proteobacteria 20200925
  if(identical(which, "Proteobacteria")) {
    taxa <- c("Alphaproteobacteria", "Betaproteobacteria", "Gammaproteobacteria", "Deltaproteobacteria", "Epsilonproteobacteria", "Zetaproteobacteria",
              "Acidithiobacillia", "Hydrogenophilalia", "Oligoflexia")
    pch <- rep(21:23, 3)
    col <- seq_along(taxa)
    if(is.null(ylim)) ylim <- c(-0.77, -0.71)
  }

  # Acidobacteria 20210529
  if(identical(which, "Acidobacteria")) {
    taxa <- c("Acidobacteriia", "Blastocatellia", "Holophagae", "Thermoanaerobaculia", "Vicinamibacteria")
    pch <- c(21, 22, 23, 21, 22)
    col <- seq_along(taxa)
  }

  # Cyanobacteria 20210529
  if(identical(which, "Cyanobacteria")) {
    taxa <- c("Synechococcales", "Nostocales", "Oscillatoriales", "Chroococcales", 
      "Pleurocapsales", "Chroococcidiopsidales", "Gloeobacterales", 
      "Spirulinales", "Gloeoemargaritales")
    pch <- rep(21:23, 3)
    col <- seq_along(taxa)
    if(is.null(xlim)) xlim <- c(-0.19, -0.13)
    if(is.null(ylim)) ylim <- c(-0.77, -0.71)
  }

  # Use semi-transparent colors for lines 20210518
  lcol <- palette()
  lcol[1] <- "#000000"  # black
  lcol[8] <- "#9e9e9e"  # gray62
  lcol <- paste0(lcol, "80")
  lcol <- rep(lcol, length.out = length(col))

  # Initialize plot
  if(is.null(xlim)) xlim <- c(-0.25, -0.05)
  if(is.null(ylim)) ylim <- c(-0.82, -0.68)
  plot(xlim, ylim, xlab = canprot::cplab$ZC, ylab = canprot::cplab$nH2O, type = "n", xaxs = "i", yaxs = "i")
  # Add horizontal lines to show range of following plot 20200925
  if(!is.null(hline)) abline(h = hline, lty = 2, col = "gray40")
  # Store all values for identify()
  ZC <- nH2O <- numeric()
  names <- character()

  # Loop 1: Get values to plot
  vals <- list()
  for(i in seq_along(taxa)) {
    thisgroup <- taxa[i]
    # Get the compositional metrics for this group
    igroup <- metrics$group == thisgroup
    if(sum(igroup) > 1) warning(paste0("found more than one ", thisgroup, " (", paste(metrics$rank[igroup], collapse = ", "), "); using the first"))
    group <- metrics[which(igroup)[1], ]
    # Get the compositional metrics for all children
    ichildren <- metrics$parent == thisgroup
    children <- metrics[ichildren, ]
    # Store the values to make the plot
    vals[[i]] <- list(group = group, children = children)

    # Keep values for identification
    ZC <- c(ZC, group$ZC, children$ZC)
    nH2O <- c(nH2O, group$nH2O, children$nH2O)
    names <- c(names, group$group, children$group)
  }

  # Loop 2: Plot lines from parents to all children 20200925
  for(i in seq_along(taxa)) {
    group <- vals[[i]]$group
    children <- vals[[i]]$children
    for(j in seq_along(children$ZC)) lines(c(group$ZC, children$ZC[j]), c(group$nH2O, children$nH2O[j]), col = lcol[i], lty = lty)
  }

  # Loop 3: Plot points for parents
  for(i in seq_along(taxa)) {
    group <- vals[[i]]$group
    children <- vals[[i]]$children
    points(group$ZC, group$nH2O, pch = pch[i], cex = 1.5, col = col[i], bg = col[i])
  }

  # Loop 4: Plot points for children
  for(i in seq_along(taxa)) {
    group <- vals[[i]]$group
    children <- vals[[i]]$children
    # Use white outline for black points 20210518
    pt.col <- ifelse((col[i] - 1) %% 8 == 0, "white", 1)
    points(children$ZC, children$nH2O, pch = pch[i], cex = 0.7, col = pt.col, bg = col[i], lwd = 0.5)
    # Label Halobacteria and Nanohaloarchaea 20200930
    thisgroup <- taxa[i]
    if(thisgroup == "Euryarchaeota") {
      ihalo <- match(c("Halobacteria", "Nanohaloarchaea", "Methanococci"), children$group)
      dy <- ifelse(which == "majorcellular", -0.0025, -0.005)
      text(children$ZC[ihalo], children$nH2O[ihalo] + dy, c(1, 2, 3))
    }
    # Label Clostridia 20200930
    if(thisgroup == "Firmicutes" & identical(which, "majorcellular")) {
      iclos <- match("Clostridia", children$group)
      text(children$ZC[iclos], children$nH2O[iclos] + 0.0025, 4)
    }
#    # Label Pisoniviricetes 20210520
#    if(thisgroup == "Pisuviricota") {
#      ipison <- match("Pisoniviricetes", children$group)
#      text(children$ZC[ipison], children$nH2O[ipison] - 0.0025, 5)
#    }
  }

  # Add legend
  len <- length(taxa)
  if(identical(which, "majorcellular")) {
    legend("bottomleft", taxa[1:8], pch = pch[1:8], col = col[1:8], pt.bg = col[1:8], cex = 0.9, bg = "white")
    legend("bottomright", taxa[9:len], pch = pch[9:len], col = col[9:len], pt.bg = col[9:len], cex = 0.9, bg = "white")
  } else if(identical(which, "Proteobacteria")) {
    taxa[taxa == "Epsilonproteobacteria"] <- "Epsilonproteobacteria *"
    legend("topright", taxa[1:6], pch = pch[1:6], col = col[1:6], pt.bg = col[1:6], cex = 0.9, bg = "white")
    legend("bottomleft", taxa[7:len], pch = pch[7:len], col = col[7:len], pt.bg = col[7:len], cex = 0.9, bg = "white")
  } else if(identical(which, "majorphyla")) {
    legend <- c("Cellular", taxa[1:6], "Viruses", taxa[7:11])
    pch <- c(NA, pch[1:6], NA, pch[7:11])
    col <- c(NA, col[1:6], NA, col[7:11])
    legend("bottomleft", legend, text.font = c(2, 1,1,1,1,1,1, 2, 1,1,1,1,1), pch = pch, col = col, pt.bg = col, cex = 0.9, bg = "white")
  } else if(!is.null(legend.x) & !identical(legend.x, NA)) legend(legend.x, taxa, pch = pch, col = seq_along(taxa), pt.bg = seq_along(taxa), cex = 0.9, bg = "white")
  if(identify) identify(ZC, nH2O, names)
}

# Plot compositional metrics for all samples in a study 20200901
plotcomp <- function(study, cn = FALSE, identify = FALSE, title = TRUE, xlim = NULL, ylim = NULL,
  plot.it = TRUE, points = TRUE, lines = FALSE, lineage = NULL, pch.up = 21, pch.down = 1, pval = TRUE) {
  # Get amino acid composition for samples
  mdat <- getmdat(study)
  RDP <- getRDP(study, cn = cn, mdat = mdat, lineage = lineage)
  metrics <- getmetrics(study, mdat = mdat, RDP = RDP, lineage = lineage)
  # Keep metadata only for samples with >= 200 counts 20201006
  mdat <- mdat[mdat$Run %in% metrics$Run, ]
  pch <- mdat$pch
  col <- mdat$col
  # Get nH2O and ZC
  nH2O <- metrics$nH2O
  ZC <- metrics$ZC

  if(plot.it) {
    # Start plot
    if(is.null(xlim)) xlim <- range(ZC)
    if(is.null(ylim)) ylim <- range(nH2O)
    plot(xlim, ylim, xlab = NA, ylab = NA, type = "n")
    if(points) {
      lmlines()
      ifill <- pch > 20
      points(ZC[ifill], nH2O[ifill], pch = pch[ifill], col = 1, bg = col[ifill])
      points(ZC[!ifill], nH2O[!ifill], pch = pch[!ifill], col = col[!ifill])
    }
    if(lines) lines(ZC, nH2O, lty = 3)
    if(isTRUE(title)) title(na.omit(mdat$name)[1], font.main = 1)
    else if(!isFALSE(title)) title(title, font.main = 1)
    # Identify points 20200903
    if(identify) {
      identify(ZC, nH2O, metrics$sample)
      ## Label points with RDP counts 20200919
      #count <- round(colSums(RDP[, -(1:3)]))
      #identify(ZC, nH2O, count)
    }
    xlab <- canprot::cplab$ZC
    ylab <- canprot::cplab$nH2O
  }

  # Calculate mean values and p-values 20201003
  iup <- pch %in% pch.up
  idn <- pch %in% pch.down
  mean <- list()
  p.nH2O <- p.ZC <- NA
  if(!is.null(pch.up) & !is.null(pch.down) & sum(iup) > 0 & sum(idn) > 0) {
    if(sum(idn) > 2 & sum(iup) > 2) {
      p.nH2O <- t.test(nH2O[idn], nH2O[iup])$p.value
      p.ZC <- t.test(ZC[idn], ZC[iup])$p.value
      print(paste("p.ZC =", round(p.ZC, 3), "p.nH2O =", round(p.nH2O, 3)))
    }
    mean <- list(ZC.dn = mean(ZC[idn]), ZC.up = mean(ZC[iup]), nH2O.dn = mean(nH2O[idn]), nH2O.up = mean(nH2O[iup]))
    if(plot.it) {
      points(mean$ZC.dn, mean$nH2O.dn, pch = 8, cex = 2, lwd = 4, col = "white")
      points(mean$ZC.dn, mean$nH2O.dn, pch = 8, cex = 2, lwd = 2, col = 1)
      points(mean$ZC.up, mean$nH2O.up, pch = 8, cex = 2, lwd = 4, col = "white")
      points(mean$ZC.up, mean$nH2O.up, pch = 8, cex = 2, lwd = 2, col = 2)
    }
  }

  # Add axis labels
  if(plot.it) {
    # Make formatted axis labels
    if(is.na(p.ZC) | !pval) xlab <- canprot::cplab$ZC else {
      # Make log10 p-value bold if p-value is less than 0.05
      log10p.ZC <- formatC(log10(p.ZC), 1, format = "f")
      if(p.ZC < 0.05) xlab <- bquote(.(canprot::cplab$ZC[[1]]) ~ "(" * bold(.(log10p.ZC)) * ")")
      else xlab <- bquote(.(canprot::cplab$ZC[[1]]) ~ "(" * .(log10p.ZC) * ")")
    }
    if(is.na(p.nH2O) | !pval) ylab <- canprot::cplab$nH2O else {
      log10p.nH2O <- formatC(log10(p.nH2O), 1, format = "f")
      if(p.nH2O < 0.05) ylab <- bquote(.(canprot::cplab$nH2O[[1]]) ~ "(" * bold(.(log10p.nH2O)) * ")")
      else ylab <- bquote(.(canprot::cplab$nH2O[[1]]) ~ "(" * .(log10p.nH2O) * ")")
    }
    mtext(xlab, side = 1, line = par("mgp")[1], cex = 0.8)
    mtext(ylab, side = 2, line = par("mgp")[1], cex = 0.8)
  }

  invisible(list(study = study, nH2O = nH2O, ZC = ZC, pch = pch, col = col, mean = mean))
}

# function to add convex hulls 20200923
addhull <- function(x, y, basecol, outline = FALSE, ...) {
  i <- chull(x, y)
  r <- as.numeric(col2rgb(basecol))
  if(outline) {
    polygon(x[i], y[i], col = NA, border = basecol, ...)
  } else {
    col <- rgb(r[1], r[2], r[3], 80, maxColorValue=255)
    polygon(x[i], y[i], col = col, border = NA, ...)
  }
}

# Composition-abundance plots for sample groups within taxonomic groups 20210520
groupcomp <- function(study = "XDZ+17", metric = "nH2O", rank = "domain", pch.up = 24, pch.down = 21, minpercent = 2,
                      xlim = NULL, ylim = NULL, xadj = NULL, yadj = NULL, study2 = NA, mdat = NULL, map = NULL, RDP = NULL) {

  # Get metadata, RDP and taxonomy mapping
  if(is.null(mdat)) mdat <- getmdat(study)
  mdat <- mdat[, c("study", "name", "Run", "sample", "pch", "col")]
  # Get data to compare two studies 20210513
  if(!is.na(study2)) {
    mdat2 <- getmdat(study2)[, c("study", "name", "Run", "sample", "pch", "col")]
    mdat$pch <- pch.up
    mdat2$pch <- pch.down
    # Get RDP classification and taxonomy mapping
    RDP <- getRDP(study, mdat = mdat)
    map <- getmap(study, RDP = RDP)
    RDP2 <- getRDP(study2, mdat = mdat2)
    map2 <- getmap(study2, RDP = RDP2)
    # Put together the data from both studies
    mdat <- rbind(mdat, mdat2)
    # Include taxonomy mapping in data frames before merging
    RDP <- cbind(RDP, map = map)
    RDP2 <- cbind(RDP2, map = map2)
    RDP <- merge(RDP, RDP2, all = TRUE)
    # Get out the mappings from the merged data
    map <- RDP$map
    RDP <- RDP[, -grep("map", colnames(RDP))]
    # Replace NA with 0
    RDP[is.na(RDP)] <- 0
  } else {
    # Drop samples with NA pch (we can't include them in the difference)
    mdat <- mdat[!is.na(mdat$pch), ]
    if(is.null(RDP)) RDP <- getRDP(study, mdat = mdat)
    if(is.null(map)) map <- getmap(study, RDP = RDP)
  }
  # Identify samples in up- and down-groups
  iup <- mdat$pch %in% pch.up
  idown <- mdat$pch %in% pch.down
  # Retrieve colors for points
  col.up <- mdat[iup, ]$col[1]
  col.down <- mdat[idown, ]$col[1]
  # Read compositional metrics for faster running
  datadir <- system.file("extdata/comp16S", package = "JMDplots")
  RefSeq_metrics <- read.csv(file.path(datadir, "RefSeq_metrics.csv"), as.is = TRUE)

  # Split the lineage text
  lsplit <- strsplit(RDP$lineage, ";")
  # Find the taxa with the specified rank in the lineage
  irank <- vapply(lsplit, function(x) match(rank, x), 0) - 1
  taxon <- mapply("[", lsplit, irank)
  # Calculate the compositional metrics for each unique taxon
  taxa <- na.omit(unique(taxon))
  Xup <- Xdown <- Pup <- Pdown <- Ptaxa <- numeric()
  Xtaxa <- character()
  for(j in seq_along(taxa)) {
    # Which organisms (by RDP classification) are in this taxon
    itaxon <- taxon == taxa[j]
    itaxon[is.na(itaxon)] <- FALSE
    thisRDP <- RDP[itaxon, ]
    thismap <- map[itaxon]
    # Calculate percent abundance of this taxon in the whole community
    thispercent <- sum(thisRDP[, -(1:3)]) / sum(RDP[, -(1:3)]) * 100
    # Skip low-abundance taxa
    if(thispercent < minpercent) next
    # If we got here, print message about taxon name and abundance
    print(paste0("groupcomp: ", rank, "_", taxa[j], " (", round(thispercent), "%)"))
    # Skip taxa with no available mappings
    if(all(is.na(thismap))) {
      print("            --- no taxonomic mappings available!")  
      next
    }
    # Keep percentages and names of used taxa
    Xtaxa <- c(Xtaxa, taxa[j])
    Ptaxa <- c(Ptaxa, thispercent)
    # Calculate the compositional metrics
    metrics <- getmetrics(study, mdat = mdat, RDP = thisRDP, map = thismap, metrics = RefSeq_metrics)
    # Get selected compositional metric
    if(metric == "nH2O") X <- metrics$nH2O
    if(metric == "ZC") X <- metrics$ZC
    # Calculate median values of compositional metrics
    Xup <- c(Xup, median(X[iup], na.rm = TRUE))
    Xdown <- c(Xdown, median(X[idown], na.rm = TRUE))
    # Calculate percent abundance within the sample groups 20210520
    upall <- RDP[, which(iup) + 3]
    upthis <- thisRDP[, which(iup) + 3]
    Pup <- c(Pup, sum(upthis) / sum(upall) * 100)
    downall <- RDP[, which(idown) + 3]
    downthis <- thisRDP[, which(idown) + 3]
    Pdown <- c(Pdown, sum(downthis) / sum(downall) * 100)
  }

  # Plot the compositions and abundances
  if(is.null(xlim)) xlim <- range(na.omit(c(Pup, Pdown)))
  if(is.null(ylim)) ylim <- range(na.omit(c(Xup, Xdown)))
  if(identical(xlim, c(0, 100))) {
    plot(extendrange(xlim), ylim, type = "n", xlab = "Abundance (%)", ylab = canprot::cplab[[metric]], xaxt = "n")
    axis(1, c(0, 50, 100))
  } else plot(xlim, ylim, type = "n", xlab = "Abundance (%)", ylab = canprot::cplab[[metric]])
  for(k in seq_along(Xtaxa)) {
    # Add points for up- and down- sample groups
    cex <- 1.5
    if(pch.up > 20) points(Pup[k], Xup[k], pch = pch.up, bg = col.up, cex = cex)
    else points(Pup[k], Xup[k], pch = pch.up, col = col.up, cex = cex)
    if(pch.down > 20) points(Pdown[k], Xdown[k], pch = pch.down, bg = col.down, cex = cex)
    else points(Pdown[k], Xdown[k], pch = pch.down, col = col.down, cex = cex)
    # Add arrow connecting the points
    arrows(Pdown[k], Xdown[k], Pup[k], Xup[k], length = 0.1)
    # Pad labels with spaces to offset from points
    label <- paste0("  ", Xtaxa[k], "  ")
    # Get adjustment from arguments if provided
    adj <- c(0, 0.5)
    if(!is.null(xadj)) if(!is.na(xadj[Xtaxa[k]])) adj[1] <- xadj[Xtaxa[k]]
    if(!is.null(yadj)) if(!is.na(yadj[Xtaxa[k]])) adj[2] <- yadj[Xtaxa[k]]
    # Add group name with spaces to offset from points
    text(Pup[k], Xup[k], label, adj = adj)
  }
  # Add lines for total composition of these groups
  OKup <- !is.na(Xup)
  Xup <- Xup[OKup]
  Pup <- Pup[OKup]
  up <- sum(Xup * Pup / sum(Pup))
  OKdown <- !is.na(Xdown)
  Xdown <- Xdown[OKdown]
  Pdown <- Pdown[OKdown]
  down <- sum(Xdown * Pdown / sum(Pdown))
  abline(h = up, lty = 2, lwd = 2, col = col.up)
  abline(h = down, lty = 3, lwd = 2, col = col.down)
  # Calculate and return total percentage of community represented by these taxa
  Ptaxa <- Ptaxa[OKup & OKdown]
  sum(Ptaxa)

}

########################
# Unexported functions #
########################

# Plot differences of nH2O and ZC 20200901
# Use iminuend and isubtrahend to identify sample pairs 20200914
diffcomp <- function(study, cn = FALSE, identify = FALSE, title = TRUE, xlim = NULL, ylim = NULL, plot.it = TRUE) {
  # Get metadata
  mdat <- getmdat(study)
  # Get compositional metrics for samples
  metrics <- getmetrics(study, cn = cn, mdat = mdat)
  # Keep metadata only for samples with >= 200 counts 20201001
  mdat <- mdat[mdat$Run %in% metrics$Run, ]
  nH2O <- metrics$nH2O
  ZC <- metrics$ZC
  if(all(is.na(mdat$minuend))) stop("minuend and subtrahend for differences are not defined")
  # Find pairs of samples for minuend and subtrahend
  pairs <- intersect(na.omit(mdat$minuend), na.omit(mdat$subtrahend))
  isubtrahend <- match(pairs, mdat$subtrahend)
  iminuend <- match(pairs, mdat$minuend)
  # Calculate differences
  DnH2O <- nH2O[iminuend] - nH2O[isubtrahend]
  DZC <- ZC[iminuend] - ZC[isubtrahend]
  # Default symbol is red circle with black outline
  pch <- rep(21, length(pairs))
  col <- rep(2, length(pairs))
  if(study == "HCH+16") {
    # Use different points for cancer and benign disease 20200914
    # This is the 'histology_cat' annotation from Biosample metadata
    col <- sapply(mdat$cohort[iminuend], switch, InvCa = 2, DCIS = 2, BBD_non_atypia = 0, Atypia = 0)
  }
  if(study == "TWB+18") {
    # Red filled circle for cancer, open circle for not cancer
    pch <- sapply(mdat$cohort[iminuend], switch, cancer = 21, "not cancer" = 1)
    col <- sapply(mdat$cohort[iminuend], switch, cancer = 2, "not cancer" = 1)
  }
  if(study == "NLZ+15") {
    # Red circle for carcinoma, blue square for adenoma
    pch <- sapply(mdat$type[iminuend], switch, tumor = 21, polyp = 22)
    col <- sapply(mdat$type[iminuend], switch, tumor = 2, polyp = 4)
  }
  if(study == "TZT+20") {
#    # Red circle for TNBC, red square for TPBC
#    pch <- sapply(mdat$cohort[iminuend], switch, TNBC = 21, TPBC = 22)
    race <- gsub("[ab]", "", sapply(strsplit(mdat$subject, "_"), "[", 1))
    pch <- sapply(race[iminuend], switch, BNH = 21, WNH = 22)
    col <- sapply(race[iminuend], switch, BNH = 2, WNH = 4)
  }
  if(study == "SKB+14_paired") {
    # Filled circle for cancer/normal or CIS/normal, filled square for dysplasia/normal, open circle for healthy normal (right/left)
    pch <- sapply(mdat$cohort[iminuend], switch, cancer = 21, CIS = 21, "pre-cancer" = 22,
                  "cancer duplicate" = 21, "pre-cancer duplicate" = 22, "healthy normal" = 1)
    # Red for cancer, blue for dysplasia, unfilled for healthy normal
    col <- sapply(mdat$cohort[iminuend], switch, "pre-cancer" = 4, "pre-cancer duplicate" = 4, "healthy normal" = 0, 2)
  }

  # Calculate mean difference and p-value 20201001
  # Use col == 2 (red) to get cancer-normal pairs
  iup <- iminuend[col == 2]
  idn <- isubtrahend[col == 2]
  D_mean_nH2O <- mean(nH2O[iup]) - mean(nH2O[idn])
  D_mean_ZC <- mean(ZC[iup]) - mean(ZC[idn])

  if(plot.it) {
    # Make plot
    if(is.null(xlim)) xlim <- range(DZC)
    if(is.null(ylim)) ylim <- range(DnH2O)
    plot(xlim, ylim, xlab = NA, ylab = NA, type = "n")
    points(DZC, DnH2O, pch = pch, col = 1, bg = col)
    abline(v = 0, h = 0, lty = 2, col = "gray60")
    points(D_mean_ZC, D_mean_nH2O, pch = 8, col = "white", lwd = 3.5, cex = 2)
    points(D_mean_ZC, D_mean_nH2O, pch = 8, lwd = 1.5, cex = 2)
    points(D_mean_ZC, D_mean_nH2O, pch = 8, lwd = 0.5, col = 2, cex = 2)
    p.nH2O <- t.test(nH2O[idn], nH2O[iup], paired = TRUE)$p.value
    log10p.nH2O <- formatC(log10(p.nH2O), 1, format = "f")
    p.ZC <- t.test(ZC[idn], ZC[iup], paired = TRUE)$p.value
    log10p.ZC <- formatC(log10(p.ZC), 1, format = "f")
    # Make log10 p-value bold if p-value is less than 0.05
    if(p.ZC < 0.05) xlab <- bquote(.(canprot::cplab$DZC[[1]]) ~ "(" * bold(.(log10p.ZC)) * ")") else xlab <- bquote(.(canprot::cplab$DZC[[1]]) ~ "(" * .(log10p.ZC) * ")")
    if(p.nH2O < 0.05) ylab <- bquote(.(canprot::cplab$DnH2O[[1]]) ~ "(" * bold(.(log10p.nH2O)) * ")") else ylab <- bquote(.(canprot::cplab$DnH2O[[1]]) ~ "(" * .(log10p.nH2O) * ")")
    mtext(xlab, side = 1, line = par("mgp")[1], cex = 0.8)
    mtext(ylab, side = 2, line = par("mgp")[1], cex = 0.8)
    # Add title
    if(isTRUE(title)) title(na.omit(mdat$name)[1], font.main = 1, cex = 0.9, xpd = NA)
    else if(!isFALSE(title)) title(title, font.main = 1, cex = 0.9, xpd = NA)
    # Identify points 20200903
    if(identify) identify(DZC, DnH2O, mdat$subject[iminuend])
  }

  # Return the mean values 20201003
  invisible(list(study = study, DZC = D_mean_ZC, DnH2O = D_mean_nH2O))
}


# Add nH2O-ZC guidelines parallel to regression for amino acids
# Modified from JMDplots::gradH2O1() and JMDplots:::lmlines() 20200901
lmlines <- function(step = 0.01) {
  if(FALSE) {
    # Calculate ZC of the amino acids
    aa <- aminoacids("")
    ZC.aa <- ZC(info(aa, "aq"))
    # Load amino acids with QCa or QEC basis 20200914
    basis(c("glutamine", "glutamic acid", "cysteine", "H2O", "O2"))
    #if(options("basis")$basis == "QCa") basis(c("glutamine", "cysteine", "acetic acid", "H2O", "O2"))
    species(aa)
    # Make linear regression
    lm <- lm(species()$H2O ~ ZC.aa)
    coef <- coef(lm)
    # Clear species!
    reset()
  } else {
    # Use previously computed intercept and slope 20200920
    coef <- c(-0.1242780, -0.3088251)
    #if(options("basis")$basis == "QCa") coef <- c(-0.4830396, 0.1579203)
  }
  x <- par("usr")[1:2]
  y <- coef[1] + coef[2] * x
  for(dy in seq(-0.48, -1.20, -step)) lines(x, y + dy, col = "gray80")
  # Add box so ends of lines don't cover plot edges 20201007
  box()
}
