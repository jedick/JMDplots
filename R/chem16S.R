# JMDplots/chem16S.R
# Calculate chemical metrics based on 16S data and RefSeq proteins 20200902
# Revised to include "unclassified" groups in RDP (i.e. classified above genera) 20200911
# Moved to JMDplots 20210416

# Utility functions
# getmdat("SVH+19")     # Metadata for this study (study, name, Run, BioSample, sample, type, cohort)
# getRDP("SVH+19")      # RDP results for all samples in this study
# getmap("SVH+19")      # Map RDP to RefSeq taxonomy (match to rows of taxon_AA.csv)
# getmetrics("SVH+19")  # Calculate chemical metrics (nH2O, ZC) for each sample

#####################
# Utility functions #
#####################

# Get metadata for a study, appending columns for pch and col 20200914
getmdat <- function(study, dropNA = TRUE) {
  # Read metadata file
  # Remove suffix after underscore 20200929
  studyfile <- gsub("_.*", "", study)
  datadir <- getOption("chem16Sdir")
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
  # Use NULL pch as flag for unavailable dataset 20210820
  pch <- NULL

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
    # For fauna, use a darkened yellow4 with transparency 20210609
    col <- sapply(type, switch, lowT = 4, highT = 2, "fauna surface" = "#757500C0", "rock/chimney" = 1, NA)
  }
  if(study == "SVH+19") { # Black Sea
    pch <- sapply(mdat$Type, switch, Oxic = 24, Suboxic = 20, Euxinic = 25, NA)
    col <- sapply(mdat$Type, switch, Oxic = 4, Suboxic = 1, Euxinic = 2, NA)
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
  if(grepl("CHM+14", study, fixed = TRUE)) {
    # Injected fluids and later produced water
    if(study == "CHM+14_injected-49") mdat <- mdat[mdat$day == 0 | mdat$day >= 49, ]
    pch <- ifelse(mdat$day >= 49, 21, 1)
    col <- ifelse(mdat$day >= 49, 2, 1)
  }
  if(grepl("HRR+18", study, fixed = TRUE)) {
    if(study == "HRR+18_injected-130") mdat <- mdat[mdat$day %in% c(0, 130, 220), ]
    pch <- ifelse(mdat$day > 10, 21, 1)
    col <- ifelse(mdat$day > 10, 2, 1)
  }
  if(grepl("ZLF+19", study, fixed = TRUE)) {
    # Source water and flowback day 18
    if(study == "ZLF+19_injected-18") mdat <- mdat[mdat$day %in% c(-1, 18), ]
    pch <- ifelse(mdat$day >= 1, 21, 1)
    col <- ifelse(mdat$day >= 1, 2, 1)
  }

  # Stratified water datasets
  if(study == "MZG+20") {
    # Identify shallowest and deepest samples from Lakes Zug and Lugano
    newdat <- lapply(c("Lake Zug", "Lake Lugano"), function(lake) {
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
    type <- rep("transition", nrow(mdat))
    type[mdat$O2 > 100] <- "oxic"
    type[mdat$O2 == 0] <- "anoxic"
    pch <- sapply(type, switch, oxic = 24, transition = 20, anoxic = 25)
    col <- sapply(type, switch, oxic = 4, transition = 1, anoxic = 2)
    type[mdat$station == "C4"] <- NA
    col[mdat$station == "C4"] <- NA
  }
  if(study == "GBL+15") {
    pch <- ifelse(mdat$depth < 100, 24, ifelse(mdat$depth > 100, 25, 20))
    col <- ifelse(mdat$depth < 100, 4, ifelse(mdat$depth > 100, 2, 1))
    pch[mdat$size != "0.2-1.6micron"] <- NA
    col[mdat$size != "0.2-1.6micron"] <- NA
  }
  if(study == "BCA+21") {
    type <- rep("transition", nrow(mdat))
    type[mdat$depth < 3] <- "oxic"
    type[mdat$depth > 4] <- "anoxic"
    pch <- sapply(type, switch, oxic = 24, transition = 20, anoxic = 25)
    col <- sapply(type, switch, oxic = 4, transition = 1, anoxic = 2)
  }
  if(study == "EH18") {
    pch <- sapply(mdat$type, switch, top = 24, transition = 20, bottom = 25)
    col <- sapply(mdat$type, switch, top = 4, transition = 1, bottom = 2)
  }

  ## Datasets for orp16S paper 20211003
  shortstudy <- study
  if(study == "RBW+14") {
    type <- rep("reducing", nrow(mdat))
    type[mdat$layer == "Top"] <- "oxidizing"
    type[mdat$layer == "SWI"] <- "transition"
    pch <- sapply(type, switch, oxidizing = 24, transition = 20, reducing = 25)
    col <- sapply(type, switch, oxidizing = 4, transition = 1, reducing = 2)
  }
  if(grepl("PCL\\+18", study)) {
    # PCL+18, PCL+18_Acidic, PCL+18_Alkaline
    Type <- sapply(strsplit(study, "_"), "[", 2)
    if(!is.na(Type)) mdat <- mdat[mdat$Type == Type, ]
    shortstudy <- "PCL+18"
  }
  if(grepl("LMBA21", study)) {
    # LMBA21, LMBA21_2013, LMBA21_2015, LMBA21_2016, LMBA21_2017, LMBA21_2018, LMBA21_2019
    Year <- sapply(strsplit(study, "_"), "[", 2)
    if(!is.na(Year)) mdat <- mdat[mdat$Year == Year, ]
    shortstudy <- "LMBA21"
  }
  if(study == "DLS21") {
    # Just look at bulk soil 20210910
    mdat <- mdat[mdat$Source == "bulk soil", ]
  }
  if(shortstudy %in% c(
    "MLL+19", "HXZ+20", "BCA+21", "RSJ+21", "RMB+17", "SBP+20", "NTB+21", "MWY+21", "SAR+13", "CTS+17",
    "SCM+18", "HDZ+19", "BOEM21", "ZHZ+19", "YHK+20", "CNA+20", "BMJ+19", "SRM+19", "HLZ+18", "XLD+20",
    "JHL+12", "PSG+20", "KSR+21", "ZCZ+21", "SKP+21", "ZZL+21", "PBU+20", "GWS+19", "KLY+20", "SRM+21",
    "MLL+18", "JDP+20", "BWD+19", "LXH+20", "LMG+20", "WHL+21", "LLL+21", "SDH+19", "GWSS21", "HSF+19",
    "ZML+17", "DTJ+20", "WFB+21", "SBW+17", "KLM+16", "LMBA21", "ZDA+20", "ZZZ+18", "BSPD17", "CWC+20",
    "BMOB18", "JVW+20", "LJC+20", "GFE+16", "ECS+18", "FAV+21", "VMB+19", "DLS21", "ZZLL21", "GWS+20",
    "CLS+19", "SMS+12", "OFY+19", "BYB+17", "MCS+21", "SVH+19", "PMM+20", "GZL21", "LLC+19", "NLE+21",
    "GSY+20", "SCH+16", "LZR+17", "GRG+20", "APV+20", "YHK+19", "WHC+19", "WHLH21", "PCL+18"
  )) {
    # General processing of metadata for orp16S datasets 20210820
    # Get Eh or ORP values (uses partial name matching, can match a column named "Eh (mV)")
    Eh <- mdat$Eh
    Ehname <- "Eh"
    if(is.null(Eh)) {
      Eh <- mdat$ORP
      Ehname <- "ORP"
    }
    # Also match a column starting with "redox" (case-insensitive)
    iEh <- grep("^redox", colnames(mdat), ignore.case = TRUE)
    if(length(iEh) > 0) {
      Eh <- mdat[, iEh[1]]
      Ehname <- strsplit(colnames(mdat)[iEh[1]], " ")[[1]][1]
    }
    if(is.null(Eh)) stop("can't find Eh or ORP column")

    # Append Ehorig column (normalized column name used for exporting data) 20210821
    mdat <- cbind(mdat, Ehorig = Eh)
    
    # Get temperature for dEh/dpH calculation 20210829
    iT <- match("T", colnames(mdat))  # matches "T"
    if(is.na(iT)) iT <- grep("^T\\ ", colnames(mdat))[1]  # matches "T (°C)" but not e.g. "Treatment"
    if(is.na(iT)) iT <- grep("^Temp", colnames(mdat))[1]  # matches "Temperature (°C)"
    if(!is.na(iT)) {
      T <- mdat[, iT]
      Ttext <- paste(round(range(na.omit(T)), 1), collapse = " to ")
    } else {
      T <- rep(25, nrow(mdat))
      Ttext <- "assumed 25"
    }

    # Adjust Eh to pH 7 20210828
    # Find pH column
    ipH <- match("pH", colnames(mdat))
    if(!is.na(ipH)) {
      # pH values
      pH <- mdat[, ipH]
      pHtext <- paste(round(range(na.omit(pH)), 1), collapse = " to ")
      ## Eh(mV)-pH slope at 25 °C
      #dEhdpH <- -59.2
      # Find dEh/dpH as a function of T 20210829
      TK <- T + 273.15
      R <- 0.0083147
      F <- 96.4935
      dEhdpH <- - (log(10) * R * TK) / F * 1000

      # Difference to pH 7
      pHdiff <- 7 - pH
      # Adjust to pH 7
      Eh7 <- round(Eh + pHdiff * dEhdpH, 1)

    } else {
      Eh7 <- Eh
      pHtext <- "assumed 7"
    }
    mdat <- cbind(mdat, Eh7 = Eh7)
    # Print message about T, pH and Eh ranges
    Ehtext <- paste(round(range(na.omit(Eh))), collapse = " to ")
    Eh7text <- paste(round(range(na.omit(Eh7))), collapse = " to ")
    print(paste0(study, ": T ", Ttext, ", pH ", pHtext, ", Eh ", Ehtext, ", Eh7 ", Eh7text))

    # Keep NA values out of clusters
    ina <- is.na(Eh7)
    # Divide data into 3 clusters
    cl <- kmeans(Eh7[!ina], 3, nstart = 100)
    # Find the clusters with the lowest and highest means
    imin <- which.min(cl$centers[, 1])
    imax <- which.max(cl$centers[, 1])
    imid <- (1:3)[-c(imin, imax)]
    # Name the clusters
    cluster <- rep(NA, nrow(mdat))
    cluster[!ina][cl$cluster == imin] <- "reducing"
    cluster[!ina][cl$cluster == imid] <- "transition"
    cluster[!ina][cl$cluster == imax] <- "oxidizing"
    # Get pch and col for plot
    # Use is.null to let previous assignments take precedence (for SVH+19 in geo16S) 20211008
    if(is.null(pch)) {
      pch <- sapply(cluster, switch, oxidizing = 24, transition = 21, reducing = 25, NA)
      col <- sapply(cluster, switch, oxidizing = 4, transition = 1, reducing = 2, NA)
    }
    # Print message about Eh7 or ORP ranges of clusters
    ro <- range(na.omit(Eh7[cluster=="oxidizing"]))
    no <- sum(na.omit(cluster=="oxidizing"))
    message(paste0("getmdat: oxidizing Eh7 ", paste(ro, collapse = " to "), " mV (n = ", no, ")"))
    rt <- range(na.omit(Eh7[cluster=="transition"]))
    nt <- sum(na.omit(cluster=="transition"))
    message(paste0("getmdat: transition Eh7 ", paste(rt, collapse = " to "), " mV (n = ", nt, ")"))
    rr <- range(na.omit(Eh7[cluster=="reducing"]))
    nr <- sum(na.omit(cluster=="reducing"))
    message(paste0("getmdat: reducing Eh7 ", paste(rr, collapse = " to "), " mV (n = ", nr, ")"))
    # Append cluster names to data frame
    mdat <- cbind(mdat, cluster)
  }

  if(is.null(pch)) stop(paste(study, "metadata file exists, but not set up for processing"))

  mdat <- cbind(mdat, pch, col)
  mdat
}

# Get RDP results for all samples in a study 20200912
getRDP <- function(study, cn = FALSE, mdat = NULL, lineage = NULL, mincount = 200) {
  # Handle missing arguments
  if(is.null(mdat)) mdat <- getmdat(study)
  # Get the sample Runs (SRR numbers)
  Run <- mdat$Run
  # Remove suffix after underscore 20200929
  studyfile <- gsub("_.*", "", study)
  # Read output of RDP classifer
  datadir <- getOption("chem16Sdir")
  file <- file.path(datadir, "RDP", paste0(studyfile, ".tab.xz"))
  # If there is no .xz file, look for a .tab file 20210607
  if(!file.exists(file)) file <- file.path(datadir, "RDP", paste0(studyfile, ".tab"))
  datorig <- dat <- read.table(file, sep = "\t", header = TRUE, check.names = FALSE)
  # Get counts for each sample
  icol <- match(Run, colnames(dat))
  # Error for samples not in RDP output 20210926
  if(any(is.na(icol))) stop(paste0("One or more runs in metadata not in RDP output for ", study, ": ", paste(Run[is.na(icol)], collapse = ", ")))
  # Keep the "lineage", "rank", "name", and counts columns
  dat <- dat[, c(2, 4, 3, icol)]

  iunclass <- grepl("unclassified_", dat$name)
  if(any(iunclass)) {
    # Find the ranks for the "unclassified" counts
    # (which means classified at higher rank than genus)
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
  if(dat$lineage[1] == "null") {
    # This is for the actual output from the RDP Classifier
    # We want to keep genus-level classifications; also:
    # Remove the counts of classifications at higher ranks (which are also counted at lower ranks)
    # *except* for the ones labeled "unclassified" (which are not counted at lower ranks)
    igenus <- dat$rank == "genus"
    out <- dat[igenus | iunclass, ]
    isRDP <- TRUE
  } else {
    # Otherwise, the taxon counts are assembled from SI or OTU tables; use all the counts
    out <- dat
    isRDP <- FALSE
  }
  # Get the total counts
  totalcounts <- colSums(out[, -(1:3), drop = FALSE])
  if(isRDP) {
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
    # Recalculate total counts 20211008
    totalcounts <- colSums(out[, -(1:3), drop = FALSE])
  }

  # Keep the rows with any counts > 0
  groupcounts <- rowSums(out[, -(1:3), drop = FALSE])
  out <- out[groupcounts > 0, ]

  # Get the number of counts classified at genus level
  igenus <- out$rank == "genus"
  genuscounts <- colSums(out[igenus, -(1:3), drop = FALSE], na.rm = TRUE)
  # Print percentage of assignments at genus level
  genuspercent <- round(100 * sum(genuscounts) / sum(totalcounts))
  print(paste0("getRDP [", study, "]: ", genuspercent, "% of classifications at genus level"))

  # Remove classifications at root and domain level (Bacteria and Archaea),
  # and Chlorophyta, Chloroplast and Bacillariophyta 20200922
  # and Eukaryota 20211012
  RDPgroups <- paste(out$rank, out$name, sep = "_")
  RDPgroups[grepl("Eukaryota", out$lineage)] <- "Eukaryota"
  rmgroups <- c("rootrank_Root", "domain_Bacteria", "domain_Archaea", "Eukaryota",
    # Cyanobacteria/Chloroplast
    "class_Chloroplast", "family_Chloroplast", "genus_Chlorophyta", "genus_Bacillariophyta")
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

  # Discard samples with < mincount total counts 20201001
  ismall <- totalcounts < mincount
  if(any(ismall)) {
    print(paste0("getRDP [", study, "]: discarding ", sum(ismall), " samples with < ", mincount, " total counts"))
    out <- out[, c(TRUE, TRUE, TRUE, !ismall)]
  }
  # Recalculate total counts
  totalcounts <- colSums(out[, -(1:3), drop = FALSE])
  # Report the median number of counts 20200917
  # Change this to range 20200924
  if(length(totalcounts) == 0) {
    print(paste0("getRDP [", study, "]: no samples contain at least ", mincount, " counts"))
  } else {
    print(paste0("getRDP [", study, "]: count range is [", paste(round(range(totalcounts)), collapse = " "), "]"))
  }

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
getmap <- function(study, RDP = NULL, lineage = NULL, mincount = 200) {
  # Handle missing arguments
  if(is.null(RDP)) {
    mdat <- getmdat(study)
    RDP <- getRDP(study, mdat = mdat, lineage = lineage, mincount = mincount)
  }
  # Make group names by combining rank and name
  RDPgroups <- paste(RDP$rank, RDP$name, sep = "_")
  # Calculate group abundances for displaying messages
  groupcounts <- as.numeric(rowSums(RDP[, -(1:3), drop = FALSE]))

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
    # 20210609 Verrucomicrobia
    "genus_Subdivision3_genera_incertae_sedis" = "family_Verrucomicrobia subdivision 3",
    # 20211215 Clostridia
    # Eubacteriales is used in NCBI; is a synonym for Clostridiales (https://lpsn.dsmz.de/order/eubacteriales)
    "order_Clostridiales" = "order_Eubacteriales",
    # Oscillospiraceae is used in NCBI; is a synonym for Ruminococcaceae (https://lpsn.dsmz.de/order/eubacteriales)
    "family_Ruminococcaceae" = "family_Oscillospiraceae",

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
  datadir <- system.file("extdata/chem16S", package = "JMDplots")
  AA <- read.csv(file.path(datadir, "taxon_AA.csv"), as.is = TRUE)
  AAgroups <- paste(AA$protein, AA$organism, sep = "_")
  # Do the mapping!
  iAA <- match(RDPgroups, AAgroups)
  # Get percentages of unmapped groups
  naAA <- is.na(iAA)
  nagroups <- RDPgroups[naAA]
  nacounts <- groupcounts[naAA]
  napercent <- nacounts / sum(groupcounts) * 100
  # Print summary of missing groups
  naorder <- order(napercent, decreasing = TRUE)
  ordergroups <- nagroups[naorder]
  orderpercent <- round(napercent[naorder], 2)
  if(sum(naAA) > 0) namsg <- paste0("getmap [", study, "]: can't map RDP group ", ordergroups[1], " (", orderpercent[1], "%)")
  if(sum(naAA) > 1) namsg <- paste0("getmap [", study, "]: can't map RDP groups ", ordergroups[1], " (", orderpercent[1], "%), ",
                                    ordergroups[2], " (", orderpercent[2], "%)")
  if(sum(naAA) > 2) namsg <- paste0("getmap [", study, "]: can't map RDP groups ", ordergroups[1], " (", orderpercent[1], "%), ",
                                    ordergroups[2], " (", orderpercent[2], "%), ", ordergroups[3], " (", orderpercent[3], "%)")
  if(sum(naAA) > 3) namsg <- paste0("getmap [", study, "]: can't map RDP groups ", ordergroups[1], " (", orderpercent[1], "%), ",
                                    ordergroups[2], " (", orderpercent[2], "%), ", sum(naAA) - 2, " others (", sum(orderpercent[-(1:2)]), "%)")
  if(sum(naAA) > 0) print(namsg)
  # Print message about total mapped percent 20200927
  mappedpercent <- formatC(100 - sum(napercent), 1, format = "f")
  print(paste0("getmap [", study, "]: mapped ", mappedpercent, "% of RDP classifications to NCBI taxonomy"))
  # Set attributes to indicate unmapped groups 20211007
  attr(iAA, "unmapped_groups") <- nagroups
  attr(iAA, "unmapped_percent") <- napercent
  # Return result
  iAA
}

# Get chemical metrics for all samples in a study 20200927
getmetrics <- function(study, cn = FALSE, mdat = NULL, RDP = NULL, map = NULL, lineage = NULL, mincount = 200, metrics = NULL, groups = NULL) {
  # Handle missing arguments
  if(is.null(mdat)) mdat <- getmdat(study)
  if(is.null(RDP)) RDP <- getRDP(study, cn = cn, mdat = mdat, lineage = lineage, mincount = mincount)
  if(is.null(map)) map <- getmap(study, RDP = RDP, lineage = lineage, mincount = mincount)
  # Keep metadata only for samples with >= mincount counts 20201001
  mdat <- mdat[mdat$Run %in% colnames(RDP), ]
  if(nrow(mdat) == 0) stop("no samples have >= mincount counts!")
  # Exclude NA mappings
  RDP <- RDP[!is.na(map), ]
  map <- na.omit(map)
  if(length(map) == 0) stop("no mappings to available RefSeq taxa!")

  # Get chemical metrics of RefSeq groups
  datadir <- system.file("extdata/chem16S", package = "JMDplots")
  if(is.null(metrics)) metrics <- read.csv(file.path(datadir, "taxon_metrics.csv"), as.is = TRUE)
  metrics <- metrics[map, ]
  # Make sure the mapping is correct
  equalrank <- RDP$rank == metrics$rank
  # Don't test particular RDP-NCBI mappings that cross ranks
  iclassCyano <- RDP$rank == "class" & RDP$name == "Cyanobacteria"
  igenusSparto <- RDP$rank == "genus" & RDP$name == "Spartobacteria_genera_incertae_sedis"
  iclassActino <- RDP$rank == "class" & RDP$name == "Actinobacteria"
  igenusMarini <- RDP$rank == "genus" & RDP$name == "Marinimicrobia_genera_incertae_sedis"
  igenusVerruco <- RDP$rank == "genus" & RDP$name == "Subdivision3_genera_incertae_sedis"
  igenusGpI <- RDP$rank == "genus" & RDP$name == "GpI"
  igenusGpIIa <- RDP$rank == "genus" & RDP$name == "GpIIa"
  equalrank <- equalrank[!(iclassCyano | igenusSparto | iclassActino | igenusMarini | igenusVerruco | igenusGpI | igenusGpIIa)]
  stopifnot(all(equalrank))

  # Get classification matrix (rows = taxa, columns = samples)
  RDPmat <- RDP[, -(1:3), drop = FALSE]
  # Calculate abundance-weighted mean nH2O for each sample
  nH2O <- colSums(RDPmat * metrics$nH2O) / colSums(RDPmat)
  # To calculate ZC, we need to compute the sum of charge (ZC * nC) and the sum of carbon atoms
  sumZ <- colSums(RDPmat * metrics$ZC * metrics$nC)
  sumC <- colSums(RDPmat * metrics$nC)
  ZC <- sumZ / sumC
  if(is.null(groups)) {
    # Create output data frame
    # Use first column name starting with "sample" or "Sample" 20210818
    sampcol <- grep("^sample", colnames(mdat), ignore.case = TRUE)[1]
    out <- data.frame(Run = colnames(RDPmat), sample = mdat[, sampcol], nH2O = nH2O, ZC = ZC)
  } else {
    # Split data into sample groups and calculate metrics for each group 20210607
    nH2O <- ZC <- numeric()
    for(i in 1:length(groups)) {
      # Use rowSums to combine all samples in each group into one meta-sample
      thisRDP <- rowSums(RDPmat[, groups[[i]], drop = FALSE])
      nH2O <- c(nH2O, sum(thisRDP * metrics$nH2O) / sum(thisRDP))
      sumZ <- sum(thisRDP * metrics$ZC * metrics$nC)
      sumC <- sum(thisRDP * metrics$nC)
      ZC <- c(ZC, sumZ / sumC)
    }
    out <- data.frame(Run = rep(NA, length(groups)), sample = 1:length(groups), nH2O = nH2O, ZC = ZC)
  }

  rownames(out) <- 1:nrow(out)
  out
}

