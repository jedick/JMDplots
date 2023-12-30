# JMDplots/orp16S_util.R
# Plotting functions used in orp16S paper 20211007

## Uncomment to source and run these functions interactively (developer mode)
#source("orp16S.R")

# Plot Zc values vs Eh7 for a single study 20210827
# Use 'groupby' (name of column with sample groups) and 'groups' (names of sample groups) to apply the pch and col to individual samples
plotEZ <- function(study, lineage = NULL, mincount = 100, pch = NULL, col = NULL, add = FALSE, type = "p", groupby = NULL, groups = NULL,
  legend.x = "topleft", show = c("lm", "points"), col.line = "gray62", lwd = 1, cex = 1, title.line = NA,
  dxlim = c(0, 0), dylim = c(0, 0), size = NULL, slope.legend = "title", ylim = NULL, ylab = cplab$Zc) {

  if(identical(lineage, "two")) {
    # Make two plots for studies that have Bacteria and Archaea 20210913
    out1 <- plotEZ(study, "Bacteria", mincount, pch, col, add, type, groupby, groups, legend.x, show, col.line, lwd, cex, title.line,
                   dxlim, dylim, size, slope.legend, ylim, ylab = ylab)
    # Don't show legend on second (Archaea) plot 20210914
    out2 <- plotEZ(study, "Archaea", mincount, pch, col, add, type, groupby, groups, legend.x = NA, show, col.line, lwd, cex, title.line,
                   dxlim, dylim, size, slope.legend, ylim, ylab = ylab)
    out <- c(out1, out2)
    return(invisible(out))
  }

  # Use try() to capture errors (with no mapped sequences for lineage = "Archaea")
  metrics.in <- try(
    suppressMessages(
      getmetrics_orp16S(study, lineage = lineage, mincount = mincount)
    ), silent = TRUE
  )
  # Print message and skip dataset with no mapped sequences
  if(inherits(metrics.in, "try-error")) {
    print(paste0(study, ": no mapped sequences for ", lineage))
    return()
  }

  # Get metadata
  mdat <- suppressMessages(getmdat_orp16S(study, metrics = metrics.in, size = size))
  metadata <- mdat$metadata
  metrics <- mdat$metrics
  # If specified lineage is Bacteria or Archaea, use only runs that are labeled as such 20210920
  if("Domain" %in% colnames(metadata)) {
    idomain <- rep(TRUE, nrow(metadata))
    if(identical(lineage, "Bacteria")) idomain <- metadata$Domain == "Bacteria"
    if(identical(lineage, "Archaea")) idomain <- metadata$Domain == "Archaea"
    metadata <- metadata[idomain, , drop = FALSE]
    metrics <- metrics[idomain, , drop = FALSE]
  }

  nsamp <- nrow(metadata)
  # Remove samples with NA Eh7 or Zc 20210822
  metadata <- metadata[!(is.na(metadata$Eh7) | is.na(metrics$Zc)), ]
  metrics <- metrics[metrics$Run %in% metadata$Run, ]
  stopifnot(all(metadata$Run == metrics$Run))
  # Print message about number of samples and Eh7 and Zc range
  Zctext <- paste(range(round(metrics$Zc, 3)), collapse = " to ")
  if(!is.null(lineage)) ltext <- paste0(lineage, ": ") else ltext <- ""
  print(paste0(study, ": ", ltext, nrow(metadata), "/", nsamp, " samples, Zc ", Zctext))

  # Assign pch and col to sample groups
  if(!is.null(groupby) & !is.null(groups)) {

    # Get default point symbols
    pchavail <- 21:25
    if(is.null(pch)) pch <- rep(pchavail, length(groups))

    # The pch and col for each sample type
    pchtype <- rep(pch, length.out = length(groups))
    coltype <- 1:length(groups)

    # The pch and col for individual samples
    pch <- col <- rep(NA, nrow(metadata))
    # The column with sample groups
    icol <- match(groupby, colnames(metadata))
    if(is.na(icol)) stop(paste(groupby, "is not a column name in metadata for", study))
    # Loop over sample groups
    for(i in seq_along(groups)) {
      # Find matching samples and set the pch and col
      itype <- metadata[, icol] == groups[i]
      pch[itype] <- pchtype[i]
      col[itype] <- orp16Scol[coltype[i]]
    }

  }

  # Defaults for pch and col if sample groups are not specified
  if(is.null(pch)) pch <- 19
  if(is.null(col)) col <- "#40404080"

  # Make data frame with Eh7 and Zc values
  # Add T and pH 20220516
  # Add O2_umol_L 20220517
  iT <- match("T", colnames(metadata))  # matches "T"
  if(is.na(iT)) iT <- grep("^T\\ ", colnames(metadata))[1]  # matches "T (°C)" but not e.g. "Treatment"
  if(is.na(iT)) iT <- grep("^Temp", colnames(metadata))[1]  # matches "Temperature (°C)"
  if(is.na(iT)) T <- NA else T <- metadata[, iT]
  ipH <- match("pH", colnames(metadata))
  if(is.na(ipH)) pH <- NA else pH <- metadata[, ipH]
  EZdat <- data.frame(T = T, pH = pH, O2_umol_L = metadata$O2_umol_L, Eh = metadata$Ehorig, Eh7 = metadata$Eh7, Zc = round(metrics$Zc, 6))
  # Get dataset name
  iname <- match("name", tolower(colnames(metadata)))
  name <- na.omit(metadata[, iname])[1]
  # Split suffix from study key 20210914
  root <- strsplit(study, "_")[[1]][1]
  suffix <- strsplit(study, "_")[[1]][2]
  # Create subtitle for environment type 20210904
  sub <- envirotype <- envirodat$group[envirodat$study == study]
  if(!add) {
    # Calculate x- and y-limits (with adjustment from arguments) 20220511
    xlim <- range(EZdat$Eh7) + dxlim
    if(is.null(ylim)) ylim <- range(EZdat$Zc) + dylim
    # Start new plot
    plot(EZdat$Eh7 / 1000, EZdat$Zc, xlab = "", ylab = ylab, type = "n", xlim = xlim / 1000, ylim = ylim)
    # Draw x-axis label with mtext to avoid getting cut off by small margin 20220517
    mtext("Eh7 (V)", side = 1, line = par("mgp")[1], cex = par("cex"))
    if(!is.null(title.line)) {
      main <- paste0(name, " (", root, ")")
      title(main = hyphen.in.pdf(main), font.main = 1, line = 2.5)
      # Include suffix in subtite 20210914
      if(!is.na(suffix)) sub <- paste(sub, "-", suffix)
      # Add lineage 20210913
      if(!is.null(lineage)) sub <- paste(sub, "-", lineage)
      title(main = sub, line = 1.5, cex.main = 1)
    }
  }
  # Add linear fit
  if("lm" %in% show) {
    adlilo <- add.linear.local(EZdat$Eh7, EZdat$Zc, col.line, lwd, slope.legend)
  }
  # Add points
  if("points" %in% show) {
    points(EZdat$Eh7 / 1000, EZdat$Zc, pch = pch, col = col, bg = col, type = type, cex = cex)
  }

  if(!is.null(groupby) & !is.null(groups)) {
    # Add legend
    legend <- as.character(groups)
    # Deal with mu character 20220522
    # Also deal with less than or equal to / greater than or equal to 20220614
    unicode <- c("\u03BC", "\u2264", "\u2265")
    legend.expr <- list()
    for(i in 1:length(legend)) {
      # Use unmodified expression unless we hit one of the special characters
      legend.expr[[i]] <- bquote(.(legend[i]))
      for(j in 1:length(unicode)) {
        if(grepl(unicode[j], legend[i])) {
          start <- strsplit(legend[i], unicode[j])[[1]][1]
          end <- strsplit(legend[i], unicode[j])[[1]][2]
          if(j == 1) legend.expr[[i]] <- bquote(.(start)*mu*.(end))
          if(j == 2) legend.expr[[i]] <- bquote(.(start) <= .(end))
          if(j == 3) legend.expr[[i]] <- bquote(.(start) >= .(end))
        } 
      }
    }

    # Make legend only for groups that are in the plot 20221001
    isgroup <- groups %in% metadata[, icol]
    legend <- as.expression(legend.expr[isgroup])
    legend(legend.x, legend, pch = pchtype[isgroup], col = orp16Scol[coltype][isgroup], pt.bg = orp16Scol[coltype][isgroup], title = groupby, cex = 0.9)
    # Add sample type (group) to output
    EZdat <- cbind(groupby = groupby, group = metadata[, icol], EZdat)

  } else {
    EZdat <- cbind(groupby = NA, group = NA, EZdat)
  }

  # Include suffix in name 20220522
  if(!is.na(suffix)) name <- paste0(name, " (", suffix, ")")
  # Return values
  # Use first column name starting with "sample" or "Sample" 20210818
  sampcol <- grep("^sample", colnames(metadata), ignore.case = TRUE)[1]
  if(is.null(lineage)) lineage <- ""
  if(length(envirotype) == 0) envirotype <- ""
  EZdat <- cbind(study = study, envirotype = envirotype, lineage = lineage, sample = metadata[, sampcol], Run = metadata$Run, EZdat)
  out <- list(study = study, name = name, envirotype = envirotype, lineage = lineage, metadata = metadata, EZdat = EZdat)
  if("lm" %in% show) out <- c(out, list(EZlm = adlilo$EZlm, Eh7lim = adlilo$Eh7lim, Zcpred = adlilo$Zcpred, pearson = adlilo$pearson))
  invisible(out)

}

# Plot Zc vs percentage of most abundant mapped taxon (MAMT) and
# show percentages of MAMT and MAUT (most abundant unmapped taxon) 20211007
plotMA <- function(study, lineage = NULL, mincount = 100, pch = NULL, col = NULL, groupby = NULL, groups = NULL, legend.x = "topright") {

  # Get RDP counts, mapping to NCBI taxonomy, and chemical metrics
  studyfile <- gsub("_.*", "", study)
  datadir <- system.file("extdata/orp16S/RDP", package = "JMDplots")
  RDPfile <- file.path(datadir, paste0(studyfile, ".tab.xz"))
  # If there is no .xz file, look for a .tab file 20210607
  if(!file.exists(RDPfile)) RDPfile <- file.path(datadir, paste0(studyfile, ".tab"))
  RDP <- read_RDP(RDPfile, lineage = lineage, mincount = mincount)
  map <- map_taxa(RDP, refdb = "RefSeq")
  metrics <- getmetrics_orp16S(study, lineage = lineage, mincount = mincount)
  mdat <- getmdat_orp16S(study, metrics)
  metadata <- mdat$metadata
  metrics <- mdat$metrics

  # Keep RDP columns for which we have data 20220508
  RDP <- cbind(RDP[, 1:4], RDP[, colnames(RDP) %in% metadata$Run])
  # Extract numeric rows
  RDPnum <- RDP[, -(1:4)]
  # Calculate sum of counts for each taxon
  taxoncounts <- rowSums(RDPnum)
  # Don't count unmapped taxa for MAMT identification
  mappedcounts <- taxoncounts
  mappedcounts[is.na(map)] <- 0
  # Find the MAMT for the entire dataset
  iMAMT <- which.max(mappedcounts)

  # Get the name and abundance of the MAMT
  MAMTname <- paste(RDP$rank[iMAMT], RDP$name[iMAMT], sep = "_")
  MAMTperc <- formatC(taxoncounts[iMAMT] / sum(taxoncounts) * 100, format = "f", digits = 1)
  # Get Zc of the MAMT
  datadir <- system.file("extdata/RefDB/RefSeq", package = "JMDplots")
  taxon_metrics <- read.csv(file.path(datadir, "taxon_metrics.csv.xz"), as.is = TRUE)
  MAMTZc <- taxon_metrics$Zc[map[iMAMT]]

  # Assign pch and col to sample groups
  if(!is.null(groupby) & !is.null(groups)) {
    # Get default point symbols
    pchavail <- 21:25
    if(is.null(pch)) pch <- rep(pchavail, length(groups))
    # The pch and col for each sample type
    pchtype <- rep(pch, length.out = length(groups))
    coltype <- 1:length(groups)
    # The pch and col for individual samples
    pch <- col <- rep(NA, nrow(metadata))
    # The column with sample groups
    icol <- match(groupby, colnames(metadata))
    if(is.na(icol)) stop(paste(groupby, "is not a column name in metadata for", study))
    # Loop over sample groups
    for(i in seq_along(groups)) {
      # Find matching samples and set the pch and col
      itype <- metadata[, icol] == groups[i]
      pch[itype] <- pchtype[i]
      col[itype] <- orp16Scol[coltype[i]]
    }
  }
  # Defaults for pch and col if sample groups are not specified
  if(is.null(pch)) pch <- 19
  if(is.null(col)) col <- "#40404080"

  # Calculate the abundance of the MAMT within each sample
  samplecounts <- colSums(RDPnum)
  perc <- as.numeric(RDPnum[iMAMT, ] / samplecounts * 100)
  # Start plot with percentage and Zc range
  plot(range(perc), range(metrics$Zc, MAMTZc), type = "n", xlab = "Abundance of MAMT (%)", ylab = cplab$Zc)
  # Add line for Zc of MAMT
  abline(h = MAMTZc, lty = 2)
  # Add points for community Zc vs abundance of MAMT
  points(perc, metrics$Zc, pch = pch, col = col, bg = col)

  # Get MAUT
  unmapped_groups <- attr(map, "unmapped_groups")
  unmapped_percent <- attr(map, "unmapped_percent")
  iMAUT <- which.max(unmapped_percent)
  MAUTname <- unmapped_groups[iMAUT]
  MAUTperc <- formatC(unmapped_percent[iMAUT], format = "f", digits = 1)

  # Add title: study, MAMT, MAUT
  main <- paste(study, lineage)
  title(main, line = 3)
  MAMTtitle <- paste0("MAMT: ", MAMTname, " (", MAMTperc, "%)")
  title(MAMTtitle, line = 2, font.main = 1)
  MAUTtitle <- paste0("MAUT: ", MAUTname, " (", MAUTperc, "%)")
  title(MAUTtitle, line = 1, font.main = 1)

  if(!is.null(groupby) & !is.null(groups)) {
    # Add legend
    legend <- as.character(groups)
    legend(legend.x, legend, pch = pchtype, col = orp16Scol[coltype], pt.bg = orp16Scol[coltype], title = groupby, cex = 0.9, bg = "white")
  }

}
