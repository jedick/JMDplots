# JMDplots/orp16S_plot.R
# Plotting functions used in orp16S paper 20211007

## Uncomment to source and run these functions interactively (developer mode)
#source("chem16S.R")
#source("orp16S.R")
#options(chem16Sdir = system.file("extdata/orp16S", package = "JMDplots"))

# Plot ZC values vs Eh7 for a single study 20210827
# Use 'groupby' (name of column with sample groups) and 'groups' (names of sample groups) to apply the pch and col to individual samples
plotEZ <- function(study, lineage = NULL, mincount = 50, pch = NULL, col = NULL, add = FALSE, type = "p", groupby = NULL, groups = NULL,
                   legend.x = "topleft", show = c("lm", "points"), col.line = "gray62", lwd = 1, cex = 1, mdat = NULL, title.line = NA) {

  if(identical(lineage, "two")) {
    # Make two plots for studies that have Bacteria and Archaea 20210913
    out1 <- plotEZ(study, "Bacteria", mincount, pch, col, add, type, groupby, groups, legend.x, show, col.line, lwd, cex)
    # Don't show legend on second (Archaea) plot 20210914
    out2 <- plotEZ(study, "Archaea", mincount, pch, col, add, type, groupby, groups, legend.x = NA, show, col.line, lwd, cex, mdat = out1$mdat)
    out <- c(out1, out2)
    return(invisible(out))
  }

  # Get metadata; use suppressMessages() to suppress messages
  if(is.null(mdat)) mdat <- suppressMessages(getmdat(study))
  mdat.orig <- mdat
  # For Bacteria or Archaea, use only runs that are labeled as such 20210920
  if("Domain" %in% colnames(mdat)) {
    if(identical(lineage, "Bacteria")) mdat <- mdat[mdat$Domain == "Bacteria", ]
    if(identical(lineage, "Archaea")) mdat <- mdat[mdat$Domain == "Archaea", ]
  }

  # Use capture.output to hide printed output
  null <- capture.output(
    # Use try() to capture errors (with no mapped sequences for lineage = "Archaea")
    met <- try(
      suppressMessages(
        getmetrics(study, mdat = mdat, lineage = lineage, mincount = mincount)
      ), silent = TRUE
    )
  )
  # Print message and skip dataset with no mapped sequences
  if(inherits(met, "try-error")) {
    print(paste0(study, ": no mapped sequences for ", lineage))
    return()
  }

  # Keep metadata only for samples with >= mincount counts 20201006
  mdat <- mdat[mdat$Run %in% met$Run, ]
  nsamp <- nrow(mdat)
  # Remove samples with NA Eh7 or ZC 20210822
  mdat <- mdat[!(is.na(mdat$Eh7) | is.na(met$ZC)), ]
  met <- met[met$Run %in% mdat$Run, ]
  stopifnot(all(mdat$Run == met$Run))
  # Print message about number of samples and Eh7 and ZC range
  ZCtext <- paste(range(round(met$ZC, 3)), collapse = " to ")
  if(!is.null(lineage)) ltext <- paste0(lineage, ": ") else ltext <- ""
  print(paste0(study, ": ", ltext, nrow(mdat), "/", nsamp, " samples, ZC ", ZCtext))

  # Assign pch and col to sample groups
  if(!is.null(groupby) & !is.null(groups)) {

    # Get default point symbols
    pchavail <- 21:25
    if(is.null(pch)) pch <- rep(pchavail, length(groups))

    # The pch and col for each sample type
    pchtype <- rep(pch, length.out = length(groups))
    coltype <- 1:length(groups)

    # The pch and col for individual samples
    pch <- col <- rep(NA, nrow(mdat))
    # The column with sample groups
    icol <- match(groupby, colnames(mdat))
    if(is.na(icol)) stop(paste(groupby, "is not a column name in metadata for", study))
    # Loop over sample groups
    for(i in seq_along(groups)) {
      # Find matching samples and set the pch and col
      itype <- mdat[, icol] == groups[i]
      pch[itype] <- pchtype[i]
      col[itype] <- orp16Scol[coltype[i]]
    }

  }

  # Defaults for pch and col if sample groups are not specified
  if(is.null(pch)) pch <- 19
  if(is.null(col)) col <- "#40404080"

  # Make data frame with Eh7 and ZC values
  EZdat <- data.frame(Eh = mdat$Ehorig, Eh7 = mdat$Eh7, ZC = round(met$ZC, 6))
  # Create subtitle for environment type 20210904
  sub <- envirotype <- envirodat$group[envirodat$study == study]
  if(!add) {
    # Start new plot
    plot(EZdat$Eh7, EZdat$ZC, xlab = "Eh7 (mV)", ylab = cplab$ZC, type = "n")
    # Take off suffix after underscore 20210914
    root <- strsplit(study, "_")[[1]][1]
    suffix <- strsplit(study, "_")[[1]][2]
    main <- paste0(na.omit(mdat$name)[1], " (", root, ")")
    title(main = main, font.main = 1, line = title.line)
    # Include suffix in subtite 20210914
    if(!is.na(suffix)) sub <- paste(sub, "-", suffix)
    # Add lineage 20210913
    if(!is.null(lineage)) sub <- paste(sub, "-", lineage)
    title(main = sub, line = 0.5, cex.main = 1)
  }
  # Add linear fit
  if("lm" %in% show) {
    EZlm <- lm(ZC ~ Eh7, EZdat)
    Eh7lim <- range(EZlm$model$Eh7)
    ZCpred <- predict.lm(EZlm, data.frame(Eh7 = Eh7lim))
    # Use solid or dashed line to indicate large or small slope 20210926
    slope <- EZlm$coefficients[2] * 1000
    if(is.na(slope)) lty <- 3 else if(abs(slope) < 0.01) lty <- 2 else lty <- 1
    lines(Eh7lim, ZCpred, col = col.line, lwd = lwd, lty = lty)
  }
  # Add points
  if("points" %in% show) {
    points(EZdat$Eh7, EZdat$ZC, pch = pch, col = col, bg = col, type = type, cex = cex)
  }

  if(!is.null(groupby) & !is.null(groups)) {
    # Add legend
    legend <- as.character(groups)
    legend(legend.x, legend, pch = pchtype, col = orp16Scol[coltype], pt.bg = orp16Scol[coltype], title = groupby, cex = 0.9)
    # Add sample type (group) to output
    EZdat <- cbind(groupby = groupby, group = mdat[, icol], EZdat)
  } else {
    EZdat <- cbind(groupby = NA, group = NA, EZdat)
  }


  # Return values
  # Use first column name starting with "sample" or "Sample" 20210818
  sampcol <- grep("^sample", colnames(mdat), ignore.case = TRUE)[1]
  if(is.null(lineage)) lineage <- ""
  EZdat <- cbind(study = study, envirotype = envirotype, lineage = lineage, sample = mdat[, sampcol], Run = mdat$Run, EZdat)
  out <- list(study = study, envirotype = envirotype, lineage = lineage, mdat = mdat.orig, EZdat = EZdat)
  if("lm" %in% show) out <- c(out, list(EZlm = EZlm, Eh7lim = Eh7lim, ZCpred = ZCpred))
  invisible(out)

}

# Plot ZC vs percentage of most abundant mapped taxon (MAMT) and
# show percentages of MAMT and MAUT (most abundant unmapped taxon) 20211007
plotMA <- function(study, lineage = NULL, mincount = 50, pch = NULL, col = NULL, groupby = NULL, groups = NULL, legend.x = "topright") {

  # Get RDP counts, mapping to NCBI taxonomy, and chemical metrics
  mdat <- getmdat(study)
  RDP <- getRDP(study, lineage = lineage, mincount = mincount, mdat = mdat)
  map <- getmap(study, lineage = lineage, mincount = mincount, RDP = RDP)
  metrics <- getmetrics(study, lineage = lineage, mincount = mincount, mdat = mdat, RDP = RDP, map = map)
  # Keep metadata only for samples with >= mincount counts 20211008
  mdat <- mdat[mdat$Run %in% metrics$Run, ]

  # Extract numeric rows
  RDPnum <- RDP[, -(1:3)]
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
  # Get ZC of the MAMT
  datadir <- system.file("extdata/chem16S", package = "JMDplots")
  taxon_metrics <- read.csv(file.path(datadir, "taxon_metrics.csv"), as.is = TRUE)
  MAMTZC <- taxon_metrics$ZC[map[iMAMT]]

  # Assign pch and col to sample groups
  if(!is.null(groupby) & !is.null(groups)) {
    # Get default point symbols
    pchavail <- 21:25
    if(is.null(pch)) pch <- rep(pchavail, length(groups))
    # The pch and col for each sample type
    pchtype <- rep(pch, length.out = length(groups))
    coltype <- 1:length(groups)
    # The pch and col for individual samples
    pch <- col <- rep(NA, nrow(mdat))
    # The column with sample groups
    icol <- match(groupby, colnames(mdat))
    if(is.na(icol)) stop(paste(groupby, "is not a column name in metadata for", study))
    # Loop over sample groups
    for(i in seq_along(groups)) {
      # Find matching samples and set the pch and col
      itype <- mdat[, icol] == groups[i]
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
  # Start plot with percentage and ZC range
  plot(range(perc), range(metrics$ZC, MAMTZC), type = "n", xlab = "Abundance of MAMT (%)", ylab = cplab$ZC)
  # Add line for ZC of MAMT
  abline(h = MAMTZC, lty = 2)
  # Add points for community ZC vs abundance of MAMT
  points(perc, metrics$ZC, pch = pch, col = col, bg = col)

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
