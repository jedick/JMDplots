# JMDplots/comp16S_plot.R
# Plotting functions for predicted RefSeq and inferred community proteomes
# Separated from comp16S.R 20210607

# taxacomp()                # nH2O-ZC plot for taxa (default: Bacteria, Archaea) and their children
# plotcomp("KGP+12")        # nH2O-ZC plot for specified study
# diffcomp("KGP+12")        # DnH2O-DZC plot for specified study

######################
# Plotting functions #
######################

# Get amino acid compositions for taxids in RefSeq, excluding
# super-sequenced species (biased to high ZC/low nH2O) 20210604
getrefseq <- function(filterspecies = TRUE) {
  # Read RefSeq amino acid compositions and taxid names
  refseq <- read.csv(system.file("extdata/refseq/protein_refseq.csv.xz", package = "JMDplots"), as.is = TRUE)
  taxa <- read.csv(system.file("extdata/refseq/taxid_names.csv.xz", package = "JMDplots"), as.is = TRUE)
  if(filterspecies) {
    # Take out species with > 20000 sequences
    ispecies <- !is.na(taxa$species)
    isuper <- refseq$chains > 20000
    isuperspecies <- ispecies & isuper
    message(paste("getrefseq: removing", sum(isuperspecies), "species with > 20000 sequences"))
    # Return both the amino acid compositions and taxid names
    # NOTE: the following genera are completely removed: Buchnera, Sorangium, Dolosigranulum, Enhygromyxa, Ruthenibacterium
    refseq <- refseq[!isuperspecies, ]
    taxa <- taxa[!isuperspecies, ]
  }
  list(refseq = refseq, taxa = taxa)
}

# Make a nH2O-ZC plot for selected taxa and all their children 20200911
taxacomp <- function(groups = c("Bacteria", "Archaea"), xlim = NULL, ylim = NULL,
  col = seq_along(groups), legend.x = "topleft", identify = FALSE, pch = NULL, hline = NULL, filterspecies = TRUE, lcol = NULL) {

  # Read compositional metrics of all taxa
  datadir <- system.file("extdata/comp16S", package = "JMDplots")
  metrics <- read.csv(file.path(datadir, "RefSeq_metrics.csv"), as.is = TRUE)
  # Default point symbols
  taxa <- groups
  if(is.null(pch)) pch <- rep(21, length(taxa))
  lty <- 1

  # For "majorphyla", get names of phyla with more than 500 representatives
  if(identical(groups, "majorphyla")) {
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
  if(identical(groups, "majorcellular")) {
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
  # How to count the representatives in each proteobacterial class:
  #> taxa <- read.csv(system.file("extdata/refseq/taxid_names.csv.xz", package = "JMDplots"), as.is = TRUE)
  #> sort(table(na.omit(taxa$class[taxa$phylum == "Proteobacteria"])), decreasing = TRUE)
  #
  #  Gammaproteobacteria   Alphaproteobacteria    Betaproteobacteria 
  #                 8269                  5667                  2456 
  #Epsilonproteobacteria   Deltaproteobacteria           Oligoflexia 
  #                  451                   441                    32 
  #    Acidithiobacillia     Hydrogenophilalia    Zetaproteobacteria 
  #                   20                    11                    11 
  if(identical(groups, "proteobacteria")) {
    taxa <- c("Alphaproteobacteria", "Betaproteobacteria", "Gammaproteobacteria", "Deltaproteobacteria", "Epsilonproteobacteria", "Zetaproteobacteria",
              "Acidithiobacillia", "Hydrogenophilalia", "Oligoflexia")
    pch <- rep(21:23, 3)
    col <- seq_along(taxa)
    if(is.null(ylim)) ylim <- c(-0.77, -0.71)
  }

  # Acidobacteria 20210529
  #> sort(table(na.omit(taxa$class[taxa$phylum == "Acidobacteria"])), decreasing = TRUE)
  #     Acidobacteriia      Blastocatellia          Holophagae Thermoanaerobaculia    Vicinamibacteria 
  #                 68                  12                   5                   1                   1 
  if(identical(groups, "acidobacteria")) {
    taxa <- c("Acidobacteriia", "Blastocatellia", "Holophagae", "Thermoanaerobaculia", "Vicinamibacteria")
    pch <- c(21, 22, 23, 21, 22)
    col <- seq_along(taxa)
  }

  # Cyanobacteria **orders** 20210529
  #> sort(table(na.omit(taxa$order[taxa$phylum == "Cyanobacteria"])), decreasing = TRUE)
  #      Synechococcales            Nostocales       Oscillatoriales 
  #                  323                   228                   121 
  #        Chroococcales        Pleurocapsales Chroococcidiopsidales 
  #                   49                    13                     9 
  #      Gloeobacterales          Spirulinales    Gloeoemargaritales 
  #                    3                     2                     1 
  if(identical(groups, "cyanobacteria")) {
    taxa <- c("Synechococcales", "Nostocales", "Oscillatoriales", "Chroococcales", 
      "Pleurocapsales", "Chroococcidiopsidales", "Gloeobacterales", 
      "Spirulinales", "Gloeoemargaritales")
    pch <- rep(21:23, 3)
    col <- seq_along(taxa)
    if(is.null(xlim)) xlim <- c(-0.19, -0.13)
    if(is.null(ylim)) ylim <- c(-0.77, -0.71)
  }

  # Use semi-transparent colors for lines 20210518
  if(is.null(lcol)) {
    lcol <- palette()
    lcol[1] <- "#000000"  # black
    lcol[8] <- "#9e9e9e"  # gray62
    lcol <- paste0(lcol, "80")
    lcol <- rep(lcol, length.out = length(col))
  }

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
  refseq_species <- NULL

  for(i in seq_along(taxa)) {
    thisgroup <- taxa[i]
    # Get the compositional metrics for this group
    igroup <- metrics$group == thisgroup
    if(sum(igroup) > 1) warning(paste0("found more than one ", thisgroup, " (", paste(metrics$rank[igroup], collapse = ", "), "); using the first"))
    group <- metrics[which(igroup)[1], ]
    if(identical(group$rank, "genus")) {
      # For a genus, look for children (species) in full RefSeq data frame 20210603
      if(is.null(refseq_species)) {
        # Read RefSeq amino acid compositions and taxon names
        gr <- getrefseq(filterspecies)
        refseq <- gr$refseq
        alltaxa <- gr$taxa
        # Keep species-level taxa that have a genus name
        irefseq <- !is.na(alltaxa$species) & !is.na(alltaxa$genus)
        refseq_species <- refseq[irefseq, ]
        # Put genus name in "abbrv" column
        refseq_species$abbrv <- alltaxa$genus[irefseq]
        # Keep species with at least 1000 sequences
        refseq_species <- refseq_species[refseq_species$chains >= 1000, ]
      }
      # Find species in this genus and calculate ZC and nH2O
      ispecies <- refseq_species$abbrv == thisgroup
      sp.refseq <- refseq_species[ispecies, ]
      sp.ZC <- ZCAA(sp.refseq)
      sp.nH2O <- H2OAA(sp.refseq)
      children <- data.frame(group = sp.refseq$ref, chains = sp.refseq$chains, ZC = sp.ZC, nH2O = sp.nH2O)
    } else {
      # Get the compositional metrics for all children
      ichildren <- metrics$parent == thisgroup
      children <- metrics[ichildren, ]
    }
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
    pt.col <- 1
    if(is.numeric(col[i])) {
      if((col[i] - 1) %% 8 == 0) pt.col <- "white"
    }
    points(children$ZC, children$nH2O, pch = pch[i], cex = 0.7, col = pt.col, bg = col[i], lwd = 0.5)
    # Label Halobacteria and Nanohaloarchaea 20200930
    thisgroup <- taxa[i]
    if(thisgroup == "Euryarchaeota") {
      ihalo <- match(c("Methanococci", "Thermococci", "Archaeoglobi", "Nanohaloarchaea", "Halobacteria"), children$group)
      dy <- ifelse(groups == "majorcellular", 0.0025, 0.005)
      dx <- c(0, 0, 0, 0.002, 0)
      text(children$ZC[ihalo] + dx, children$nH2O[ihalo] + dy, c(1, 2, 3, 4, 5))
    }
    # Label Clostridia 20200930
    if(thisgroup == "Firmicutes" & identical(groups, "majorcellular")) {
      iclos <- match("Clostridia", children$group)
      text(children$ZC[iclos], children$nH2O[iclos] + 0.0025, 6)
    }
#    # Label Pisoniviricetes 20210520
#    if(thisgroup == "Pisuviricota") {
#      ipison <- match("Pisoniviricetes", children$group)
#      text(children$ZC[ipison], children$nH2O[ipison] - 0.0025, 5)
#    }
  }

  # Add legend
  len <- length(taxa)
  if(identical(groups, "majorcellular")) {
    legend("bottomleft", taxa[1:8], pch = pch[1:8], col = col[1:8], pt.bg = col[1:8], cex = 0.9, bg = "white")
    legend("bottomright", taxa[9:len], pch = pch[9:len], col = col[9:len], pt.bg = col[9:len], cex = 0.9, bg = "white")
  } else if(identical(groups, "proteobacteria")) {
    taxa[taxa == "Epsilonproteobacteria"] <- "Epsilonproteobacteria *"
    legend("topright", taxa[1:6], pch = pch[1:6], col = col[1:6], pt.bg = col[1:6], cex = 0.9, bg = "white")
    legend("bottomleft", taxa[7:len], pch = pch[7:len], col = col[7:len], pt.bg = col[7:len], cex = 0.9, bg = "white")
  } else if(identical(groups, "majorphyla")) {
    legend <- c("Cellular", taxa[1:6], "Viruses", taxa[7:11])
    pch <- c(NA, pch[1:6], NA, pch[7:11])
    col <- c(NA, col[1:6], NA, col[7:11])
    legend("bottomleft", legend, text.font = c(2, 1,1,1,1,1,1, 2, 1,1,1,1,1), pch = pch, col = col, pt.bg = col, cex = 0.9, bg = "white")
  } else if(!is.null(legend.x) & !identical(legend.x, NA)) legend(legend.x, taxa, pch = pch, col = seq_along(taxa), pt.bg = seq_along(taxa), cex = 0.9, bg = "white")
  if(identify) identify(ZC, nH2O, names)
  # Return values invisibly 20210603
  invisible(vals)
}

# Plot compositional metrics for all samples in a study 20200901
plotcomp <- function(study, cn = FALSE, identify = FALSE, title = TRUE, xlim = NULL, ylim = NULL,
  plot.it = TRUE, points = TRUE, lines = FALSE, lineage = NULL, pch1 = 1, pch2 = 21) {
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
    plot(xlim, ylim, xlab = canprot::cplab$ZC, ylab = canprot::cplab$nH2O, type = "n")
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
  }

  # Calculate sample group values 20210607
  i1 <- pch %in% pch1
  i2 <- pch %in% pch2
  metrics <- getmetrics(study, mdat = mdat, RDP = RDP, lineage = lineage, groups = list(i1, i2))
  group <- list()
  if(!is.null(pch2) & !is.null(pch1) & sum(i2) > 0 & sum(i1) > 0) {
    group <- list(ZC1 = metrics$ZC[1], ZC2 = metrics$ZC[2], nH2O1 = metrics$nH2O[1], nH2O2 = metrics$nH2O[2])
    if(plot.it) {
      col1 <- mdat$col[mdat$pch == pch1]
      col2 <- mdat$col[mdat$pch == pch2]
      points(group$ZC1, group$nH2O1, pch = 8, cex = 2, lwd = 4, col = "white")
      points(group$ZC1, group$nH2O1, pch = 8, cex = 2, lwd = 2, col = col1)
      points(group$ZC2, group$nH2O2, pch = 8, cex = 2, lwd = 4, col = "white")
      points(group$ZC2, group$nH2O2, pch = 8, cex = 2, lwd = 2, col = col2)
    }
  }

  invisible(list(study = study, nH2O = nH2O, ZC = ZC, pch = pch, col = col, group = group))
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

# Get abundances and compositional metrics for taxonomic groups
# to compare samples (within a study or between studies) 20210606
getgroup <- function(study = "XDZ+17", metric = "nH2O", rank = "domain", pch1 = 21, pch2 = 24,
  minpercent = 2, study2 = NA, mdat = NULL, map = NULL, RDP = NULL) {

  # Get metadata, RDP and taxonomy mapping
  if(is.null(mdat)) mdat <- getmdat(study)
  mdat <- mdat[, c("study", "name", "Run", "sample", "pch", "col")]
  # Get data to compare two studies 20210513
  if(!is.na(study2)) {
    mdat2 <- getmdat(study2)[, c("study", "name", "Run", "sample", "pch", "col")]
    mdat$pch <- pch2
    mdat2$pch <- pch1
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
  # Drop samples that are excluded by getRDP (because of low counts or NA name) 20210606
  mdat <- mdat[mdat$Run %in% colnames(RDP), ]
  # Identify samples in each group
  i1 <- mdat$pch %in% pch1
  i2 <- mdat$pch %in% pch2
  # Stop if no samples match 20210608
  pchavail <- paste(unique(mdat$pch), collapse = " ")
  if(!any(i1)) stop(paste0("pch1 = ", pch1, " doesn't match any samples (available values are ", pchavail, ")"))
  if(!any(i2)) stop(paste0("pch2 = ", pch2, " doesn't match any samples (available values are ", pchavail, ")"))
  # Retrieve colors for points
  col1 <- mdat[i1, ]$col[1]
  col2 <- mdat[i2, ]$col[1]
  # Read compositional metrics for faster running
  datadir <- system.file("extdata/comp16S", package = "JMDplots")
  RefSeq_metrics <- read.csv(file.path(datadir, "RefSeq_metrics.csv"), as.is = TRUE)
  # Calculate the change in compositional metric for the whole community
  metrics <- getmetrics(study, mdat = mdat, RDP = RDP, map = map, metrics = RefSeq_metrics, groups = list(i1, i2))
  # Get selected compositional metric
  if(metric == "nH2O") X <- metrics$nH2O
  if(metric == "ZC") X <- metrics$ZC
  Xwhole <- X

  # Split the lineage text
  lsplit <- strsplit(RDP$lineage, ";")
  # Find the taxa with the specified rank in the lineage
  irank <- vapply(lsplit, function(x) match(rank, x), 0) - 1
  RDPtaxa <- mapply("[", lsplit, irank)
  # Calculate the compositional metrics for each unique taxon
  taxa <- na.omit(unique(RDPtaxa))
#  # Include "Other" 20210608
#  taxa <- c(taxa, "Other")
  X2 <- X1 <- P2 <- P1 <- Pboth <- numeric()
  taxon <- character()
  allI <- logical(length(RDPtaxa))
  for(j in seq_along(taxa)) {
    # Which organisms (by RDP classification) are in this taxon
    if(taxa[j] == "Other") iRDP <- !allI else {
      iRDP <- RDPtaxa == taxa[j]
      iRDP[is.na(iRDP)] <- FALSE
      allI <- allI | iRDP
    }
    thisRDP <- RDP[iRDP, ]
    thismap <- map[iRDP]
    # Calculate percent abundance of this taxon in the whole community
    thispercent <- sum(thisRDP[, -(1:3)]) / sum(RDP[, -(1:3)]) * 100
    # Skip low-abundance taxa
    if(taxa[j] != "Other") if(thispercent < minpercent) next
    # If we got here, print message about taxon name and abundance
    print(paste0("getgroup: ", rank, "_", taxa[j], " (", round(thispercent), "%)"))
    # Skip taxa with no available mappings
    if(all(is.na(thismap))) {
      print("            --- no taxonomic mappings available!")  
      next
    }
    # Keep percentages and names of used taxa
    taxon <- c(taxon, taxa[j])
    Pboth <- c(Pboth, thispercent)

    # Calculate the compositional metrics
    metrics <- getmetrics(study, mdat = mdat, RDP = thisRDP, map = thismap, metrics = RefSeq_metrics, groups = list(i1, i2))
    # Get selected compositional metric
    if(metric == "nH2O") X <- metrics$nH2O
    if(metric == "ZC") X <- metrics$ZC
    X1 <- c(X1, X[1])
    X2 <- c(X2, X[2])
    # Calculate percent abundance within the sample groups 20210520
    all1 <- RDP[, which(i1) + 3]
    this1 <- thisRDP[, which(i1) + 3]
    P1 <- c(P1, sum(this1) / sum(all1) * 100)
    all2 <- RDP[, which(i2) + 3]
    this2 <- thisRDP[, which(i2) + 3]
    P2 <- c(P2, sum(this2) / sum(all2) * 100)
  }

  # Calculate change in metric contributed by each taxon 20210606
  # Replace NA values for metric with 0 (where a taxon has zero abundance in one of the sample groups)
  x1 <- X1
  x1[is.na(x1)] <- 0
  x2 <- X2
  x2[is.na(x2)] <- 0
  DX <- (x2 - Xwhole[1]) * P2 / 100 - (x1 - Xwhole[1]) * P1 / 100
  names(DX) <- taxon
  # Report total contribution 
  DXpercent <- round(sum(DX / diff(Xwhole)) * 100, 2)
  message(paste0("getgroup: total contribution to \u0394", metric, " by individual taxa is ", DXpercent, "% of whole"))
  # Assemble output
  out <- list(study, metric, pch1, pch2, col1, col2, X1, X2, P1, P2, Pboth, DX, Xwhole, taxon)
  names(out) <- c("study", "metric", "pch1", "pch2", "col1", "col2", "X1", "X2", "P1", "P2", "Pboth", "DX", "Xwhole", "taxon")
  out

}

# Composition-abundance plots for taxonomic groups within sample groups 20210520
groupcomp <- function(..., xlim = NULL, ylim = NULL, xadj = NULL, yadj = NULL) {

  gg <- getgroup(...)

  # Plot the compositions and abundances
  if(is.null(xlim)) xlim <- range(na.omit(c(gg$P2, gg$P1)))
  if(is.null(ylim)) ylim <- range(na.omit(c(gg$X2, gg$X1)))
  if(identical(xlim, c(0, 100))) {
    plot(extendrange(xlim), ylim, type = "n", xlab = "Abundance (%)", ylab = canprot::cplab[[gg$metric]], xaxt = "n")
    axis(1, c(0, 50, 100))
  } else plot(xlim, ylim, type = "n", xlab = "Abundance (%)", ylab = canprot::cplab[[gg$metric]])
  for(k in seq_along(gg$taxon)) {
    # Add points for sample groups
    cex <- 1.5
    if(gg$pch1 > 20) points(gg$P1[k], gg$X1[k], pch = gg$pch1, bg = gg$col1, cex = cex)
    else points(gg$P1[k], gg$X1[k], pch = gg$pch1, col = gg$col1, cex = cex)
    if(gg$pch2 > 20) points(gg$P2[k], gg$X2[k], pch = gg$pch2, bg = gg$col2, cex = cex)
    else points(gg$P2[k], gg$X2[k], pch = gg$pch2, col = gg$col2, cex = cex)
    # Add arrow connecting the points
    arrows(gg$P1[k], gg$X1[k], gg$P2[k], gg$X2[k], length = 0.1)
    # Pad labels with spaces to offset from points
    label <- paste0("  ", gg$taxon[k], "  ")
    # Get adjustment from arguments if provided
    adj <- c(0, 0.5)
    if(!is.null(xadj)) if(!is.na(xadj[gg$taxon[k]])) adj[1] <- xadj[gg$taxon[k]]
    if(!is.null(yadj)) if(!is.na(yadj[gg$taxon[k]])) adj[2] <- yadj[gg$taxon[k]]
    # Add group name with spaces to offset from points
    text(gg$P2[k], gg$X2[k], label, adj = adj)
  }
  # Add lines for total composition of the plotted taxonomic groups
  OK1 <- !is.na(gg$X1)
  X1 <- gg$X1[OK1]
  P1 <- gg$P1[OK1]
  total1 <- sum(X1 * P1 / sum(P1))
  OK2 <- !is.na(gg$X2)
  X2 <- gg$X2[OK2]
  P2 <- gg$P2[OK2]
  total2 <- sum(X2 * P2 / sum(P2))
  abline(h = total1, lty = 3, lwd = 2, col = gg$col1)
  abline(h = total2, lty = 2, lwd = 2, col = gg$col2)
  # Calculate and return total percentage of community represented by these taxa
  Pboth <- gg$Pboth[OK2 & OK1]
  sum(Pboth)
}

# Plot contribution to difference of ZC or nH2O (percent of total) vs abundance for taxonomic groups 20210606
groupperc <- function(..., xlim = NULL, ylim = NULL, xadj = NULL, yadj = NULL) {
  # Calculate per-taxon contributions to total change of ZC
  gg <- getgroup(...)
  DXpercent <- gg$DX / diff(gg$Xwhole) * 100
  # Start plot
  if(is.null(xlim)) xlim <- range(gg$P1, gg$P2)
  if(is.null(ylim)) ylim <- range(DXpercent)
  plot(xlim, ylim, xlab = "Abundance (%)", ylab = paste("Contribution to", gg$metric, "change (%)"), type = "n")
  # Loop over taxa
  for(k in seq_along(gg$taxon)) {
    # Add points for sample groups
    cex <- 1.5
    if(gg$pch2 > 20) points(gg$P2[k], DXpercent[k], pch = gg$pch2, bg = gg$col2, cex = cex)
    else points(gg$P2[k], DXpercent[k], pch = gg$pch2, col = gg$col2, cex = cex)
    if(gg$pch1 > 20) points(gg$P1[k], DXpercent[k], pch = gg$pch1, bg = gg$col1, cex = cex)
    else points(gg$P1[k], DXpercent[k], pch = gg$pch1, col = gg$col1, cex = cex)
    # Add arrow connecting the points
    arrows(gg$P1[k], DXpercent[k], gg$P2[k], DXpercent[k], length = 0.1)
    # Add labels
    label <- gg$taxon[k]
  #  # Get adjustment from arguments if provided
  #  adj <- c(05, 0.5)
  #  if(!is.null(xadj)) if(!is.na(xadj[gg$Xtaxa[k]])) adj[1] <- xadj[gg$Xtaxa[k]]
  #  if(!is.null(yadj)) if(!is.na(yadj[gg$Xtaxa[k]])) adj[2] <- yadj[gg$Xtaxa[k]]
  #  # Add group name with spaces to offset from points
  #  text(gg$P2[k], gg$Xup[k], label, adj = adj)
    labx <- mean(c(gg$P1[k], gg$P2[k]))
    dy <- diff(ylim) / 30
    if(DXpercent[k] > 50) dy <- -dy
    text(labx, DXpercent[k] + dy, label)
  }
  invisible(gg)
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
  # Keep metadata only for samples with sufficient counts 20201001
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
  i2 <- iminuend[col == 2]
  i1 <- isubtrahend[col == 2]
  D_mean_nH2O <- mean(nH2O[i2]) - mean(nH2O[i1])
  D_mean_ZC <- mean(ZC[i2]) - mean(ZC[i1])

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
    p.nH2O <- t.test(nH2O[i1], nH2O[i2], paired = TRUE)$p.value
    log10p.nH2O <- formatC(log10(p.nH2O), 1, format = "f")
    p.ZC <- t.test(ZC[i1], ZC[i2], paired = TRUE)$p.value
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
