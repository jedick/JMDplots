# JMDplots/chem16S_plot.R
# Plotting functions for reference (RefSeq) and estimated community proteomes
# Separated from chem16S.R 20210607

# taxacomp()               # nH2O-ZC plot for taxa (default: Bacteria, Archaea) and their children
# plotmet("SVH+19")        # nH2O-ZC plot for specified study

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

  # Read chemical metrics of all taxa
  datadir <- system.file("extdata/chem16S", package = "JMDplots")
  metrics <- read.csv(file.path(datadir, "taxon_metrics.csv"), as.is = TRUE)
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
    # Get the chemical metrics for this group
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
      # Get the chemical metrics for all children
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
      ihalo <- match(c("Methanococci", "Archaeoglobi", "Thermococci", "Halobacteria", "Nanohaloarchaea"), children$group)
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

# Plot chemical metrics for all samples in a study 20200901
plotmet <- function(study, cn = FALSE, identify = FALSE, title = TRUE, xlim = NULL, ylim = NULL,
  plot.it = TRUE, points = TRUE, lines = FALSE, lineage = NULL, mincount = 200, pch1 = 1, pch2 = 21, dropNA = TRUE,
  return = "data", extracolumn = NULL, add = FALSE, plot.bg = TRUE) {
  # Get amino acid composition for samples
  mdat <- getmdat(study, dropNA = dropNA)
  RDP <- getRDP(study, cn = cn, mdat = mdat, lineage = lineage, mincount = mincount)
  metrics <- getmetrics(study, mdat = mdat, RDP = RDP, lineage = lineage, mincount = mincount)
  # Keep metadata only for samples with >= mincount counts 20201006
  mdat <- mdat[mdat$Run %in% metrics$Run, ]
  pch <- mdat$pch
  col <- mdat$col
  # Get nH2O and ZC
  nH2O <- metrics$nH2O
  ZC <- metrics$ZC

  if(plot.it) {
    # Get axis limits, excluding values of non-plotted points 20210820
    # Also exclude NA values (for Bison Pool site Q with lineage = "Archaea") 20210916
    if(is.null(xlim)) xlim <- range(na.omit(ZC[!is.na(pch)]))
    if(is.null(ylim)) ylim <- range(na.omit(nH2O[!is.na(pch)]))
    # Start plot
    if(!add) plot(xlim, ylim, xlab = canprot::cplab$ZC, ylab = canprot::cplab$nH2O, type = "n")
    if(points) {
      # Add background nH2O-ZC correlation (from basis species)
      if(plot.bg) lmlines()
      # Plot points for samples
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

  i1 <- pch %in% pch1
  i2 <- pch %in% pch2
  means <- list()
  if(!is.null(pch2) & !is.null(pch1) & sum(i2) > 0 & sum(i1) > 0) {
    # Calculate mean of sample values 20201003
    means <- list(ZC1 = mean(ZC[i1]), ZC2 = mean(ZC[i2]), nH2O1 = mean(nH2O[i1]), nH2O2 = mean(nH2O[i2]))
#    # Calculate values for aggregated samples 20210607
#    metrics <- getmetrics(study, mdat = mdat, RDP = RDP, lineage = lineage, meanss = list(i1, i2))
#    means <- list(ZC1 = metrics$ZC[1], ZC2 = metrics$ZC[2], nH2O1 = metrics$nH2O[1], nH2O2 = metrics$nH2O[2])
    if(plot.it) {
      col1 <- na.omit(mdat$col[mdat$pch == pch1])[1]
      col2 <- na.omit(mdat$col[mdat$pch == pch2])[1]
      points(means$ZC1, means$nH2O1, pch = 8, cex = 2, lwd = 4, col = "white")
      points(means$ZC1, means$nH2O1, pch = 8, cex = 2, lwd = 2, col = col1)
      points(means$ZC2, means$nH2O2, pch = 8, cex = 2, lwd = 4, col = "white")
      points(means$ZC2, means$nH2O2, pch = 8, cex = 2, lwd = 2, col = col2)
    }
  }

  # Return either the means means or individual values 20210831
  if(return == "means") out <- means
  if(return == "data") {
    name <- na.omit(mdat$name)[1]
    out <- data.frame(study = study, name = name, metrics, pch = pch, col = col)
    if(!is.null(extracolumn)) {
      # Add an extra column (e.g. 'type') to the output 20210901
      extracols <- mdat[, extracolumn, drop = FALSE]
      out <- cbind(out, extracols)
    }
  }
  invisible(out)
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

########################
# Unexported functions #
########################

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
