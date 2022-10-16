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
  refseq <- read.csv(system.file("extdata/RefSeq/protein_refseq.csv.xz", package = "JMDplots"), as.is = TRUE)
  taxa <- read.csv(system.file("extdata/RefSeq/taxid_names.csv.xz", package = "JMDplots"), as.is = TRUE)
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
  datadir <- system.file("extdata/RefSeq", package = "chem16S")
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
  #> taxa <- read.csv(system.file("extdata/RefSeq/taxid_names.csv.xz", package = "JMDplots"), as.is = TRUE)
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
      ihalo <- match(c("Methanococci", "Archaeoglobi", "Thermococci", "Halobacteria"), children$group)
      dy <- ifelse(groups == "majorcellular", 0.0025, 0.005)
      dx <- c(0, 0, 0, 0.002)
      text(children$ZC[ihalo] + dx, children$nH2O[ihalo] + dy, c(1, 2, 3, 4))
    }
    # Label Clostridia 20200930
    if(thisgroup == "Firmicutes" & identical(groups, "majorcellular")) {
      iclos <- match("Clostridia", children$group)
      text(children$ZC[iclos], children$nH2O[iclos] + 0.0025, 5)
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

