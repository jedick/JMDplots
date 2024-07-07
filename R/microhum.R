# JMDplots/microhum.R
# Plots for chem16S paper 20220202
# Renamed to sars16S 20220907
# Moved to JMDplots 20230211
# Renamed to hum16S 20230530
# Renamed to microhum 20230723

# Plot symbols and colors for body sites 20221125
pch_Oral <- 21
pch_Nasal <- 22
pch_Skin <- 23
pch_Gut <- 24
pch_IBD <- 25
pch_UG <- 25
col_Oral <- "#D62728"
col_Nasal <- "#56B4E9"
col_Skin <- "#9467BD"
col_Gut <- "#E69F00"
col_IBD <- "#009E73"
col_UG <- "#9E9E9E"
# For semi-transparent symbol outline 20220122
black50 <- adjustcolor(1, alpha.f = 0.5)
black25 <- adjustcolor(1, alpha.f = 0.25)

# Location of data files
getdatadir <- function() {
  system.file("extdata/microhum", package = "JMDplots")
}

##########################
### Plotting Functions ###
##########################

# Consistency between shotgun metagenomes and community reference proteomes
# 20211218 First version: comparison of Zc used in geo16S paper
#          Based on samples used by Aßhauer et al. (2015) (Tax4Fun paper) with additions by Dick and Tan (2023)
# 20231218 Added human DNA screening and nH2O-nO2 plots for microhum paper
microhum_1 <- function(pdf = FALSE) {

  if(pdf) pdf("Figure_1.pdf", width = 10, height = 6)
  mat <- matrix(c(1,1,2,2,3,3,4,4, 0,0,0,0,0,0,0,0, 0,5,5,5,6,6,6,0), nrow = 3, byrow = TRUE)
  layout(mat, heights = c(1, 0.08, 1))
  par(mar = c(3.5, 3.5, 2, 1), mgp = c(2.5, 1, 0))
  par(cex = 0.8)

  # Get chemical metrics for community reference proteomes
  metrics <- getmetrics_microhum("HMP12")
  # Get sample metadata
  metadata <- getmdat_microhum("HMP12")
  # Define colors and symbols
  bg <- sapply(metadata$"Body site", switch, "Skin" = col_Skin, "Nasal cavity" = col_Nasal, "Oral cavity" = col_Oral, "GI tract" = col_Gut, "UG tract" = col_UG)
  pch <- sapply(metadata$"Body site", switch, "Skin" = pch_Skin, "Nasal cavity" = pch_Nasal, "Oral cavity" = pch_Oral, "GI tract" = pch_Gut, "UG tract" = pch_UG)

  # Define cutoff for number of proteins (% of MG reads)
  cutoff <- 40

  ## Top row: nO2 and nH2O of community reference proteome vs shotgun metagenome
  # Loop over metrics
  for(metric in c("nO2", "nH2O")) {
    
    if(metric == "Zc") {
      xlim <- c(-0.2, -0.08)
      ylim <- c(-0.2, -0.12)
      main <- "Carbon oxidation state (Zc)"
    }

    if(metric == "nO2") {
      xlim <- c(-0.8, -0.55)
      ylim <- c(-0.8, -0.6)
      main <- quote("Stoichiometric oxidation state (" * italic(n)[O[2]] * ")")
    }

    if(metric == "nH2O") {
      xlim <- c(-1.05, -0.72)
      ylim <- c(-0.8, -0.72)
      main <- quote("Stoichiometric hydration state (" * italic(n)[H[2]*O] * ")")
    }

    # Compute min/max limits for 1:1 line
    xylim <- c(min(xlim, ylim), max(xlim, ylim))

    # Loop over MG pipelines
    for(pipeline in c("no_screening", "screening")) {

      # Read MG datasets processed without or with contamination removal
      if(pipeline=="no_screening") {
        MG_aafile <- file.path(getdatadir(), "ARAST/HMP12_no_screening_aa.csv")
        xlab <- "Metagenomes - NO screening"
      }
      if(pipeline=="screening") {
        MG_aafile <- file.path(getdatadir(), "ARAST/HMP12_aa.csv")
        xlab <- "Metagenomes - WITH screening"
      }
      MG_aa <- read.csv(MG_aafile)

      # Put MG and 16S samples in same order
      i16S <- match(metadata$Run, metrics$Run)
      metrics <- metrics[i16S, ]
      iMG <- match(metadata$Metagenome, MG_aa$protein)
      MG_aa <- MG_aa[iMG, ]
      # Make sure the 16S and metagenomes are paired correctly
      stopifnot(all(na.omit(metrics$Run == metadata$Run)))
      stopifnot(all(MG_aa$protein == metadata$Metagenome))

      # Get chemical metric for community reference proteomes
      metric_16S <- metrics[, metric]
      # Get chemical metric for metagenomes
      # This executes canprot::Zc(), canprot::nO2(), or canprot::nH2O()
      metric_MG <- get(metric)(MG_aa)

      # Start plot
      if(metric == "nO2" & pipeline == "no_screening") ylab <- "Community reference proteomes" else ylab <- ""
      plot(xlim, ylim, type = "n", xlab = xlab, ylab = ylab)
      lines(xylim, xylim, lty = 2, col = "gray40")
      # Use open symbols (no color) for MG with low numbers of protein fragments 20231218
      ilow <- MG_aa$chains / MG_aa$ref * 100 < cutoff
      ilow[is.na(ilow)] <- TRUE
      mybg <- bg
      mybg[ilow] <- "transparent"
      col <- ifelse(ilow, black25, black50)
      points(metric_MG, metric_16S, pch = pch, bg = mybg, col = col)

      # Show R-squared values excluding MG with low numbers of protein fragments 20231219
      mylm <- lm(metric_16S[!ilow] ~ metric_MG[!ilow])
      R2 <- round(summary(mylm)$r.squared, 2)
      R2txt <- bquote(italic(R)^2 == .(R2))
      legend("bottomleft", legend = R2txt, bty = "n", cex = 0.85)

      if(pipeline == "no_screening") {
        if(metric == "nO2") {
          legend("topleft", c("Skin", "Nasal cavity", "Oral cavity", "GI tract", "UG tract"),
            pch = c(pch_Skin, pch_Nasal, pch_Oral, pch_Gut, pch_UG),
            pt.bg = c(col_Skin, col_Nasal, col_Oral, col_Gut, col_UG), col = black50, bty = "n", cex = 0.7)
          # Use mtext to put title over two plots
          mtext(main, at = -0.5, line = 0.5)
          label.figure("A", font = 2, cex = 1.8, yfrac = 0.97)
        }
        if(metric == "nH2O") {
          mtext(main, at = -0.65, line = 0.5)
          label.figure("B", font = 2, cex = 1.8, yfrac = 0.97)
        }
      }

      if(pipeline == "screening" & metric == "nO2") {
        cutoff_txt <- c(paste("<", cutoff, "% of reads"), paste(">", cutoff, "% of reads"))
        legend("bottomright", cutoff_txt, pch = pch_UG, pt.bg = c("transparent", col_UG), col = c(black25, black50),
          bty = "n", cex = 0.7, title = "Protein prediction rate")
        legend("bottomright", c("", ""), pch = pch_Nasal, pt.bg = c("transparent", col_Nasal), col = c(black25, black50),
          bty = "n", cex = 0.7, inset = c(0.46, 0))
      }

    }
  }

  ## Bottom row: nH2O vs nO2 for community reference proteomes and metagenomes
  xlim <- c(-0.77, -0.63)
  ylim <- c(-0.82, -0.72)
  par(mar = c(4, 4, 3, 1))
  par(cex.lab = 1.2)
  plotmet_microhum("HMP12", title = FALSE, pt.open.col = black50, xlim = xlim, ylim = ylim)
  legend("bottomright", c("Skin", "Nasal cavity", "Oral cavity", "GI tract", "UG tract"),
    pch = c(pch_Skin, pch_Nasal, pch_Oral, pch_Gut, pch_UG),
    pt.bg = c(col_Skin, col_Nasal, col_Oral, col_Gut, col_UG), col = black50, bty = "n", cex = 0.9)
  title("Community reference proteomes", font.main = 1)
  label.figure("C", font = 2, cex = 1.8, xfrac = 0.12, yfrac = 0.92)

  # Read amino acid composition of proteins inferred from metagenome
  MG_aafile <- file.path(getdatadir(), "ARAST/HMP12_aa.csv")
  MG_aa <- read.csv(MG_aafile)
  # Exclude samples with low numbers of proteins
  ilow <- MG_aa$chains / MG_aa$ref * 100 < cutoff
  ilow[is.na(ilow)] <- TRUE
  MG_aa[ilow, ] <- NA
  # Put into order of metadata table
  iMG <- match(metadata$Metagenome, MG_aa$protein)
  MG_aa <- MG_aa[iMG, ]
  plot(nO2(MG_aa), nH2O(MG_aa), xlim = xlim, ylim = ylim, xlab = chemlab("nO2"), ylab = chemlab("nH2O"), pch = pch, bg = bg, col = col)
  title("Metagenomes - WITH screening", font.main = 1)

  if(pdf) dev.off()

}

# Chemical metrics are broadly different among genera and are similar between GTDB and low-contamination genomes from UHGG
# 20231230
microhum_2 <- function(pdf = FALSE) {

  if(pdf) pdf("Figure_2.pdf", width = 12, height = 7)
  mat <- matrix(c(1,1,2,2, 3,4,5,6), byrow = TRUE, nrow = 2)
  layout(mat, heights = c(2, 1))
  par(mgp = c(2.9, 1, 0))
  par(cex.lab = 1.2)
  par(mar = c(5.1, 4.1, 3.1, 2.1))
  # Get genera from Figure 4
  genera <- readLines(system.file("extdata/microhum/Figure_5_genera.txt", package = "JMDplots"))
  # Get colors
  col <- hcl.colors("Dynamic", n = length(genera))
  # Loop over reference databases
  refdb <- c("GTDB_220", "UHGG_2.0.1")
  main <- c("GTDB (used for community reference proteomes)", "UHGG (contamination < 2% and completeness > 95%)")
  genus_vals <- list()
  metrics <- c("nO2", "nH2O")
  for(i in seq_along(refdb)) {
    # Make plot
    vals <- plot_starburst(genera, metrics, xlim = c(-0.81, -0.58), ylim = c(-0.83, -0.69), refdb = refdb[i], pch = NA, col = col)
    # Add points for genera
    gvals <- do.call(rbind, sapply(vals, "[", 1))
    X <- gvals$Xvals
    Y <- gvals$Yvals
    points(X, Y, pch = ".", cex = 2, col = "#00000080")
    # Get labels for genera
    labels <- paste0(" ", gvals$taxon)
    xadj <- rep(0, length(labels))
    yadj <- rep(0.5, length(labels))
    if(grepl("GTDB", refdb[i])) {
      yadj[gvals$taxon == "Bacteroides"] <- 0.2
      yadj[gvals$taxon == "Phocaeicola"] <- 1.2
      xadj[gvals$taxon == "Phocaeicola"] <- 0.1
      yadj[gvals$taxon == "Blautia_A"] <- 1.2
      xadj[gvals$taxon == "Blautia_A"] <- 0.5
      xadj[gvals$taxon == "Escherichia"] <- 0.2
      yadj[gvals$taxon == "Escherichia"] <- 1.4
      yadj[gvals$taxon == "Anaerostipes"] <- -0.3
      xadj[gvals$taxon == "Anaerostipes"] <- 0.08
      yadj[gvals$taxon == "Enterococcus_B"] <- 0.4
      yadj[gvals$taxon == "Haemophilus_D"] <- 0.8
      yadj[gvals$taxon == "Latilactobacillus"] <- -0.1
      xadj[gvals$taxon == "Lactobacillus"] <- 1.03
    }
    if(grepl("UHGG", refdb[i])) {
      yadj[gvals$taxon == "Finegoldia"] <- 0
      yadj[gvals$taxon == "Latilactobacillus"] <- 0
      yadj[gvals$taxon == "Campylobacter_B"] <- 0.8
      xadj[gvals$taxon == "Enterococcus_B"] <- 1.03
      yadj[gvals$taxon == "Lactobacillus"] <- 0.3
      yadj[gvals$taxon == "Bacteroides"] <- 0
      yadj[gvals$taxon == "Phocaeicola"] <- 1
    }
    for(j in seq_along(labels)) text(X[j], Y[j], labels[j], adj = c(xadj[j], yadj[j]), cex = 0.8)
    title(main = main[i], font.main = 1)
    # Keep points for correlation plot
    colnames(gvals)[colnames(gvals) == "Xvals"] <- metrics[1]
    colnames(gvals)[colnames(gvals) == "Yvals"] <- metrics[2]
    genus_vals[[i]] <- gvals
    label.figure(LETTERS[i], font = 2, cex = 2, xfrac = 0.025)
  }
  names(genus_vals) <- refdb

  # Plot nO2 and nH2O for GTDB vs UHGG
  par(mar = c(4, 4.1, 2.1, 2.1))
  for(i in seq_along(metrics)) {
    x <- genus_vals$UHGG[, metrics[i]]
    y <- genus_vals$GTDB[, metrics[i]]
    # Start plot
    plot(x, y, xlab = "UHGG", ylab = "GTDB", type = "n")
    # Compute min/max limits for 1:1 line
    xylim <- extendrange(c(min(x, y, na.rm = TRUE), max(x, y, na.rm = TRUE)))
    lines(xylim, xylim, lty = 2, col = "gray40")
    # Plot points
    points(x, y, pch = 16, col = col)
    # Show R-squared values
    mylm <- lm(y ~ x)
    R2 <- round(summary(mylm)$r.squared, 2)
    R2txt <- bquote(italic(R)^2 == .(R2))
    legend("bottomright", legend = R2txt, bty = "n", cex = 0.9)
    title(chemlab(metrics[i]))
    if(i == 1) {
      mtext("Genus reference proteomes", line = 1, adj = 3.8, cex = 0.8, xpd = NA)
      label.figure("C", font = 2, cex = 2)
    }
  }

  # Plot nO2 and nH2O for GTDB vs RDP
  # Classifications made using GTDB 16S rRNA training set (this study)
  met_microhum <- getmetrics_microhum("HMP12")
  mdat_microhum <- getmdat_microhum("HMP12", metrics = met_microhum)
  pch <- mdat_microhum$metadata$pch
  col <- adjustcolor(mdat_microhum$metadata$col, alpha.f = 0.6)
  # Classifications made using rDP 16S rRNA training set (geo16S paper)
  met_geo16S <- getmetrics_geo16S("HMP12")
  mdat_geo16S <- getmdat_geo16S("HMP12", metrics = met_geo16S)
  stopifnot(all.equal(mdat_microhum$metrics$Run, mdat_geo16S$metrics$Run))
  for(i in seq_along(metrics)) {
    x <- mdat_geo16S$metrics[, metrics[i]]
    y <- mdat_microhum$metrics[, metrics[i]]
    # Start plot
    plot(x, y, xlab = "RDP training set", ylab = "GTDB training set", type = "n")
    # Compute min/max limits for 1:1 line
    xylim <- extendrange(c(min(x, y), max(x, y)))
    lines(xylim, xylim, lty = 2, col = "gray40")
    # Plot points
    points(x, y, pch = pch, bg = col)
    # Show R-squared values
    mylm <- lm(y ~ x)
    R2 <- round(summary(mylm)$r.squared, 2)
    R2txt <- bquote(italic(R)^2 == .(R2))
    legend("bottomright", legend = R2txt, bty = "n", cex = 0.9)
    title(chemlab(metrics[i]))
    if(i == 1) {
      legend("topleft", c("Skin", "Nasal cavity", "Oral cavity", "GI tract", "UG tract"),
        pch = c(pch_Skin, pch_Nasal, pch_Oral, pch_Gut, pch_UG),
        pt.bg = c(col_Skin, col_Nasal, col_Oral, col_Gut, col_UG), col = black50, bty = "n", cex = 0.7)
      mtext("HMP 16S rRNA samples", line = 1, adj = 2.8, cex = 0.8, xpd = NA)
      label.figure("D", font = 2, cex = 2)
    }
  }

  if(pdf) dev.off()

}

# Chemical variation of microbial proteins across body sites (multi-omics comparison) 20221125
microhum_3 <- function(pdf = FALSE) {

  if(pdf) pdf("Figure_3.pdf", width = 7, height = 6)
  par(mfrow = c(2, 2))
  par(mgp = c(2.5, 1, 0))

  ## Panel A: Community reference proteomes for body sites
  ## based on 16S rRNA gene sequences from Boix-Amoros et al. (2021) 20220814
  xlim <- c(-0.77, -0.59)
  ylim <- c(-0.8, -0.72)

  # nH2O-nO2 plot for all samples
  par(mar = c(4, 4, 3, 1))
  # Plot data for "no treatment" samples
  BodySites <- plotmet_microhum("BPB+21_NoTreatment", extracolumn = c("Subject", "Site"),
    title = FALSE, pt.open.col = black50, xlim = xlim, ylim = ylim)
  # Keep samples for subjects (not controls)
  BodySites <- BodySites[BodySites$Site != "Control", ]
  # Draw convex hull around all samples 20221125
  BShull <- chull(BodySites$nO2, BodySites$nH2O)
  polygon(BodySites$nO2[BShull], BodySites$nH2O[BShull], border = 8, lty = 2)

  # Plot p-values 20230204
  # Use metrics for untreated samples
  met <- getmetrics_microhum("BPB+21_NoTreatment")
  met <- met[match(BodySites$Run, met$Run), ]
  # Use metrics for gut and oral samples
  igut <- BodySites$Site == "feces"
  ioral <- BodySites$Site == "Oral cavity"
  plot.p.values(met$nO2[igut], met$nO2[ioral], met$nH2O[igut], met$nH2O[ioral], ypos = "bottom")

  # Make legend
  legend("topright", legend = c("Skin", "Nasal", "Oral", "Gut"), pch = c(pch_Skin, pch_Nasal, pch_Oral, pch_Gut),
    pt.bg = c(col_Skin, col_Nasal, col_Oral, col_Gut), col = black50, bty = "n", cex = 0.8)
  title(hyphen.in.pdf("Community reference proteomes\n(data from Boix-Amor\u00f3s et al., 2021)"), font.main = 1)
  label.figure("A", font = 2, cex = 1.8, yfrac = 0.97)

  ## Panel B: Community reference proteomes for controls in COVID-19 datasets 20220822
  # Setup plot
  plot(xlim, ylim, xlab = cplab$nO2, ylab = cplab$nH2O, type = "n")
  # Colors and point symbols for sample types
  col <- list(oro = col_Oral, naso = col_Nasal, gut = col_Gut)
  col <- lapply(col, adjustcolor, alpha.f = 0.8)
  pch <- list(oro = pch_Oral, naso = pch_Nasal, gut = pch_Gut)
  # Read precomputed mean values 20220823
  means <- read.csv(file.path(getdatadir(), "16S/dataset_metrics_all.csv"))
  # Loop over body sites
  for(type in c("gut", "oro", "naso")) {
    itype <- means$type == type
    # Plot points for control samples
    points(means$nO2_dn[itype], means$nH2O_dn[itype], pch = pch[[type]], bg = col[[type]], col = black50)
  }
  # Plot convex hull from Panel A
  polygon(BodySites$nO2[BShull], BodySites$nH2O[BShull], border = 8, lty = 2)
  # Plot p-values 20230204
  gut <- subset(means, type == "gut")
  oro <- subset(means, type == "oro")
  plot.p.values(gut$nO2_dn, oro$nO2_dn, gut$nH2O_dn, oro$nH2O_dn, ypos = "bottom")
  legend("topright", legend = c("Nasopharyngeal", "Oropharyngeal", "Gut"), pch = c(pch_Nasal, pch_Oral, pch_Gut),
    pt.bg = c(col_Nasal, col_Oral, col_Gut), col = black50, bty = "n", cex = 0.8)
  title(hyphen.in.pdf("Community reference proteomes\n(controls in COVID-19 studies)"), font.main = 1)
  label.figure("B", font = 2, cex = 1.8, yfrac = 0.96)

  ## Panel C: Metagenomes from different body sites 20221124
  ylim <- c(-0.84, -0.60)

  # Setup plot
  plot(xlim, ylim, xlab = cplab$nO2, ylab = cplab$nH2O, type = "n")
  # Studies are for gut, oral, nasal
  studies <- c("ZZL+20", "CZH+22", "LLZ+21")
  pchs <- c(pch_Gut, pch_Oral, pch_Nasal)
  cols <- c(col_Gut, col_Oral, col_Nasal)
  cols <- sapply(cols, adjustcolor, alpha.f = 0.8)
  for(i in 1:length(studies)) {
    file <- file.path(getdatadir(), "ARAST", paste0(studies[i], "_aa.csv"))
    aa <- read.csv(file)
    # Calculate chemical metrics
    nO2 <- nO2(aa)
    nH2O <- nH2O(aa)
    # Identify samples with low numbers of predicted proteins 20240102
    cutoff <- 40
    ilow <- aa$chains / aa$ref * 100 < cutoff
    mybg <- rep(cols[i], nrow(aa))
    mybg[ilow] <- "transparent"
    col <- ifelse(ilow, black25, black50)
    points(nO2, nH2O, pch = pchs[i], bg = mybg, col = col)
    if(studies[i] == "ZZL+20") gut <- list(nO2 = nO2[!ilow], nH2O = nH2O[!ilow])
    if(studies[i] == "CZH+22") oro <- list(nO2 = nO2[!ilow], nH2O = nH2O[!ilow])
  }
  # Plot convex hull from Panel A
  polygon(BodySites$nO2[BShull], BodySites$nH2O[BShull], border = 8, lty = 2)
  # Plot p-values 20230204
  plot.p.values(gut$nO2, oro$nO2, gut$nH2O, oro$nH2O, ypos = "bottom")
  # Add legend
  legend <- c("Nasopharyngeal (Liu'21)", "Oropharyngeal (de Castilhos'22)", "Gut (Zuo'20)")
  legend("topright", legend, pch = rev(pchs), pt.bg = c(col_Nasal, col_Oral, col_Gut), col = black50, bty = "n", cex = 0.8)
  legend("topright", legend = character(3), pch = rev(pchs), pt.bg = "transparent", col = black25, bty = "n", cex = 0.8, inset = c(0.68, 0))
  legend("topright", legend = c("", "", "", "<40 >40  Protein prediction rate (%)"), bty = "n", cex = 0.8, inset = c(0.115, 0))
  title(hyphen.in.pdf("Proteins from metagenomes\n(controls and COVID-19 patients)"), font.main = 1)
  label.figure("C", font = 2, cex = 1.8, yfrac = 0.96)

  ## Panel D: Metaproteomes from various body sites 20221114
  # Setup plot
  plot(xlim, ylim, xlab = cplab$nO2, ylab = cplab$nH2O, type = "n")
  # Define studies and point symbols
  studies_MP <- c(
    "TWC+22", "MLL+17", # Gut
    "GNT+21_cells", "JZW+22" # Oral
  )
  pchs <- c(
    pch_Gut, 25,
    pch_Oral, pch_Oral
  )
  cols <- c(
    col_Gut, col_Gut,
    col_Oral, col_Oral
  )
  cols <- sapply(cols, adjustcolor, alpha.f = 0.8)
  cexs <- ifelse(studies_MP %in% c("JZW+22"), 0.7, 1)
  # Store all data to calculate p-values 20230203
  gut.nO2 <- gut.nH2O <- oral.nO2 <- oral.nH2O <- numeric()
  # Loop over studies
  for(i in 1:length(studies_MP)) {
    # Get amino acid composition from metaproteome
    studydir <- strsplit(studies_MP[i], "_")[[1]][1]
    aa <- read.csv(file.path(getdatadir(), "metaproteome", studydir, paste0(studies_MP[i], "_aa.csv")))
    # Calculate nO2 and nH2O
    nO2 <- nO2(aa)
    nH2O <- nH2O(aa)
    # Add points
    points(nO2, nH2O, pch = pchs[i], bg = cols[i], cex = cexs[i], col = black50)
    if(studies_MP[i] %in% c("TWC+22", "MLL+17")) gut.nO2 <- c(gut.nO2, nO2)
    if(studies_MP[i] %in% c("TWC+22", "MLL+17")) gut.nH2O <- c(gut.nH2O, nH2O)
    if(studies_MP[i] %in% c("GNT+21_cells", "JZW+22")) oral.nO2 <- c(oral.nO2, nO2)
    if(studies_MP[i] %in% c("GNT+21_cells", "JZW+22")) oral.nH2O <- c(oral.nH2O, nH2O)
  }
  # Plot convex hull from Panel A
  polygon(BodySites$nO2[BShull], BodySites$nH2O[BShull], border = 8, lty = 2)
  # Plot p-values 20230204
  plot.p.values(gut.nO2, oral.nO2, gut.nH2O, oral.nH2O, ypos = "bottom")
  # Add legend
  legend("topright", c("Jiang'22", "Granato'21"), title = "Oral",
    pch = c(pch_Oral, pch_Oral), pt.bg = col_Oral, col = black50,
    pt.cex = c(0.7, 1), cex = 0.8, bg = "transparent")
  legend("topleft", c(hyphen.in.pdf("Thuy-Boun'22"), "Maier'17"), title = "Gut", title.adj = 0.4,
    pch = c(pch_Gut, 25), pt.bg = col_Gut, col = black50,
    cex = 0.8, bg = "transparent", inset = c(0, 0.2))
  title(hyphen.in.pdf("Metaproteomes\n(controls and patients, not COVID-19)"), font.main = 1)
  label.figure("D", font = 2, cex = 1.8, yfrac = 0.96)

  if(pdf) dev.off()

}

# Differences of chemical metrics between controls and COVID-19 or IBD patients 20220806
microhum_4 <- function(pdf = FALSE) {

  # Start plot
  if(pdf) pdf("Figure_4.pdf", width = 13, height = 12)
  mat <- matrix(c(
    1,1,1,1,1,1,1,1,     2,2,2,2,2,2,2,2,     3,3,3,3,3,3,3,3,
    4,4,4, 5,5,5, 6,6,6, 7,7,7, 8,8,8, 9,9,9, 10,10,10, 11,11,11,
    12,12,12,12,12,12,12,12, 13,13,13,13, 14,14,14,14, 15,15,15,15, 16,16,16,16
    ), nrow = 3, byrow = TRUE
  )
  layout(mat)

  # Define plot settings
  par(cex = 1.2)
  par(mgp = c(2.5, 1, 0))
  startplot <- function(xlim = c(-0.01, 0.01), ylim = c(-0.015, 0.020)) {
    plot(xlim, ylim, type  = "n", pch = ".", xlab = quote(Delta*italic(n)[O[2]]), ylab = quote(Delta*italic(n)[H[2]*O]))
    abline(h = 0, v = 0, lty = 2, col = 8)
  }
  # Colors and point symbols for sample types
  col <- list(naso = col_Nasal, oro = col_Oral, gut = col_Gut, IBD = col_IBD)
  col <- lapply(col, adjustcolor, alpha.f = 0.8)
  pch <- list(naso = pch_Nasal, oro = pch_Oral, gut = pch_Gut, IBD = pch_IBD)
  # Read precomputed mean values 20220823
  means <- read.csv(file.path(getdatadir(), "16S/dataset_metrics_all.csv"))

  ## Panel A: nH2O-nO2 plots for nasopharyngeal, oral/oropharyngeal, and gut communities
  par(mar = c(4, 4, 3, 1))
  for(type in c("naso", "oro", "gut")) {

    if(type == "naso") startplot(c(-0.01, 0.017)) else if(type == "oro") startplot(c(-0.017, 0.01)) else startplot()
    itype <- means$type == type
    label <- 1:sum(itype)

    # Add points
    points(means$D_nO2[itype], means$D_nH2O[itype], pch = pch[[type]], col = black50, bg = col[[type]])
    dx <- rep(0, sum(itype))
    dy <- 0.002
    if(type == "gut") dy <- 0.0022
    dy <- rep(dy, sum(itype))
    if(type == "naso") {
      dx[5] <- 0.0002
    }
    if(type == "gut") {
      dy[c(2, 5, 6, 9)] <- -0.0018
      dx[5] <- -0.0002
      dx[6] <- 0.0003
      dy[7] <- 0
      dx[7] <- -0.0007
      dx[8] <- 0.0003
    }
    # Label points
    text(means$D_nO2[itype] + dx, means$D_nH2O[itype] + dy, label, cex = 0.8)
    # Plot p-values 20230204
    plot.p.values(means$nO2_dn[itype], means$nO2_up[itype], means$nH2O_dn[itype], means$nH2O_up[itype], paired = TRUE)
    # Add plot title
    titles <- hyphen.in.pdf(c(naso = "Nasopharyngeal (COVID-19)", oro = "Oropharyngeal (COVID-19)", gut = "Gut (COVID-19)"))
    title(titles[type], font.main = 1, line = 0.5)
    # Add panel title
    start <- hyphen.in.pdf("A. Community reference proteomes (")
    covid <- hyphen.in.pdf("COVID-19")
    label <- bquote(bold(.(start) * Delta == .(covid) ~ "minus control)"))
    if(type == "naso") label.figure(label, font = 2, yfrac = 0.95, adj = 0.02)
  }

  # Common nH2O and nO2 limits for boxplots
  ylims <- list(
    nH2O = c(-0.87, -0.68),
    nO2 = c(-0.75, -0.65)
  )

  # Function to make boxplots
  boxplotfun <- function(metric, x_list, ylim, squeeze = FALSE) {
    if(squeeze) opar <- par(mar = c(4, 2.5, 3.5, 0.5)) else opar <- par(mar = c(4, 4, 3, 1))
    # Add number of samples to group names
    len <- sapply(x_list, length)
    labels <- paste0(names(x_list), " (", len, ")")
    names(x_list) <- ""
    # Make boxplot
    boxplot(x_list, ylim = ylim, ylab = cplab[[metric]], xlab = "")
    # Make rotated labels (modified from https://www.r-bloggers.com/rotated-axis-labels-in-r-plots/)
    text(x = (1:2)+0.5, y = par()$usr[3] - 1.5 * strheight("A"), labels = labels, srt = 25, adj = 1, xpd = NA)
    # Add p-value
    x_pvalue <- wilcox.test(x_list[[1]], x_list[[2]])$p.value
    legend <- bquote(italic(p) == .(format(signif(x_pvalue, 1), scientific = 2)))
    if(x_pvalue < 0.05) legend <- bquote(bolditalic(p) == bold(.(format(signif(x_pvalue, 1), scientific = 2))))
    if(squeeze) inset <- c(-0.25, -0.05) else inset <- c(-0.2, -0.05)
    legend("bottomleft", legend = legend, bty = "n", inset = inset, cex = 0.9)
    # Show median difference
    x_diff <- median(x_list[[2]]) - median(x_list[[1]])
    diffval <- signif(x_diff, 2)
    if(metric == "nO2") difftxt <- bquote(Delta*italic(n)[O[2]] == .(diffval))
    if(metric == "nH2O") difftxt <- bquote(Delta*italic(n)[H[2]*O] == .(diffval))
    if(squeeze) cex.main <- 0.9 else cex.main <- 1
    if(squeeze) adj <- 0.6 else adj <- 0.5
    title(difftxt, line = 0.7, cex.main = cex.main, adj = adj, xpd = NA)
    par(opar)
  }

  ## Panel B: Boxplots for nH2O and nO2 in fecal MAGs 20221029
  # Read amino acid compositions of proteins predicted from MAGs
  aa <- read.csv(file.path(getdatadir(), "KWL22/KWL22_MAGs_prodigal_aa.csv.xz"))
  # BioSample metadata
  dat <- read.csv(file.path(getdatadir(), "KWL22/BioSample_metadata.csv"))
  # Loop over SRA run prefix to choose study:
  # SRR1307: Yeoh et al. (PRJNA650244)
  # SRR1232: Zuo et al. (PRJNA624223)
  for(SRAprefix in c("SRR1307", "SRR1232")) {
    # Loop over metrics
    for(metric in c("nH2O", "nO2")) {
      # Get amino acid compositions for this BioProject
      iaa <- grep(SRAprefix, aa$protein)
      thisaa <- aa[iaa, ]
      # Calculate nO2 or nH2O
      x <- get(metric)(thisaa)
      ylim <- ylims[[metric]]
      if(metric == "nO2") ylim <- c(-0.85, -0.55)
      # Get names of groups
      idat <- match(thisaa$protein, dat$Run)
      group <- dat$Group[idat]
      # Make list of values in each group
      if(SRAprefix == "SRR1232") x_list <- list(Control = x[group == "Healthy_controls"], "COVID-19" = x[group == "COVID19"])
      if(SRAprefix == "SRR1307") x_list <- list(Control = x[group == "non-COVID-19"], "COVID-19" = x[group == "COVID-19"])
      # Make boxplot
      boxplotfun(metric, x_list, ylim, squeez = TRUE)
      if(metric == "nO2" & SRAprefix == "SRR1307") label.figure(hyphen.in.pdf("B. Gut MAGs (COVID-19)"), font = 2, adj = 0, yfrac = 1.02)
      if(metric == "nO2" & SRAprefix == "SRR1307") label.figure("(Ke et al., 2022; Yeoh et al., 2021)", font = 2, adj = 0)
      if(metric == "nO2" & SRAprefix == "SRR1232") label.figure("(Ke et al., 2022; Zuo et al., 2020)", font = 2, adj = 0)
    }
  }

  ## Panel C: boxplots of nO2 and nH2O of bacterial metaproteomes in control and COVID-19 patients 20220830
  # Limit to bacterial taxonomy
  taxonomy <- "Bacteria"
  # Loop over datasets
  for(study in c("HZX+21", "GPM+22")) {
    if(study == "HZX+21") {
      if(taxonomy == "All") aa <- read.csv(file.path(getdatadir(), "metaproteome/HZX+21/HZX+21_aa.csv"))
      if(taxonomy == "Bacteria") aa <- read.csv(file.path(getdatadir(), "metaproteome/HZX+21/HZX+21_bacteria_aa.csv"))
      # Identify control and COVID-19 patients
      icontrol <- grep("Ctrl", aa$organism)
      icovid <- grep("P", aa$organism)
    }
    if(study == "GPM+22") {
      if(taxonomy == "All") aa <- read.csv(file.path(getdatadir(), "metaproteome/GPM+22/GPM+22_aa.csv"))
      if(taxonomy == "Bacteria") aa <- read.csv(file.path(getdatadir(), "metaproteome/GPM+22/GPM+22_bacteria_aa.csv"))
      # Identify control and COVID-19 patients
      icontrol <- grep("negative", aa$abbrv)
      icovid <- grep("positive", aa$abbrv)
    }
    # Loop over nH2O and nO2
    for(metric in c("nH2O", "nO2")) {
      # Calculate nO2 or nH2O
      x <- get(metric)(aa)
      ylim <- ylims[[metric]]
      if(metric == "nO2" & study == "GPM+22") ylim <- c(-0.8, -0.6)
      x_list <- list(Control = x[icontrol], "COVID-19" = x[icovid])
      names(x_list)[2] <- hyphen.in.pdf(names(x_list)[2])
      boxplotfun(metric, x_list, ylim, squeeze = TRUE)
      if(metric == "nO2" & study == "HZX+21") label.figure(hyphen.in.pdf("C. Gut Metaproteomes (COVID-19)"), font = 2, yfrac = 1.02, xfrac = 0.5)
      if(metric == "nO2" & study == "HZX+21") label.figure("(He et al., 2021)", font = 2, xfrac = 0.25)
      if(metric == "nO2" & study == "GPM+22") label.figure("(Grenga et al., 2022)", font = 2, xfrac = 0.3)
    }
  }

  ## Panel D: nH2O-nO2 plots for community reference proteomes in IBD 20230723
  par(mar = c(4, 4, 3, 1))
  type <- "IBD"
  startplot(c(-0.030, 0.005), c(-0.01, 0.025))
  itype <- means$type == type
  label <- 1:sum(itype)
  # Add points
  points(means$D_nO2[itype], means$D_nH2O[itype], pch = pch[[type]], col = black50, bg = col[[type]])
  # Label points
  dx <- rep(0, sum(itype))
  dy <- rep(0.0018, sum(itype))
  dx[12] <- -0.0015
  dy[12] <- 0.001
  dx[14] <- 0.0002
  dx[8] <- 0.0012
  dy[8] <- -0.0005
  dx[4] <- 0.0014
  dy[4] <- 0.0005
  dx[c(10, 11, 15)] <- 0.0015
  dy[c(10, 11, 15)] <- -0.0015
  text(means$D_nO2[itype] + dx, means$D_nH2O[itype] + dy, label, cex = 0.8)
  # Plot p-values 20230204
  plot.p.values(means$nO2_dn[itype], means$nO2_up[itype], means$nH2O_dn[itype], means$nH2O_up[itype], paired = TRUE)
  # Add plot title
  title("Gut (IBD)", font.main = 1, line = 0.5)
  # Add panel title
  start <- hyphen.in.pdf("D. CRPs (")
  ibd <- hyphen.in.pdf("IBD")
  label <- bquote(bold(.(start) * Delta == .(ibd) ~ "minus control)"))
  label.figure(label, font = 2, yfrac = 0.95, adj = 0.04)

  ## Panel E: Boxplots for IBD metagenomes 20240102
  file <- file.path(getdatadir(), "ARAST/LAA+19_aa.csv")
  aa <- read.csv(file)
  # Exclude samples with low numbers of predicted proteins 20240102
  cutoff <- 40
  ilow <- aa$chains / aa$ref * 100 < cutoff
  aa <- aa[!ilow, ]
  for(disease in c("UC", "CD")) {
    for(metric in c("nH2O", "nO2")) {
      # Calculate nO2 or nH2O
      x <- get(metric)(aa)
      ylim <- ylims[[metric]]
      # Make list of values in control and patient groups
      x_list <- list(Control = x[aa$abbrv == "Control"], IBD = x[aa$abbrv == disease])
      names(x_list)[2] <- disease
      boxplotfun(metric, x_list, ylim)
      figlab <- hyphen.in.pdf("E. Gut metagenomes (Lloyd-Price et al., 2019)")
      if(metric == "nH2O" & disease == "UC") label.figure(figlab, font = 2, adj = 0)
    }
  }

  if(pdf) dev.off()
}

# Differences of relative abundances of genera between controls and patients 20231227
microhum_5 <- function(pdf = FALSE) {

  if(pdf) pdf("Figure_5.pdf", width = 10, height = 8.5)

  # First, get the genera with largest abundance differences in individual COVID-19 gut and IBD datasets
  covid <- get_abundance("gut")
  ibd <- get_abundance("IBD")
  # Get the unique genus names
  genus_names <- unique(c(colnames(covid), colnames(ibd)))
  # Then, get the abundance differences for these genera in COVID-19 gut and IBD datasets
  COVID <- get_abundance("gut", genus_names)
  IBD <- get_abundance("IBD", genus_names)

  # Setup plot
  layout(matrix(c(3, 1, 2)), heights = c(0.3, 1, 1))
  par(mar = c(7, 5, 0, 4))
  par(tcl = -0.3)
  par(las = 1)
  par(cex.lab = 1.2)
  par(mgp = c(3.2, 1, 0))
  # Define colors and breaks for heatmap
  col <- function(n) hcl.colors(n, "RdYlBu")
  breaks <- c(-.25, -.10, -.05, 0, .05, .10, .25)

  for(disease in c("COVID", "IBD")) {

    # Get differential abundance data
    D_abundance <- get(disease)
    # Convert to matrix to use plot.matrix
    D_abundance <- as.matrix(D_abundance)
    # Get amino acid compositions of reference proteomes for genera and calculate nO2
    AAcomp <- taxon_AA[["GTDB_220"]]
    AAcomp <- AAcomp[match(colnames(D_abundance), AAcomp$organism), ]
    nO2 <- calc_metrics(AAcomp, "nO2")[, 1]
    # Reorder genera from most reduced to most oxidized
    onO2 <- order(nO2)
    D_abundance <- D_abundance[, onO2]
    nO2 <- nO2[onO2]

    # Truncate values to the range for the color scale
    islo <- D_abundance < min(breaks)
    D_abundance[islo] <- min(breaks)
    ishi <- D_abundance > max(breaks)
    D_abundance[ishi] <- max(breaks)

    # Plot heatmap
    ## This isn't needed because import("plot.matrix") has been added to NAMESPACE
    #requireNamespace("plot.matrix")
    # Temporarily suppress the x-axis
    opar <- par(xaxt = "n")
    ylab <- paste(disease, "dataset")
    ylab <- gsub("COVID", hyphen.in.pdf("COVID-19"), ylab)
    plot(D_abundance, col = col, breaks = breaks, main = "", xlab = "", ylab = ylab)
    par(opar)

    # Add triangles to indicate values beyond the color scale
    pch <- c(6, 2)
    islh <- list(islo, ishi)
    for(i in 1:2) {
      isi <- islh[[i]]
      if(any(isi)) {
        xy <- which(isi, arr.ind = TRUE)
        x <- xy[, 2]
        y <- nrow(D_abundance) + 1 - xy[, 1]
        points(x, y, pch = pch[i], col = "white")
      }
    }

    # Find oxygen tolerance of genera
    genus <- colnames(D_abundance)
    oxygen.tolerance <- get.oxytol(genus)
    # Add symbols to indicate obligate anaerobes
    ianaerobe <- oxygen.tolerance == "obligate anaerobe"
    labels <- genus
    labels[ianaerobe] <- paste(labels[ianaerobe], "*")
    labels[labels == "Bifidobacterium"] <- "Bifidobacterium +"
    # Make rotated labels (modified from https://www.r-bloggers.com/rotated-axis-labels-in-r-plots/)
    text(x = seq_along(labels), y = par()$usr[3] - strheight("A"), labels = labels, srt = 40, adj = 1, xpd = TRUE)
    # Add tick marks
    axis(1, at = seq_along(labels), labels = FALSE)
    if(disease == "COVID") {
      # Add legend title
      text(ncol(D_abundance) + 1.5, -2, "Change in\nrelative abundance\n(patient - control)", xpd = NA, cex = 1.2)
    }

  }

  # Plot nO2 of genera at top
  par(mar = c(1, 5, 1, 4))
  # Calculate x-axis limits to account for width of boxes and legend in heatmap
  xlim <- c(1 - 0.5, ncol(D_abundance) + 2)
  plot(xlim, range(nO2), xaxs = "i", xaxt = "n", xlab = "", yaxs = "i", ylab = quote(italic(n)[O[2]]~"of genus RP"), type = "n", bty = "n")
  for(i in 1:ncol(D_abundance)) lines(c(i-0.5, i+0.5), rep(nO2[i], 2), lwd = 2)

  if(pdf) dev.off()

  # Return list of genera 20231230
  #writeLines(genus, "Figure_5_genera.txt")
  invisible(genus)

}

# Oxygen tolerance of genera in body sites, COVID-19, and IBD 20230726
microhum_6 <- function(pdf = FALSE) {

  # Start plot
  if(pdf) pdf("Figure_6.pdf", width = 12, height = 9)
  mat <- cbind(matrix(1:12, nrow = 3, byrow = TRUE), c(13, 0, 0))
  layout(mat, widths = c(2, 2, 2, 2, 1))
  par(mar = c(4, 4, 2.8, 1), mgp = c(2.5, 1, 0), cex.lab = 1.2)

  # Panels A and B: oxygen tolerance of genera in body sites / in selected COVID-19 and IBD studies
  for(site in c("Nasal", "Oral", "Skin", "Feces", "COVID_control", "COVID", "IBD_control", "IBD")) {
    dat <- calc.oxytol(site)
    plot.oxytol(dat)
    main <- hyphen.in.pdf(gsub("COVID", "COVID-19", site))
    title(main, font.main = 1, line = 0.5)
    # Add panel title
    atitle <- bquote(bold("A. Genus abundance vs"~bolditalic(n)[O[2]]~.(
      hyphen.in.pdf("for body sites (data from Boix-Amor\u00f3s et al., 2021)"))))
    if(site == "Nasal") label.figure(atitle, font = 2, xfrac = 1.17, yfrac = 0.965, cex = 1.5)
    btitle <- bquote(bold("B. Genus abundance vs"~bolditalic(n)[O[2]]~.(
      hyphen.in.pdf("for COVID-19 (data from Schult et al., 2022) and IBD (data from Lloyd-Price et al., 2019)"))))
    if(site == "COVID_control") label.figure(btitle, font = 2, xfrac = 1.7, yfrac = 0.965, cex = 1.5)
  }

  # Panel C: Percent of aerotolerant genera in COVID-19 or IBD vs controls 20230726
  startplot <- function(ylab, xymax = 100) {
    plot(c(0, xymax), c(0, xymax), type = "n", xlab = "Controls (% aerotolerant)", ylab = ylab)
    lines(c(0, xymax), c(0, xymax), lty = 2, col = 8)
  }
  # Colors and point symbols for sample types
  col <- list(naso = col_Nasal, oro = col_Oral, gut = col_Gut, IBD = col_IBD)
  col <- lapply(col, adjustcolor, alpha.f = 0.8)
  pch <- list(naso = pch_Nasal, oro = pch_Oral, gut = pch_Gut, IBD = pch_IBD)
  # Read precomputed values
  metrics <- read.csv(file.path(getdatadir(), "16S/dataset_metrics_all.csv"))

  # Calculate percentage of aerotolerant genera among those with known oxygen tolerance 20230726
  control <- metrics$control_aerotolerant / (metrics$control_aerotolerant + metrics$control_anaerobe) * 100
  disease <- metrics$disease_aerotolerant / (metrics$disease_aerotolerant + metrics$disease_anaerobe) * 100

  # Loop over sample groups
  for(type in c("naso", "oro", "gut", "IBD")) {
    if(type == "IBD") ylab <- "IBD (% aerotolerant)" else ylab <- hyphen.in.pdf("COVID-19 (% aerotolerant)")
    if(type %in% c("gut", "IBD")) xymax <- 60 else xymax <- 100
    startplot(ylab, xymax)
    itype <- metrics$type == type
    label <- 1:sum(itype)

    # Add points
    points(control[itype], disease[itype], pch = pch[[type]], col = black50, bg = col[[type]])
    dx <- rep(0, sum(itype))
    dy <- 4
    if(type == "gut") dy <- 3
    if(type == "IBD") dy <- 2.2
    dy <- rep(dy, sum(itype))
    if(type == "naso") {
      dy[8] <- -4
    }
    if(type == "oro") {
      dy[9] <- 0
      dx[9] <- -4
    }
    if(type == "gut") {
      dy[c(3, 8)] <- -2.5
      dy[1] <- 0
      dx[1] <- -2
    }
    if(type == "IBD") {
      dy[c(7, 10)] <- -3
      dy[c(12, 15)] <- 0
      dx[c(12, 15)] <- 3
      dx[4] <- 2
      dy[4] <- 1
      dx[c(5, 6)] <- -2
      dx[6] <- -1.5
      dy[5] <- 0
      dy[6] <- -2
      dx[9] <- -1
    }
    # Label points
    text(control[itype] + dx, disease[itype] + dy, label)
    # Add plot title
    titles <- hyphen.in.pdf(c(naso = "Nasopharyngeal (COVID-19)", oro = "Oropharyngeal (COVID-19)", gut = "Gut (COVID-19)", IBD = "Gut (IBD)"))
    title(titles[type], font.main = 1, line = 0.5)

    if(type == "naso") label.figure(
      hyphen.in.pdf("C. Cumulative abundance of aerotolerant genera in COVID-19 and IBD (data sources listed in Table 1)"),
      font = 2, xfrac = 1.505, yfrac = 0.965, cex = 1.5
    )
  }

  # Add legend for colors 20231228
  plot.new()
  cols <- adjustcolor(rev(c(2, 4, 8)), alpha.f = 0.3)
  legend("topright", c("Unassigned", "Aerotolerant", "Obligate anaerobe"), title = "Oxygen\ntolerance",
    pch = 15, col = cols, pt.cex = 2, bty = "n", xpd = NA, cex = 1.2, inset = c(-0.2, 0))


  if(pdf) dev.off()

}

# Amount of putative human DNA removed from HMP metagenomes in screening step 20231222
microhum_1_1 <- function(pdf = FALSE) {
  if(pdf) pdf("Figure_1-1.pdf", width = 10, height = 8)
  # Get sequence processing statistics
  statsfile <- file.path(getdatadir(), "ARAST/HMP12_stats.csv")
  stats <- read.csv(statsfile)
  # Get sample metadata and put in same order as processed sequences
  mdat <- getmdat_microhum("HMP12")
  imdat <- match(stats$Run, mdat$Metagenome)
  mdat <- mdat[imdat, ]
  # Calculate number of sequences removed by screening
  n_removed <- stats$scrubbed_sequences - stats$screened_sequences
  # Calculate number of removed sequences as percentage of input sequences
  perc_removed <- n_removed / stats$input_sequences * 100
  # Setup plot
  par(mfrow = c(3, 2))
  par(cex = 0.8)
  # Loop over body sites
  sites <- c("Skin", "Nasal cavity", "Oral cavity", "GI tract", "UG tract")
  col <- c(col_Skin, col_Nasal, col_Oral, col_Gut, col_UG)
  for(i in seq_along(sites)) {
    hist(perc_removed[mdat$"Body site" == sites[i]], breaks = seq(0, 100, 20), col = col[i],
      xlab = "Percentage of sequences removed by human DNA screening", main = sites[i])
  }
  if(pdf) dev.off()
}

# Differences of nO2 and nH2O between untreated and viral-inactivated samples 20221125
microhum_3_1 <- function(pdf = FALSE) {

  if(pdf) pdf("Figure_3-1.pdf", width = 6, height = 4)
  par(mfrow = c(2, 3))
  par(mgp = c(2.5, 1, 0))
  par(mar = c(4, 4, 2, 1))

  ## Community reference proteomes for body sites and viral inactivation
  ## based on 16S rRNA gene sequences from Boix-Amoros et al. (2021) 20220814

  # Get data for any treatment
  Any <- plotmet_microhum("BPB+21_AnyTreatment", extracolumn = c("Subject", "Site", "Treatment"), plot.it = FALSE)
  # Get data for no treatment
  No <- plotmet_microhum("BPB+21_NoTreatment", extracolumn = c("Subject", "Site"), plot.it = FALSE)

  # Keep samples for subjects (not controls)
  Any <- Any[grep("^S", Any$sample), ]
  No <- No[grep("^S", No$sample), ]

  plottreated <- function(Treatment) {
    # Get chemical metrics for treated samples
    Treated <- Any[Any$Treatment == Treatment, ]
    # Put in same order as untreated samples
    No_subject_site <- paste(No$Subject, No$Site)
    Treated_subject_site <- paste(Treated$Subject, Treated$Site)
    iTreated <- match(No_subject_site, Treated_subject_site)
    Treated <- Treated[iTreated, ]
    # Get colors and symbols
    pch <- sapply(Treated$Site, switch, "Oral cavity" = pch_Oral, "Nasal cavity" = pch_Nasal, "Skin of forearm" = pch_Skin, "feces" = pch_Gut, NA)
    col <- sapply(Treated$Site, switch, "Oral cavity" = col_Oral, "Nasal cavity" = col_Nasal, "Skin of forearm" = col_Skin, "feces" = col_Gut, NA)
    col <- sapply(col, adjustcolor, alpha.f = 0.8)
    # Calculate D_nH2O and D_nO2
    D_nH2O <- Treated$nH2O - No$nH2O
    D_nO2 <- Treated$nO2 - No$nO2
    # Plot D_nH2O and D_nO2
    plot(D_nO2, D_nH2O, xlab = quote(Delta*italic(n)[O[2]]), ylab = quote(Delta*italic(n)[H[2]*O]), xlim = c(-0.03, 0.03), ylim = c(-0.03, 0.015), type = "n")
    abline(h = 0, v = 0, lty = 2, col = 8)
    points(D_nO2, D_nH2O, pch = pch, bg = col, col = NA)
    # Plot p-values 20230204
    plot.p.values(No$nO2, Treated$nO2, No$nH2O, Treated$nH2O, paired = TRUE)
    title(Treatment, font.main = 1)
  }

  plottreated("Ethanol")
  plottreated("Formaldehyde")
  plottreated("Heat")
  # Make legend
  plot.new()
  legend("center", legend = c("Skin", "Nasal", "Oral", "Gut"), pch = c(pch_Skin, pch_Nasal, pch_Oral, pch_Gut),
    pt.bg = c(col_Skin, col_Nasal, col_Oral, col_Gut), col = NA, bty = "n", cex = 1.5)
  plottreated("Psoralen")
  plottreated("Trizol")

  if(pdf) dev.off()

}

# Differences of Zc and nH2O between obligate anaerobic and aerotolerant genera 20221017
# Changed Zc to nO2 20230729
microhum_5_1 <- function(pdf = FALSE) {

  # Setup plot
  if(pdf) pdf("Figure_5-1.pdf", width = 7, height = 4)
  par(mfrow = c(1, 2))
  par(mar = c(3, 4, 1, 1))
  par(mgp = c(2.5, 1, 0))

  # Read table from Million and Raoult, 2018
  dat <- read.csv(system.file("extdata/microhum/MR18_Table_S1_modified.csv", package = "JMDplots"))
  # Drop organisms with unknown oxygen tolerance
  dat <- dat[dat$Obligate.anerobic.prokaryote %in% c(0, 1, 2), ]
  # Find obligate anaerobes and aerotolerant genera
  isAnaerobe <- dat$Obligate.anerobic.prokaryote %in% c(1, 2)
  isAerotolerant <- dat$Obligate.anerobic.prokaryote == 0
  # Make sure no genus has conflicting assignments
  stopifnot(length(intersect(dat$Genus.name[isAnaerobe], dat$Genus.name[isAerotolerant])) == 0)

  # Use GTDB-based reference proteomes
  refdb <- "GTDB_220"
  # Get amino acid compositions for genera
  aa <- taxon_AA[[refdb]]
  aa <- aa[aa$protein == "genus", ]
  ngenera_GTDB_all <- nrow(aa)
  # Remove alphabetical suffixes
  aa$organism <- sapply(strsplit(aa$organism, "_"), "[", 1)
  # Keep genera with known oxygen tolerance
  aa <- aa[aa$organism %in% dat$Genus.name, ]
  ngenera_GTDB_known <- nrow(aa)
  # Make data frames for obligate anaerobes and aerotolerant genera
  aaAnaerobe <- aa[aa$organism %in% dat$Genus.name[isAnaerobe], ]
  aaAerotolerant <- aa[aa$organism %in% dat$Genus.name[isAerotolerant], ]

  # Print numbers of genera
  print(paste("Reference proteomes are available for", sum(dat$Genus.name[isAnaerobe] %in% aaAnaerobe$organism),
    "of", sum(isAnaerobe), "obligately anaerobic genera in the list"))
  print(paste("Reference proteomes are available for", sum(dat$Genus.name[isAerotolerant] %in% aaAerotolerant$organism),
    "of", sum(isAerotolerant), "aerotolerant genera in the list"))
  print(paste("Oxygen tolerance is listed for", ngenera_GTDB_known, "of", ngenera_GTDB_all, "genera in GTDB"))

  # Loop over nO2 and nH2O
  for(metric in c("nO2", "nH2O")) {

    # Calculate nO2 or nH2O
    if(metric == "nO2") {
      # Get values for obligate anaerobes and aerotolerant genera
      values <- list(
        Anaerobe = nO2(aaAnaerobe),
        Aerotolerant = nO2(aaAerotolerant)
      )
    }
    if(metric == "nH2O") {
      values <- list(
        Anaerobe = nH2O(aaAnaerobe),
        Aerotolerant = nH2O(aaAerotolerant)
      )
    }
    boxplot(values, ylab = cplab[[metric]], ylim = c(-0.83, -0.59), col = c("palevioletred", "lightblue"))
    # Calculate p-value
    pval <- wilcox.test(values$Anaerobe, values$Aerotolerant)$p.value
    ptext <- bquote(italic(p) == .(format(signif(pval, 1), scientific = 2)))
    legend("bottomright", legend = ptext, bty = "n", inset = c(0, -0.03))
    # Show mean difference
    md <- format(round(diff(sapply(values, mean, na.rm = TRUE)), 4), scientific = 1)
    if(metric == "nO2") difftxt <- bquote(Delta*italic(n)[O[2]] == .(md))
    if(metric == "nH2O") difftxt <- bquote(Delta*italic(n)[H[2]*O] == .(md))
    legend("topleft", legend = difftxt, bty = "n", inset = c(-0.05, -0.02))

  }

  if(pdf) dev.off()

}

# Differences of nO2 and nH2O between subcommunities of obligate anaerobes and aerotolerant genera in controls and patients 20240212
microhum_6_1 <- function(pdf = FALSE) {

  if(pdf) pdf("Figure_6-1.pdf", width = 10, height = 3.5)
  par(mfrow = c(1, 3))

  # Colors and point symbols for sample types
  col <- list(naso = col_Nasal, oro = col_Oral, gut = col_Gut, IBD = col_IBD)
  col <- lapply(col, adjustcolor, alpha.f = 0.8)
  pch <- list(naso = pch_Nasal, oro = pch_Oral, gut = pch_Gut, IBD = pch_IBD)

  # Read precomputed values
  anaerobe <- read.csv(file.path(getdatadir(), "16S/dataset_metrics_anaerobe.csv"))
  aerotolerant <- read.csv(file.path(getdatadir(), "16S/dataset_metrics_aerotolerant.csv"))

  # nO2 differences
  plot(aerotolerant$nO2_dn - anaerobe$nO2_dn, aerotolerant$nO2_up - anaerobe$nO2_up,
    xlab = "Control", ylab = "Patient", pch = unlist(pch[anaerobe$type]), bg = unlist(col[anaerobe$type]))
  abline(h = 0, v = 0, lty = 2, col = 8)
  legend("bottomright", legend = c("Nasopharyngeal", "Oropharyngeal", "Gut"), pch = c(pch_Nasal, pch_Oral, pch_Gut),
    pt.bg = c(col_Nasal, col_Oral, col_Gut), col = black50, bty = "n", cex = 0.8, title = hyphen.in.pdf("COVID-19"), inset = c(0, 0.12))
  legend("bottomright", legend = "Gut", pch = pch_IBD, pt.bg = col_IBD, col = black50, bty = "n", cex = 0.8, title = "                  IBD", inset = c(0.135, 0))
  title(quote(Delta*italic(n)[O[2]] ~ "(aerotolerant - obligate anaerobe)"), font.main = 1)

  # nH2O differences
  plot(aerotolerant$nH2O_dn - anaerobe$nH2O_dn, aerotolerant$nH2O_up - anaerobe$nH2O_up,
    xlab = "Control", ylab = "Patient", pch = unlist(pch[anaerobe$type]), bg = unlist(col[anaerobe$type]))
  abline(h = 0, v = 0, lty = 2, col = 8)
  title(quote(Delta*italic(n)[H[2]*O] ~ "(aerotolerant - obligate anaerobe)"), font.main = 1)

  # Table of p-values
  plot.new()
  xs <- seq(0.4, 1, length.out = 4)
  text(xs, rep(0.8, 4), c("Control", "Patient", "Control", "Patient"), adj = 1)
  text(c(0.46, 0.88), rep(0.85, 2), c(expression(italic(n)[O[2]]), expression(italic(n)[H[2]*O])), adj = 1)
  text(0.63, 0.92, "P-values")
  ys <- seq(0.7, 0.3, length.out = 4)
  text(rep(-0.2, 4), ys, hyphen.in.pdf(c("Nasopharyngeal", "Oropharyngeal", "Gut (COVID-19)", "Gut (IBD)")), adj = 0, xpd = NA)
  # Loop over metrics
  metrics <- c("nO2", "nH2O")
  for(i in 1:length(metrics)) {
    # Loop over control and patient
    groups <- c("dn", "up")
    for(j in 1:length(groups)) {
      # Name of column with metrics
      colname <- paste(metrics[i], groups[j], sep = "_")
      # Loop over dataset types
      types <- c("naso", "oro", "gut", "IBD")
      for(k in 1:length(types)) {
        itype <- anaerobe$type == types[k]
        # Get metrics for anaerobic and aerotolerant subcommunities
        anaer <- anaerobe[itype, colname]
        aero <- aerotolerant[itype, colname]
        # Calculate p-value
        p.value <- wilcox.test(anaer, aero, paired = TRUE)$p.value
        ptext <- format(signif(p.value, 1), scientific = 2)
        # Use bold text for p < 0.05
        if(p.value < 0.05) font <- 2 else font <- 1
        text(xs[j*2 - 2 + i], ys[k], ptext, adj = 1, font = font)
      }
    }
  }

  if(pdf) dev.off()

}

##################################
### Data Processing Functions  ###
##################################

# Calculate mean values of chemical metrics for patients and controls 20220823
# Include abundances of genera in oxygen tolerance groups 20230725
dataset_metrics <- function() {

  # Define oxygen tolerance groups
  oxytols <- c("all", "anaerobe", "aerotolerant", "unknown")

  # Function to calculate mean values of metrics for patients and controls
  getmeans <- function(study) {

    print(study)
    # Loop over subsets by oxygen tolerance 20240211
    means <- sapply(oxytols, function(oxytol) {

      # Calculate chemical metrics for each run
      metrics <- getmetrics_microhum(study, oxytol)
      # Remove runs with NA metrics (for oxytol != "all")
      metrics <- metrics[!is.na(metrics$Zc), ]
      # Get metadata for the remaining runs
      mdat <- getmdat_microhum(study, metrics)
      # "up" for disease/positive, "down" for control/negative
      is.up <- sapply(mdat$metadata$pch == 25, isTRUE)
      is.dn <- sapply(mdat$metadata$pch == 24, isTRUE)

      # Calculate means of chemical metrics
      Zc_dn <- mean(mdat$metrics$Zc[is.dn])
      Zc_up <- mean(mdat$metrics$Zc[is.up])
      nO2_dn <- mean(mdat$metrics$nO2[is.dn])
      nO2_up <- mean(mdat$metrics$nO2[is.up])
      nH2O_dn <- mean(mdat$metrics$nH2O[is.dn])
      nH2O_up <- mean(mdat$metrics$nH2O[is.up])
      # Calculate p-values 20220905
      Zc_pvalue <- wilcox.test(mdat$metrics$Zc[is.dn], mdat$metrics$Zc[is.up])$p.value
      nO2_pvalue <- wilcox.test(mdat$metrics$nO2[is.dn], mdat$metrics$nO2[is.up])$p.value
      nH2O_pvalue <- wilcox.test(mdat$metrics$nH2O[is.dn], mdat$metrics$nH2O[is.up])$p.value
      # Include number of samples 20220905
      mymeans <- data.frame(n_dn = sum(is.dn), n_up = sum(is.up),
           Zc_dn = Zc_dn, Zc_up = Zc_up, Zc_pvalue = Zc_pvalue, 
           nO2_dn = nO2_dn, nO2_up = nO2_up, nO2_pvalue = nO2_pvalue, 
           nH2O_dn = nH2O_dn, nH2O_up = nH2O_up, nH2O_pvalue = nH2O_pvalue
      )
      if(oxytol == "all") {
        # Sum abundances of genera in oxygen tolerance groups 20230725
        disease <- calc.oxytol(study = study)
        disease_anaerobe <- sum(disease$abundance[disease$oxygen.tolerance == "obligate anaerobe"])
        disease_aerotolerant <- sum(disease$abundance[disease$oxygen.tolerance == "aerotolerant"])
        disease_unknown <- sum(disease$abundance[disease$oxygen.tolerance == "unknown"])
        control <- calc.oxytol("control", study = study)
        control_anaerobe <- sum(control$abundance[control$oxygen.tolerance == "obligate anaerobe"])
        control_aerotolerant <- sum(control$abundance[control$oxygen.tolerance == "aerotolerant"])
        control_unknown <- sum(control$abundance[control$oxygen.tolerance == "unknown"])
        # Include number of samples 20220905
        mymeans <- cbind(mymeans, data.frame(
             control_anaerobe = control_anaerobe, disease_anaerobe = disease_anaerobe, 
             control_aerotolerant = control_aerotolerant, disease_aerotolerant = disease_aerotolerant,
             control_unknown = control_unknown, disease_unknown = disease_unknown
        ))
      }
      mymeans

    })

    names(means) <- oxytols
    means

  }

  out <- list(all = list(), anaerobe = list(), aerotolerant = list(), unknown = list())

  # List COVID-19 and IBD datasets
  microhum_studies <- list(
    # COVID-19 nasopharyngeal
    naso = c("PMM+22", "SGC+21", "HMH+21", "VCV+21", "SRS+22", "CSC+22", "GKJ+22", "MLW+21_Nasopharyngeal"),
    # COVID-19 oral/oropharyngeal
    oro = c("RFH+22_Oral", "IZC+21", "GBS+22", "WCJ+21_Oral", "XLZ+21", "MAC+21", "MLW+21_Oropharyngeal", "GWL+21", "RWC+21_Oral"),
    # COVID-19 gut
    gut = c("ZZZ+21", "RFH+22_Gut", "KMG+21", "WCJ+21_Gut", "CGC+22", "GCW+20", "NGH+21",
            "RDM+22", "MIK+22", "WZL+23", "AHM+21", "FBD+22", "RWC+21_Gut", "SRK+22"),
    # IBD gut
    IBD = c("ZTG+21", "ASM+23", "DKK+23", "AAM+20", "PYL+23", "MLL+16", "LZD+19",
            "BKK+17", "MDV+22", "LAA+19", "REP+23", "HBL+17", "GKD+14", "WGL+19", "RAF+20")
  )

  # Loop over groups of datasets
  for(i in 1:4) {
    # Calculate means for each dataset
    means <- lapply(microhum_studies[[i]], getmeans)
    # Loop over oxygen tolerance groups 20240211
    for(j in 1:4) {
      mymeans <- sapply(means, "[", j)
      mymeans <- do.call(rbind, mymeans)
      mymeans <- cbind(type = names(microhum_studies)[i], study = microhum_studies[[i]], mymeans)
      out[[j]][[i]] <- mymeans
    }
  }

  # Loop over oxygen tolerance groups
  for(j in 1:4) {
    # Put together groups of datasets
    myout <- do.call(rbind, out[[j]])
    # Calculate mean differences of chemical metrics
    D_Zc <- myout$Zc_up - myout$Zc_dn
    D_nO2 <- myout$nO2_up - myout$nO2_dn
    D_nH2O <- myout$nH2O_up - myout$nH2O_dn
    myout <- cbind(myout, D_Zc, D_nO2, D_nH2O)
    # Round values 20230212
    if(j == 1) myout[, 5:22] <- signif(myout[, 5:22], 6)
    if(j > 1) myout[, 5:16] <- signif(myout[, 5:16], 6)
    file <- paste0("16S/dataset_metrics_", oxytols[j], ".csv")
    write.csv(myout, file, row.names = FALSE, quote = FALSE)
  }

}

#############################
### Unexported Functions  ###
#############################

# Function to add p-values to x and y axes 20230204
plot.p.values <- function(X.1, X.2, Y.1, Y.2, paired = FALSE, ypos = "top") {
  p.value.X <- wilcox.test(X.1, X.2, paired = paired)$p.value
  p.value.Y <- wilcox.test(Y.1, Y.2, paired = paired)$p.value
  p.X <- format(signif(p.value.X, 1), scientific = 2)
  p.Y <- format(signif(p.value.Y, 1), scientific = 2)
  p.X.txt <- bquote(italic(p) == .(p.X))
  p.Y.txt <- bquote(italic(p) == .(p.Y))
  # Use bold text for significant p-values 20230723
  if(p.value.X < 0.05) p.X.txt <- bquote(bolditalic(p) == bold(.(p.X)))
  if(p.value.Y < 0.05) p.Y.txt <- bquote(bolditalic(p) == bold(.(p.Y)))
  pu <- par("usr")
  dx <- (pu[2] - pu[1]) / 30
  dy <- (pu[4] - pu[3]) / 30
  text(pu[2] - dx, pu[3] + dy/2, p.X.txt, adj = c(1, 0), cex = 0.9)
  if(ypos == "bottom") text(pu[1] + dx/2, pu[3] + dy, p.Y.txt, srt = 90, adj = c(0, 1), cex = 0.9)
  if(ypos == "top") text(pu[1] + dx/2, pu[4] - dy, p.Y.txt, srt = 90, adj = c(1, 1), cex = 0.9)
}

# Get oxygen tolerance for specified genera 20231227
get.oxytol <- function(genus) {

  # Read the table of oxygen tolerance for genera
  dat <- read.csv(system.file("extdata/microhum/MR18_Table_S1_modified.csv", package = "JMDplots"))
  obligate.anaerobe <- dat$Genus.name[dat$Obligate.anerobic.prokaryote %in% c(1, 2)]
  aerotolerant <- dat$Genus.name[dat$Obligate.anerobic.prokaryote == 0]
  # Remove suffixes used in GTDB (_A, _B, etc.)
  genus <- sapply(strsplit(genus, "_"), "[", 1)
  # Assign oxygen tolerance
  oxygen.tolerance <- rep("unknown", length(genus))
  oxygen.tolerance[genus %in% obligate.anaerobe] <- "obligate anaerobe"
  oxygen.tolerance[genus %in% aerotolerant] <- "aerotolerant"
  oxygen.tolerance

}

# Summarize abundances of genera in segments (body sites, disease, or control) and list their oxygen tolerance 20230725
calc.oxytol <- function(segment = "Feces", study = NULL) {
  
  if(tolower(segment) %in% c("feces", "nasal", "oral", "skin")) {
    # Identify RDP file
    RDPfile <- file.path(getdatadir(), "16S/RDP-GTDB/BPB+21.tab.xz")
    # Read metadata
    mdat <- getmdat_microhum("BPB+21_NoTreatment")
    # Get run IDs for this body segment
    Run <- mdat$Run[grep(tolower(segment), tolower(mdat$Site))]
  }

  if(tolower(segment) %in% c("covid", "covid_control")) {
    RDPfile <- file.path(getdatadir(), "16S/RDP-GTDB/SRK+22.tab.xz")
    mdat <- getmdat_microhum("SRK+22")
    if(tolower(segment) == "covid") Run <- mdat$Run[mdat$pch == 25]
    if(tolower(segment) == "covid_control") Run <- mdat$Run[mdat$pch == 24]
  }

  if(tolower(segment) %in% c("ibd", "ibd_control")) {
    RDPfile <- file.path(getdatadir(), "16S/RDP-GTDB/LAA+19.tab.xz")
    mdat <- getmdat_microhum("LAA+19")
    if(tolower(segment) == "ibd") Run <- mdat$Run[mdat$pch == 25]
    if(tolower(segment) == "ibd_control") Run <- mdat$Run[mdat$pch == 24]
  }

  # Get classifications and metadata for any COVID-19 or IBD study
  if(!is.null(study)) {
    # Remove suffix after underscore 20200929
    studyfile <- gsub("_.*", "", study)
    RDPfile <- file.path(getdatadir(), "16S/RDP-GTDB", paste0(studyfile, ".tab.xz"))
    mdat <- getmdat_microhum(study)
    if(tolower(segment) == "control") Run <- mdat$Run[mdat$pch == 24] else Run <- mdat$Run[mdat$pch == 25]
  }

  # Read RDP file
  RDP <- read_RDP(RDPfile, quiet = TRUE)
  # Keep only genus-level classifications
  RDP <- RDP[RDP$rank == "genus", ]
  # Use genus names for row names
  rownames(RDP) <- RDP$name
  # Keep only selected samples
  Run <- Run[Run %in% colnames(RDP)]
  RDP <- RDP[, Run]
  # Remove genera with zero counts
  RDP <- RDP[rowSums(RDP) > 0, ]
  # Sum abundances across samples and normalize
  abundance <- rowSums(RDP)
  abundance <- abundance / sum(abundance) * 100

  # Get amino acid compositions of reference proteomes for genera and calculate nO2
  AAcomp <- taxon_AA[["GTDB_220"]]
  AAcomp <- AAcomp[match(names(abundance), AAcomp$organism), ]
  nO2 <- calc_metrics(AAcomp, "nO2")[, 1]
  # Get oxygen tolerance
  genus <- rownames(RDP)
  oxygen.tolerance <- get.oxytol(genus)

  # Make summary table
  data.frame(oxygen.tolerance, abundance, nO2)

}

# Line plot: genus abundance vs nO2, grouped by oxygen tolerance 20230725
plot.oxytol <- function(dat) {
  # Start with 0 total abundance
  total.abundance <- 0
  plot(c(-0.81, -0.595), c(0, 100), xlab = chemlab("nO2"), ylab = "Abundance (%)", type = "n")
  # List classifications and colors
  oxytols <- c("obligate anaerobe", "aerotolerant", "unknown")
  cols <- c(2, 4, 8)
  # Loop over classifications
  for(i in 1:3) {
    # Get genera with this oxygen tolerance
    thisdat <- dat[dat$oxygen.tolerance == oxytols[i], ]
    # Order by nO2
    thisdat <- thisdat[order(thisdat$nO2), ]
    # Get y values
    y2 <- total.abundance + cumsum(thisdat$abundance)
    y1 <- head(c(total.abundance, y2), -1)
    # Plot a rectangle to encompass entire range
    rect(head(thisdat$nO2, 1), head(y1, 1), tail(thisdat$nO2, 1), tail(y2, 1), col = adjustcolor(cols[i], alpha.f = 0.3), border = NA)
    # Plot a line for the mean
    nO2.mean <- sum(thisdat$nO2 * thisdat$abundance) / sum(thisdat$abundance)
    lines(c(nO2.mean, nO2.mean), c(head(y1, 1), tail(y2, 1)), lwd = 3, col = "white")
    # Loop over genera
    for(j in 1:nrow(thisdat)) {
      # Plot a line with length corresponding to abundance of this genus
      lines(c(thisdat$nO2[j], thisdat$nO2[j]), c(y1[j], y2[j]), col = cols[i], lwd = 2)
      # Label lines for most-abundant genera
      genus.abundance <- y2[j] - y1[j]
      if(genus.abundance > 3) {
        ymid <- (y1[j] + y2[j]) / 2
        if(thisdat$nO2[j] > -0.7) adj <- 1 else adj <- 0
        text(thisdat$nO2[j], ymid, paste0(" ", rownames(thisdat)[j], " "), adj = adj)
      }
    }
    # Update the total abundance
    total.abundance <- tail(y2, 1)
 }
}

# Get differences of relative abundance of genera in COVID and IBD
# 20231225 first version jmd
get_abundance <- function(study_type, genus_names = NULL) {

  # List studies with datasets
  metrics <- read.csv(file.path(getdatadir(), "16S/dataset_metrics_all.csv"))
  studies <- metrics$study[metrics$type == study_type]

  # Loop over studies
  for(istudy in seq_along(studies)) {

    # Remove suffix after underscore
    studyname <- strsplit(studies[istudy], "_")[[1]][1]
    RDPfile <- file.path(getdatadir(), "16S/RDP-GTDB", paste0(studyname, ".tab.xz"))
    mdat <- getmdat_microhum(studies[istudy])

    # Read RDP file
    RDP <- read_RDP(RDPfile, quiet = TRUE)
    # Keep only genus-level classifications
    RDP <- RDP[RDP$rank == "genus", ]
    # Use genus names for row names
    rownames(RDP) <- RDP$name

    # Loop over control and disease
    for(segment in c("control", "disease")) {
      if(segment == "control") in.segment <- sapply(mdat$pch == 24, isTRUE)
      if(segment == "disease") in.segment <- sapply(mdat$pch == 25, isTRUE)
      myRuns <- mdat$Run[in.segment]
      # Remove runs without taxonomic classifications
      # (possibly due to low-count samples discarded by read_RDP)
      myRuns <- myRuns[myRuns %in% colnames(RDP)]
      # Get genus-level classifications for these runs
      myRDP <- RDP[, myRuns]
      # Sum genus abundances for all samples in this segment
      myabundance <- t(rowSums(myRDP))
      if(istudy == 1) abundance <- myabundance else {
        abundance <- get(segment)
        abundance <- merge(abundance, myabundance, all = TRUE, sort = FALSE)
      }
      assign(segment, abundance)
    }

  }

  # Make sure genus names are the same for control and disease
  stopifnot(all(colnames(control) == colnames(disease)))
  # Set NA counts to 0 and normalize counts to sum to 1
  control[is.na(control)] <- 0
  control <- control / rowSums(control)
  disease[is.na(disease)] <- 0
  disease <- disease / rowSums(disease)
  # Calculate differences of abundance of genera between disease and control
  D_abundance_all <- disease - control

  if(!is.null(genus_names)) {
    return(D_abundance_all[, genus_names])
  } else {
    # Find the genera with the largest change
    cutoff <- 0.05
    is_top <- apply(abs(D_abundance_all) > cutoff, 2, any)
    return(D_abundance_all[, is_top])
  }

}
