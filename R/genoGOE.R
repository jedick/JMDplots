# JMDplots/genoGEO.R
# Plots for the paper: Evolutionary oxidation of proteins in Earth's history
# 20240328 First JMDplots commit: eukaryotic gene age groups
# 20240409 Add Rubisco plots
# 20240528 Analyze methanogen genomes
# 20240803 Compare Rubisco proteins and unrelated genomes (stability_comparison)
# 20241224 Add Zc and stability diagram for S-cycling genomes
# 20250325 Add ancestral nitrogenase
# 20250625 Put all ancestral proteins in one plot and add thioredoxin and IPMDH
# 20250626 Put stability diagrams in one plot
# 20250627 Use average affinity instead of average rank of affinity for groupwise relative stability
# 20250903 Add Zc range diagram

# Figure 1: Ranges of carbon oxidation state for organic compounds, amino acids, and proteins
genoGOE_1 <- function(pdf = FALSE) {

  if(pdf) pdf("Figure_1.pdf", width = 11, height = 6)
  par(mar = c(0, 0, 0, 0))

  # Start with an empty plot
  plot(c(0, 10), c(0, 10), type = "n", axes = FALSE, xlab = "", ylab = "")
  # y-limits
  ybottom <- 0.5
  ytop <- 9.5
  # Get gradient colors
  col <- smoothColors("2", 100, "4")
  # Loop over x-left values
  xs <- c(0, 3, 6)
  for(xleft in xs) {
    # Draw a smooth gradient across the bar
    gradient.rect(xleft, ybottom, xleft + 0.5, ytop, col = col, border = NA, gradient = "y")
  }

  # Read protein data
  aa <- read.csv(system.file("RefDB/organisms/UP000000625_83333.csv.xz", package = "JMDplots"))
  Zc <- Zc(aa)

  # Get 90% interval
  q90 <- quantile(Zc, probs = c(0.05, 0.95))
  # Get elemental formulas for reduced and oxidized proteins
  ilow <- which.min(abs(Zc - q90[1]))
  pf_low <- as.chemical.formula(protein.formula(aa[ilow, ]))
  ihigh <- which.min(abs(Zc - q90[2]))
  pf_high <- as.chemical.formula(protein.formula(aa[ihigh, ]))

  # Add labels

  text(xs[1], ybottom, "-4", cex = 2, font = 2, adj = c(0, 1.5))
  text(xs[1], ytop, "+4", cex = 2, font = 2, adj = c(0, -0.5))
  text(xs[1] + 0.6, ybottom, expr.species("CH4"), cex = 2, adj = 0)
  text(xs[1] + 0.6, ytop, expr.species("CO2"), cex = 2, adj = 0)

  text(xs[2], ybottom, "-1", cex = 2, font = 2, adj = c(0, 1.5))
  text(xs[2], ytop, "+1", cex = 2, font = 2, adj = c(0, -0.5))
  text(xs[2] + 0.6, ybottom, expr.species(info(info("leucine"))$formula), cex = 2, adj = 0)
  text(xs[2] + 0.6, ybottom + 0.8, "Leucine", cex = 2, adj = 0)
  text(xs[2] + 0.6, ytop, expr.species(info(info("aspartic acid"))$formula), cex = 2, adj = 0)
  text(xs[2] + 0.6, ytop - 1, "Aspartic\nacid", cex = 2, adj = 0)

  text(xs[3], ybottom, round(q90[1], 3), cex = 2, font = 2, adj = c(0, 1.5))
  text(xs[3], ytop, round(q90[2], 3), cex = 2, font = 2, adj = c(0, -0.5))
  text(xs[3] + 1.1, ybottom, expr.species(pf_low), cex = 2, adj = 0)
  text(xs[3] + 1.1, ytop, expr.species(pf_high), cex = 2, adj = 0)

  # Add lines

  # Figure out y positions for amino acids range
  aa_range <- predict(lm(data.frame(y = c(ybottom, ytop), x = c(-4, 4))), data.frame(x = c(-1, 1)))
  lines(c(xs[1] + 0.53, xs[2] - 0.03), c(aa_range[1], ybottom + 0.02), lwd = 3, col = 2)
  lines(c(xs[1] + 0.53, xs[2] - 0.03), c(aa_range[2], ytop - 0.02), lwd = 3, col = 4)

  # Figure out y positions for protein range
  p_range <- predict(lm(data.frame(y = c(ybottom, ytop), x = c(-1, 1))), data.frame(x = c(q90[1], q90[2])))
  lines(c(xs[2] + 0.53, xs[3] - 0.03), c(p_range[1], ybottom + 0.02), lwd = 3, col = 2)
  lines(c(xs[2] + 0.53, xs[3] - 0.03), c(p_range[2], ytop - 0.02), lwd = 3, col = 4)

  # Add arrows and text for E. coli
  arrows(xs[3] + 0.7, ybottom + 0.02, xs[3] + 0.7, ytop - 0.02, code = 3, length = 0.2, angle = 35, lwd = 3)
  range_txt <- "90% of\nproteins\nin\nare in this\nrange"
  text(xs[3] + 1.1, mean(c(ybottom, ytop)), range_txt, cex = 2, adj = 0)
  E_coli_txt <- "    E. coli"
  text(xs[3] + 1.1, mean(c(ybottom, ytop)), E_coli_txt, cex = 2, adj = 0, font = 3)

  if(pdf) dev.off()
}

# Figure 2: Genome-wide differences of oxidation state between two lineages of methanogens
genoGOE_2 <- function(pdf = FALSE, panel = NULL) {

  if(is.null(panel)) {
    if(pdf) pdf("Figure_2.pdf", width = 8, height = 6)
    mat <- matrix(c(1,2,3, 1,2,4, 5,5,5), nrow = 3, byrow = TRUE)
    layout(mat, heights = c(1, 1, 2))
  }
  panels <- if(is.null(panel)) LETTERS[1:4] else panel
  opar <- par(mgp = c(2.8, 1, 0), mar = c(5.1, 4.1, 2.1, 2.1))

  # Read methanogen genomes information
  mgfile <- system.file("extdata/genoGOE/GTDB/methanogen_genomes.csv", package = "JMDplots")
  mg <- read.csv(mgfile)
  # Genomes in Halobacteriota
  Halo <- mg$Genome[mg$Methanogen_class == "II"]
  # Genomes in Methanobacteriota
  Methano <- mg$Genome[mg$Methanogen_class == "I"]

  # Panels A-B: Zc and GC of marker genes 20240528

  # Read data for GTDB marker genes
  markerfile <- system.file("extdata/genoGOE/GTDB/ar53_msa_marker_info_r220_XHZ+06.csv", package = "JMDplots")
  markerdat <- read.csv(markerfile)
  markerid <- sapply(strsplit(markerdat$Marker.Id, "_"), "[", 2)
  methanogendir <- system.file("extdata/genoGOE/methanogen", package = "JMDplots")

  # NULL values for variables used in subset()
  protein <- gene <- NULL

  get_Zc <- function(phylum = "Halo") {
    files <- file.path(methanogendir, "marker", "faa", paste0(markerid, ".faa.xz"))
    aalist <- lapply(files, function(file) {
      aa <- suppressMessages(read_fasta(file))
      fivenum(Zc(subset(aa, protein %in% get(phylum))))
    })
    do.call(rbind, aalist)
  }

  get_GC <- function(phylum = "Halo") {
    files <- file.path(methanogendir, "marker", "fna", paste0(markerid, ".fna.xz"))
    nalist <- lapply(files, function(file) {
      na <- suppressMessages(read_fasta(file, molecule = "DNA"))
      na <- subset(na, gene %in% get(phylum))
      GC <- (na$G + na$C) / (na$G + na$C + na$A + na$T)
      fivenum(GC)
    })
    do.call(rbind, nalist)
  }

  # Get Zc for species in each phylum
  Zc_Halo <- get_Zc("Halo")
  Zc_Methano <- get_Zc("Methano")

  # Order by median Zc of Methanobacteriota
  iord <- order(Zc_Methano[, 3])
  Zc_Halo <- Zc_Halo[iord, ]
  Zc_Methano <- Zc_Methano[iord, ]

  if("A" %in% panels) {
    # Plot IQR of Zc
    plot(c(1, 53), c(-0.28, -0.04), xlab = "Marker gene", ylab = quote("Protein"~italic(Z)[C]), type = "n")
    for(i in 1:53) {
      lines(c(i, i) - 0.1, Zc_Methano[i, c(2, 4)], col = 2)
      lines(c(i, i) + 0.1, Zc_Halo[i, c(2, 4)], col = 4)
    }
    # Add legend for Class I and II methanogens
    legend("bottomright", "Class I", lty = 1, col = 2, bty = "n")
    legend("topleft", "Class II", lty = 1, col = 4, bty = "n")
    if(is.null(panel)) label.figure("A", font = 2, cex = 1.6)
    # Calculate p-value 20250304
    # Use median value in each group (3rd column) and paired observations
    p <- t.test(Zc_Halo[, 3], Zc_Methano[, 3], paired = TRUE)$p.value
    ptext <- bquote(italic(p) == .(signif(p, 2)))
    text(5, par("usr")[3], ptext, adj = c(0, -0.5))
  }

  # Get GC for species in each phylum
  GC_Halo <- get_GC("Halo")
  GC_Methano <- get_GC("Methano")
  GC_Halo <- GC_Halo[iord, ]
  GC_Methano <- GC_Methano[iord, ]

  if("B" %in% panels) {
    # Plot IQR of GC
    plot(c(1, 53), c(0.25, 0.65), xlab = "Marker gene", ylab = "GC content", type = "n")
    for(i in 1:53) {
      lines(c(i, i) - 0.1, GC_Methano[i, c(2, 4)], col = 2)
      lines(c(i, i) + 0.1, GC_Halo[i, c(2, 4)], col = 4)
    }
    if(is.null(panel)) label.figure("B", font = 2, cex = 1.6)
    # Calculate p-value 20250304
    # Use median value in each group (3rd column) and paired observations
    p <- t.test(GC_Halo[, 3], GC_Methano[, 3], paired = TRUE)$p.value
    ptext <- bquote(italic(p) == .(signif(p, 2)))
    text(5, par("usr")[3], ptext, adj = c(0, -0.5))
  }

  # Panel C: Delta Zc for marker genes

  if("C" %in% panels) {

    # If plotting only this panel, only make the abundance plot
    if(is.null(panel)) {
      # Plot Delta Zc vs Delta GC
      par(mar = c(4.1, 4.1, 1.1, 2.1))
      Delta_Zc <- na.omit(Zc_Halo[, 3] - Zc_Methano[, 3])
      Delta_GC <- na.omit(GC_Halo[, 3] - GC_Methano[, 3])
      plot(Delta_GC, Delta_Zc, xlab = quote(Delta*"GC"),
        ylab = quote(Delta*italic(Z)[C]~"(Class II - Class I)                                                         "),
        pch = 19, col = adjustcolor(1, alpha.f = 0.5), xpd = NA)
      # Calculate linear fit
      mylm <- lm(Delta_Zc ~ Delta_GC)
      x <- range(Delta_GC)
      y <- predict.lm(mylm, data.frame(Delta_GC = x))
      # Plot linear fit and show R2
      lines(x, y, lty = 2, lwd = 1.5, col = 8)
      R2 <- summary(mylm)$r.squared
      R2_txt <- bquote(italic(R)^2 == .(formatC(R2, digits = 2, format = "f")))
      legend("topleft", legend = R2_txt, bty = "n", inset = c(-0.05, 0))
      if(is.null(panel)) label.figure("C", font = 2, cex = 1.6, yfrac = 0.9)
    }

    ylab <- if(is.null(panel)) "" else quote(Delta*italic(Z)[C]~"(Class II - Class I)")
    # Plot Delta Zc vs log10 protein abundance in M. maripaludis 20240531
    Delta_Zc <- Zc_Halo[, 3] - Zc_Methano[, 3]
    abundance <- markerdat$Redundant.Peptides / markerdat$MW
    log10a <- log10(abundance)
    plot(log10a, Delta_Zc, xlab = quote(log[10]~"protein abundance in"~italic("M. maripaludis")), ylab = ylab, pch = 19, col = adjustcolor(1, alpha.f = 0.5))
    # Calculate linear fit
    mylm <- lm(Delta_Zc ~ log10a)
    x <- range(log10a)
    y <- predict.lm(mylm, data.frame(log10a = x))
    # Plot linear fit and show R2
    lines(x, y, lty = 2, lwd = 1.5, col = 8)
    # Add horizontal line at Delta ZC = 0
    abline(h = 0, lty = 3)
    R2 <- summary(mylm)$r.squared
    R2_txt <- bquote(italic(R)^2 == .(formatC(R2, digits = 2, format = "f")))
    legend("topleft", legend = R2_txt, bty = "n", inset = c(-0.05, 0))

  }

  par(opar)

  # Panel D: Zc controlled for various factors 20240529

  if("D" %in% panels) {
   
    # Get values of Zc, GC, and Cost
    genomes <- mg$Genome
    values <- lapply(genomes, function(genome) {
      aa <- read.csv(file.path(methanogendir, "aa", paste0(genome, "_aa.csv.xz")))
      data.frame(
        Zc = Zc(aa),
        GC = aa$abbrv,
        Cost = Cost(aa)
      )
    })
    names(values) <- genomes

    # NULL values for variables used in subset()
    GC <- Cost <- NULL
    # Get mean Zc for segment (phylum x condition)
    get_mean_Zc <- function(phylum = "Halo", condition = "all") {
      genomes <- get(phylum)
      myval <- values[genomes]
      if(condition == "low_GC") myval <- lapply(myval, subset, GC < 0.34)
      if(condition == "mid_GC") myval <- lapply(myval, subset, GC >= 0.34 & GC <= 0.36)
      if(condition == "high_GC") myval <- lapply(myval, subset, GC > 0.36)
      if(condition == "low_Cost") myval <- lapply(myval, subset, Cost < 23)
      if(condition == "mid_Cost") myval <- lapply(myval, subset, Cost >= 23 & Cost <= 24)
      if(condition == "high_Cost") myval <- lapply(myval, subset, Cost > 25)
      print(paste("median", median(sapply(myval, nrow)), "proteins for", phylum, condition))
      sapply(sapply(myval, "[", "Zc"), "mean")
    }

    Zc <- data.frame(
      Methano_all = get_mean_Zc("Methano", "all"),
      Halo_all = get_mean_Zc("Halo", "all"),
      Methano_low_GC = get_mean_Zc("Methano", "low_GC"),
      Halo_low_GC = get_mean_Zc("Halo", "low_GC"),
      Methano_mid_GC = get_mean_Zc("Methano", "mid_GC"),
      Halo_mid_GC = get_mean_Zc("Halo", "mid_GC"),
      Methano_high_GC = get_mean_Zc("Methano", "high_GC"),
      Halo_high_GC = get_mean_Zc("Halo", "high_GC"),
      Methano_low_Cost = get_mean_Zc("Methano", "low_Cost"),
      Halo_low_Cost = get_mean_Zc("Halo", "low_Cost"),
      Methano_mid_Cost = get_mean_Zc("Methano", "mid_Cost"),
      Halo_mid_Cost = get_mean_Zc("Halo", "mid_Cost"),
      Methano_high_Cost = get_mean_Zc("Methano", "high_Cost"),
      Halo_high_Cost = get_mean_Zc("Halo", "high_Cost")
    )

    # Don't plot overall line
    what <- c(FALSE, TRUE, TRUE, TRUE)
    # Start beanplot with all proteins
    par(mar = c(3.1, 4.1, 4.1, 2.1))
    bp <- beanplot(Zc[, 1:2], side = "both", col = list(c(2, 7, 2, 2), c(4, 3, 4, 4)), xlim = c(0.5, 7.5), what = what, names = "")
    # Add means for species in each phylum
    abline(h = bp$stats[1], col = 2, lty = 2)
    abline(h = bp$stats[2], col = 4, lty = 2)
    # Add labels for methanogen classes
    text(0.6, bp$stats[1] + 0.008, "Class I")
    text(1.6, bp$stats[2] + 0.008, "Class II")
    # Add beans for GC and Cost
    beanplot(Zc[, 3:14], side = "both", col = list(c(2, 7, 2, 2), c(4, 3, 4, 4)), xlim = c(0.5, 7.5), what = what, names = character(6), add = TRUE, at = 2:7)
    mtext(quote("Protein"~italic(Z)[C]), 2, line = 2.8, cex = par("cex"))

    # Add group names
    axis(1, at = 2:4, labels = c("GC < 0.34", "0.34 < GC < 0.36", "GC > 0.36"), gap.axis = 0)
    axis(1, at = 5:7, labels = c("Cost < 23", "23 < Cost < 25", "Cost > 25"), gap.axis = 0)
    axis(3, at = c(1, 3, 6), labels = c("Entire genomes", "Control for GC content", "Control for metabolic cost"), tick = FALSE, font = 2)

    if(is.null(panel)) label.figure("D", font = 2, cex = 1.6, xfrac = 0.018)

  }

  if(pdf & is.null(panel)) dev.off()

}

# Figure 3: Carbon oxidation state of proteins in eukaryotic gene age groups
genoGOE_3 <- function(pdf = FALSE, metric = "Zc") {

  if(pdf) pdf("Figure_3.pdf", width = 7, height = 6)
  layout(matrix(1:2), heights = c(1.2, 1.7))

  # Gene ages from Liebeskind et al. (2016)
  datadir <- system.file("extdata/evdevH2O/LMM16", package = "JMDplots")
  modeAges <- read.csv(file.path(datadir, "modeAges_names.csv"))
  # Read list of reference proteomes
  refprot <- read.csv(file.path(datadir, "reference_proteomes.csv"))
  # Read summed amino acid compositions for proteins in each modeAge in each organism 20231218
  aa <- read.csv(file.path(datadir, "modeAges_aa.csv"))

  if(metric == "Zc") {
    ylab <- "Carbon oxidation state of proteins"
    ylim <- c(-0.18, -0.06)
    yOE <- -0.07
    ytick.at <- seq(-0.18, -0.06, 0.04)
    ylabels.at <- c(-0.18, -0.06)
  }
  if(metric == "nO2") {
    ylab <- "Stoichiometric oxidation state of proteins"
    ylim <- c(-0.75, -0.55)
    yOE <- -0.01
    ytick.at <- seq(-0.75, -0.55, 0.05)
    ylabels.at <- c(-0.75, -0.55)
  }
  if(metric == "nH2O") {
    ylab <- "Stoichiometric hydration state of proteins"
    ylim <- c(-0.8, -0.65)
    yOE <- -0.67
    ytick.at <- seq(-0.8, -0.65, 0.05)
    ylabels.at <- c(-0.8, -0.65)
  }
  if(metric == "Cost") {
    ylab <- "Metabolic cost"
    ylim <- c(22, 24)
    yOE <- 11
    ytick.at <- seq(22, 24, 0.5)
    ylabels.at <- c(22, 24, 1)
  }

  # Loop over lineages
  col <- c(adjustcolor("blue1", 0.6), adjustcolor("green4", 0.6))
  coltext <- c("blue1", "green4")
  lineages <- c("Mammalia", "Saccharomyceta")
  for(j in seq_along(lineages)) {
    # Adjust margins and colors
    if(j == 1) {
      par(mar = c(2.5, 4.4, 2, 3))
    }
    if(j == 2) {
      par(mar = c(7.3, 4.4, 2.5, 3))
    }
    # Start plot
    plot(c(1, 7), ylim, xaxt = "n", xlab = "", yaxt = "n", ylab = "", yaxs = "i", type = "n", font.lab = 2)
    axis(2, ytick.at, labels = FALSE)
    axis(2, ylabels.at, tick = FALSE)
    # Loop over organisms in this lineage
    ilineage <- which(modeAges$X8 %in% lineages[j])
    for(i in ilineage) {
      # Get modeAge and amino acid composition for this organism
      OSCODE <- refprot$OSCODE[i]
      myaa <- aa[aa$organism == OSCODE, ]
      modeAge <- myaa$protein
      # Get chemical metric for each modeAge
      vals <- get(metric)(myaa)
      # Remove Euk+Bacteria (non-phylogenetic age category) 20231210
      vals <- vals[-3]
      modeAge <- head(modeAge, -1)
      # Add lines to plot
      lines(modeAge, vals, col = col[j])
    }
    if(j == 1) inset <- c(-0.02, 0.37) else inset <- c(-0.02, 0.22)
    legend("topleft", paste(lineages[j], "genomes", sep = "\n"), bty = "n", inset = inset)
    # Lines for GOE
    lines(c(2.1, 2.5), c(yOE, yOE), lwd = 4, col = 2)
    if(j == 1) {
      # Top axis: GOE and NOE
      mtext("GOE", side = 3, at = 2.3, line = 0.5, font = 2)
      mtext("NOE", side = 3, at = 5, line = 0.5, font = 2)
      # Lines for NOE
      lines(c(4.5, 6.5), c(yOE, yOE), lwd = 4, col = "red3")
      # Divergence times (TimeTree 5)
      Mya <- c(4250, "3500-2600", 1598, 1275, 743, 563, 180)
      # Put Mya labels on bottom of first panel ...
      Mya_axis <- 1
    } else {
      # ... or top of second panel
      Mya_axis <- 3
      Mya <- c(4250, "3500-2600", 1598, 1275, 642, 528, 523)
      # Lines for NOE
      lines(c(4.5, 5.5), c(yOE, yOE), lwd = 4, col = 2)
      # Bottom axis: tick marks and names for shared ancestry
      axis(1, at = 1:7, labels = FALSE)
      labels <- c("Cellular organisms", "Eukaryotes + Archaea", "Eukaryota", "Opisthokonta")
      text(x = 1:4, y = par()$usr[3] - 1.5 * strheight("A"), labels = labels, srt = 45, adj = 1, xpd = TRUE)
      # Bottom axis: names for Mammalia and Saccharomyceta lineages
      dx <- 0.2
      labels_slash <- c("/", "/", "/")
      text(x = (5:7) - dx, y = par()$usr[3] - 1.5 * strheight("A"), labels = labels_slash, srt = 45, adj = 1, xpd = TRUE)
      labels_Mammalia <- c("Eumetazoa  ", "Vertebrata  ", "Mammalia  ")
      text(x = (5:7) - dx, y = par()$usr[3] - 1.5 * strheight("A"), labels = labels_Mammalia, srt = 45, adj = 1, xpd = TRUE, col = coltext[1])
      labels_Saccharomyceta <- c("Dikarya", "Ascomycota", "Saccharomyceta")
      text(x = (5:7) + dx, y = par()$usr[3] - 1.5 * strheight("A"), labels = labels_Saccharomyceta, srt = 45, adj = 1, xpd = TRUE, col = coltext[2])
    }
    # Add labels for divergence times
    at <- seq_along(Mya)
    axis(Mya_axis, at, Mya)
    # Add plot label 20240803
    label.plot(LETTERS[j], font = 2, cex = 1.2, xfrac = 0.04)
  }
  # Outer axis labels
  mtext(ylab, side = 2, line = 3, adj = -0.68, font = 2)
  mtext("Gene\nage\n(Ma)", side = 4, line = 1.5, font = 2, las = 1, at = -0.0215, adj = 0.5)

  if(pdf) dev.off()

}

# Figure 4: Carbon oxidation state of reconstructed ancestral sequences and extant proteins 20250625
genoGOE_4 <- function(pdf = FALSE) {

  # Function to add Zc labels with red/blue colors 20250625
  label_y_axis <- function() {
    # Add tick labels with red/blue colors 20250625
    axis(2, at = seq(-0.18, -0.12, 0.02), col.axis = 4)
    axis(2, at = -0.20)
    axis(2, at = seq(-0.28, -0.22, 0.02), col.axis = 2)
    # Add axis label
    mtext(cplab$Zc, side = 2, line = 3.5, las = 0, cex = par("cex"))
  }

  # Plot Zc of rubisco from Kacar et al. (2017)  20240407
  plot_Rubisco <- function(ylim = c(-0.20, -0.12)) {
    # Read amino acid compositions
    fasta_file <- system.file("extdata/fasta/KHAB17.fasta", package = "canprot")
    aa <- read_fasta(fasta_file)
    # Assign protein names
    aa$protein <- sapply(strsplit(aa$protein, "_"), "[", 2)

    # Get point locations
    xs <- 1:6
    ys <- Zc(aa)
    # Start plot
    xlab <- "Ancestral sequences (older to younger)"
    plot(xs, ys, type = "n", xaxt = "n", xlab = xlab, ylab = "", ylim = ylim, yaxt = "n")
    label_y_axis()
    # Plot main branch (excluding Anc I/III')
    lines(xs[-3], ys[-3], type = "b", pch = 19)
    # Add point for Anc I/III'
    points(xs[3], ys[3], pch = 19, col = 8)
    axis(1, at = 1:6, aa$protein)
    abline(v = 3.5, lty = 2, lwd = 2)
    text(2.9, -0.13, "GOE\n(estimated)")
  }

  # Plot Zc of ancestral and extant nitrogenases from Garcia et al. (2020)  20250325
  plot_nitrogenase <- function(ylim = c(-0.20, -0.12)) {

    ## Read FASTA file of ancient and extant sequences,
    ## downloaded from https://github.com/kacarlab/AncientNitrogenase.git
    #aa <- canprot::read_fasta("GMKK20/Extant-MLAnc_Align.fasta")
    # Get amino acid sequences precomputed from Extant-MLAnc_Align.fasta
    aa <- read.csv(system.file("extdata/genoGOE/GMKK20/nitrogenase_aa.csv", package = "JMDplots"))

    # List forms of nitrogenase and their ancestors
    form_to_anc <- list(
      "Clfx" = "D",
      "F-Mc" = "C",
      "Mb-Mc" = "B",
      "Anf" = "A",
      "Vnf" = "A",
      "Nif-II" = "E",
      "Nif-I" = "E"
    )

    # Start plot
    plot(extendrange(c(1, 7)), c(-0.20, -0.12), xlab = "Form of nitrogenase", ylab = "", type = "n", axes = FALSE, ylim = ylim)
    axis(side = 1, at = seq_along(form_to_anc), labels = hyphen.in.pdf(names(form_to_anc)), gap.axis = 0)
    label_y_axis()
    box()

    # Loop over nitrogenase forms
    set.seed(42)
    for(iform in seq_along(form_to_anc)) {

      # Calculate Zc of the ancestral proteins
      node <- form_to_anc[[iform]]
      ianc <- grepl(paste0("^Anc", node), aa$protein)
      Zc_anc <- Zc(aa[ianc, ])
      # Plot points with jitter
      xvals <- jitter(rep(iform, length(Zc_anc)), amount = 0.1)
      points(xvals, Zc_anc, pch = 19)

      # Calculate Zc of the extant proteins
      form <- names(form_to_anc[iform])
      iext <- which(aa$ref == form)
      Zc_ext <- Zc(aa[iext, ])
      xvals <- jitter(rep(iform, length(Zc_ext)), amount = 0.1)
      points(xvals, Zc_ext)

    }

    # Add legend and title
    legend("bottomright", c("Extant", "Ancestral"), pch = c(1, 19), bty = "n")

  }

  # Plot Zc of IPMDH from Cui et al. (2025)  20250407
  plot_IPMDH <- function(ylim = c(-0.20, -0.12)) {
    aa <- read_fasta(system.file("extdata/genoGOE/CDY+25/IPMDH.fasta", package = "JMDplots"))
    Zc <- Zc(aa)
    # Ages from Table 1 of Cui et al., 2025
    ages <- c(
      2980, 2960, 2910, 2590,
      2360, 2160, 2140, 1570,
      1200, 932, 624, 0
    )
    pch <- rep(19, length(Zc))
    pch[length(pch)] <- 1
    plot(ages / 1000, Zc, xlim = c(3, 0), xlab = "Age (Ga)", ylab = "", type = "b", ylim = ylim, yaxt = "n", pch = pch)
    label_y_axis()
  }

  # Plot Zc for ancestral thioredoxins from Perez-Jimenez et al. (2011)  20250625
  plot_thioredoxin <- function(ylim = c(-0.28, -0.18)) {
    # Read data file with ages and PDB IDs from Del Galdo et al. (2019)
    dat <- read.csv(system.file("extdata/genoGOE/PIZ+11/DAAD19.csv", package = "JMDplots"))
    # Read amino acid compositions
    aa <- read_fasta(system.file("extdata/genoGOE/PIZ+11/thioredoxin.fasta", package = "JMDplots"))
    # Calculate Zc
    Zc <- Zc(aa)
    # Setup plot
    xlim <- c(4.2, 0)
    plot(xlim, range(Zc), xlim = xlim, xlab = "Age (Ga)", ylab = "", type = "n", ylim = ylim, yaxt = "n")
    label_y_axis()
    # Add separate lines for each lineage
    iBac <- dat$Lineage == "Bacteria"
    pch <- rep(19, sum(iBac))
    pch[length(pch)] <- 1
    lines(dat$Age[iBac], Zc[iBac], type = "b", pch = pch)
    text(3, -0.22, "Bacteria")
    iArcEuk <- dat$Lineage == "Arc-Euk"
    pch <- rep(15, sum(iArcEuk))
    pch[length(pch)] <- 0
    lines(dat$Age[iArcEuk], Zc[iArcEuk], type = "b", pch = pch)
    text(3, -0.265, "Archaea+Eukaryota")
  }

  if(pdf) pdf("Figure_4.pdf", width = 9, height = 6)
  par(mfrow = c(2, 2))
  par(las = 1)
  par(mar = c(4.0, 5.0, 2.5, 1.0), mgp = c(2.5, 1, 0))
  plot_thioredoxin()
  title("Thioredoxin")
  label.figure("A", font = 2, cex = 1.5)
  plot_IPMDH()
  title("IPMDH")
  label.figure("B", font = 2, cex = 1.5)
  plot_Rubisco()
  title("Rubisco")
  label.figure("C", font = 2, cex = 1.5)
  plot_nitrogenase()
  title("Nitrogenase")
  label.figure("D", font = 2, cex = 1.5)
  if(pdf) dev.off()

}

# Figure 5: From carbon oxidation state to relative stability diagrams
genoGOE_5 <- function(pdf = FALSE, panel = NULL) {

  if(is.null(panel)) {
    if(pdf) pdf("Figure_5.pdf", width = 12, height = 8)
    layout(matrix(c(1,1,1,1,1, 2,2,2,2,2, 3,3,3,3,3, 4,4,4,
                    rep(5, 9), rep(6, 9)), nrow = 2, byrow = TRUE))
    par(cex = 1)
  }
  panels <- if(is.null(panel)) LETTERS[1:7] else panel

  # Read amino acid compositions
  fasta_file <- system.file("extdata/fasta/KHAB17.fasta", package = "canprot")
  aa <- read_fasta(fasta_file)
  # Assign protein names
  aa$protein <- sapply(strsplit(aa$protein, "_"), "[", 2)

  # Add proteins to CHNOSZ
  ip <- add.protein(aa, as.residue = TRUE)
  # Setup basis species and swap O2 for e- to make Eh-pH diagram
  basis("QEC+")
  swap.basis("O2", "e-")
  # Set plot resolution
  res <- 300
  
  # Panel A: Pairwise stability boundaries for Rubisco

  if("A" %in% panels) {

    # Loop over individual pairs
    for(pre in 1:3) {
      for(post in 4:6) {
        add <- TRUE
        if(pre == 1 & post == 4) add <- FALSE
        a <- affinity(pH = c(0, 14, res), Eh = c(-0.5, 0.8, res), iprotein = ip[c(pre, post)])
        d <- diagram(a, names = "", lty = 2, col = "#000000b0", add = add, limit.water = !add, fill.NA = "gray80", xlab = "pH", ylab = axis.label("Eh"))
        # Only label lines for reaction with Anc. I
        if(post == 4) {
          # Sort x values and get x and y values of boundary line
          order <- order(d$linesout[[1]])
          xs <- d$linesout[[1]][order]
          ys <- d$linesout[[2]][order]
          # Get a single value along the length of the line
          ilab <- floor(length(xs) * pre * 3 / 10)
          x <- xs[ilab]
          y <- ys[ilab]
          text(x, y - 0.03, aa$protein[pre], cex = 0.6)
          text(x, y + 0.03, aa$protein[post], cex = 0.6)
        }
      }
    }

    text(5.5, 0.67, hyphen.in.pdf("Higher affinity\nfor post-GOE protein\nin each pair"), cex = 0.8)
    text(4.8, -0.15, hyphen.in.pdf("Higher affinity for\npre-GOE protein in each pair"), cex = 0.8, srt = -37)
    title("Pairwise Rubiscos", font.main = 1)
    if(is.null(panel)) label.figure("A", cex = 1.5, font = 2, yfrac = 0.936)

  }

  if("B" %in% panels) {

    # Panel B: Groupwise stability boundaries for Rubisco

    # Calculate affinity of composition reactions for all proteins
    aout <- affinity(pH = c(0, 14, res), Eh = c(-0.5, 0.8, res), iprotein = ip)
    # Set up groups for affinity ranking:
    # 3 pre-GOE and 3 post-GOE proteins
    groups <- list(pre = c(TRUE, TRUE, TRUE, FALSE, FALSE, FALSE), post = c(FALSE, FALSE, FALSE, TRUE, TRUE, TRUE))
    amean <- agg.affinity(aout, groups = groups)
    diagram(amean, lwd = 2, col = 4, names = "", xlab = "pH", ylab = axis.label("Eh"), balance = 1)
    text(6, -0.17, hyphen.in.pdf("Higher mean affinity\nfor pre-GOE proteins"), col = 4, font = 2, cex = 0.8, srt = -33)
    text(6.5, 0.1, hyphen.in.pdf("Higher mean affinity\nfor post-GOE proteins"), col = 4, font = 2, cex = 0.8, srt = -33)
    title("Groupwise Rubiscos", font.main = 1)
    if(is.null(panel)) label.figure("B", cex = 1.5, font = 2, yfrac = 0.936)

  }

  if("C" %in% panels) {

    # Panel D: Comparison between Rubiscos and methanogen and Nitrososphaeria genomes
    stability_comparison(res = res)

    # Label lines
    text(6.8, -61.4, "Rubiscos", cex = 0.8, adj = 0)
    text(6.5, -61.8, hyphen.in.pdf("Pre-GOE"), cex = 0.75, adj = 1, srt = 30)
    text(6.5, -60.7, hyphen.in.pdf("Post-GOE"), cex = 0.75, adj = 1, srt = 30)

    text(6.8, -67, "Methanogens", cex = 0.8, adj = 0)
    text(6, -67.6, "Class I", cex = 0.75)
    text(6, -66.7, "Class II", cex = 0.75)

    text(6.8, -69.4, "Nitrososphaeria", cex = 0.8, adj = 0, font = 3)
    text(6, -70.2, "Basal", cex = 0.75)
    text(6, -69.4, "Terrestrial", cex = 0.75)

    if(is.null(panel)) {
      title("   Rubiscos and all proteins in genomes", font.main = 1, xpd = NA)
      label.figure("C", cex = 1.5, font = 2, yfrac = 0.936)
    }

  }

  if("D" %in% panels) {

    # Add more arrows 20240812
    plot.new()
    opar <- par(xpd = NA)
    arrows(-0.3, 0.08, -0.3, 0.28, length = 0.2, lwd = 2, col = 2)
    arrows(-0.1, 0.22, -0.1, 0.42, length = 0.2, lwd = 2, col = 7)
    arrows(-0.3, 0.65, -0.3, 0.85, length = 0.2, lwd = 2, col = 4)
    text(0.05, 0.25, "Oxidation in\nmany lineages\naround GOE", adj =0)
    text(-0.15, 0.75, hyphen.in.pdf("Rubisco transitions\nat more oxidizing\nconditions"), adj = 0)
    par(opar)

  }

  # Evolutionary oxidation and relative stabilities for genomes with S-cycling genes 20241211
  # genoGOE/sulfur_genomes.R
  # 20241211 add Eh-pH affinity ranking
  # 20241223 convert to logfO2-pH

  # List genomes with single sulfur-cycling genes
  genomes <- list(
    dsrAB = c("GCA_002782605.1", "GCF_000517565.1", "GCA_002878135.1"),
    soxC = c("GCF_000153205.1", "GCF_000024725.1", "GCA_001914955.1", "GCA_002731275.1", 
      "GCA_002007425.1", "GCA_001780165.1", "GCA_002721445.1", "GCF_900129635.1", 
      "GCA_002162915.1", "GCA_002712885.1", "GCF_000484535.1", "GCF_000969705.1", 
      "GCF_002148795.1", "GCF_002514725.1", "GCF_900106035.1", "GCA_002687025.1", 
      "GCA_003222815.1", "GCF_900187885.1", "GCA_002712165.1", "GCA_002705185.1", 
      "GCA_003228115.1"),
    soxABXYZ = c("GCF_000021565.1", "GCF_900142435.1", "GCA_000830255.1", "GCF_000227215.1"),
    aprAB = c("GCA_001800245.1", "GCA_002898195.1", "GCA_002717185.1", "GCF_000328625.1", 
      "GCF_002252565.1", "GCA_001805205.1", "GCA_001784555.1", "GCA_001443375.1"),
    dmsA = c("GCF_000384115.1", "GCA_001593855.1", "GCF_000020005.1", "GCA_003242675.1", 
      "GCA_001771285.1", "GCA_002717245.1", "GCA_003223635.1", "GCF_001051235.1", 
      "GCA_001304035.1", "GCF_000772535.1", "GCA_002898895.1", "GCF_001049895.1", 
      "GCF_001860525.1", "GCA_001775395.1", "GCA_001775995.1", "GCA_001830835.1", 
      "GCF_000487995.1", "GCA_002747435.1", "GCA_002839495.1", "GCA_001515205.2", 
      "GCA_001742785.1", "GCA_001775755.1", "GCA_001768675.1"),
    mddA = c("GCA_003223145.1", "GCA_001872725.1", "GCF_000018105.1", "GCA_001447805.1", 
      "GCA_002400775.1", "GCA_001563325.1", "GCA_002746235.1", "GCA_001780825.1", 
      "GCF_000970205.1", "GCA_002746185.1", "GCA_002790835.1", "GCF_001886815.1", 
      "GCF_000192575.1", "GCA_002256595.1", "GCF_900111015.1", "GCA_001664505.1", 
      "GCA_002841995.1", "GCA_002699105.1", "GCF_002563855.1", "GCA_003219195.1"),
    dmdA = c("GCA_002707655.1", "GCF_900102465.1", "GCA_001800745.1", "GCF_001029505.1", 
      "GCA_002717565.1", "GCA_002701885.1", "GCA_002722565.1")
  )

  # Get amino acid compositions and Zc for genomes
  aa <- read.csv(system.file("extdata/genoGOE/MCK+23/genome_aa.csv", package = "JMDplots"))
  Zcvals <- Zc(aa)
  # List Zc for each genome in list
  Zclist <- lapply(genomes, function(genome) Zcvals[aa$organism %in% genome])
  # Use colors from Mateos et al., 2023
  dsr <- "#9c92ae" # "#bcb2ce"
  sox <- "#45b78d"
  mdd <- "#c24a96"
  # Colors for protein groups
  col <- c(dsr, sox, sox, dsr, mdd, mdd, mdd)

  # Plot Zc of genomes with S-cycling genes from Mateos et al. (2023)  20240916
  sulfur_Zc <- function() {
    # Setup plot
    # Need to make some adjustments after plotting with CHNOSZ::diagram()
    opar <- par(mar = c(4.1, 4.1, 2.1, 2.1), mgp = c(3.1, 1, 0), tcl = -0.5, yaxs = "r")
    n <- length(Zclist)
    boxplot(Zclist, col = col, names = character(n), xlab = "Age of earliest gene event (Ga)", ylab = "Zc of all proteins in genome", ylim = c(-0.25, -0.08))
    text(2.1, -0.24, hyphen.in.pdf("Sulfate-sulfite-sulfide"), col = dsr)
    text(2.8, -0.115, hyphen.in.pdf("Sulfate-\nthiosulfate"), col = sox)
    text(6.2, -0.22, "Organic sulfur", col = mdd)
    # Ages from Table 2 of Mateos et al.
    ages <- c("3.3-3.35  ", "  2.65-2.88", "2.6", "2.33-2.47", "2.28", "1.77", "0-2.38")
    # Make rotated labels (modified from https://www.r-bloggers.com/rotated-axis-labels-in-r-plots/)
    text(x = 1:n, y = par()$usr[3] - 1.5 * strheight("A"), labels = ages, srt = 30, adj = 1, xpd = TRUE)
    genelab <- names(Zclist)
    genelab[genelab == "soxABXYZ"] <- "ABXYZ"
    axis(3, at = 1:n, labels = genelab, line = -0.8, lwd = 0, font = 3, gap.axis = -1)
    axis(3, at = 3, labels = "sox", line = 0.1, lwd = 0, font = 3)
    # Add number of genomes to labels
    n_genomes <- paste0("(", sapply(Zclist, length), ")")
    axis(3, at = 1:n, labels = n_genomes, line = -2, lwd = 0)
    #title(hyphen.in.pdf("Sulfur-cycling gene or gene cluster"), font.main = 1, line = 3)
    par(opar)
  }

  # Affinity ranking for genomes with different S-cycling genes
  sulfur_affinity <- function(panel) {
    if(is.null(panel)) par(mar = c(4.1, 4.1, 2.1, 4.1), mgp = c(2.5, 1, 0))
    basis("QEC+")
    # Keep genomes with single sulfur-cycling genes listed above
    myaa <- aa[aa$organism %in% unlist(genomes), ]
    # Load proteins for each genome
    ip <- add.protein(myaa, as.residue = TRUE)
    # Calculate affinity of composition reactions as a function of Eh and pH
    a <- affinity(pH = c(3, 10, res), O2 = c(-72.5, -58, res), iprotein = ip)
    # Group genomes according to presence of sulfur-cycling genes
    groups <- sapply(genomes, function(genome) match(genome, myaa$organism))
    # Calculate normalized sum of ranks for each group and make diagram
    amean <- agg.affinity(a, groups)
    # Lighten colors
    fill <- adjustcolor(col, alpha.f = 0.3)
    # Adjust labels
    names <- names(genomes)
    dx <- dy <- rep(0, length(names))
    dx[names == "soxABXYZ"] <- 1
    dy[names == "soxABXYZ"] <- -0.2
    # We need balance = 1 here to balance on residues 20250627
    diagram(amean, fill = fill, lty = 1, lwd = 2, font = 3, names = names, cex.names = 0.8, dx = dx, dy = dy, col = "gray20", balance = 1)
  }

  if("D" %in% panels) {
    sulfur_Zc()
    if(is.null(panel)) label.figure("D", font = 2, cex = 1.6)
  }
  if("E" %in% panels) {
    sulfur_affinity(panel)
    if(is.null(panel)) {
      ## Overlay stability boundaries for other genomes
      #stability_comparison(res = res, add = TRUE, pHlim = c(3, 10), alpha.f = 0.8, Eh7_las = 0)
      title(main = hyphen.in.pdf("All proteins in genomes with specific S-cycling genes"), font.main = 1)
      label.figure("E", font = 2, cex = 1.6)
    } else {
      # Only add Eh7 axis
      stability_comparison(res = res, add = TRUE, Eh7_las = 0, datasets = numeric())
    }
  }

  if(pdf & is.null(panel)) dev.off()

}

# Comparison of Rubiscos and methanogen and Nitrososphaeria genomes
stability_comparison <- function(res = 400, add = FALSE, lwd = 2, lty = 1, pHlim = c(4, 10), O2lim = c(-72.5, -58), alpha.f = 1, Eh7_las = 1, datasets = 1:3) {

  # Setup basis species
  basis("QEC+")

  for(i in datasets) {

    if(i == 1) {
      # Methanogen proteomes 20220424
      # Read amino acid composition and compute Zc
      aa <- read.csv(system.file("extdata/utogig/methanogen_AA.csv", package = "JMDplots"))
      # Indices of Class I and Class II methanogens
      iI <- 20:36
      iII <- 1:19
      # Get the species in each group
      groups <- list("Class I" = iI, "Class II" = iII)
      col <- 7
      add <- add
    }

    if(i == 2) {
      # Thaumarchaeota (now Nitrososphaeria) 20220414
      # Amino acid compositions of predicted (Glimmer) and database (NCBI or IMG) proteomes
      predicted <- read.csv(system.file("extdata/utogig/Thaumarchaeota_predicted_AA.csv", package = "JMDplots"))
      database <- read.csv(system.file("extdata/utogig/Thaumarchaeota_database_AA.csv", package = "JMDplots"))
      # If both are available, use predicted instead of database
      aa <- rbind(predicted, database)
      aa <- aa[!duplicated(aa$organism), ]
      groupnames <- c("Basal", "Terrestrial", "Shallow", "Deep")
      # Get the species in each group
      groups <- sapply(groupnames, function(group) aa$protein == group, simplify = FALSE)
      # Compare Basal to Terrestrial 20240802
      groups <- groups[1:2]
      col <- 2
      add <- TRUE
    }

    if(i == 3) {
      # Rubisco 20240802
      # Read amino acid compositions
      fasta_file <- system.file("extdata/fasta/KHAB17.fasta", package = "canprot")
      aa <- read_fasta(fasta_file)
      # Assign protein names
      aa$protein <- sapply(strsplit(aa$protein, "_"), "[", 2)
      # Set up groups for affinity ranking:
      # 3 pre-GOE and 3 post-GOE proteins
      groups = list(pre = c(TRUE, TRUE, TRUE, FALSE, FALSE, FALSE), post = c(FALSE, FALSE, FALSE, TRUE, TRUE, TRUE))
      col <- 4
      add <- TRUE
    }

    # Make affinity ranking plot 20220602
    # Load proteins and calculate affinity
    ip <- add.protein(aa, as.residue = TRUE)
    aout <- affinity(pH = c(pHlim, res), O2 = c(O2lim, res), iprotein = ip)
    # Calculate average ranking for each group and make diagram
    amean <- agg.affinity(aout, groups)
    diagram(amean, col = adjustcolor(col, alpha.f = alpha.f), lwd = lwd, lty = lty, add = add, names = "", balance = 1)

  }

  # Add Eh7 axis 20241218
  # Calculate range of pe at pH = 7
  # H2O = 0.5 O2(gas) + 2 H+ + 2 e- 
  # logK = -41.55
  # --> pe7 = 0.25 * logfO2 + 13.775
  # --> logfO2 = 4 * pe7 - 55.1
  Eh7_to_logfO2 <- function(Eh7) {
    pe7 <- convert(Eh7, "pe")
    logfO2 <- 4 * pe7 - 55.1
    logfO2
  }
  Eh7ticks <- seq(-0.25, 0.05, 0.05)
  logfO2ticks <- Eh7_to_logfO2(Eh7ticks)
  axis(4, at = logfO2ticks, labels = Eh7ticks, tcl = -0.3, mgp = c(2, 0.5, 0))
  mtext("Eh7 (V)", 4, line = 3, las = Eh7_las)

}

# Calculate average affinities for species in different groups
# 20220416 jmd first version (rank.affinity)
# 20250626 use aggregate function (mean or median) instead of average rank
agg.affinity <- function(aout, groups, fun = "mean") {

  # Put the affinities into matrix form
  amat <- sapply(aout$values, as.numeric)
  # Keep track of empty groups
  is_empty_group <- logical()

  # Get the average affinity for species in each group
  agg_values <- sapply(groups, function(group) {

    # Get number of species in this group
    if(inherits(group, "logical")) n <- sum(group)
    if(inherits(group, "integer")) n <- length(group)
    # Also handle indices classed as numeric 20250522
    if(inherits(group, "numeric")) n <- length(group)
    # Aggregate affinities
    group_values <- apply(amat[, group, drop = FALSE], 1, fun)

    # Remember empty group 20250527
    if(n == 0) {
      is_empty_group <<- c(is_empty_group, TRUE)
    } else {
      is_empty_group <<- c(is_empty_group, FALSE)
    }
    group_values

  })

  # Remove empty groups 20250527
  if(any(is_empty_group)) {
    agg_values <- agg_values[, !is_empty_group, drop = FALSE]
    empty_groups <- names(groups)[is_empty_group]
    message(paste("aggregate.affinity: removing empty groups:", paste(empty_groups, collapse = ", ")))
    groups <- groups[!is_empty_group, drop = FALSE]
  }

  # Restore dims
  dims <- dim(aout$values[[1]])
  if(getRversion() < "4.1.0") {
    # Using 'simplify = FALSE' in R < 4.1.0 caused error: 3 arguments passed to 'dim<-' which requires 2
    alist <- lapply(lapply(apply(agg_values, 2, list), "[[", 1), "dim<-", dims)
  } else {
    # apply() got 'simplify' argument in R 4.1.0 20230313
    alist <- apply(agg_values, 2, "dim<-", dims, simplify = FALSE)
  }
  aout$values <- alist

  # Rename species to group names (for use by diagram())
  aout$species <- aout$species[1:length(groups), ]
  aout$species$name <- names(groups)
  aout

}
