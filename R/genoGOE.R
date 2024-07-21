# JMDplots/genoGEO.R
# Make plots for paper:
# Genomes record the Great Oxidation Event
# 20231206 jmd first version
# 20240328 Moved to JMDplots
# 20240409 Add Rubisco plots
# 20240528 Analyze methanogen genomes

# Figure 1: Comparison of methanogen genomes
genoGOE_1 <- function(pdf = FALSE) {

  if(pdf) cairo_pdf("Figure_1.pdf", width = 8, height = 6)
  mat <- matrix(c(1,2,3, 1,2,4, 5,5,5), nrow = 3, byrow = TRUE)
  layout(mat, heights = c(1, 1, 2))
  opar <- par(mgp = c(2.8, 1, 0), mar = c(5.1, 4.1, 2.1, 2.1))

  # Read methanogen genomes information
  mgfile <- system.file("extdata/genoGOE/methanogen_genomes.csv", package = "JMDplots")
  mg <- read.csv(mgfile)
  # Genomes in Halobacteriota
  Halo <- mg$Genome[mg$Methanogen_class == "II"]
  # Genomes in Methanobacteriota
  Methano <- mg$Genome[mg$Methanogen_class == "I"]

  # Panel A: Zc and GC of marker genes 20240528

  # Read data for GTDB marker genes
  markerfile <- system.file("extdata/genoGOE/ar53_msa_marker_info_r220_XHZ+06.csv", package = "JMDplots")
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

  # Plot IQR of Zc
  plot(c(1, 53), c(-0.28, -0.04), xlab = "Marker gene", ylab = quote("Protein"~italic(Z)[C]), type = "n")
  for(i in 1:53) {
    lines(c(i, i) - 0.1, Zc_Methano[i, c(2, 4)], col = 2)
    lines(c(i, i) + 0.1, Zc_Halo[i, c(2, 4)], col = 4)
  }
  # Add legend for methanogen Class I and II
  legend("bottomright", "Class I", lty = 1, col = 2, bty = "n")
  legend("topleft", "Class II", lty = 1, col = 4, bty = "n")
  label.figure("A", font = 2, cex = 1.5)

  # Get GC for species in each phylum
  GC_Halo <- get_GC("Halo")
  GC_Methano <- get_GC("Methano")
  GC_Halo <- GC_Halo[iord, ]
  GC_Methano <- GC_Methano[iord, ]

  # Plot IQR of GC
  plot(c(1, 53), c(0.25, 0.65), xlab = "Marker gene", ylab = "GC content", type = "n")
  for(i in 1:53) {
    lines(c(i, i) - 0.1, GC_Methano[i, c(2, 4)], col = 2)
    lines(c(i, i) + 0.1, GC_Halo[i, c(2, 4)], col = 4)
  }

  # Plot Delta Zc vs Delta GC
  par(mar = c(4.1, 4.1, 1.1, 2.1))
  Delta_Zc <- na.omit(Zc_Halo[, 3] - Zc_Methano[, 3])
  Delta_GC <- na.omit(GC_Halo[, 3] - GC_Methano[, 3])
  plot(Delta_GC, Delta_Zc, xlab = quote(Delta*"GC"),
    ylab = quote(Delta*italic(Z)[C]~"(Class II \u2212 Class I)                                                         "),
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

  # Plot Delta Zc vs log10 protein abundance in M. maripaludis 20240531
  Delta_Zc <- Zc_Halo[, 3] - Zc_Methano[, 3]
  abundance <- markerdat$Redundant.Peptides / markerdat$MW
  log10a <- log10(abundance)
  plot(log10a, Delta_Zc, xlab = quote(log[10]~"protein abundance in"~italic("M. maripaludis")), ylab = "", pch = 19, col = adjustcolor(1, alpha.f = 0.5))
  # Calculate linear fit
  mylm <- lm(Delta_Zc ~ log10a)
  x <- range(log10a)
  y <- predict.lm(mylm, data.frame(log10a = x))
  # Plot linear fit and show R2
  lines(x, y, lty = 2, lwd = 1.5, col = 8)
  R2 <- summary(mylm)$r.squared
  R2_txt <- bquote(italic(R)^2 == .(formatC(R2, digits = 2, format = "f")))
  legend("topleft", legend = R2_txt, bty = "n", inset = c(-0.05, 0))

  par(opar)

  # Panel B: Zc controlled for various factors 20240529
   
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
  # Add beans for GC and Cost
  beanplot(Zc[, 3:14], side = "both", col = list(c(2, 7, 2, 2), c(4, 3, 4, 4)), xlim = c(0.5, 7.5), what = what, names = character(6), add = TRUE, at = 2:7)
  mtext(quote("Protein"~italic(Z)[C]), 2, line = 2.8, cex = par("cex"))

  # Add group names
  axis(1, at = 2:4, labels = c("GC < 0.34", "0.34 \u2264 GC \u2264 0.36", "GC > 0.36"))
  axis(1, at = 5:7, labels = c("Cost < 23", "23 \u2264 Cost \u2264 25", "Cost > 25"))
  axis(3, at = c(1, 3, 6), labels = c("All Proteins", "Control for GC content", "Control for metabolic cost"), tick = FALSE, font = 2)

  label.figure("B", font = 2, cex = 1.5, xfrac = 0.018)

  if(pdf) dev.off()

}

# Carbon oxidation state of proteins as a function of gene age in two lineages
genoGOE_2 <- function(pdf = FALSE, metric = "Zc") {

  if(pdf) pdf("Figure_2.pdf", width = 7, height = 6)
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
    if(j == 1) inset <- c(-0.02, 0.5) else inset <- c(-0.02, 0.28)
    legend("topleft", lineages[j], bty = "n", inset = inset)
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
  }
  # Outer axis labels
  mtext(ylab, side = 2, line = 3, adj = -0.68, font = 2)
  mtext("Gene\nage\n(Ma)", side = 4, line = 1.5, font = 2, las = 1, at = -0.0215, adj = 0.5)

  if(pdf) dev.off()

}

# Chemical analysis and thermodynamic calculations for ancestral Rubiscos
# Can be used to make logaH2O-logfO2, logfO2-pH, or Eh-pH diagram (x = O2 or pH, y = H2O, O2, or Eh)
# 'basis' can be QEC or CHNOS
genoGOE_3 <- function(pdf = FALSE, x = "pH", y = "Eh", basis = "QEC") {

  if(pdf) cairo_pdf("Figure_3.pdf", width = 9, height = 4.5)
  par(mfrow = c(1, 2))

  # Read amino acid compositions
  fasta_file <- system.file("extdata/fasta/KHAB17.fasta", package = "canprot")
  aa <- read_fasta(fasta_file)
  # Assign protein names
  aa$protein <- sapply(strsplit(aa$protein, "_"), "[", 2)

  # Panel A: Zc vs ancestry

  xlab <- "Ancestral sequences (older to younger)"
  par(mar = c(4.0, 4.0, 2.5, 1.0), mgp = c(2.5, 1, 0))
  # Get point locations
  xs <- 1:6
  ys <- Zc(aa)
  # Start plot
  plot(xs, ys, type = "n", xaxt = "n", xlab = xlab, ylab = cplab$Zc)
  # Plot main branch (excluding Anc I/III')
  lines(xs[-3], ys[-3], type = "b", pch = NA)
  # Add points
  points(xs[-3], ys[-3], pch = 19)
  points(xs[3], ys[3], pch = 19, col = 8)
  axis(1, at = 1:6, aa$protein)
  abline(v = 3.5, lty = 2, col = "darkgreen", lwd = 2)
  axis(3, at = 3.5, "GOE (proposed)")
  label.figure("A", cex = 1.5, font = 2)

  # Panel B: Relative stability diagram

  # Add proteins to CHNOSZ
  ip <- add.protein(aa)
  # Set plot resolution
  res <- 500

  # Setup basis species for charged proteins
  if(x == "pH" | y == "Eh") basis <- paste0(basis, "+")
  basis(basis)
  # Setup basis species for Eh diagram
  if(y == "Eh") swap.basis("O2", "e-")
  
  # Get axis limits
  lims <- list(pH = c(0, 14), Eh = c(-0.5, 0.8), O2 = c(-90, 0), H2O = c(-10, 10))
  # Create argument list for affinity()
  aff_args0 <- list(c(lims[[x]], res), c(lims[[y]], res), iprotein = ip)
  names(aff_args0)[1:2] <- c(x, y)

  # Calculate maximum affinity among all proteins
  a0 <- do.call(affinity, aff_args0)
  d0 <- diagram(a0, plot.it = FALSE)

  # The affinity range
  aff_range <- range(d0$predominant.values)
  # The fraction of the balanced range (equal negative and positive endpoints) that is missing from the actual range
  miss_frac <- ( max(abs(aff_range)) - max(aff_range) ) / max(abs(aff_range)) / 2
  # Total number of colors in the diverging scale
  ntot <- 200
  # How many colors to take away to center the scale on zero affinity
  nout <- round(miss_frac * ntot)
  nuse <- ntot - nout
  col <- hcl.colors(ntot, palette = "Blue-Yellow 2")[1:nuse]
  # Start plot with colors for affinity
  thermo.plot.new(lims[[x]], lims[[y]], xlab = axis.label(x), ylab = axis.label(y))
  image(d0$vals[[1]], d0$vals[[2]], d0$predominant.values, add = TRUE, col = col, useRaster = TRUE)

  # Plot stability lines
  # Anc_I is metastable; remove IB and IA/B to see it
  aff_args1 <- aff_args0
  aff_args1$iprotein <- ip[c(2,4)]
  a1 <- do.call(affinity, aff_args1)
  dx <- c(0.5, 2.8)
  dy <- c(0.36, -0.03)
  diagram(a1, lty = 2, lwd = 2, font = 2, col = "darkgreen", col.names = "darkgreen",
    names = names, add = TRUE, limit.water = TRUE, fill.NA = "gray80", dx = dx, dy = dy)
  # Visualize Anc_I - Anc_IA/B boundary 20240603
  aff_args2 <- aff_args0
  aff_args2$iprotein <- ip[c(4,5)]
  a2 <- do.call(affinity, aff_args2)
  dx <- c(3.8, -0.3)
  dy <- c(-0.05, -0.44)
  diagram(a2, lty = 2, lwd = 2, font = 2, col = "darkorange2", col.names = "darkorange2",
    names = names, add = TRUE, limit.water = FALSE, fill.NA = "gray80", dx = dx, dy = dy)

  # Overlay lines for all proteins
  dx <- c(0.2, 2.2, 0, 0, -1.4, 3.0)
  dy <- c(-0.02, -0.10, 0, 0, -0.17, -0.3)
  d <- diagram(a0, font = 2, lwd = 2, add = TRUE, limit.water = TRUE, dx = dx, dy = dy)
  # Add water stability lines
  water.lines(d, lty = 1, col = 8)
  # Add contour line at zero affinity
  contour(d0$vals[[1]], d0$vals[[2]], d0$predominant.values, levels = 0, col = 4, lty = 4, lwd = 2, add = TRUE, drawlabels = FALSE)
  # Replot frame and axis ticks
  box()
  thermo.axis()

  label.figure("B", cex = 1.5, font = 2)

  if(pdf) dev.off()

}
