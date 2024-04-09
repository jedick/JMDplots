# JMDplots/genoGEO.R
# Plot for the manuscript by Dick (2024) about the genomic record of Earth's surface oxygenation
# 20231206 jmd first version
# 20240328 moved to JMDplots
# 20240409 add genoGOE_2()

# Carbon oxidation state of proteins as a function of gene age in two lineages
genoGOE_1 <- function(pdf = FALSE, metric = "Zc") {

  if(pdf) pdf("Figure_1.pdf", width = 7, height = 6)
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
genoGOE_2 <- function(pdf = FALSE, x = "pH", y = "Eh", basis = "QEC") {

  if(pdf) pdf("Figure_2.pdf", width = 10, height = 5)
  par(mfrow = c(1, 2))

  # Read amino acid compositions
  fasta_file <- system.file("extdata/fasta/KHAB17.fasta", package = "canprot")
  aa <- read_fasta(fasta_file)
  # Assign protein names
  aa$protein <- sapply(strsplit(aa$protein, "_"), "[", 2)

  # Panel A: Zc vs ancestry

  xlab <- "Ancestral sequences (older to younger)"
  plot(Zc(aa), type = "b", xaxt = "n", xlab = xlab, ylab = cplab$Zc)
  axis(1, at = 1:6, aa$protein)
  abline(v = 3.5, lty = 2, col = 6, lwd = 2)
  axis(3, at = 3.5, "GOE (proposed)")
  label.figure("A", cex = 1.5)

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

  # Calculate maximum affinity among all protein
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
  col <- hcl.colors(ntot, palette = "Blue-Red 3")[1:nuse]
  # Start plot with colors for affinity
  thermo.plot.new(lims[[x]], lims[[y]], xlab = axis.label(x), ylab = axis.label(y))
  image(d0$vals[[1]], d0$vals[[2]], d0$predominant.values, add = TRUE, col = col, useRaster = TRUE)

  # Plot stability lines
  # Anc_I is metastable; remove IB and IA/B to see it
  aff_args1 <- aff_args0
  aff_args1$iprotein <- ip[1:4]
  a1 <- do.call(affinity, aff_args1)
  #diagram(a1, lty = 2, lwd = 2, col = 6, col.names = 6, names = names, add = TRUE)
  diagram(a1, lty = 2, lwd = 2, col = 6, col.names = 6, names = names, add = TRUE, limit.water = TRUE, fill.NA = "gray80")
  # Overlay lines for all proteins
  d <- diagram(a0, font = 2, lwd = 3, add = TRUE, limit.water = TRUE)
  # Add water stability lines
  water.lines(d, lty = 1, col = 8)
  # Add contour line at zero affinity
  contour(d0$vals[[1]], d0$vals[[2]], d0$predominant.values, levels = 0, col = 4, lty = 4, lwd = 2, add = TRUE, drawlabels = FALSE)
  # Replot frame and axis ticks
  box()
  thermo.axis()

  label.figure("B", cex = 1.5)

  if(pdf) dev.off()

}
