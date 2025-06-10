# JMDplots/genoGEO.R
# Make plots for paper:
# Genomes record the Great Oxidation Event
# 20231206 jmd first version
# 20240328 Moved to JMDplots
# 20240409 Add Rubisco plots
# 20240528 Analyze methanogen genomes
# 20240803 Compare Rubisco proteins and unrelated genomes (Fig. 3C)

# Figure 1: Genome-wide differences of oxidation state between two lineages of methanogens
genoGOE_1 <- function(pdf = FALSE, panel = NULL) {

  if(is.null(panel)) {
    if(pdf) pdf("Figure_1.pdf", width = 8, height = 6)
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

# Figure 2: Carbon oxidation state of proteins in eukaryotic gene age groups
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

# Figure 3: Evolutionary oxidation of ancestral Rubiscos and thermodynamic prediction of redox boundaries around the GOE
genoGOE_3 <- function(pdf = FALSE, panel = NULL) {

  if(is.null(panel)) {
    if(pdf) pdf("Figure_3.pdf", width = 9, height = 11)
    layout(matrix(c(0,1,1,0, 2,2,3,3, 4,4,5,5), nrow = 3, byrow = TRUE))
    par(cex = 1)
  }
  panels <- if(is.null(panel)) LETTERS[1:5] else panel

  # Read amino acid compositions
  fasta_file <- system.file("extdata/fasta/KHAB17.fasta", package = "canprot")
  aa <- read_fasta(fasta_file)
  # Assign protein names
  aa$protein <- sapply(strsplit(aa$protein, "_"), "[", 2)

  if("A" %in% panels) {

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
    abline(v = 3.5, lty = 2, col = 4, lwd = 2)
    text(2.7, -0.13, "GOE\n(estimated)")
    if(is.null(panel)) {
      title("Evolutionary oxidation of Rubiscos", font.main = 1)
      label.figure("A", cex = 1.5, font = 2, yfrac = 0.936)
    }

  }

  # Add proteins to CHNOSZ
  ip <- add.protein(aa, as.residue = TRUE)
  # Setup basis species and swap O2 for e- to make Eh-pH diagram
  basis("QEC+")
  swap.basis("O2", "e-")
  # Set resolution
  res <- 300
  
  # Panel B: Pairwise stability boundaries for Rubisco

  if("B" %in% panels) {

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

    text(5.5, 0.64, hyphen.in.pdf("Higher AFFINITY\nfor post-GOE protein\nin each pair"), cex = 0.8)
    text(4.5, -0.15, hyphen.in.pdf("Higher AFFINITY for\npre-GOE protein in each pair"), cex = 0.8, srt = -24)
    title("Pairwise relative stabilities of Rubiscos", font.main = 1)
    if(is.null(panel)) label.figure("B", cex = 1.5, font = 2, yfrac = 0.936)

  }

  if("C" %in% panels) {

    # Panel C: Groupwise stability boundaries for Rubisco

    # Calculate affinity of composition reactions for all proteins
    aout <- affinity(pH = c(0, 14, res), Eh = c(-0.5, 0.8, res), iprotein = ip)
    # Set up groups for affinity ranking:
    # 3 pre-GOE and 3 post-GOE proteins
    groups <- list(pre = c(TRUE, TRUE, TRUE, FALSE, FALSE, FALSE), post = c(FALSE, FALSE, FALSE, TRUE, TRUE, TRUE))
    arank <- rank.affinity(aout, groups = groups)
    diagram(arank, lwd = 2, col = 4, names = "", xlab = "pH", ylab = axis.label("Eh"))
    text(3.5, -0.07, hyphen.in.pdf("Higher affinity RANKING\nfor pre-GOE proteins"), col = 4, font = 2, cex = 0.8, srt = -24)
    text(3.5, 0.25, hyphen.in.pdf("Higher affinity RANKING\nfor post-GOE proteins"), col = 4, font = 2, cex = 0.8, srt = -24)
    title("Groupwise relative stabilities of Rubiscos", font.main = 1)
    # Add arrow
    arrows(7, -0.2, 7, 0.1, length = 0.2, lwd = 2, col = 4)
    text(9.5, 0.15, "Rubiscos became\nmore oxidized\nover the GOE", cex = 0.9)
    if(is.null(panel)) label.figure("C", cex = 1.5, font = 2, yfrac = 0.936)

  }

  if("D" %in% panels) {

    # Panel D: Comparison between Rubiscos and methanogen and Nitrososphaeria genomes
    yvar <- "O2"
    genoGOE_3D(yvar, res)

    # Label lines
    if(yvar == "Eh") {
      text(8.1, -0.07, "Rubiscos", cex = 0.8)
      text(7.4, -0.081, hyphen.in.pdf("Post-GOE"), srt = -27, cex = 0.75)
      text(7.33, -0.105, hyphen.in.pdf("Pre-GOE"), srt = -27, cex = 0.75)

      text(8.3, -0.19, "Methanogen\ngenomes", cex = 0.8)
      text(7.35, -0.215, "Class I", srt = -33, cex = 0.75)
      text(7.5, -0.2, "Class II", srt = -33, cex = 0.75)

      text(4.9, -0.18, "Nitrososphaeria\ngenomes", cex = 0.8)
      text(6.03, -0.17, "Basal", srt = -33, cex = 0.75)
      text(6.2, -0.155, "Terrestrial", srt = -33, cex = 0.75)
    }
    if(yvar == "O2") {
      text(7, -61, " Rubiscos", cex = 0.8, adj = 0)
      text(7, -61, hyphen.in.pdf("Pre-GOE "), cex = 0.75, adj = 1, srt = 20)
      text(7, -59.8, hyphen.in.pdf("Post-GOE "), cex = 0.75, adj = 1, srt = 20)

      text(7, -67, " Methanogen genomes", cex = 0.8, adj = 0)
      text(6.5, -67.8, "Class I ", cex = 0.75)
      text(6.5, -66.7, "Class II ", cex = 0.75)

      text(7, -69.2, " Nitrososphaeria genomes", cex = 0.8, adj = 0)
      text(6.5, -70.2, "Basal ", cex = 0.75)
      text(6.5, -69.2, "Terrestrial ", cex = 0.75)
    }

    if(is.null(panel)) {
      title("Groupwise relative stabilities of Rubiscos\nand unrelated genomes", font.main = 1, xpd = NA)
      label.figure("D", cex = 1.5, font = 2, yfrac = 0.936)
    }

  }

  if("E" %in% panels) {

    # Add more arrows 20240812
    plot.new()
    par(xpd = NA)
    if(yvar == "Eh") {
      arrows(0, 0.1, 0, 0.3, length = 0.2, lwd = 2, col = 7)
      arrows(0.1, 0.1, 0.1, 0.3, length = 0.2, lwd = 2, col = 2)
      arrows(0.05, 0.4, 0.05, 0.6, length = 0.2, lwd = 2, col = 4)
      text(0.5, 0.35, "1: Proteins in many organisms\nbecame more oxidized\nover the GOE")
      text(0.05, 0.85, hyphen.in.pdf("2: Rubisco records\nmore oxidizing conditions\ncompared to genomes of\nnon-photosynthetic organisms\nat the time of the GOE"))
    }
    if(yvar == "O2") {
      arrows(-0.08, 0.12, -0.08, 0.32, length = 0.2, lwd = 2, col = 7)
      arrows(-0.02, 0.18, -0.02, 0.38, length = 0.2, lwd = 2, col = 2)
      arrows(-0.05, 0.7, -0.05, 0.9, length = 0.2, lwd = 2, col = 4)
      text(0.05, 0.25, "1: Protein sequences in many lineages\nwere oxidized over the GOE", adj =0)
      text(0.05, 0.8, hyphen.in.pdf("2: Rubisco records more oxidizing conditions\ncompared to genomes of non-photosynthetic\norganisms at the time of the GOE"), adj = 0)
    }

  }

  if(pdf & is.null(panel)) dev.off()

}


# Code for Figure 3D
# Comparison of Rubiscos and methanogen and Nitrososphaeria genomes
genoGOE_3D <- function(yvar = "Eh", res = 400, add = FALSE, lwd = 2, pHlim = c(4, 10), Ehlim = c(-0.3, 0.1), O2lim = c(-72.5, -58), alpha.f = 1, Eh7_las = 1, datasets = 1:3) {

  # Setup basis species
  basis("QEC+")
  # Swap O2 for e- to make Eh-pH diagram
  if(yvar == "Eh") swap.basis("O2", "e-")

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
      col <- 2
      add <- add
    }

    if(i == 2) {
      # Thaumarchaeota 20220414
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
      col <- 7
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
    if(yvar == "Eh") aout <- affinity(pH = c(pHlim, res), Eh = c(Ehlim, res), iprotein = ip)
    if(yvar == "O2") aout <- affinity(pH = c(pHlim, res), O2 = c(O2lim, res), iprotein = ip)
    # Calculate average ranking for each group and make diagram
    arank <- rank.affinity(aout, groups)
    diagram(arank, col = adjustcolor(col, alpha.f = alpha.f), lwd = lwd, add = add, names = "")

  }

  # Add Eh7 axis 20241218
  if(yvar == "O2") {
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

}

# Evolutionary oxidation and relative stabilities for genomes with S-cycling genes 20241211
genoGOE_4 <- function(pdf = FALSE, panel = NULL) {
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
  dsr <- "#bcb2ce"
  sox <- "#45b78d"
  mdd <- "#c24a96"
  # Colors for protein groups
  col <- c(dsr, sox, sox, dsr, mdd, mdd, mdd)

  # Plot Zc of genomes with S-cycling genes from Mateos et al. (2023)  20240916
  sulfur_Zc <- function() {
    # Setup plot
    par(mar = c(4.1, 4.1, 4.1, 2.1))
    n <- length(Zclist)
    boxplot(Zclist, col = col, names = character(n), xlab = "Age of earliest gene event (Ga)", ylab = "Zc of all proteins in genome")
    text(2.5, -0.24, hyphen.in.pdf("Dissimilatory sulfate-sulfite-sulfide oxidation/reduction"), col = dsr)
    text(2.8, -0.12, hyphen.in.pdf("Sulfate-thiosulfate\noxidation/reduction"), col = sox)
    text(6, -0.22, "Organic sulfur cycling", col = mdd)
    # Ages from Table 2 of Mateos et al.
    ages <- c("3.3-3.35", "2.65-2.88", "2.6", "2.33-2.47", "2.28", "1.77", "0-2.38")
    axis(1, at = 1:n, labels = ages, lwd = 0)
    axis(3, at = 1:n, labels = names(Zclist), line = 0.5, lwd = 0, font = 3)
    # Add number of genomes to labels
    n_genomes <- paste0("(", sapply(Zclist, length), ")")
    axis(3, at = 1:n, labels = n_genomes, line = -0.5, lwd = 0)
    title(hyphen.in.pdf("Sulfur-cycling gene or gene cluster (# of exclusive genomes)"), font.main = 1, line = 3)
  }

  # Affinity ranking for genomes with different S-cycling genes
  sulfur_affinity <- function(panel) {
    if(is.null(panel)) par(mar = c(4.1, 4.1, 4.1, 4.1))
    basis("QEC+")
    # Keep genomes with single sulfur-cycling genes listed above
    myaa <- aa[aa$organism %in% unlist(genomes), ]
    # Load proteins for each genome
    ip <- add.protein(myaa, as.residue = TRUE)
    # Set resolution
    res <- 150
    # Calculate affinity of composition reactions as a function of Eh and pH
    a <- affinity(pH = c(3, 10, res), O2 = c(-72.5, -58, res), iprotein = ip)
    # Group genomes according to presence of sulfur-cycling genes
    groups <- sapply(genomes, function(genome) match(genome, myaa$organism))
    # Calculate normalized sum of ranks for each group and make diagram
    arank <- rank.affinity(a, groups)
    # Lighten colors
    fill <- adjustcolor(col, alpha.f = 0.3)
    # Adjust labels
    names <- names(genomes)
    dy <- rep(0, length(names))
    #dy[names == "mddA"] <- 1
    #dy[names == "soxC"] <- 1
    #names[5] <- ""
    diagram(arank, fill = fill, lty = 1, lwd = 1.5, font = 3, names = names, dy = dy, col = "gray40")
    #text(4.3, -68.6, "dmsA", font = 3)
    #lines(c(4.04, 3.68), c(-68.60, -68.76))
    #lines(c(3.28, 3.17), c(-65.30, -66.25))
  }

  if(is.null(panel)) {
    if(pdf) pdf("Figure_4.pdf", width = 8, height = 12)
    layout(matrix(1:2), heights = c(5, 7))
  }
  panels <- if(is.null(panel)) LETTERS[1:4] else panel

  if("A" %in% panels) {
    sulfur_Zc()
    if(is.null(panel)) label.figure("A", font = 2, cex = 1.6)
  }
  if("B" %in% panels) {
    sulfur_affinity(panel)
    if(is.null(panel)) {
      # Overlay stability boundaries for other genomes
      genoGOE_3D("O2", add = TRUE, lwd = 4, pHlim = c(3, 10), alpha.f = 0.7, Eh7_las = 0)
      title(main = hyphen.in.pdf("Groupwise relative stabilities of proteins in\ngenomes with different S-cycling genes"), font.main = 1)
      label.figure("B", font = 2, cex = 1.6)
    } else {
      # Only add Eh7 axis
      genoGOE_3D("O2", add = TRUE, Eh7_las = 0, datasets = numeric())
    }
  }
  if(pdf & is.null(panel)) dev.off()

}

# Carbon oxidation state of ancestral and extant nitrogenases
# 20240216 first version
# 20250325 added to JMDplots
# 20250327 use more representative extant nitrogenases (Garcia et al. Fig. 6)
genoGOE_5 <- function(pdf = FALSE, panel = NULL) {

  if(pdf & is.null(panel)) pdf("Figure_5.pdf", width = 7, height = 5.5)

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
  # Use colors from Garcia et al. (2020) Fig. 2
  cols <- c(A = "#df674f", B = "#b779dc", C = "#27a08d", D = "#e8c44c", E = "#759dce")
  pt_cols <- adjustcolor(cols, alpha.f = 0.75)
  names(pt_cols) <- names(cols)

  # Start plot
  plot(extendrange(c(1, 7)), c(-0.20, -0.12), xlab = "Form of nitrogenase", ylab = axis.label("ZC"), type = "n", axes = FALSE)
  axis(side = 1, at = seq_along(form_to_anc), labels = hyphen.in.pdf(names(form_to_anc)), gap.axis = 0)
  axis(side = 2)
  box()

  # Loop over nitrogenase forms
  for(iform in seq_along(form_to_anc)) {

    # Calculate Zc of the ancestral proteins
    node <- form_to_anc[[iform]]
    ianc <- grepl(paste0("^Anc", node), aa$protein)
    Zc_anc <- Zc(aa[ianc, ])

    # Plot lines for ancestral proteins
    dx <- 0.25
    for(Zc in Zc_anc) lines(c(iform-dx, iform+dx), c(Zc, Zc), col = cols[node])

    # Calculate Zc of the extant proteins
    form <- names(form_to_anc[iform])
    iext <- which(aa$ref == form)
    Zc_ext <- Zc(aa[iext, ])
    points(rep(iform, length(Zc_ext)), Zc_ext, pch = 19, col = pt_cols[node])

  }

  # Add legend and title
  legend("bottomright", c("Extant", "Ancestral"), lty = c(NA, 1), pch = c(19, NA), bty = "n")
  if(is.null(panel)) title("Carbon oxidation state of nitrogenase sequences", font.main = 1)

  if(pdf & is.null(panel)) dev.off()

}
