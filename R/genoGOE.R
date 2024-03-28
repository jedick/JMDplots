# JMDplots/genoGEO.R
# Plot for the manuscript by Dick (2024) about the genomic record of Earth's surface oxygenation
# 20231206 jmd first version
# 20240328 moved to JMDplots

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
