# JMDplots/utegig.R
# Plots for perspective paper on oxidation state of biomolecules
# 20200418 jmd first version
# 20210516 Add methanogen tree
# 20210727 Add nH2O-ZC overview plot [deleted]
# 20220220 Make hierarchy diagram
# 20220420 Move to JMDplots

# Create bold axis labels
ZClab <- quote(bolditalic(Z)[bold(C)])
Tlab <- quote(bolditalic(T)~bold("("*degree*C*")"))
Alab <- quote(bold("Affinity"~(kJ~"(mol C)"^{-1})))
Aresiduelab <- quote(bold("Affinity"~(kJ~"(mol residue)"^{-1})))
Toptlab <- quote(bolditalic(T)[bold(opt)]~bold("("*degree*C*")"))
logaH2maxlab <- quote(bold(log)~bolditalic(a)[bold(H[2])]~bold("for maximum activity"))
logalab <- quote(bold(log)~bolditalic(a))
logaH2lab <- quote(bold(log)~bolditalic(a)[bold(H[2])])
logaNH4lab <- quote(bold(log)~bolditalic(a)[bold(NH[4]^"+")])

# Identify species used in Shock and Canovas (2010)
# Used in Energy() and calc_logaH2_intermediate()
specieslist <- list(
  "C1 and C2" = c("CH4", "methanol", "formaldehyde", "CO", "formic acid",
           "ethane", "ethylene", "ethyne", "ethanol", "acetaldehyde", "acetic acid", "glycolic acid", "oxalic acid"),
  Acid = c("formic acid", "acetic acid", "propanoic acid", "n-hexanoic acid", "n-nonanoic acid", "n-dodecanoic acid"),
  "Amino acid" = c("glycine", "alanine", "valine", "leucine", "aspartic acid", "asparagine", "glutamic acid", "glutamine", "serine", "proline", "phenylalanine", "tryptophan"),
  Sugar = c("ribose", "ribulose", "deoxyribose", "glucose"),
  Nucleobase = c("adenine", "guanine", "thymine", "cytosine"),
  # TCA cycle metabolites added 20211110
  "TCA cycle" = c("pyruvate", "oxaloacetate-2", "citrate-3", "cis-aconitate-3", "isocitrate-3", "a-ketoglutarate-2", "succinate-2", "fumarate-2", "malate-2")
)

# Get seawater from Amend and Shock (1998)
Seawater.AS98 <- data.frame(T = 18, CO2 = log10(1e-4), H2 = log10(2e-9), pH = -log10(5e-9), "NH4+" = log10(5e-8), H2S = log10(1e-15), check.names = FALSE)

# Make semi-transparent colors 20220223
addalpha <- function(col, alpha) {
  x <- col2rgb(col)
  newcol <- rgb(x[1], x[2], x[3], maxColorValue = 255)
  newcol <- paste0(newcol, alpha)
  newcol
}
col4 <- addalpha(4, "b0")
col8 <- addalpha(8, "b0")
col2 <- addalpha(2, "b0")

# Read methanogen tree
methanogen_tree <- read.tree(system.file("extdata/utegig/methanogen_tree.txt", package = "JMDplots"))
# Match species names without underscore used for labeling the tree
methanogens <- gsub("_", " ", methanogen_tree$tip.label)
# Indices of Class I and Class II methanogens
iI <- 20:36
iII <- 1:19
# Assign point symbols and colors for methanogens
mcol <- mbg <- mpch <- rep(NA, length(methanogens))
# Open red circle for Class I
mpch[iI] <- 21
mbg[iI] <- "transparent"
mcol[iI] <- col2
# Filled red circle for Methanocaldococcus 20220220
mbg[grepl("Methanocaldococcus", methanogens)] <- col2
# Grey-filled red circle for outliers
iout <- methanogens %in% c("Methanothermobacter marburgensis", "Methanothermobacter thermautotrophicus", "Methanopyrus kandleri")
mbg[iout] <- col8
mcol[iout] <- col2
# Filled blue square for Class II
mpch[iII] <- 22
mbg[iII] <- col4
mcol[iII] <- col4

## To make methanogen_AA.csv
## Read RefSeq amino acid compositions
## NOTE: taxids are in 'organism' column, species names are in 'ref' column
#refseq <- read.csv(system.file("extdata/refseq/protein_refseq.csv.xz", package = "chem16S"), as.is = TRUE)
#irefseq <- match(methanogens, refseq$ref)
## Amino acid compositions of reference proteomes of methanogens
#methanogen_AA <- refseq[irefseq, ]
#write.csv(methanogen_AA, "methanogen_AA.csv", row.names = FALSE, quote = FALSE)

# ZC of proteins and lipids in hot springs 20210516
utegig1 <- function(pdf = FALSE) {

  if(pdf) pdf("Figure_1.pdf", width = 8.5, height = 4)
  layout(matrix(1:3, nrow = 1), widths = c(1.07, 1, 0.4))
  par(cex = 1)

  ## Plot A: Proteins
  par(mar = c(4, 4.2, 0.5, 1))
  plot(c(50, 100), c(-0.22, -0.12), type = "n", xlab = Tlab, ylab = ZClab, xaxs = "i", axes = FALSE)
  # Remove leading zeros for better visual parallel with lipids
  at <- seq(-0.22, -0.12, 0.02)
  labels <- gsub("-0", "-", formatC(at, digits = 2, format = "f"))
  axis(2, at, labels)
  axis(1)
  box()
  # Get ZC of Bison Pool proteins aggregated by phylum
  aa.phyla <- read.csv(system.file("extdata/bison/DS13.csv", package = "JMDplots"))
  ZC <- ZC(protein.formula(aa.phyla))
  # Get T values
  site <- gsub("bison", "", aa.phyla$protein)
  T <- sapply(site, switch, N = 93.3, S = 79.4, R = 67.5, Q = 65.3, P = 57.1)
  # Plot lines for each phlum
  for(phylum in unique(aa.phyla$organism)) {
    iphylum <- aa.phyla$organism == phylum
    lines(T[iphylum], ZC[iphylum], col = "gray")
  }

  # Make different color groups based on temperature
  ilo <- T < 60
  imid <- T >= 60 & T <= 70
  ihi <- T > 70
  points(T[ilo], ZC[ilo], pch = 22, bg = col4)
  points(T[imid], ZC[imid], pch = 23, bg = col8)
  points(T[ihi], ZC[ihi], pch = 21, bg = col2)
  # Add labels
  text(83, -0.14, "Proteins grouped\nby major phyla")
  label.figure("A", cex = 1.5, font = 2, xfrac = 0.06)

  ## Plot B: Lipids
  par(mar = c(4, 3, 0.5, 1))
  plot(c(25, 100), c(-2.2, -1.2), type = "n", xlab = Tlab, ylab = "", xaxs = "i")
  # Add "100" because it doesn't show up 20220408
  axis(1, 100, tick = FALSE)
  # IPL data from Boyer et al., 2020 Fig. 7
  T <- c(29.0, 35.1, 38.1, 40.5, 51.6, 53, 59.8, 60.7, 63.1, 64.8, 70.5, 73.3, 77.3, 80.9, 82.2, 85.4, 89.0, 91.0)
  ZC.alkyl <- c(-1.73, -1.71, -1.77, -1.80, -1.83, -1.82, -1.89, -1.90, -1.88, -1.92, -1.92, -1.84, -1.93, -1.89, -1.96, -1.95, -1.92, -1.94)
  ZC.IPL <- c(-1.36, -1.33, -1.37, -1.38, -1.44, -1.43, -1.48, -1.54, -1.46, -1.52, -1.54, -1.45, -1.61, -1.53, -1.60, -1.56, -1.56, -1.68)
  # Make different color groups based on temperature
  ilo <- T < 60
  imid <- T >= 60 & T <= 70
  ihi <- T > 70
  # Plot points
  points(T[ilo], ZC.alkyl[ilo], pch = 22, bg = col4)
  points(T[imid], ZC.alkyl[imid], pch = 23, bg = col8)
  points(T[ihi], ZC.alkyl[ihi], pch = 21, bg = col2)
  points(T[ilo], ZC.IPL[ilo], pch = 22, bg = col4)
  points(T[imid], ZC.IPL[imid], pch = 23, bg = col8)
  points(T[ihi], ZC.IPL[ihi], pch = 21, bg = col2)
  # Label molecular groups
  text(70, -1.32, "Intact polar lipids")
  text(50, -2, "Alkyl chains")
  label.figure("B", cex = 1.5, font = 2, xfrac = 0.01)

  # Add legend 20220221
  # Move to right panel 20220408
  plot.new()
  par(mar = c(0, 0, 0, 0))
  par(xpd = NA)
  legend <- as.expression(c(quote(" "*italic(T) < 60~degree*C*";"), " less reducing", "",
                            quote(" Transition zone"), "",
                            quote(" "*italic(T) > 70~degree*C*";"), " more reducing"))
  legend("left", legend = legend, pch = c(22, NA, NA, 23, NA, 21, NA), pt.bg = c(col4, NA, NA, col8, NA, col2, NA),
    bty = "n", inset = -1, y.intersp = 1)
  par(xpd = FALSE)

  if(pdf) dev.off()

}

# Overlay ZC on phylogenetic tree of methanogens 20210516
utegig2 <- function(pdf = FALSE) {

  if(pdf) pdf("Figure_2.pdf", width = 12, height = 8)
  layout(matrix(1:2, nrow = 1), widths = c(1.5, 0.7))
  par(mar = c(2, 0, 1, 32))

  # Use colors to distinguish Class I and Class II
  edge.color <- c(rep(4, 37), rep(2, 37))
  methanogen_tree <- read.tree(system.file("extdata/utegig/methanogen_tree.txt", package = "JMDplots"))
  pp <- plot.phylo(methanogen_tree, root.edge = TRUE, edge.color = edge.color, show.tip.label = FALSE)
  # Add tip labels 20220531
  labels <- methanogen_tree$tip.label
  labels <- gsub("_", " ", labels)
  FStxt <- "sp. FS406-22"
  iFS <- grep(FStxt, labels, fixed = TRUE)
  labels <- gsub(FStxt, "", labels, fixed = TRUE)
  dy <- -0.25
  text(3100, 1:36 + dy, labels, font = 3, adj = c(0, 0), xpd = NA)
  text(4170, iFS + dy, FStxt, adj = c(0, 0), xpd = NA)

  # Calculate ZC from protein formulas
  methanogen_AA <- read.csv(system.file("extdata/utegig/methanogen_AA.csv", package = "JMDplots"))
  ZC <- ZC(protein.formula(methanogen_AA))

  # Add labels
  text(0, 29.5, "Class I", font = 2, adj = 0, cex = 1.2)
  text(0, 5, "Class II", font = 2, adj = 0, cex = 1.2)
  text(3200, -2.4, ZClab, adj = c(0, 0), xpd = NA)
  label.figure("A", cex = 1.5, font = 2, yfrac = 0.97)
  # Setup a new plot in the empty space between the tree and names
  par(new = TRUE)
  par(mar = c(2, 9, 1, 18))
  plot(range(ZC), c(1, length(ZC)), type = "n", ylab = "", axes = FALSE)
  axis(1, at = seq(-0.25, -0.15, 0.01), labels = FALSE, tcl = -0.3)
  axis(1, at = c(-0.25, -0.2, -0.15))
  axis(3, at = seq(-0.25, -0.15, 0.01), labels = FALSE, tcl = -0.3)
  # Plot guidelines 20220405
  abline(v = c(-0.25, -0.20, -0.15), lty = 2, col = 8)
#  abline(h = 1:36, lty = 3, col = "gray40")
  # Plot lines from ZC to organism name
  for(i in 1:36) lines(c(ZC[i], -0.145), c(i, i), lty = 3, col = "gray40")

  # Plot the ZC
  points(ZC, c(iII, iI), pch = mpch, col = mcol, bg = mbg, cex = 1.2)

  # ZC-Topt plot 20220224
  par(mar = c(9, 2, 7, 1))
  # Read Topt
  dat <- read.csv(system.file("extdata/utegig/Topt.csv", package = "JMDplots"))
  # Append ZC column
  methanogen_AA <- read.csv(system.file("extdata/utegig/methanogen_AA.csv", package = "JMDplots"))
  ZC <- ZC(protein.formula(methanogen_AA))
  dat$ZC <- ZC[match(dat$species, methanogens)]
  # Use par(xpd = NA) to show the y-axis label 20220401
  par(xpd = NA)
  plot(dat$Topt, dat$ZC, pch = mpch, col = mcol, bg = mbg, xlab = Toptlab, ylab = ZClab, ylim = c(-0.26, -0.14))
  par(xpd = FALSE)
  text(50, -0.24, "Class I", font = 2)
  text(40, -0.145, "Class II", font = 2)
  # Uncomment this to check the alignment of separate text items below
  #text(85, -0.168, "High-ZC\nthermophiles\n(Class I)", col = 2)
  text(85, -0.163, quote("High-"*italic(Z)[C]))
  text(85, -0.1705, "thermophiles\n(Class I)")
  text(80, -0.257, "Methanocaldococcus", font = 3)
  # Add convex hulls 20220401
  dat2 <- dat[iII, ]
  i2 <- chull(dat2$Topt, dat2$ZC)
  polygon(dat2$Topt[i2], dat2$ZC[i2], col = "#90CAF980", border = NA)
  dat1 <- dat[setdiff(iI, which(iout)), ]
  i1 <- chull(dat1$Topt, dat1$ZC)
  polygon(dat1$Topt[i1], dat1$ZC[i1], col = "#E5737380", border = NA)
  label.figure("B", cex = 1.5, font = 2, yfrac = 0.85, xfrac = -0.05)

  if(pdf) dev.off()

}

# Affinities of organic synthesis reactions 20211109
utegig3 <- function(pdf = FALSE) {

  if(pdf) pdf("Figure_3.pdf", width = 7, height = 6)
  mat <- matrix(c(1,1,2,2,3,3,4, 0,5,5,5,5,0,0), byrow = TRUE, nrow = 2)
  layout(mat, heights = c(1, 0.7), widths = c(1.1,1,1,1,1,1,1))
  par(mar = c(4, 4, 0.5, 1), mgp = c(2.7, 1, 0))

  # Set temperature and logaH2 activities
  T <- 25
  logaH2s <- c(-5, -10.27, -15)
  # Use seawater composition from Amend and Shock (1998)
  # NOTE: H2S is present but is not in any formation reactions 20220403
  basis(c("CO2", "H2", "NH4+", "H2O", "H2S", "H+"), c(Seawater.AS98$CO2, Seawater.AS98$H2, Seawater.AS98$"NH4+", 0, Seawater.AS98$H2S, -Seawater.AS98$pH))
  # Position of slope text
  slopey <- c(-112, -160, -70)

  # Loop over logaH2 activities
  for(ilogaH2 in c(1, 2, 3)) {

    logaH2 <- logaH2s[ilogaH2]
    basis("H2", logaH2)

    # Start plot
    xlim <- c(-4, 3)
    ylim <- c(-310, 100)
    if(ilogaH2 == 1) ylab <- Alab
    if(ilogaH2 > 1) {
      ylab <- ""
      par(mar = c(4, 3, 0.5, 1))
    }

    plot(xlim, ylim, type = "n", xlab = ZClab, ylab = ylab)
    abline(h = 0, lty = 4)

    pchs <- c(21:25, 20)
    cols <- sapply(c(7, 5, 8, 5, 3, 1), addalpha, alpha = "b0")

    # Loop over groups of species
    # Start with n-carboxylic acids (i = 2)
    AC_all <- ZC_all <- numeric()
    for(i in c(2, 1, 3, 4, 5, 6)) {

      # Load species
      myspecies <- specieslist[[i]]
      species(myspecies, -3)
      # Calculate affinities
      a <- affinity(T = T)
      A <- unlist(a$values)
      # Convert to Gibbs energy (J/mol)
      TK <- convert(T, "K")
      A <- -convert(A, "G", TK)
      if(packageVersion("CHNOSZ") <= "1.4.3") A <- convert(A, "J")
      # Convert to kJ/mol
      A <- A / 1000
      # Divide A by number of carbon
      nC <- sapply(makeup(info(myspecies)), "[", "C")
      AC <- A / nC
      # Calculate carbon oxidation state
      ZC <- ZC(info(myspecies))
      # Adjust point size for n-carboxylic acids
      cex <- 1
      if(i == 2) cex <- 1.5
      if(i < 6) points(ZC, AC, pch = pchs[i], bg = cols[i], cex = cex)
      else points(ZC, AC, pch = pchs[i], col = cols[i], cex = cex)
      AC_all <- c(AC_all, AC)
      ZC_all <- c(ZC_all, ZC)

    }

    # Calculate linear fit
    thislm <- lm(AC_all ~ ZC_all)
    x <- range(ZC_all)
    y <- predict.lm(thislm, data.frame(ZC_all = x))
    # Plot linear fit
    lines(x, y, lty = 2, lwd = 2, col = 8)
    # Add legend with slope
    slope <- thislm$coefficients[[2]]
    word <- "Slope"
    if(ilogaH2 == 1) {
      word <- "linear fit"
      text(-1.53, -95, "Slope of")
    }
    slopetxt <- bquote(.(word) == .(formatC(slope, digits = 1, format = "f")))
    if(round(slope, 1) == 0) slopetxt <- "Slope = 0"
    text(-0.5, slopey[ilogaH2], slopetxt)

    # Add title: logaH2
    Htxt <- bquote(bold(log*bolditalic(a)[H[2]] == .(logaH2)))
    title(Htxt, line = -1, cex.main = 1)
    if(ilogaH2 == 3) {
      # Add oxidation-intermediate-reduction arrow 20220221
      par(xpd = NA)
      x1 <- -21.0
      x2 <- -0.5
      dx <- 3.2
      rect(x1 - dx, -317, x2 + dx, -227, col = "white", lty = 3)
      x1.1 <- x1 + 6
      x2.1 <- x2 - 6
      arrows(x1.1, -300, x1, -300, col = "#D50000", length = 0.1, lwd = 2)
      arrows(x2.1, -300, x2, -300, col = "#366BFB", length = 0.1, lwd = 2)
      lines(c(x1.1, x2.1), c(-300, -300), col = 8, lwd = 2)
      # Uncomment this to check the alignment of separate text items below
      #text(x1, -288, "Reducing conditions:\nReduced compounds\nare more stable", adj = c(0.5, 0), col = 2)
      text(x1, -248, "Reducing conditions:", adj = c(0.5, 0), font = 3)
      text(x1, -288, "Reduced compounds\nare more stable", adj = c(0.5, 0))
      #text(x2, -288, "Oxidizing conditions:\nOxidized compounds\nare more stable", adj = c(0.5, 0), col = 2)
      text(x2, -248, "Oxidizing conditions:", adj = c(0.5, 0), font = 3)
      text(x2, -288, "Oxidized compounds\nare more stable", adj = c(0.5, 0))
      #text(mean(c(x1, x2)), -288, "Intermediate conditions:\nReduced and oxidized compounds\nare approximately equally stable", adj = c(0.5, 0), col = 2)
      text(mean(c(x1, x2)), -248, "Intermediate conditions:", adj = c(0.5, 0), font = 3)
      text(mean(c(x1, x2)), -288, "Reduced and oxidized compounds\nare approximately equally stable", adj = c(0.5, 0))
      par(xpd = FALSE)
    }
    if(ilogaH2 == 1) {
      text(-3.1, 72, quote(CH[4]))
      label.figure("A", cex = 1.5, font = 2, yfrac = 0.97)
    }

  }

  # Add legend 20220407
  plot.new()
  par(mar = c(0, 0, 0, 0))
  par(xpd = NA)
  legend <- paste0(" ", names(specieslist))
  legend(-0.8, 0.7, legend, pch = pchs, pt.cex = c(1, 1.5, 1, 1, 1, 1),
         pt.bg = cols, col = c(1, 1, 1, 1, 1, cols[6]), bty = "n", y.intersp = 1.2)
  dT <- describe.property("T", T)
  dbasis <- describe.basis(ibasis = c(1, 3, 6))
  legend("topright", legend = c(dT, dbasis), bty = "n", y.intersp = 1.2)
  par(xpd = FALSE)

  # Plot intermediate logaH2 20220401
  intermediate_logaH2()
  label.figure("B", cex = 1.5, font = 2)

  if(pdf) dev.off()

}

# Affinities for methanogen proteomes and logaH2 of habitable niches 20220402
utegig4 <- function(pdf = FALSE) {

  if(pdf) pdf("Figure_4.pdf", width = 7, height = 5)
  mat <- matrix(c(1,1,2,2,3,3,4, 5,5,5,6,6,6,0), nrow = 2, byrow = TRUE)
  layout(mat, widths = c(1.2,1,1,1,1,1,1.2))
  par(mar = c(4, 4, 0.5, 1), mgp = c(2.7, 1, 0))

  # Define basis species and temperature
  basis(c("CO2", "H2", "NH4+", "H2O", "H2S", "H+"), c(Seawater.AS98$CO2, Seawater.AS98$H2, Seawater.AS98$"NH4+", 0, Seawater.AS98$H2S, -Seawater.AS98$pH))
  T <- 25

  # Calculate ZC of methanogen proteomes excluding high-ZC outliers
  methanogen_AA <- read.csv(system.file("extdata/utegig/methanogen_AA.csv", package = "JMDplots"))
  ZC <- ZC(protein.formula(methanogen_AA[!iout, ]))
  # Add proteins to CHNOSZ
  ip <- add.protein(methanogen_AA[!iout, ])
  pl <- protein.length(ip)

  # Create combinations between Class I and II methanogens
  # Class I methanogens excluding high-ZC outliers
  iIin <- setdiff(1:length(ip), iII)
  ncomb <- length(iIin) * length(iII)
  i1 <- rep(iIin, length.out = ncomb)
  i2 <- rep(iII, each = length(iIin))

  ## Plot A: affinities at specified logaH2
  logaH2s <- c(-5, -7.467, -10)
  for(i in 1:3) {
    if(i == 1) ylab <- Aresiduelab
    if(i > 1) {
      ylab <- ""
      par(mar = c(4, 2.5, 0.5, 1))
    }
    basis("H2", logaH2s[i])
    # Calculate affinities
    a <- affinity(iprotein = ip, T = T)
    A <- unlist(a$values)
    # Convert to Gibbs energy (J/mol)
    TK <- convert(T, "K")
    A <- -convert(A, "G", TK)
    if(packageVersion("CHNOSZ") <= "1.4.3") A <- convert(A, "J")
    # Convert to kJ/mol
    A <- A / 1000
    # Calculate affinity per residue
    Ares <- A / pl
    plot(ZC, Ares, pch = mpch[!iout], col = mcol[!iout], bg = mbg[!iout], xlab = ZClab, ylab = ylab)
    # Percentage of combinations that favor Class I over Class II
    p1 <- round(sum(sign(Ares[i1] - Ares[i2]) == 1) / ncomb * 100)
    legend <- bquote(.(p1)*"%"~A[I] > A[II])
    if(i < 3) inset <- c(-0.1, 0) else inset <- c(0.05, 0)
    legend("bottomleft", legend = legend, bty = "n", inset = inset)
    if(i == 1) legend("bottomleft", c("Pairwise affinity", "differences:", ""), inset = c(-0.1, 0), bty = "n")
    # Add legend: logaH2
    Htxt <- bquote(bold(log*bolditalic(a)[H[2]] == .(logaH2s[i])))
    if(i > 1) legend.x <- "topleft" else legend.x <- "top"
    legend(legend.x, legend = Htxt, bty = "n")
    if(i==1) label.figure("A", cex = 1.5, font = 2, yfrac = 0.97, xfrac = 0.06)
  }

  # Add legend 20220403
  plot.new()
  par(mar = c(0, 0, 0, 0))
  par(xpd = NA)
  legend <- c("Class I", "", "Class II")
  legend(-0.6, 0.6, legend = legend, pch = c(21, NA, 22), col = c(col2, NA, col4), pt.bg = c("white", NA, col4), bty = "n", y.intersp = 1.2)
  legend(-0.45, 0.52, legend = c("Methanocaldococcus"), pch = 21, col = col2, pt.bg = col2, text.font = 3, bty = "n", cex = 0.9, pt.cex = 1)
  dT <- describe.property("T", T)
  dbasis <- describe.basis(ibasis = c(1, 3, 5, 6))
  legend("topright", legend = c(dT, dbasis), bty = "n", y.intersp = 1.1)
  par(xpd = FALSE)

  ## Plot B: Contour plot for affinity differences between Class I and Class II proteomes 20220402
  # Calculate affinities as a function of logaH2 and T
  par(mar = c(4, 4, 1, 2.2))
  H2 <- c(-10, -1.5, 200)
  T <- c(20, 100, 200)
  a <- affinity(T = T, H2 = H2, iprotein = ip)
  # Calculate affinity per residue
  Ares <- mapply("/", a$values, pl, SIMPLIFY = FALSE)
  # Calculate sign of difference between all combinations
  Ares.diff <- mapply("-", Ares[i1], Ares[i2], SIMPLIFY = FALSE)
  Ares.diff.sign <- sapply(Ares.diff, sign, simplify = FALSE)
  Ares.diff.positive <- sapply(Ares.diff.sign, "==", y = 1, simplify = FALSE)
  # Sum to get percentage of combinations that favor Class I over Class I
  Ares.diff.positive.sum <- Reduce("+", Ares.diff.positive)
  p1 <- Ares.diff.positive.sum / ncomb * 100
  levels <- c(10, 50, 100)
  contour(a$vals$T, a$vals$H2, p1, xlab = Tlab, ylab = logaH2lab, xaxs = "i", yaxs = "i",
          levels = c(10, 50, 100), drawlabels = FALSE,
          col = c("#366BFB", 1, "#D50000"))
  # Plot intermediate logaH2 background
  intermediate_logaH2(add = TRUE, lines = FALSE, pH = FALSE, NH4 = FALSE, redox = FALSE)
  # Replot contour lines to make them stand out
  contour(a$vals$T, a$vals$H2, p1, add = TRUE,
          levels = c(10, 50, 100), drawlabels = FALSE,
          col = c("#366BFB", 1, "#D50000"))
  arrows(60, -7.8, 60, -4.6, length = 0.1, code = 3)
  text(55, -5.6, "More reducing", adj = 1)
  text(55, -7.9, "Less reducing", adj = 1)
  text(80, -2.6, "Class I proteomes\nare most stable", cex = 0.8)
  text(80, -4.5, "Class I proteomes\nare more stable", cex = 0.8)
  text(80, -6.8, "Class II proteomes\nare more stable", cex = 0.8)
  # Add contour labels in margin 20220407
  CL <- contourLines(a$vals$T, a$vals$H2, p1, levels = c(10, 50, 100))
  par(xpd = NA)
  text(101, tail(CL[[1]]$y, 1), "10%", adj = 0, col = "#366BFB")
  text(101, tail(CL[[2]]$y, 1), "50%", adj = 0, col = 1)
  text(101, tail(CL[[3]]$y, 1), "100%", adj = 0, col = "#D50000")
  par(xpd = FALSE)
  label.figure("B", cex = 1.5, font = 2, yfrac = 0.97)

  ## Plot C: logaH2-T plot for habitable niches 20220303
  par(mar = c(4, 2.5, 1, 2.5))
  plot(T[1:2], H2[1:2], xlab = Tlab, ylab = "", type = "n", xaxs = "i", yaxs = "i")

  # Fill area for habitable region
  pu <- par("usr")
  rect(pu[1], pu[3], pu[2], pu[4], col = "#C5E1A580")

  # Add line for mixed fluids at Rainbow 20220303
  srt <- 11
  Rainbow <- read.csv(system.file("extdata/mjenergy/SC10_Rainbow.csv", package = "JMDplots"))
  Tvals <- Rainbow$T
  logaH2.max <- Rainbow$H2
  polygon(c(Tvals, rev(Tvals)), c(logaH2.max, rep(par("usr")[4], length(Tvals))), col = "gray80", lty = 0)
  lines(Tvals, logaH2.max, lty = 2)
  text(Tvals[14], logaH2.max[14], "Rainbow", srt = 0.7 * srt, adj = c(0, 1.3), col = "green4")

  # Add logaCO2 = logaCH4 line 20220223
  Tvals <- seq(par("usr")[1], par("usr")[2], length.out = 100)
  logaH2.min <- subcrt(c("CO2", "H2", "CH4", "H2O"), c(-1, -4, 1, 2), T = Tvals)$out$logK / -4
  polygon(c(Tvals, rev(Tvals)), c(logaH2.min, rep(par("usr")[3], length(Tvals))), col = "gray80", lty = 0)
  lines(Tvals, logaH2.min, lty = 3)
  box()
  text(83, -6.4, quote(italic(a)[CH[4]] / italic(a)[CO[2]] == 10^0), srt = srt)
  # Add lines for aCO2/aCH4 = 10^2 and 10^4 20220303
  lines(Tvals[33:100], logaH2.min[33:100] + 2, lty = 3)
  text(94, -4.1, quote(10^2), srt = srt)
  lines(Tvals, logaH2.min + 4, lty = 3)
  text(94, -2.0, quote(10^4), srt = srt)

  # Add box for Lost City 20220303
  # 40 to 90 °C (Kelley et al., 2005)
  x1 <- 40
  x2 <- 90
  # 1 to 15 mmol/kg H2 (Kelley et al., 2005)
  y1 <- log10(1e-3)
  y2 <- log10(15e-3)
  rect(x1, y1, x2, y2, col = "white", lty = 0)
  rect(x1, y1, x2, y2, col = "#C5E1A540", border = "green4")
  text(mean(c(x1, x2)), mean(c(y1, y2)), "Lost City", col = "green4")

  # Add lines for sediments (Lovley and Goodwin, 1988) 20220303
  sediment <- read.csv(system.file("extdata/utegig/LG88_Fig1.csv", package = "JMDplots"))
  for(i in 1:nrow(sediment)) {
    if(sediment$type[i] == "Methane") col <- "green4" else col <- "gray40"
    if(sediment$type[i] == "Methane") lty <- 1 else lty <- 2
    lines(c(15, 25), rep(log10(sediment$H2[i]*1e-9), 2), col = col, lty = lty)
  }
  
  # Add affinity contours 20220403
  lines(as.data.frame(do.call(cbind, CL[[1]][c("x", "y")])), col = "#366BFB")
  lines(as.data.frame(do.call(cbind, CL[[2]][c("x", "y")])), col = 1)
  lines(as.data.frame(do.call(cbind, CL[[3]][c("x", "y")])), col = "#D50000")

  # Add labels 20220221
  text(21, -6.4, "Methanogenic\nsediments", adj = 0, cex = 0.9, col = "green4")
  lines(c(26, 22), c(-7, -7.9), col = "green4")
  text(40, -9.5, "Non-methanogenic sediments", adj = 0, cex = 0.9, col = "gray40")
  lines(c(26, 39), c(-9, -9.5), col = "gray40")
  par(xpd = NA)
  lines(c(100, 130), c(-2.34, -2.34), lty = 2)
  lines(c(100, 130), rep(tail(logaH2.min, 1), 2), lty = 3)
  text(132, -4, "Habitable\nniches", font = 3, col = "green4", adj = 0)
  arrows(130,-2.34, 130, tail(logaH2.min, 1), code = 3, col = "green4", length = 0.1)
  text(101, tail(CL[[1]]$y, 1), quote("10%"~A[I] > A[II]), adj = 0, col = "#366BFB")
  text(101, tail(CL[[2]]$y, 1), quote("50%"~A[I] > A[II]), adj = 0, col = 1)
  text(101, tail(CL[[3]]$y, 1), quote("100%"~A[I] > A[II]), adj = 0, col = "#D50000")
  par(xpd = FALSE)
  label.figure("C", cex = 1.5, font = 2, yfrac = 0.97, xfrac = -0.02)

  if(pdf) dev.off()

}

# Thaumarchaeota ZC and thermodynamic analysis 20220601
utegig5 <- function(pdf = FALSE) {

  if(pdf) pdf("Figure_5.pdf", width = 17, height = 12)
  layout(matrix(1:6, nrow = 2, byrow = TRUE), widths = c(2.5, 4, 4))
  par(cex = 1.5, mar = c(4, 4, 3, 1))
  ylim <- c(-0.25, -0.10)
  
  # Function to add significant difference letters 20220531
  cldfun <- function(ZClist, bp) {
    # Ugly one-liner to turn a list into a data frame with "group" column taken from names of the list elements
    ZCdat <- do.call(rbind, sapply(1:length(ZClist), function(i) data.frame(group = names(ZClist)[i], ZC = ZClist[[i]]), simplify = FALSE))
    # One-way ANOVA and Tukey's Honest Significant Differences
    # Adapted from https://statdoe.com/one-way-anova-and-box-plot-in-r/
    anova <- aov(ZC ~ group, data = ZCdat)
    tukey <- TukeyHSD(anova)
    # Compact letter display
    cld <- multcompLetters4(anova, tukey, reversed = TRUE)$group$Letters
    # Get into same order as data
    cld <- cld[match(names(ZClist), names(cld))]
    # Add to plot
    n <- length(ZClist)
    text((1:n) + 0.35, bp$stats[4, ] + 0.0045, cld)
  }

  # Methanogen proteomes 20220424
  # Read amino acid composition and compute ZC
  methanogen_AA <- read.csv(system.file("extdata/utegig/methanogen_AA.csv", package = "JMDplots"))
  ZC <- ZC(protein.formula(methanogen_AA))
  # Indices of Class I and Class II methanogens
  iI <- 20:36
  iII <- 1:19
  # Faded colors
  col4 <- addalpha(4, "b0")
  col2 <- addalpha(2, "b0")
  # Make plot
  ZClist <- list("Class I" = ZC[iI], "Class II" = ZC[iII])
  names(ZClist) <- paste0(names(ZClist), "\n(", sapply(ZClist, length), ")")
  bp <- boxplot(ZClist, ylab = ZClab, col = c(col2, col4), ylim = ylim, names = character(2))
  axis(1, at = 1:2, labels = names(ZClist), line = 1, lwd = 0)
  cldfun(ZClist, bp)
  text(1, -0.12, "Anoxic\nhabitats", font = 2, cex = 0.8)
  text(2, -0.12, "Anoxic\nand oxic\nhabitats", font = 2, cex = 0.8)
  abline(v = 1.5, lty = 2, lwd = 2, col = 8)
  title("Methanogens\n(Lyu and Lu, 2018)", font.main = 1, cex.main = 1)
  label.figure("A", cex = 1.5, font = 2, xfrac = 0.06)

  # Nif-bearing organisms 20220531
  ZClist <- NifProteomes()$ZClist
  bp <- boxplot(ZClist, ylab = ZClab, col = c(col2, col2, col4, col4), ylim = ylim, names = character(4))
  names(ZClist) <- paste0(names(ZClist), "\n(", sapply(ZClist, length), ")")
  axis(1, at = 1:4, labels = names(ZClist), line = 1, lwd = 0)
  # Names with "-" confuse multcompLetters4()
  names(ZClist) <- gsub("-", "", names(ZClist))
  cldfun(ZClist, bp)
  text(1.5, -0.12, "Anaerobic", font = 2, cex = 0.8)
  text(3.2, -0.11, "Anaerobic\nand aerobic", font = 2, cex = 0.8)
  abline(v = 2.5, lty = 2, lwd = 2, col = 8)
  title("Nif-bearing organisms\n(Poudel et al., 2018)", font.main = 1, cex.main = 1)

  # Thaumarchaeota 20220414
  aa <- read.csv(system.file("extdata/utegig/RFH+19_aa.csv", package = "JMDplots"))
  group <- c("Basal", "Terrestrial", "Shallow", "Deep")
  ZC <- ZCAA(aa)
  ZClist <- lapply(group, function(g) ZC[aa$protein == g])
  names(ZClist) <- group
  names(ZClist) <- paste0(names(ZClist), "\n(", sapply(ZClist, length), ")")
  bp <- boxplot(ZClist, ylab = ZClab, col = c(col2, col4, col4, col4), ylim = ylim, names = character(4))
  cldfun(ZClist, bp)
  axis(1, at = 1:4, labels = names(ZClist), line = 1, lwd = 0, gap.axis = 0)
  text(0.9, -0.12, "Pre-GOE\nemergence", font = 2, cex = 0.8)
  text(2.1, -0.12, "Post-GOE\nemergence", font = 2, cex = 0.8)
  abline(v = 1.5, lty = 2, lwd = 2, col = 8)
  title("Thaumarchaeota\n(Ren et al., 2019)", font.main = 1, cex.main = 1)

  # Define basis species and temperature
  basis(c("CO2", "H2", "NH4+", "H2O", "H2S", "H+"), c(Seawater.AS98$CO2, Seawater.AS98$H2, Seawater.AS98$"NH4+", 0, Seawater.AS98$H2S, -Seawater.AS98$pH))
  T <- 25

  # Add legend
  plot.new()
  par(mar = c(0, 0, 0, 0))
  par(xpd = NA)
  dT <- describe.property("T", T)
  dbasis <- describe.basis(ibasis = c(1, 5, 6))
  legend("topright", legend = c(dT, dbasis), bty = "n", y.intersp = 1.5, title = "Thaumarchaeota", title.cex = 1.1)
  par(xpd = FALSE)
  par(cex = 1.5, mar = c(3, 4, 3, 1))
  label.figure("B", cex = 1.5, font = 2, xfrac = 0.06)

  if(!packageVersion("CHNOSZ") > "1.4.3") {
    warning("Not making average affinity ranking plots because of insufficient CHNOSZ version")
  } else {

    # Affinity ranking diagram for Thaumarchaeota 20220415
    # Load protein residues
    ip <- add.protein(aa, as.residue = TRUE)
    # Names of groups and colors from Ren et al. (2019)
    groupnames <- c("Basal", "Terrestrial", "Shallow", "Deep")
    col <- c("#b2427e", "#c78d55", "#00a06f", "#4085c3")
    # Get the species in each group
    groups <- sapply(groupnames, function(group) aa$protein == group, simplify = FALSE)

    # Calculate logK for H2 = 2H+ + 2e- for conversion to Eh
    logK <- subcrt(c("H2", "H+", "e-"), c(-1, 2, 2), T = T)$out$logK
    pH <- Seawater.AS98$pH

    # Make two plots with different NH4+ activity
    dx <- list(c(0, 0, -2, -1), c(0, -1, 0, -1))
    dy <- list(c(1, 2, 5, -3), c(1, 2, 2, -3))
    for(i in 1:2) {
      if(i == 2) basis("NH4+", 0)
      a <- affinity(H2 = c(0, -15), iprotein = ip, T = T)
      # Calculate normalized sum of ranks for each group and make diagram
      arank <- CHNOSZ::affinity_rank(a, groups)
      diagram(arank, xlim = c(0, -15), col = col, lty = 1, lwd = 2, dx = dx[[i]], dy = dy[[i]],
              xlab = logaH2lab, ylab = "Average affinity ranking", ylim = c(0, 50))
      # Add logaNH4 legend
      dbasis <- describe.basis(ibasis = 3)
      legend("bottomright", legend = dbasis, bty = "n")

      for(j in 1:2) {
        # Don't draw Basal-Shallow on first plot
        if(j > i) next
        # Draw line at Basal-Terrestrial or Basal-Shallow transition
        if(j == 1) itrans <- which.min(abs(arank$values$Basal - arank$values$Terrestrial))
        if(j == 2) itrans <- which.min(abs(arank$values$Basal - arank$values$Shallow))
        logaH2 <- arank$vals$H2[itrans]
        A <- arank$values$Basal[itrans]
        lines(c(logaH2, logaH2), c(A, 53), lty = 2, lwd = 2, col = 8, xpd = NA)

        # Calculate pe and Eh (mV)
        pe <- -pH - 0.5*logaH2 - 0.5*logK
        Eh <- 1000 * convert(pe, "Eh", pH = pH)
        Ehtext <- paste(round(Eh), "mV")
        text(logaH2, 54, Ehtext, xpd = NA)
      }
    }

  }

  if(pdf) dev.off()

}


# Plot intermediate logaH2 20220220
# Add class argument 20220418
intermediate_logaH2 <- function(class = NULL, add = FALSE, parargs = list(mar = c(4, 4, 1, 1), mgp = c(2.7, 1, 0)),
  xlim = c(0, 300), ylim = NULL, lines = TRUE, redox = TRUE, pH = TRUE, NH4 = TRUE, pH.y = -9.4, NH4.y = -6.6) {
  if(!add) do.call(par, parargs)
  file <- "H2_intermediate.csv"
  if(!is.null(class)) file <- paste0("H2_intermediate_", gsub(" ", "", class), ".csv")
  dat <- read.csv(file.path(system.file("extdata/utegig", package = "JMDplots"), file))
  x <- dat[, "T", drop = FALSE]
  y <- dat[, 2:5]
  if(!add) matplot(x, y, type = "n", xlab = Tlab, ylab = logaH2lab, xlim = xlim, ylim = ylim, xaxs = "i")
  # Fill area between lines
  polygon(c(dat$T, rev(dat$T)), c(dat$loN_lopH, rev(dat$hiN_hipH)), border = NA, col = "#9E9E9E80")
  polygon(c(dat$T, rev(dat$T)), c(dat$hiN_lopH, rev(dat$loN_hipH)), border = NA, col = "#9E9E9EB0")
  # Fill reducing and oxidizing regions 20220401
  top <- rep(par("usr")[4], length(dat$T))
  polygon(c(dat$T, rev(dat$T)), c(top, rev(dat$hiN_hipH)), border = NA, col = "#E5737380")
  bot <- rep(par("usr")[3], length(dat$T))
  polygon(c(dat$T, rev(dat$T)), c(bot, rev(dat$loN_lopH)), border = NA, col = "#90CAF980")

  # Don't plot ends of lines 20220401
  Tlim <- c(10, 290)
  if(!is.null(class)) Tlim <- c(17, 270)
  if(pH) dat <- dat[dat$T > Tlim[1], ]
  if(NH4) dat <- dat[dat$T < Tlim[2], ]

  x <- dat[, "T", drop = FALSE]
  y <- dat[, 2:5]
  # dashed/solid lines for lo/hi N
  lty <- c(2, 2, 1, 1)
  ## red/blue lines for lo/hi pH
  #col <- c(2, 4, 2, 4)
  col <- c(1, 1, 1, 1)

  if(lines) {
    matlines(x, y, lty = lty, col = col)
  }
  if(pH) {
    # Add pH labels at low T
    pHdat <- dat[1, ]
    pHlab <- character(ncol(pHdat))
    pHlab[grep("lopH", colnames(pHdat))] <- 3
    pHlab[grep("hipH", colnames(pHdat))] <- 9
    pHdat.x <- 10
    pH.x <- 8
    dy <- -0.1
    if(!is.null(class)) {
      pHdat.x <- c(14, 7, 14, 7)
      pH.x <- 12
      dy <- c(-0.1, -0.6, -0.1, -0.6)
    }
    text(pHdat.x, pHdat + dy, pHlab)
    text(pH.x, pH.y, "pH", font = 2)
  }
  if(NH4) {
    # Add logaNH4+ labels at high T
    NH4dat <- dat[nrow(dat), ]
    NH4lab <- character(ncol(NH4dat))
    NH4lab[grep("loN", colnames(NH4dat))] <- -8
    NH4lab[grep("hiN", colnames(NH4dat))] <- -3
    NH4dat.x <- 293
    NH4.x <- 278
    dy <- 0
    if(!is.null(class)) {
      NH4dat.x <- c(278, 290, 290, 278)
      NH4.x <- 270
      dy <- c(0, 0.5, 0.5, 0)
    }
    text(NH4dat.x, NH4dat + dy, NH4lab)
    text(NH4.x, NH4.y, logaNH4lab, font = 2)
  }
  if(redox) {
    # Add more/less oxidizing/reducing labels 20220401
    arrows(100, -12.2, 100, -9.8, length = 0.1, code = 3)
    text(110, -12, "More oxidizing", adj = 0)
    text(110, -10, "Less oxidizing", adj = 0)
    arrows(100, -6.6, 100, -3.6, length = 0.1, code = 3)
    text(90, -3.8, "More reducing", adj = 1)
    text(90, -6.4, "Less reducing", adj = 1)
  }
}

# Find intermediate logaH2: where affinity vs ZC has a slope of 0  20220219
# Add class argument (calculate only for single compound class) 20220418
calc_logaH2_intermediate <- function(class = NULL) {

  # Put species names into vector
  allspecies <- unlist(specieslist)
  # Use single compound class if given 20220418
  if(!is.null(class)) allspecies <- unlist(specieslist[class])

  # Calculate carbon number and oxidation state
  ispecies <- suppressMessages(info(allspecies))
  nC <- sapply(makeup(ispecies), "[", "C")
  ZC <- ZC(ispecies)

  # Load basis species with default activities
  reset()
  basis(c("CO2", "H2", "NH4+", "H2O", "H+"), c(-3, -6, -4, 0, -7))
  # Load formed species
  species(allspecies, -3)
  # Species indices of basis and formed species
  basis_and_species <- c(basis()$ispecies, species()$ispecies)

  # Function to calculate slope of affinity - ZC linear fit
  slopefun <- function(H2, T, sout, plot.it = FALSE) {
    # Set H2 activity
    basis("H2", H2)
    # Calculate affinities
    a <- suppressMessages(affinity(T = T, sout = sout))
    A <- unlist(a$values)
    # Convert to Gibbs energy (J/mol)
    TK <- convert(T, "K")
    A <- -convert(A, "G", TK)
    if(packageVersion("CHNOSZ") <= "1.4.3") A <- convert(A, "J")
    # Convert to kJ/mol
    A <- convert(A, "J") / 1000
    # Divide A by number of carbon
    AC <- A / nC

    # Calculate linear fit
    thislm <- lm(AC ~ ZC)
    # Get slope
    slope <- thislm$coefficients[[2]]
    slope
  }

  # Function to find root at given temperature(s)
  unifun <- function(T = 25) {
    sapply(T, function(thisT) {
      ## Calculate affinity to get subcrt output at this temperature
      #a <- affinity()
      #sout <- a$sout
      # Get subcrt output at this temperature
      sout <- suppressMessages(subcrt(basis_and_species, T = thisT, property = "logK"))
      H2 <- uniroot(slopefun, c(-50, 50), T = thisT, sout = sout)
      H2$root
    })
  }

  # Find root for different combinations of activities of non-H2 basis species
  T = seq(0, 300, 5)
  basis(c("NH4+", "pH"), c(-8, 3))
  loN_lopH <- round(unifun(T = T), 6)
  basis(c("NH4+", "pH"), c(-8, 9))
  loN_hipH <- round(unifun(T = T), 6)
  basis(c("NH4+", "pH"), c(-3, 3))
  hiN_lopH <- round(unifun(T = T), 6)
  basis(c("NH4+", "pH"), c(-3, 9))
  hiN_hipH <- round(unifun(T = T), 6)
  out <- data.frame(T, loN_lopH, loN_hipH, hiN_lopH, hiN_hipH)
  file <- "H2_intermediate.csv"
  if(!is.null(class)) file <- paste0("H2_intermediate_", gsub(" ", "", class), ".csv")
  write.csv(out, file, row.names = FALSE, quote = FALSE)

}

# Intermediate logaH2 calculated for different compound classes 20220418
utegigS1 <- function(pdf = FALSE) {
  # To generate files:
  # lapply(names(specieslist), calc_logaH2_intermediate)
  if(pdf) pdf("Figure_S1.pdf", width = 6, height = 8)
  par(mfrow = c(3, 2))
  parargs <- list(mar = c(4, 4, 3, 1), mgp = c(2.7, 1, 0))
  ylim <- c(-25, 0)
  titlefun <- function(class) title(paste0(class, " (", length(specieslist[[class]]), ")"))
  intermediate_logaH2("C1 and C2",  parargs = parargs, ylim = ylim, redox = FALSE, pH = FALSE, NH4 = FALSE)
  titlefun("C1 and C2")
  intermediate_logaH2("Acid",       parargs = parargs, ylim = ylim, redox = FALSE, pH = FALSE, NH4 = FALSE)
  titlefun("Acid")
  intermediate_logaH2("Amino acid", parargs = parargs, ylim = ylim, redox = FALSE, pH.y = -8.5, NH4.y = -10)
  titlefun("Amino acid")
  intermediate_logaH2("Sugar",      parargs = parargs, ylim = ylim, redox = FALSE, pH = FALSE, NH4 = FALSE)
  titlefun("Sugar")
  intermediate_logaH2("Nucleobase", parargs = parargs, ylim = ylim, redox = FALSE, pH.y = -12, NH4.y = -17)
  titlefun("Nucleobase")
  intermediate_logaH2("TCA cycle",  parargs = parargs, ylim = ylim, redox = FALSE, NH4 = FALSE, pH.y = -10)
  titlefun("TCA cycle")
  if(pdf) dev.off()
}
