# JMDplots/aaaq.R
# plots for the paper by Dick et al. (2006) (doi:10.5194/bg-3-311-2006)

# pKa of sidechain and amino acid backbone groups
aaaq4 <- function(pdf = FALSE) {
  # 20050718 jmd first version for paper (aa/plots/pKa/plot.R)
  # 20060212 make labels prettier; only consider one set of sidechain groups
  # 20060721 change 'pK' label to 'pH', add water
  # 20191020 modified for inclusion in JMDplots package

  # Plot the pKa of sidechain and backbone groups as a function of temperature.
  # Also plot the pH of neutrality of pure water.

  ## set charges of groups in order to activate the g-function
  ## (inactive by default)
  #mod.OBIGT("[AABB+]", z = 1)
  #mod.OBIGT("[AABB-]", z = -1)
  #mod.OBIGT("[Asp-]", z = -1)
  #mod.OBIGT("[Glu-]", z = -1)
  #mod.OBIGT("[His+]", z = 1)
  #mod.OBIGT("[Cys-]", z = -1)
  #mod.OBIGT("[Tyr-]", z = -1)
  #mod.OBIGT("[Lys+]", z = 1)
  #mod.OBIGT("[Arg+]", z = 1)

  # define temperature values
  T <- 0:300

  # start the plot
  if(pdf) pdf(width = 8, height = 6.5, file = "aaaq4.pdf")
  par(cex = 2, mar = c(3, 3, 1, 1), mgp = c(1.8, 0.3, 0), yaxs = "i", xaxs = "i", lwd = 3)
  ylab <- axis.label("T")
  plot(c(0, 13), range(T), xlab = "pH", ylab = ylab, las = 1, tcl = 0, type = "n", xaxt = "n")

  # calculate pKa for amino acid backbone
  pK_AABB_plus <- subcrt(c("[AABB]", "H+", "[AABB+]"), c(-1, -1, 1), T = T)$out$logK
  pK_AABB_minus <- subcrt(c("[AABB-]", "H+", "[AABB]"), c(-1, -1, 1), T = T)$out$logK

  # draw filled area for between AABB ionizations
  polygon(c(pK_AABB_plus, rev(pK_AABB_minus)), c(T, rev(T)), col = "snow3")

  # calculate pKa for sidechains
  pK_AA <- list(
    D = subcrt(c("[Asp-]", "H+", "[Asp]"), c(-1, -1, 1), T = T)$out$logK,
    E = subcrt(c("[Glu-]", "H+", "[Glu]"), c(-1, -1, 1), T = T)$out$logK,
    H = subcrt(c("[His]", "H+", "[His+]"), c(-1, -1, 1), T = T)$out$logK,
    C = subcrt(c("[Cys-]", "H+", "[Cys]"), c(-1, -1, 1), T = T)$out$logK,
    Y = subcrt(c("[Tyr-]", "H+", "[Tyr]"), c(-1, -1, 1), T = T)$out$logK,
    K = subcrt(c("[Lys]", "H+", "[Lys+]"), c(-1, -1, 1), T = T)$out$logK,
    R = subcrt(c("[Arg]", "H+", "[Arg+]"), c(-1, -1, 1), T = T)$out$logK
  )

  for(i in 1:length(pK_AA)) {
    # draw the lines
    lines(pK_AA[[i]], T, lwd = 1)
    # add the labels
    # solid white circle, black border, letter inside
    points(pK_AA[[i]][26], 25, cex = 1.8, pch = 21, bg = "white", lwd = 1)
    text(pK_AA[[i]][26], 25, names(pK_AA)[i], cex = 0.7)
  }

  # add tick marks
  # (we do it here so the ticks show up above the shading)
  axis(1, seq(1, 13, 2), tcl = 0.3, lwd = 2)
  axis(2, labels = FALSE, tcl = 0.3, lwd = 2)
  axis(3, seq(1, 13, 2), labels= FALSE, tcl = 0.3, lwd = 2)
  axis(4, labels = FALSE, tcl = 0.3, lwd = 2)

  # label the AABB fields
  # first, color the edges white so the characters stand out above the bg
  label.aabb <- function(cex, xadj = 0, yadj = 0) {
    text(1.5 + xadj, 275 + yadj, expression(group('[', AABB^{bold('+')}, ']')), cex = cex)
    text(4.5 + xadj, 275 + yadj, expression(group('[', AABB, ']')), cex = cex)
    text(7.6 + xadj, 275 + yadj, expression(group('[', AABB^{bold('-')}, ']')), cex = cex)
  }
  par(fg = 'white')
  for(xadj in seq(-0.15, 0.15, 0.05)) {
    for(yadj in seq(-5, 5, 0.8)) {
       label.aabb(0.8, xadj, yadj)
    }
  }
  par(fg = 'black')
  label.aabb(0.8)

  # add water
  pK_H2O <- subcrt(c("OH-", "H+", "H2O"), c(-1, -1, 1), T = T)$out$logK / 2
  lines(pK_H2O, T, col = 4, lwd = 3)
  text(7, 50, expr.species("H2O"), srt = -65, cex = 0.8)

  ## add a title
  #title("Ionizable Groups in Proteins", font.main = 1, cex.main = 0.8)

  if(pdf) {
    dev.off()
    addexif("aaaq4", "pKa of sidechain and amino acid backbone groups", "https://doi.org/10.5194/bg-3-311-2006")
  }
}

# net charge of proteins as a function of pH and T
# moved from CHNOSZ/ionize.aa.Rd (LYSC_CHICK) and expanded 20191020
aaaq10 <- function(pdf = FALSE) {
  # set up plot
  if(pdf) pdf("aaaq10.pdf", width = 6, height = 5)
  layout(matrix(c(1, 2, 4, 3), nrow = 2))
  par(mar = c(3, 3, 1.5, 1))
  par(mgp = c(2, 0.5, 0))
  par(xaxs = "i", yaxs = "i")
  par(tcl = 0.4, las = 1)
  # list with per-plot metadata
  pldat <- list(
    protein = c("LYSC_CHICK", "RNAS1_BOVIN", "NUC_STAAU"),
    ymin = c(-25, -25, -30),
    ymax = c(25, 25, 35),
    expt = c("RT71", "NT67", "WG00"),
    ref = c("Roxby and Tanford (1971)", "Nozaki and Tanford (1967)", "Whitten and Garc\u00EDa-\nMoreno E. (2000)")
  )

  for(i in 1:3) {
    # start plot
    plot(c(0, 14), c(pldat$ymin[i], pldat$ymax[i]), type = "n", xlab = "pH", ylab = "net charge (Z)")

    # the rownumber of the protein in thermo()$protein
    ip <- pinfo(pldat$protein[i])
    # its amino acid composition
    aa <- pinfo(ip)

    # additive charges of unfolded protein at 25, 100, 150 degrees C, as a function of pH
    pH <- seq(0, 14, 0.1)
    Z.25 <- ionize.aa(aa, T = 25, pH = pH)[, 1]
    Z.100 <- ionize.aa(aa, T = 100, pH = pH)[, 1]
    Z.150 <- ionize.aa(aa, T = 150, pH = pH)[, 1]
    lines(pH, Z.25)
    lines(pH, Z.100, col="red")
    lines(pH, Z.150, col="orange")

    # suppress ionization of cysteine as if it was oxidized to form non-ionizable cystine disulfide bonds
    Z.25.ox <- ionize.aa(aa, T = 25, pH = pH, suppress.Cys = TRUE)[, 1]
    lines(pH, Z.25.ox, lty = 3, lwd = 2)

    # add experimental points
    file <- system.file(paste0("extdata/aaaq/", pldat$expt[i], ".csv"), package = "JMDplots")
    expt <- read.csv(file)
    points(expt$pH, expt$Z)
    legend("topright", pch = 1, legend = pldat$ref[i], bty = "n")

    title(main = pldat$protein[i], font.main = 1)
    label.figure(letters[i], paren = TRUE, italic = TRUE)
  }

  # add figure legend
  plot.new()
  ltxt <- c("25 \u00B0C, oxidized Cys", "25 \u00B0C", "100 \u00B0C", "150 \u00B0C")
  legend("center", ltxt, lty = c(3, 1, 1, 1), col = c(1, 1, 2, "orange"), cex = 1.2, lwd = 2, bty = "n")
  if(pdf) {
    dev.off()
    addexif("aaaq10", "Net charge of proteins as a function of pH and T", "https://doi.org/10.5194/bg-3-311-2006")
  }
}

# Eh-pH diagram for extracellular alpha-amylases
aaaq13 <- function(pdf = FALSE) {
  if(pdf) pdf("aaaq13.pdf", width = 5, height = 4)
  par(mar = c(3, 3, 1, 1))
  # use old data for methionine sidechain
  OldAAfile <- "extdata/OBIGT/OldAA.csv"
  add.OBIGT(system.file(OldAAfile, package = "JMDplots"))
  basis("CHNOSe")
  basis(c("NH3", "H2S"), c(-6, -3))
  species(c("AMY_BACSU", "O08452_PYRFU"))
  a <- affinity(pH = c(2, 10, 500), Eh = c(-1, 1, 500))
  diagram(a, balance = "CO2", lwd = 2, names = "")
  water.lines(a)
  # plot the points
  BKMdat <- read.csv(system.file("extdata/aaaq/BKM60_Fig7.csv", package = "JMDplots"))
  points(BKMdat$pH, BKMdat$Eh, pch = 20)
  # plus signs: Spear et al., 2005
  points(c(8.5, 6.2, 4.2, 6.7), c(0.018, 0.223, 0.022, 0.067), pch = 3, col = "red", lwd = 2)
  # triangles: Stefánsson and Arnórsson, 2002
  points(c(8.24, 7.17, 8.11), c(-0.44, -0.41, -0.49), pch = 17, col = "red")
  # overlay the next plot (higher T)
  a <- affinity(pH = c(2, 10, 500), Eh = c(-1, 1, 500), T = 100)
  d <- diagram(a, balance = "CO2", lwd = 2, col = "red", add = TRUE, names = "")
  water.lines(a, col = "red")
  # make space for names
  for(xadj in seq(-0.05, 0.05, 0.025)) {
    for(yadj in seq(-0.02, 0.02, 0.01)) {
       text(d$namesx[1] + xadj, d$namesy[1] + yadj, "AMY_BACSU", col = "white", cex = 1.5)
       text(d$namesx[2] + xadj, d$namesy[2] + yadj, "O08452_PYRFU", col = "white", cex = 1.5)
    }
  }
  # plot names
  text(d$namesx[1], d$namesy[1], "AMY_BACSU", col = "slateblue1", cex = 1.5)
  text(d$namesx[2], d$namesy[2], "O08452_PYRFU", col = "slateblue1", cex = 1.5)
  # add legends
  ltext <- c(describe.property("T", 25), describe.property("T", 100))
  legend("topright", legend = ltext, lty = 1, col = c("black", "red"), bg = "white")
  ltext <- c("Soils (Bass Becking et al., 1960)", "Yellowstone (Spear et al., 2005)", "Iceland (Stef\u00E1nsson and Arn\u00F3rsson, 2002)")
  legend("bottomleft", legend = ltext, pch = c(20, 3, 17), col = c(1, 2, 2), bty = "n")
  # reset the database to use current parameters for later vignettes
  reset()
  if(pdf) {
    dev.off()
    addexif("aaaq13", "Eh-pH diagram for extracellular alpha-amylases", "https://doi.org/10.5194/bg-3-311-2006")
  }
}
