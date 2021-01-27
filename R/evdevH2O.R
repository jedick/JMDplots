# JMDplots/evdevH2O.R
# Make plots for the paper:
# Redox potential linked to water loss from proteins in evolution and development
# 20201216 First version
# 20210127 Added to JMDplots

# Requires:
# CHNOSZ > 1.4.0 to get species' activities in predominant.values

# Create bold axis labels
ZClab <- quote(bolditalic(Z)[bold(C)])
nH2Olab <- quote(bolditalic(n)[bold(H[2]*O)])
nAAlab <- quote(bolditalic(n)[bold(AA)])
DnH2Olab <- quote(bold(Delta)*bolditalic(n)[bold(H[2]*O)])
DPSlab <- quote(bold(Delta*PS))
DPSHPAlab <- quote(bold(Delta*PS~"(HPA)"))
DPSTCGAlab <- quote(bold(Delta*PS~"(TCGA)"))
logH2Olab <- quote(bold(log)*bolditalic(a)[bold(H[2]*O)])
logO2lab <- quote(bold(log)*bolditalic(f)[bold(O[2])])

# Compositional analysis of Trigos and Liebeskind datasets 20201216
evdevH2O1 <- function(pdf = FALSE) {
  if(pdf) pdf("evdevH2O1.pdf", width = 10, height = 5)
  par(mar = c(4, 4, 1, 1), mgp = c(2.5, 1, 0))
  mat <- matrix(1:8, nrow = 2, byrow = TRUE)
  layout(mat, widths = c(1.2, 1, 1, 1))
  par(cex = 0.8)

  # Make legend for Trigos phylostrata
  opar <- par(mar = c(0, 0, 1.5, 0))
  plot.new()
  par(xpd = NA)
  ys <- seq(0.95, 0.2, length.out = 8)
  text1 <- c("Cellular organisms", "Eukaryota", "Opisthokonta", "Metazoa", "Eumetazoa", "Bilateria", "Chordata", "Euteleostomi")
  text2 <- c("Ammiota", "Mammalia", "Theria", "Eutheria", "Euarchontoglires", "Catarrhini", "Homininae", "")
  text(0.06, ys, 1:8, adj = c(1, 0.5))
  text(0.10, ys, text1, adj = c(0, 0.5))
  text(0.59, ys, 9:16, adj = c(1, 0.5))
  text(0.63, ys, text2, adj = c(0, 0.5))
  text(0.63, ys[8], "Homo sapiens", font = 3, adj = c(0, 0.5))
  title("Trigos phylostrata", font.main = 1)
  par(xpd = FALSE)
  par(opar)

  # Plots 1-3: nAA, ZC, nH2O for Trigos phylostrata
  memo <- plotphylo("nAA", PS_source = "TPPG17")
  label.figure("A", cex = 1.7, yfrac = 0.96, xfrac = 0.05)
  plotphylo("ZC", PS_source = "TPPG17", memo = memo)
  plotphylo("nH2O", PS_source = "TPPG17", memo = memo)

  # Make legend for Liebeskind phylostrata
  opar <- par(mar = c(0, 0, 1.5, 0))
  plot.new()
  par(xpd = NA)
  ys <- seq(0.95, 0.2, length.out = 8)
  text1 <- c("Cellular_organisms", "Euk_Archaea", "Euk+Bac", "Eukaryota", "Opisthokonta", "Eumetazoa", "Vertebrata", "Mammalia")
  text(0.29, ys, 1:8, adj = c(1, 0.5))
  text(0.33, ys, text1, adj = c(0, 0.5))
  title("Liebeskind gene ages", font.main = 1)
  par(xpd = FALSE)
  par(opar)

  # Plots 4-6: nAA, ZC, nH2O for Liebeskind phylostrata
  memo <- plotphylo("nAA", PS_source = "LMM16", xlab = "GA")
  label.figure("B", cex = 1.7, yfrac = 0.96, xfrac = 0.05)
  plotphylo("ZC", PS_source = "LMM16", memo = memo, xlab = "GA")
  plotphylo("nH2O", PS_source = "LMM16", memo = memo, xlab = "GA")

  if(pdf) {
    dev.off()
    addexif("evdevH2O1", "Compositional analysis of Trigos and Liebeskind datasets", "Dick (2021) (preprint)")
  }
}

# Thermodynamic analysis of optimal logaH2O and logfO2 for target proteins 20201219
evdevH2O2 <- function(pdf = FALSE) {

  if(pdf) pdf("evdevH2O2.pdf", width = 7, height = 3)
  par(mar = c(3, 3.1, 3, 1), mgp = c(2, 0.5, 0))

  # Get mean amino acid compositions
  gpa <- getphyloaa("TPPG17")
  PS <- gpa$aa$protein
  # Load PS model proteins 
  ipPS <- add.protein(gpa$aa)
  # Set up system
  basis("QEC")
  O2 <- c(-72, -67)
  H2O <- c(-5, 2)

  # Plot A: predominance diagram for all PS
  a <- affinity(O2 = c(O2, 128), H2O = c(H2O, 128), iprotein = ipPS)
  e <- equilibrate(a, as.residue = TRUE, loga.balance = 0)
  d <- diagram(e, plot.it = FALSE)
  # Make color image for activities
  par.orig <- my.filled.contour(e$vals$O2, e$vals$H2O, d$predominant.values, xlab = logO2lab, ylab = logH2Olab,
    nlevels = 50, col = hcl.colors(75, "YlGnBu")[20:75], frame.plot = FALSE,
    # use plot.axes to label the contour plot (see ?filled.contour)
    plot.axes = {
      box()
      #diagram(e, add = TRUE)
      # Use a higher resolution for making the lines
      a <- affinity(O2 = c(O2, 400), H2O = c(H2O, 400), iprotein = ipPS)
      diagram(a, as.residue = TRUE, add = TRUE)
      opar <- par(tcl = 0.3)
      thermo.axis()
      axis(1)
      axis(2)
      par(opar)
      # Show location of maximum activity for each PS model protein
      par(xpd = TRUE)
      optO2 <- optH2O <- numeric()
      for(i in 1:16) {
        imax <- arrayInd(which.max(e$loga.equil[[i]]), dim(e$loga.equil[[i]]))
        optO2 <- c(optO2, e$vals$O2[imax[1]])
        optH2O <- c(optH2O, e$vals$H2O[imax[2]])
      }
      set.seed(3)
      points(jitter(optO2, 0.4), jitter(optH2O, 0.05), pch = 21, bg = 7, cex = 1.5)
      par(xpd = FALSE)
      title("16 PS model proteins", font.main = 1)
      label.figure("A", cex = 1.5, yfrac = 0.937)
    },
    key.axes = {
      opar <- par(tcl = 0)
      axis(4, at = par("usr")[3:4], labels = round(par("usr")[3:4], 2))
      title(quote(bold(log)*bolditalic(a)[bold(protein)]), cex.main = 1, line = 1)
      par(opar)
    },
    add2 = FALSE
  )

  # Plot B: 16 PS model proteins and n = 200 sample of human proteins
  set.seed(3)
  iind <- sample(1:nrow(gpa$pcomp$aa), 200)
  ipind <- add.protein(gpa$pcomp$aa[iind, ])
  a <- affinity(O2 = c(O2, 128), H2O = c(H2O, 128), iprotein = c(ipPS, ipind))
  e <- equilibrate(a, as.residue = TRUE, loga.balance = 0)
  d <- diagram(e, plot.it = FALSE)
  my.filled.contour(e$vals$O2, e$vals$H2O, d$predominant.values, xlab = logO2lab, ylab = logH2Olab,
    nlevels = 50,
    col = hcl.colors(75, "YlGnBu")[20:75],
    # use plot.axes to label the contour plot (see ?filled.contour)
    plot.axes = {
      names <- sapply(strsplit(d$species$name, "\\|"), "[", 2)
      #diagram(e, add = TRUE, names = names, format.names = FALSE)
      # Use a higher resolution for making the lines
      a <- affinity(O2 = c(O2, 400), H2O = c(H2O, 400), iprotein = c(ipPS, ipind))
      diagram(a, as.residue = TRUE, add = TRUE, names = names, format.names = FALSE)
      opar <- par(tcl = 0.3)
      thermo.axis()
      axis(1)
      axis(2)
      par(opar)
      # Show location of maximum activity for each PS model protein
      par(xpd = TRUE)
      for(i in 1:16) {
        imax <- arrayInd(which.max(e$loga.equil[[i]]), dim(e$loga.equil[[i]]))
        optO2 <- e$vals$O2[imax[1]]
        optH2O <- e$vals$H2O[imax[2]]
        points(optO2, optH2O, pch = 21, bg = 7, cex = 1.5)
      }
      par(xpd = FALSE)
      title("16 PS model + 200 human proteins", font.main = 1)
      label.figure("B", cex = 1.5, yfrac = 0.937)
    },
    key.axes = {
      opar <- par(tcl = 0)
      axis(4, at = par("usr")[3:4], labels = round(par("usr")[3:4], 2))
      title(quote(bold(log)*bolditalic(a)[bold(protein)]), cex.main = 1, line = 1)
      par(opar)
    },
    add2 = TRUE
  )

  # Reset plot parameters
  par(par.orig)

  if(pdf) {
    dev.off()
    addexif("evdevH2O2", "Thermodynamic analysis of optimal logaH2O and logfO2 for target proteins", "Dick (2021) (preprint)")
  }

}

# Optimal logaH2O and logfO2 and effective Eh for target proteins 20201218
evdevH2O3 <- function(pdf = FALSE) {

  # Setup plot
  if(pdf) pdf("evdevH2O3.pdf", width = 8, height = 5)
  par(mgp = c(2.5, 1, 0), las = 1, font.lab = 2)
  layout(matrix(1:6, nrow = 2))
  par(mar = c(3.8, 4.5, 1.2, 1))
  # Set basis species to make axis labels (logaH2O, logfO2)
  basis("QEC")

  # Function to make plots
  plotfun <- function(PS_source, fig.lab) {

    # Read results
    datadir <- system.file("extdata/evdevH2O", package = "JMDplots")
    H2Ofile <- file.path(datadir, paste0(PS_source, "_H2O.csv"))
    O2file <- file.path(datadir, paste0(PS_source, "_O2.csv"))
    H2O <- read.csv(H2Ofile, as.is = TRUE, check.names = FALSE)
    O2 <- read.csv(O2file, as.is = TRUE, check.names = FALSE)
    # Get phylostrata
    iPS <- 2:ncol(H2O)
    PS <- as.numeric(colnames(H2O)[iPS])
    # Get mean values of logaH2O and logfO2
    meanH2O <- colMeans(H2O[, iPS])
    meanO2 <- colMeans(O2[, iPS])

    if(!is.na(fig.lab[1])) {
      # Make logaH2O plot
      plot(range(PS), c(-2.5, 2.5), xlab = NA, ylab = logH2Olab, type = "n", xaxt = "n", xaxs = "i", yaxs = "i")
      mtext("PS", 1, 2.2, font = 2, cex = par("cex"))
      axis(1, at = 1:16, labels = c(1,2,NA,NA,5,NA,NA,NA,NA,10,NA,NA,NA,NA,NA,16))
      # FIXME: redraw "2" because previous command doesn't plot it (spacing too tight)
      axis(1, at = 2, labels = 2)
      for(i in 1:nrow(H2O)) lines(PS, H2O[i, iPS], lwd = 0.5, col = "gray")
      lines(PS, meanH2O, lwd = 2, col = 2) 
      label.figure(fig.lab[1], cex = 1.6)
    }

    if(!is.na(fig.lab[2])) {
      # Make logfO2 plot
      plot(range(PS), c(-72, -67), xlab = NA, ylab = logO2lab, type = "n", xaxt = "n", xaxs = "i", yaxs = "i")
      mtext("PS", 1, 2.2, font = 2, cex = par("cex"))
      axis(1, at = 1:16, labels = c(1,2,NA,NA,5,NA,NA,NA,NA,10,NA,NA,NA,NA,NA,16))
      axis(1, at = 2, labels = 2)
      for(i in 1:nrow(O2)) lines(PS, O2[i, iPS], lwd = 0.5, col = "gray")
      lines(PS, meanO2, lwd = 2, col = 2) 
      label.figure(fig.lab[2], cex = 1.6)
    }

    if(!is.na(fig.lab[3])) {
      # Make logaH2O plot (mean and logaH2O = 0)
      plot(range(PS), c(-2.5, 2.5), xlab = NA, ylab = logH2Olab, type = "n", xaxt = "n", xaxs = "i", yaxs = "i")
      xlab <- switch(PS_source, TPPG17 = "PS", LMM16 = "GA")
      mtext(xlab, 1, 2.2, font = 2, cex = par("cex"))
      if(PS_source == "TPPG17") {
        axis(1, at = 1:16, labels = c(1,2,NA,NA,5,NA,NA,NA,NA,10,NA,NA,NA,NA,NA,16))
        axis(1, at = 2, labels = 2)
        # Shaded area starting at PS 10 (Mammalia)
        rect(10, par("usr")[3], 16, par("usr")[4], col = "lightgray")
        # PS 1 (Cellular organisms), 2 (Eukaryota), 5 (Eumetazoa), 10 (Mammalia)
        abline(v = c(1.04, 2, 5, 10), col = 5, lwd = 2)
      }
      if(PS_source == "LMM16") {
        axis(1, at = 1:8, labels = c(1,NA,NA,4,NA,6,NA,8))
        # PS 1 (Cellular_organisms), 4 (Eukaryota), 6 (Eumetazoa), 8 (Mammalia)
        abline(v = c(1.02, 4, 6, 7.98), col = 5, lwd = 2)
      }
      abline(h = 0, lty = 2, lwd = 1.5, col = "slategray4")
      lines(PS, meanH2O, lwd = 2, col = 2) 
      label.figure(fig.lab[3], cex = 1.6)
    }

    if(!is.na(fig.lab[4])) {
      # Calculate Eh
      # H2O = 0.5 O2 + 2 H+ + 2 e- 
      # logK = 0.5 logfO2 - 2 pH - 2 pe - logaH2O
      # pe = 0.25 logfO2 - pH - 0.5 logaH2O - 0.5 logK
      logK <- subcrt(c("H2O", "oxygen", "H+", "e-"), c(-1, 0.5, 2, 2), T = 25)$out$logK
      pH <- 7
      pe <- 0.25 * meanO2 - pH - 0.5 * meanH2O - 0.5 * logK
      Eh <- convert(pe, "Eh")
      # Another way to calculate using CHNOSZ::convert
      Eh2 <- convert(meanO2, "E0", pH = pH, logaH2O = meanH2O)
      stopifnot(all.equal(Eh, Eh2))
      mV <- Eh * 1000

      # Make Eh plot
      plot(range(PS), c(-325, -150), xlab = NA, ylab = NA, type = "n", xaxt = "n", xaxs = "i", yaxt = "n")
      xlab <- switch(PS_source, TPPG17 = "PS", LMM16 = "GA")
      mtext(xlab, 1, 2.2, font = 2, cex = par("cex"))
      axis(2, at = seq(-300, -150, 50))
      mtext("Eh (mV)", side = 2, line = 3, las = 0, cex = par("cex"), font = 2)
      if(PS_source == "TPPG17") {
        axis(1, at = 1:16, labels = c(1,2,NA,NA,5,NA,NA,NA,NA,10,NA,NA,NA,NA,NA,16))
        # FIXME: redraw "2" because previous command doesn't plot it (spacing too tight)
        axis(1, at = 2, labels = 2)
        # Shaded area starting at PS 10 (Mammalia)
        rect(10, par("usr")[3], 16, par("usr")[4], col = "lightgray")
        # PS 1 (Cellular organisms), 2 (Eukaryota), 5 (Eumetazoa), 10 (Mammalia)
        abline(v = c(1.04, 2, 5, 10), col = 5, lwd = 2)
        xtext <- 6.2
      }
      if(PS_source == "LMM16") {
        axis(1, at = 1:8, labels = c(1,NA,NA,4,NA,6,NA,8))
        # PS 1 (Cellular_organisms), 4 (Eukaryota), 6 (Eumetazoa), 8 (Mammalia)
        abline(v = c(1.02, 4, 6, 7.98), col = 5, lwd = 2)
        xtext <- 5.5
      }
      # Eh = -150 mV (plasma GSH/GSSG) Jones and Sies, 2015
      # Eh = -199 mV (erythrocyte GSH/GSSG) van 't Erve et al., 2013
      # Eh = -241 mV (cytosolic NADH/NAD+) Jones and Sies, 2015
      # Eh = -318 mV (mitochondrial NADH/NAD+) Jones and Sies, 2015
      abline(h = c(-150, -199, -241, -318), lty = 2, lwd = 1.5, col = "slategray4")
      lines(PS, mV, lwd = 2, col = 2) 
      text(xtext, -150, "Plasma GSH/GSSG", adj = c(0.5, 1.3))
      text(xtext, -199, "Erythrocyte GSH/GSSG", adj = c(0.5, -0.3))
      text(xtext, -241, "Cytosolic NADH/NAD+", adj = c(0.5, 1.3))
      text(xtext, -318, "Mitochondrial NADH/NAD+", adj = c(0.5, -0.3))
      label.figure(fig.lab[4], cex = 1.6)
    }

  }

  plotfun("TPPG17", c("A", "B", "C", "D"))
  plotfun("LMM16", c(NA, NA, "E", "F"))

  if(pdf) {
    dev.off()
    addexif("evdevH2O3", "Optimal logaH2O and logfO2 and effective Eh for target proteins", "Dick (2021) (preprint)")
  }
}

# Compositional and thermodynamic analysis of B. subtilis biofilm transcriptome and proteome 20201221
evdevH2O4 <- function(pdf = FALSE) {

  # Setup plot
  if(pdf) pdf("evdevH2O4.pdf", width = 7, height = 4.5)
  par(mfrow = c(2, 3))
  par(mar = c(4, 4, 1, 1), las = 1, mgp = c(3, 0.8, 0))

  # Read the overall amino acid compositions calculated from Futo et al., 2020 data
  aa <- read.csv("FOK+20_overall.csv", as.is = TRUE)
  # Identify rows with transcriptome and proteome data
  isT <- aa$organism == "transcriptome"
  isP <- aa$organism == "proteome"
  # Identify stages with proteomic data
  iP <- match(aa$protein[isP], aa$protein[isT])

  # Plot A: protein length
  pl <- protein.length(aa)
  plot(c(1, 11), range(pl), type = "n", xlab = NA, ylab = "Protein length", xaxt = "n", font.lab = 2)
  # Make rotated labels (modified from https://www.r-bloggers.com/rotated-axis-labels-in-r-plots/)
  text(x = 1:11, y = par()$usr[3] - 1.5 * strheight("A"), labels = aa$protein[isT], srt = 45, adj = 1, xpd = TRUE)
  axis(1, at = 1:11, labels = NA)
  abline(v = c(5, 9), lty = 3, lwd = 1.5, col = "gray40")
  lines(1:11, pl[isT], type = "b", lty = 2, col = 3, pch = 19, cex = 1.3)
  lines(iP, pl[isP], type = "b", col = 4, pch = 15, cex = 1.3)
  legend("bottomleft", c("Proteome", "Transcriptome"), lty = c(1, 2), pch = c(15, 19), col = c(4, 3), pt.cex = 1.3, bty = "n")
  label.figure("A", cex = 1.5, xfrac = 0.025)

  # Plot B: nH2O
  nH2O <- H2OAA(aa)
  plot(c(1, 11), range(nH2O), type = "n", xlab = NA, ylab = nH2Olab, xaxt = "n", font.lab = 2)
  text(x = 1:11, y = par()$usr[3] - 1.5 * strheight("A"), labels = aa$protein[isT], srt = 45, adj = 1, xpd = TRUE)
  axis(1, at = 1:11, labels = NA)
  abline(v = c(5, 9), lty = 3, lwd = 1.5, col = "gray40")
  lines(1:11, nH2O[isT], type = "b", lty = 2, col = 3, pch = 19, cex = 1.3)
  lines(iP, nH2O[isP], type = "b", col = 4, pch = 15, cex = 1.3)
  label.figure("B", cex = 1.5, xfrac = 0.025)

  # Plot C: ZC
  ZC <- ZCAA(aa)
  plot(c(1, 11), range(ZC), type = "n", xlab = NA, ylab = ZClab, xaxt = "n", font.lab = 2)
  text(x = 1:11, y = par()$usr[3] - 1.5 * strheight("A"), labels = aa$protein[isT], srt = 45, adj = 1, xpd = TRUE)
  axis(1, at = 1:11, labels = NA)
  abline(v = c(5, 9), lty = 3, lwd = 1.5, col = "gray40")
  lines(1:11, ZC[isT], type = "b", lty = 2, col = 3, pch = 19, cex = 1.3)
  lines(iP, ZC[isP], type = "b", col = 4, pch = 15, cex = 1.3)
  label.figure("C", cex = 1.5, xfrac = 0.025)

  # Plot D: logaH2O
  par(mgp = c(2.8, 0.8, 0))
  T_H2O <- read.csv("transcriptome_H2O.csv", as.is = TRUE)[, -1]
  P_H2O <- read.csv("proteome_H2O.csv", as.is = TRUE)[, -1]
  # Get mean values of logaH2O
  T_meanH2O <- colMeans(T_H2O)
  P_meanH2O <- colMeans(P_H2O)
  plot(c(1, 11), range(c(T_meanH2O, P_meanH2O)), type = "n", xlab = "Biofilm stage", ylab = logH2Olab, xaxt = "n", font.lab = 2)
  text(x = 1:11, y = par()$usr[3] - 1.5 * strheight("A"), labels = aa$protein[isT], srt = 45, adj = 1, xpd = TRUE)
  axis(1, at = 1:11, labels = NA)
  abline(v = c(5, 9), lty = 3, lwd = 1.5, col = "gray40")
  abline(h = 0, lty = 2, lwd = 1.5, col = "slategray4")
  lines(1:11, T_meanH2O, type = "b", lty = 2, col = 3, pch = 19, cex = 1.3)
  lines(iP, P_meanH2O, type = "b", col = 4, pch = 15, cex = 1.3)
  label.figure("D", cex = 1.5, xfrac = 0.025)

  # Plot E: logfO2
  T_O2 <- read.csv("transcriptome_O2.csv", as.is = TRUE)[, -1]
  P_O2 <- read.csv("proteome_O2.csv", as.is = TRUE)[, -1]
  # Get mean values of logfO2
  T_meanO2 <- colMeans(T_O2)
  P_meanO2 <- colMeans(P_O2)
  plot(c(1, 11), range(c(T_meanO2, P_meanO2)), type = "n", xlab = "Biofilm stage", ylab = logO2lab, xaxt = "n", font.lab = 2)
  text(x = 1:11, y = par()$usr[3] - 1.5 * strheight("A"), labels = aa$protein[isT], srt = 45, adj = 1, xpd = TRUE)
  axis(1, at = 1:11, labels = NA)
  abline(v = c(5, 9), lty = 3, lwd = 1.5, col = "gray40")
  abline(h = 0, lty = 2, lwd = 1.5, col = "slategray4")
  lines(1:11, T_meanO2, type = "b", lty = 2, col = 3, pch = 19, cex = 1.3)
  lines(iP, P_meanO2, type = "b", col = 4, pch = 15, cex = 1.3)
  label.figure("E", cex = 1.5, xfrac = 0.025)

  # Plot F: Eh
  logK <- subcrt(c("H2O", "oxygen", "H+", "e-"), c(-1, 0.5, 2, 2), T = 25)$out$logK
  pH <- 7
  T_pe <- 0.25 * T_meanO2 - pH - 0.5 * T_meanH2O - 0.5 * logK
  P_pe <- 0.25 * P_meanO2 - pH - 0.5 * P_meanH2O - 0.5 * logK
  T_Eh <- convert(T_pe, "Eh") * 1000
  P_Eh <- convert(P_pe, "Eh") * 1000
  plot(c(1, 11), range(c(T_Eh, P_Eh)), type = "n", xlab = "Biofilm stage", ylab = "Eh (mV)", xaxt = "n", font.lab = 2)
  text(x = 1:11, y = par()$usr[3] - 1.5 * strheight("A"), labels = aa$protein[isT], srt = 45, adj = 1, xpd = TRUE)
  axis(1, at = 1:11, labels = NA)
  abline(v = c(5, 9), lty = 3, lwd = 1.5, col = "gray40")
  lines(1:11, T_Eh, type = "b", lty = 2, col = 3, pch = 19, cex = 1.3)
  lines(iP, P_Eh, type = "b", col = 4, pch = 15, cex = 1.3)
  label.figure("F", cex = 1.5, xfrac = 0.025)

  if(pdf) {
    dev.off()
    addexif("evdevH2O4", "Compositional and thermodynamic analysis of B. subtilis biofilm transcriptome and proteome", "Dick (2021) (preprint)")
  }
}

# Organismal water content and proteomic nH2O for fruit fly development 20210116
evdevH2O5 <- function(pdf = FALSE) {

  # Setup plot
  if(pdf) pdf("evdevH2O5.pdf", width = 12, height = 5)
  par(mfrow = c(1, 2), font.lab = 2, las = 1)

  # Developmental time points and water content values from Fig. 4 of Church and Robertson, 1966
  # Hours 110, 120, 130 are for P1, P2, Adult
  Hours <- c(13.47, 17.4, 23.6, 30.5, 36.17, 48.05, 53.82, 60.02, 72.64, 84.06, 95.17, 103.95, 110, 120, 130)
  Water <- c(79.57, 75.23, 75.33, 81.69, 82.4, 82.81, 76.81, 76.77, 79.58, 76.19, 76.1, 76.01, 66.72, 66.38, 72.04)
  plot(Hours, Water, type = "b", pch = 19, xlim = c(0, 140), xlab = "Hours after hatching", ylab = "Water content (%)",
       xaxs = "i", xaxt = "n", cex = 1.5)
  axis(1, at = seq(0, 100, 20))
  axis(1, at = 110, "P1")
  axis(1, at = 120, "P2")
  axis(1, at = 130, "Adult")
  title("Organismal water content\n(Church and Robertson, 1966)", font.main = 1)
  label.figure("A", cex = 1.6)

  # Read mean amino acid compositions for developmental time points
  aa <- read.csv("CBS+17_mean_aa.csv", as.is = TRUE)
  # Plot nH2O
  plot(H2OAA(aa), type = "b", xaxt = "n", xlab = "Developmental stage", ylab = nH2Olab, cex = 1.5)
  labels <- gsub("p", "P", gsub("e", "E", aa$protein))
  # Make rotated labels (modified from https://www.r-bloggers.com/rotated-axis-labels-in-r-plots/)
  text(x = 1:17, y = par()$usr[3] - 1.5 * strheight("A"), labels = labels, srt = 45, adj = 1, xpd = TRUE)
  axis(1, at = 1:17, labels = NA)
  title("Stoichiometric hydration state of proteins\n(This study)", font.main = 1)
  label.figure("B", cex = 1.6)

  if(pdf) {
    dev.off()
    addexif("evdevH2O5", "Organismal water content and proteomic nH2O for fruit fly development", "Dick (2021) (preprint)")
  }
}

# Calculate optimal logaH2O and logfO2 for phylostrata 20201218
# Make it work for B. subtilis biofilm dataset (Futo et al., 2020) 20201221
optimal_activity <- function(dataset = "TPPG17", seed = 1:100) {

  # Process 'dataset' argument
  if(dataset %in% c("TPPG17", "LMM16")) {
    if(dataset == "TPPG17") xlab <- "Trigos phylostrata"
    if(dataset == "LMM16") xlab <- "Liebeskind gene ages"
    # Get mean amino acid compositions and corresponding phylostrata (e.g. 1..16)
    gpa <- getphyloaa(dataset)
    names <- gpa$aa$protein
    # Load PS model proteins 
    iptarget <- add.protein(gpa$aa)
    O2 <- c(-72, -67)
    H2O <- c(-2, 3)
  } else if(dataset %in% c("transcriptome", "proteome")) {
    xlab <- paste("Biofilm", dataset)
    # Read amino acid compositions of overall proteins in each biofilm stage
    aa <- read.csv("FOK+20_overall.csv")
    aa <- aa[aa$organism == dataset, ]
    names <- aa$protein
    # Load overall proteins
    iptarget <- add.protein(aa)
    # We need this to get all the human proteins 20201221
    gpa <- getphyloaa("TPPG17")
    O2 <- c(-72, -65)
    H2O <- c(-2, 5)
  }

  # Set up system
  basis("QEC")
  # Initialize output values and plot
  outO2 <- outH2O <- list()
  split.screen(c(2, 1))
  screen(1)
  par(mar = c(4, 4, 1, 1), mgp = c(2.5, 1, 0), las = 1)
  plot(range(1:length(names)), O2, xlab = xlab, ylab = logO2lab, type = "n", xaxt = "n", xaxs = "i", yaxs = "i", font.lab = 2)
  axis(1, 1:length(names), names)
  screen(2)
  par(mar = c(4, 4, 1, 1), mgp = c(2.5, 1, 0), las = 1)
  plot(range(1:length(names)), H2O, xlab = xlab, ylab = logH2Olab, type = "n", xaxt = "n", xaxs = "i", yaxs = "i", font.lab = 2)
  axis(1, 1:length(names), names)

  # Prevent multicore usage (runs out of memory)
  thermo("opt$paramin" = 10000)
  # Loop over random seeds
  for(iseed in seq_along(seed)) {
    # Calculate affinities for target proteins and a sample of human proteins (background)
    set.seed(seed[iseed])
    iback <- sample(1:nrow(gpa$pcomp$aa), 2000)
    ipback <- add.protein(gpa$pcomp$aa[iback, ])
    a <- affinity(O2 = O2, H2O = H2O, iprotein = c(iptarget, ipback))
    # Equilibrate and find maximum activity for each target protein
    e <- equilibrate(a, as.residue = TRUE, loga.balance = 0)
    optO2 <- optH2O <- numeric()
    for(i in seq_along(names)) {
      imax <- arrayInd(which.max(e$loga.equil[[i]]), dim(e$loga.equil[[i]]))
      optO2 <- c(optO2, e$vals$O2[imax[1]])
      optH2O <- c(optH2O, e$vals$H2O[imax[2]])
    }
    # Plot logfO2 and logaH2O values
    # FIXME: need two screen() calls to make this work 20201218
    screen(1, FALSE); screen(1, FALSE)
    lines(1:length(names), optO2, lwd = 0.5, col = "gray")
    screen(2, FALSE); screen(2, FALSE)
    lines(1:length(names), optH2O, lwd = 0.5, col = "gray")
    # Store results
    outO2[[iseed]] <- optO2
    outH2O[[iseed]] <- optH2O
  }
  # Plot mean values
  outO2 <- round(do.call(rbind, outO2), 3)
  outH2O <- round(do.call(rbind, outH2O), 3)
  meanO2 <- colMeans(outO2)
  meanH2O <- colMeans(outH2O)
  screen(1, FALSE); screen(1, FALSE)
  lines(1:length(names), meanO2, col = 2, lwd = 2)
  screen(2, FALSE); screen(2, FALSE)
  lines(1:length(names), meanH2O, col = 2, lwd = 2)
  close.screen(all.screens = TRUE)
  # Save results
  outO2 <- data.frame(outO2)
  colnames(outO2) <- names
  outO2 <- cbind(seed = seed, outO2)
  outH2O <- data.frame(outH2O)
  colnames(outH2O) <- names
  outH2O <- cbind(seed = seed, outH2O)
  write.csv(outO2, paste0(dataset, "_O2.csv"), row.names = FALSE, quote = FALSE)
  write.csv(outH2O, paste0(dataset, "_H2O.csv"), row.names = FALSE, quote = FALSE)
  savePlot(paste0(gsub(" ", "_", xlab), ".png"))
}


############################
### UNEXPORTED FUNCTIONS ###
############################

# Mean ZC and nH2O of phylostrata 20191122
plotphylo <- function(vars = c("ZC", "nH2O"), PS_source = "TPPG17", memo = NULL, xlab = "PS") {
  if(is.null(memo)) {
    dat <- read.csv(system.file("extdata/phylostrata/TPPG17.csv.xz", package = "canprot"), as.is = TRUE)
    if(PS_source == "LMM16") {
      dat <- read.csv(system.file("extdata/phylostrata/LMM16.csv.xz", package = "canprot"), as.is = TRUE)
      colnames(dat)[c(1,3)] <- c("Entry", "Phylostrata")
      # remove entries that have ENSP instead of UniProt IDs
      dat <- dat[!grepl("^ENSP", dat$Entry), ]
    }
    dat <- check_IDs(dat, "Entry")
    dat <- cleanup(dat, "Entry")
    pcomp <- protcomp(dat$Entry)
  } else {
    dat <- memo$dat
    pcomp <- memo$pcomp
  }
  # Use metric functions (H2OAA, ZCAA) because protcomp no longer returns these values 20201216
  nH2O <- H2OAA(pcomp$aa)
  ZC <- ZCAA(pcomp$aa)
  nAA <- protein.length(pcomp$aa)
  # get mean ZC and nH2O for each phylostratum
  PS <- sort(unique(dat$Phylostrata))
  cum.ZC <- cum.nH2O <- cum.nAA <- mean.ZC <- mean.nH2O <- mean.nAA <- numeric()
  for(p in PS) {
    # point mean
    mean.ZC <- c(mean.ZC, mean(ZC[dat$Phylostrata == p]))
    mean.nH2O <- c(mean.nH2O, mean(nH2O[dat$Phylostrata == p]))
    mean.nAA <- c(mean.nAA, mean(nAA[dat$Phylostrata == p]))
    # cumulative mean
    cum.ZC <- c(cum.ZC, mean(ZC[dat$Phylostrata <= p]))
    cum.nH2O <- c(cum.nH2O, mean(nH2O[dat$Phylostrata <= p]))
    cum.nAA <- c(cum.nAA, mean(nAA[dat$Phylostrata <= p]))
  }
  if("ZC" %in% vars) {
    plot(PS, mean.ZC, type = "b", xlab = xlab, ylab = ZClab, font.lab = 2)
    lines(PS, cum.ZC, col = 2, lty = 2)
  }
  if("nH2O" %in% vars) {
    plot(PS, mean.nH2O, type = "b", xlab = xlab, ylab = nH2Olab, font.lab = 2)
    lines(PS, cum.nH2O, col = 2, lty = 2)
  }
  if("nAA" %in% vars) {
    plot(PS, mean.nAA, type = "b", xlab = xlab, ylab = nAAlab, font.lab = 2)
    lines(PS, cum.nAA, col = 2, lty = 2)
  }
  # return the dat and pcomp for memoization 20191211
  invisible(list(dat = dat, pcomp = pcomp))
}

# Get mean amino acid composition for each phylostratum 20201219
getphyloaa <- function(PS_source) {
  dat <- read.csv(system.file(paste0("extdata/phylostrata/", PS_source, ".csv.xz"), package = "canprot"), as.is = TRUE)
  if(PS_source == "LMM16") {
    colnames(dat)[c(1,3)] <- c("Entry", "Phylostrata")
    # remove entries that have ENSP instead of UniProt IDs
    dat <- dat[!grepl("^ENSP", dat$Entry), ]
  }
  dat <- check_IDs(dat, "Entry")
  dat <- cleanup(dat, "Entry")
  pcomp <- protcomp(dat$Entry)
  # Set up blank amino acid data frame
  PS <- sort(unique(dat$Phylostrata))
  aa <- thermo()$protein[rep(1, length(PS)), ]
  aa$protein <- PS
  aa$organism <- PS_source
  aa$ref <- aa$abbrv <- NA
  aa$chains <- 1
  aa[, 6:25] <- 0
  # Loop over phylostrata
  for(i in seq_along(PS)) {
    iPS <- dat$Phylostrata == PS[i]
    aaPS <- pcomp$aa[iPS, ]
    aamean <- colMeans(aaPS[, 6:25])
    aa[i, 6:25] <- aamean
  }
  # Return both the mean compositions (aa) and all proteins (pcomp)
  list(aa = aa, pcomp = pcomp)
}


### Modification of filled.contour.R by Jeffrey M. Dick 20201219
### - Add border = NA to rect()
### - Add 'add2' argument (FALSE to initialize 2x2 plot layout; TRUE to add to 2x2 layout)
###   - With non-NULL 'add2', don't do on.exit(par(par.orig)) but return par.orig invisibly

#  File src/library/graphics/R/filled.contour.R
#  Part of the R package, https://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  https://www.R-project.org/Licenses/

my.filled.contour <-
function (x = seq(0, 1, length.out = nrow(z)),
          y = seq(0, 1, length.out = ncol(z)),
          z,
          xlim = range(x, finite=TRUE),
          ylim = range(y, finite=TRUE),
          zlim = range(z, finite=TRUE),
          levels = pretty(zlim, nlevels), nlevels = 20,
          color.palette = function(n) hcl.colors(n, "YlOrRd", rev = TRUE),
          col = color.palette(length(levels) - 1),
          plot.title, plot.axes, key.title, key.axes,
          asp = NA, xaxs = "i", yaxs = "i", las = 1, axes = TRUE,
          frame.plot = axes, add2 = NULL, ...)
{
    if (missing(z)) {
        if (!missing(x)) {
            if (is.list(x)) {
                z <- x$z
                y <- x$y
                x <- x$x
            }
            else {
                z <- x
                x <- seq.int(0, 1, length.out = nrow(z))
            }
        }
        else stop("no 'z' matrix specified")
    }
    else if (is.list(x)) {
        y <- x$y
        x <- x$x
    }
    if (any(diff(x) <= 0) || any(diff(y) <= 0))
        stop("increasing 'x' and 'y' values expected")

    mar.orig <- (par.orig <- par(c("mar","las","mfrow")))$mar
    if(is.null(add2)) on.exit(par(par.orig))

    w <- (3 + mar.orig[2L]) * par("csi") * 2.54
    if(is.null(add2)) layout(matrix(c(2, 1), ncol = 2L), widths = c(1, lcm(w)))
    else if(!add2) layout(matrix(c(2,1,4,3), ncol = 4L), widths = c(1, lcm(w/1.5), 1, lcm(w/1.5)))
    par(las = las)

    ## Plot the 'plot key' (scale):
    mar <- mar.orig
    mar[4L] <- mar[2L]
    mar[2L] <- 1
    par(mar = mar)
    plot.new()
    plot.window(xlim = c(0,1), ylim = range(levels), xaxs = "i", yaxs = "i")
    rect(0, levels[-length(levels)], 1, levels[-1L], col = col, border = NA)
    if (missing(key.axes)) {
        if (axes)
            axis(4)
    }
    else key.axes
    box()
    if (!missing(key.title))
	key.title

    ## Plot contour-image::
    mar <- mar.orig
    mar[4L] <- 1
    par(mar = mar)
    plot.new()
    plot.window(xlim, ylim, "", xaxs = xaxs, yaxs = yaxs, asp = asp)

    .filled.contour(x, y, z, levels, col)
    if (missing(plot.axes)) {
        if (axes) {
            title(main = "", xlab = "", ylab = "")
            Axis(x, side = 1)
            Axis(y, side = 2)
        }
    }
    else plot.axes
    if (frame.plot) box()
    if (missing(plot.title))
        title(...)
    else
	plot.title
    invisible(par.orig)
}
