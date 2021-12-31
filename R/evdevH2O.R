# JMDplots/evdevH2O.R
# Make plots for the paper:
# Redox potential linked to water loss from proteins in evolution and development
# 20201216 First version
# 20210127 Added to JMDplots

# Requires:
# CHNOSZ > 1.4.0 to get species' activities in predominant.values

# Create bold axis labels
ZClab <- expression(bolditalic(Z)[bold(C)])
nH2Olab <- expression(bolditalic(n)[bold(H[2]*O)])
DnH2Olab <- expression(bold(Delta)*bolditalic(n)[bold(H[2]*O)])
DPSlab <- expression(bold(Delta*PS))
DPSHPAlab <- expression(bold(Delta*PS~"(HPA)"))
DPSTCGAlab <- expression(bold(Delta*PS~"(TCGA)"))
logaH2Olab <- expression(bold(log)*bolditalic(a)[bold(H[2]*O)])
logfO2lab <- expression(bold(log)*bolditalic(f)[bold(O[2])])

# Chemical analysis of phylostrata and gene age datasets 20201216
evdevH2O1 <- function(pdf = FALSE, boot.R = 99) {
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
  text2 <- c("Amniota", "Mammalia", "Theria", "Eutheria", "Euarchontoglires", "Catarrhini", "Homininae", "")
  text(0.06, ys, 1:8, adj = c(1, 0.5))
  text(0.10, ys, text1, adj = c(0, 0.5))
  text(0.59, ys, 9:16, adj = c(1, 0.5))
  text(0.63, ys, text2, adj = c(0, 0.5))
  text(0.63, ys[8], "Homo sapiens", font = 3, adj = c(0, 0.5))
  title("Phylostrata (Trigos et al., 2017)", font.main = 1)
  par(xpd = FALSE)
  par(opar)

  # Plots 1-3: Protein length and number, ZC, nH2O for Trigos phylostrata
  memo <- plotphylo("nAA", PS_source = "TPPG17", boot.R = boot.R)
  plotphylo("n", PS_source = "TPPG17", memo = memo)
  text(10.5, 500, "length", cex = 0.8)
  text(6, 320, "number / 10", cex = 0.8)
  label.figure("a", cex = 1.7, yfrac = 0.96, xfrac = 0.05, font = 2)
  plotphylo("ZC", PS_source = "TPPG17", memo = memo, boot.R = boot.R)
  plotphylo("nH2O", PS_source = "TPPG17", memo = memo, boot.R = boot.R)

  # Make legend for Liebeskind phylostrata
  opar <- par(mar = c(0, 0, 1.5, 0))
  plot.new()
  par(xpd = NA)
  ys <- seq(0.95, 0.2, length.out = 8)
  text1 <- c("Cellular_organisms", "Euk_Archaea", "Euk+Bac", "Eukaryota", "Opisthokonta", "Eumetazoa", "Vertebrata", "Mammalia")
  text(0.29, ys, 1:8, adj = c(1, 0.5))
  text(0.33, ys, text1, adj = c(0, 0.5))
  title("Gene ages (Liebeskind et al., 2016)", font.main = 1)
  par(xpd = FALSE)
  par(opar)

  # Plots 4-6: Protein length and number, ZC, nH2O for Liebeskind gene ages
  memo <- plotphylo("nAA", PS_source = "LMM16", xlab = "GA", boot.R = boot.R)
  plotphylo("n", PS_source = "LMM16", memo = memo)
  text(2.5, 620, "length", cex = 0.8)
  text(4.5, 70, "number / 10", cex = 0.8)
  label.figure("b", cex = 1.7, yfrac = 0.96, xfrac = 0.05, font = 2)
  plotphylo("ZC", PS_source = "LMM16", memo = memo, xlab = "GA", boot.R = boot.R)
  plotphylo("nH2O", PS_source = "LMM16", memo = memo, xlab = "GA", boot.R = boot.R)

  if(pdf) {
    dev.off()
    addexif("evdevH2O1", "Chemical analysis of Trigos and Liebeskind datasets", "Dick (2021) (preprint)")
  }
}

# Thermodynamic analysis of optimal logaH2O and logfO2 for target proteins 20201219
evdevH2O2 <- function(pdf = FALSE) {

  if(pdf) pdf("evdevH2O2.pdf", width = 7, height = 3)
  par(mar = c(3, 3.1, 3, 1), mgp = c(2, 0.5, 0))

  # Get mean amino acid compositions
  gpa <- getphyloaa("TPPG17")
  PS <- gpa$aa$protein
  # Load target proteins 
  ipPS <- add.protein(gpa$aa)
  # Set up system
  basis("QEC")
  res <- 128
  O2 <- c(-72, -67)
  H2O <- c(-3, 3)

  # Plot A: predominance diagram for all PS
  a <- affinity(O2 = c(O2, res), H2O = c(H2O, res), iprotein = ipPS)
  e <- equilibrate(a, as.residue = TRUE, loga.balance = 0)
  d <- diagram(e, plot.it = FALSE)
  # Make color image for activities
  par.orig <- my.filled.contour(e$vals$O2, e$vals$H2O, d$predominant.values, xlab = logfO2lab, ylab = logaH2Olab,
    nlevels = 50, col = hcl.colors(75, "YlGnBu")[20:75], frame.plot = FALSE,
    # Use plot.axes to label the contour plot (see ?filled.contour)
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
      # Show location of maximum activity for each target protein
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
      title("16 target proteins", font.main = 1)
      label.figure("a", cex = 1.5, yfrac = 0.937, font = 2)
    },
    key.axes = {
      opar <- par(tcl = 0)
      axis(4, at = par("usr")[3:4], labels = round(par("usr")[3:4], 2))
      title(quote(bold(log)*bolditalic(a)[bold(protein)]), cex.main = 1, line = 1)
      par(opar)
    },
    add2 = FALSE
  )

  # Plot B: 16 target proteins and n = 200 sample of human proteins
  seed <- 24
  nbackground <- 200
  # Use same background proteins as in MaximAct()
  TPPG17 <- read.csv(system.file(paste0("extdata/phylostrata/TPPG17.csv.xz"), package = "canprot"), as.is = TRUE)
  LMM16 <- read.csv(system.file(paste0("extdata/phylostrata/LMM16.csv.xz"), package = "canprot"), as.is = TRUE)
  Entry <- na.omit(intersect(TPPG17$Entry, LMM16$UniProt))
  aaback <- protcomp(Entry)$aa
  set.seed(seed)
  iback <- sample(1:nrow(aaback), nbackground)
  ipback <- add.protein(aaback[iback, ])
  # Calculate affinities
  a <- affinity(O2 = c(O2, res), H2O = c(H2O, res), iprotein = c(ipPS, ipback))
  e <- equilibrate(a, as.residue = TRUE, loga.balance = 0)
  d <- diagram(e, plot.it = FALSE)
  # Get optimal logfO2 and logaH2O from MaximAct() for cross-checking 20210404
  MA <- MaximAct(gpa$aa, seed = seed, nbackground = nbackground, O2 = c(O2, res), H2O = c(H2O, res), plot.it = FALSE)
  # Make color-scale diagram for predominant protein activities
  my.filled.contour(e$vals$O2, e$vals$H2O, d$predominant.values, xlab = logfO2lab, ylab = logaH2Olab,
    nlevels = 50,
    col = hcl.colors(75, "YlGnBu")[20:75],
    # Use plot.axes to label the contour plot (see ?filled.contour)
    plot.axes = {
      names <- sapply(strsplit(d$species$name, "\\|"), "[", 2)
      #diagram(e, add = TRUE, names = names, format.names = FALSE)
      # Use a higher resolution for making the lines
      a <- affinity(O2 = c(O2, 400), H2O = c(H2O, 400), iprotein = c(ipPS, ipback))
      diagram(a, as.residue = TRUE, add = TRUE, names = names, format.names = FALSE)
      opar <- par(tcl = 0.3)
      thermo.axis()
      axis(1)
      axis(2)
      par(opar)
      # Show location of maximum activity for each target protein
      par(xpd = TRUE)
      for(i in 1:16) {
        imax <- arrayInd(which.max(e$loga.equil[[i]]), dim(e$loga.equil[[i]]))
        optO2 <- e$vals$O2[imax[1]]
        optH2O <- e$vals$H2O[imax[2]]
        points(optO2, optH2O, pch = 21, bg = 7, cex = 1.5)
        # Cross-check values
        stopifnot(optO2 == MA$O2[i+1])
        stopifnot(optH2O == MA$H2O[i+1])
      }
      par(xpd = FALSE)
      title("16 target + 200 background proteins", font.main = 1)
      label.figure("b", cex = 1.5, yfrac = 0.937, font = 2)
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

# Optimal logaH2O and logfO2 and virtual Eh for target proteins 20201218
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
    datadir <- system.file("extdata/evdevH2O/MaximAct", package = "JMDplots")
    H2Ofile <- file.path(datadir, paste0(PS_source, "_H2O_Hsa.csv"))
    O2file <- file.path(datadir, paste0(PS_source, "_O2_Hsa.csv"))
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
      plot(range(PS), c(-2.5, 2.5), xlab = NA, ylab = logaH2Olab, type = "n", xaxt = "n", xaxs = "i", yaxs = "i")
      mtext("PS", 1, 2.2, font = 2, cex = par("cex"))
      axis(1, at = 1:16, labels = c(1,2,NA,NA,5,NA,NA,NA,NA,10,NA,NA,NA,NA,NA,16))
      # FIXME: redraw "2" because previous command doesn't plot it (spacing too tight)
      axis(1, at = 2, labels = 2)
      # Use transparent gray 20211223
      for(i in 1:nrow(H2O)) lines(PS, H2O[i, iPS], lwd = 0.5, col = "#bababa80")
      lines(PS, meanH2O, lwd = 2, col = 2) 
      label.figure(fig.lab[1], cex = 1.6, font = 2)
    }

    if(!is.na(fig.lab[2])) {
      # Make logfO2 plot
      plot(range(PS), c(-72, -67), xlab = NA, ylab = logfO2lab, type = "n", xaxt = "n", xaxs = "i", yaxs = "i")
      mtext("PS", 1, 2.2, font = 2, cex = par("cex"))
      axis(1, at = 1:16, labels = c(1,2,NA,NA,5,NA,NA,NA,NA,10,NA,NA,NA,NA,NA,16))
      axis(1, at = 2, labels = 2)
      for(i in 1:nrow(O2)) lines(PS, O2[i, iPS], lwd = 0.5, col = "#bababa80")
      lines(PS, meanO2, lwd = 2, col = 2) 
      label.figure(fig.lab[2], cex = 1.6, font = 2)
    }

    if(!is.na(fig.lab[3])) {
      # Make logaH2O plot (mean and logaH2O = 0)
      plot(range(PS), c(-2.5, 2.5), xlab = NA, ylab = logaH2Olab, type = "n", xaxt = "n", xaxs = "i", yaxs = "i")
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
      abline(h = 0, lty = 4, lwd = 1.5, col = "slategray4")
      lines(PS, meanH2O, lwd = 2, col = 2) 
      label.figure(fig.lab[3], cex = 1.6, font = 2)
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
      # Also calculate Eh assuming logaH2O = 0  20210406
      Eh_noH2O <- convert(meanO2, "E0", pH = pH, logaH2O = 0)
      mV_noH2O <- Eh_noH2O * 1000

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
        xtext <- 6.5
      }
      if(PS_source == "LMM16") {
        axis(1, at = 1:8, labels = c(1,NA,NA,4,NA,6,NA,8))
        # PS 1 (Cellular_organisms), 4 (Eukaryota), 6 (Eumetazoa), 8 (Mammalia)
        abline(v = c(1.02, 4, 6, 7.98), col = 5, lwd = 2)
        xtext <- 5.5
      }
      # Plot virtual Eh
      lines(PS, mV, lwd = 2, col = 2) 
      lines(PS, mV_noH2O, lty = 2, lwd = 1, col = 2) 
      # Add lines for measured Eh
      abline(h = c(-150, -199, -241, -318), lty = 3, lwd = 1.5, col = "slategray4")
      # Eh = -150 mV (plasma GSH/GSSG) Jones and Sies, 2015
      text(xtext, -150, "Plasma GSH/GSSG", adj = c(0.5, 1.3))
      # Eh = -199 mV (erythrocyte GSH/GSSG) van 't Erve et al., 2013
      if(PS_source=="TPPG17") text(xtext, -199, "Intracellular GSH/GSSG", adj = c(0.49, -0.3))
      if(PS_source=="LMM16") text(xtext, -199, "Intracellular         GSH/GSSG     ", adj = c(0.5, -0.3))
      # Eh = -241 mV (cytosolic NADH/NAD+) Jones and Sies, 2015
      text(xtext, -241, "Cytosolic NADH/NAD+", adj = c(0.5, 1.3))
      # Eh = -318 mV (mitochondrial NADH/NAD+) Jones and Sies, 2015
      text(xtext, -318, "Mitochondrial NADH/NAD+", adj = c(0.5, -0.3))
      label.figure(fig.lab[4], cex = 1.6, font = 2)
    }

  }

  plotfun("TPPG17", c("a", "b", "c", "d"))
  plotfun("LMM16", c(NA, NA, "e", "f"))

  if(pdf) {
    dev.off()
    addexif("evdevH2O3", "Optimal logaH2O and logfO2 and virtual Eh for target proteins", "Dick (2021) (preprint)")
  }
}

# Chemical metrics for and thermodynamic parameters with different background proteomes 20210711
evdevH2O4 <- function(pdf = FALSE) {

  if(pdf) pdf("evdevH2O4.pdf", width = 8, height = 5)
  mat <- matrix(c(7,1,1,9,4, 8,1,1,0,5, 2,2,3,3,5, 2,2,3,3,6), nrow = 4, byrow = 4)
  layout(mat, heights = c(2,1,1,2), widths = c(1,1,1,1,2))
  par(mar = c(3.1, 3.1, 2, 1), mgp = c(2.1, 0.7, 0))

  # Background human proteome as used in MaximAct()
  TPPG17 <- read.csv(system.file(paste0("extdata/phylostrata/TPPG17.csv.xz"), package = "canprot"), as.is = TRUE)
  LMM16 <- read.csv(system.file(paste0("extdata/phylostrata/LMM16.csv.xz"), package = "canprot"), as.is = TRUE)
  Entry <- na.omit(intersect(TPPG17$Entry, LMM16$UniProt))
  Hsa <- protcomp(Entry)$aa
  Hsa_ZC <- ZCAA(Hsa)
  Hsa_nH2O <- H2OAA(Hsa)
  # Fly proteome
  Dme <- read.csv(system.file("extdata/organisms/UP000000803_7227.csv.xz", package = "JMDplots"), as.is = TRUE)
  Dme_ZC <- ZCAA(Dme)
  Dme_nH2O <- H2OAA(Dme)
  # Bacillus subtilis proteome
  Bsu <- read.csv(system.file("extdata/organisms/UP000001570_224308.csv.xz", package = "JMDplots"), as.is = TRUE)
  Bsu_ZC <- ZCAA(Bsu)
  Bsu_nH2O <- H2OAA(Bsu)

  # Phylostrata target proteins
  gpa <- getphyloaa("TPPG17")
  PS_ZC <- ZCAA(gpa$aa)
  PS_nH2O <- H2OAA(gpa$aa)
  # Biofilm target proteins
  devodir <- system.file("extdata/evdevH2O/devodata", package = "JMDplots")
  biofilm <- read.csv(file.path(devodir, "FOK+21_mean_aa.csv"), as.is = TRUE)
  biofilm_ZC <- ZCAA(biofilm)
  biofilm_nH2O <- H2OAA(biofilm)
  # Fly target proteins
  devodir <- system.file("extdata/evdevH2O/devodata", package = "JMDplots")
  fly <- read.csv(file.path(devodir, "CBS+17_mean_aa.csv"), as.is = TRUE)
  fly_ZC <- ZCAA(fly)
  fly_nH2O <- H2OAA(fly)

  # Get total range for all proteomes
  ZClim <- range(Hsa_ZC, Dme_ZC, Bsu_ZC)
  nH2Olim <- range(Hsa_nH2O, Dme_nH2O, Bsu_nH2O)

  # Human background
  smoothScatter(Hsa_ZC, Hsa_nH2O, xlim = ZClim, ylim = nH2Olim, xlab = ZClab, ylab = nH2Olab)
  abline(v = median(Hsa_ZC), h = median(Hsa_nH2O), lty = 2, col = "darkgray")
  title(paste(length(Hsa_ZC), "human proteins"), font.main = 1)
  # PS target
  points(PS_ZC, PS_nH2O, pch = 15, col = 3, cex = 0.5)
  # biofilm target
  points(biofilm_ZC, biofilm_nH2O, pch = 16, col = 5, cex = 0.5)
  # fly target
  points(fly_ZC, fly_nH2O, pch = 17, col = 6, cex = 0.5)

  # Fly background
  smoothScatter(Dme_ZC, Dme_nH2O, xlim = ZClim, ylim = nH2Olim, xlab = ZClab, ylab = nH2Olab)
  abline(v = median(Dme_ZC), h = median(Dme_nH2O), lty = 2, col = "darkgray")
  title(bquote(.(length(Dme_ZC))~italic("D. melanogaster")~"proteins"), font.main = 1)
  # PS target
  points(PS_ZC, PS_nH2O, pch = 15, col = 3, cex = 0.5)
  # biofilm target
  points(biofilm_ZC, biofilm_nH2O, pch = 16, col = 5, cex = 0.5)
  # fly target
  points(fly_ZC, fly_nH2O, pch = 17, col = 6, cex = 0.5)

  # B. subtilis background
  smoothScatter(Bsu_ZC, Bsu_nH2O, xlim = ZClim, ylim = nH2Olim, xlab = ZClab, ylab = nH2Olab)
  abline(v = median(Bsu_ZC), h = median(Bsu_nH2O), lty = 2, col = "darkgray")
  title(bquote(.(length(Bsu_ZC))~italic("B. subtilis")~"proteins"), font.main = 1)
  # PS target
  points(PS_ZC, PS_nH2O, pch = 15, col = 3, cex = 0.5)
  # biofilm target
  points(biofilm_ZC, biofilm_nH2O, pch = 16, col = 5, cex = 0.5)
  # fly target
  points(fly_ZC, fly_nH2O, pch = 17, col = 6, cex = 0.5)

  # Plot logaH2O, logfO2, and Eh calculated for phylostrata target proteins
  # with background proteins from different organisms  20210712

  # Read results
  datadir <- system.file("extdata/evdevH2O/MaximAct", package = "JMDplots")
  PS_source <- "TPPG17"
  H2O_Hsa <- read.csv(file.path(datadir, paste0(PS_source, "_H2O_Hsa.csv")), as.is = TRUE, check.names = FALSE)
  O2_Hsa <- read.csv(file.path(datadir, paste0(PS_source, "_O2_Hsa.csv")), as.is = TRUE, check.names = FALSE)
  H2O_Dme <- read.csv(file.path(datadir, paste0(PS_source, "_H2O_Dme.csv")), as.is = TRUE, check.names = FALSE)
  O2_Dme <- read.csv(file.path(datadir, paste0(PS_source, "_O2_Dme.csv")), as.is = TRUE, check.names = FALSE)
  H2O_Bsu <- read.csv(file.path(datadir, paste0(PS_source, "_H2O_Bsu.csv")), as.is = TRUE, check.names = FALSE)
  O2_Bsu <- read.csv(file.path(datadir, paste0(PS_source, "_O2_Bsu.csv")), as.is = TRUE, check.names = FALSE)
  # Get phylostrata
  iPS <- 2:ncol(H2O_Hsa)
  PS <- as.numeric(colnames(H2O_Hsa)[iPS])
  # Get mean values of logaH2O and logfO2
  meanH2O_Hsa <- colMeans(H2O_Hsa[, iPS])
  meanO2_Hsa <- colMeans(O2_Hsa[, iPS])
  meanH2O_Dme <- colMeans(H2O_Dme[, iPS])
  meanO2_Dme <- colMeans(O2_Dme[, iPS])
  meanH2O_Bsu <- colMeans(H2O_Bsu[, iPS])
  meanO2_Bsu <- colMeans(O2_Bsu[, iPS])

  # Make logaH2O plot
  par(mgp = c(2, 0.7, 0))
  plot(range(PS), c(-1.5, 1.5), xlab = NA, ylab = logaH2Olab, type = "n", xaxt = "n", xaxs = "i", yaxs = "i")
  abline(h = 0, lty = 4, lwd = 1.5, col = "slategray4")
  mtext("PS", 1, 2.1, font = 2, cex = par("cex"))
  axis(1, at = 1:16, labels = c(1,NA,NA,NA,5,NA,NA,NA,NA,10,NA,NA,NA,NA,15,NA))
  lines(PS, meanH2O_Hsa, lwd = 1.5, col = 4)
  lines(PS, meanH2O_Dme, lwd = 1.5, col = 4, lty = 2)
  lines(PS, meanH2O_Bsu, lwd = 1.5, col = 4, lty = 3)

  # Make logfO2 plot
  plot(range(PS), c(-72, -65), xlab = NA, ylab = logfO2lab, type = "n", xaxt = "n", xaxs = "i", yaxs = "i")
  mtext("PS", 1, 2.1, font = 2, cex = par("cex"))
  axis(1, at = 1:16, labels = c(1,NA,NA,NA,5,NA,NA,NA,NA,10,NA,NA,NA,NA,15,NA))
  lines(PS, meanO2_Hsa, lwd = 1.5, col = 4)
  lines(PS, meanO2_Dme, lwd = 1.5, col = 4, lty = 2)
  lines(PS, meanO2_Bsu, lwd = 1.5, col = 4, lty = 3)

  # Calculate Eh
  Eh_Hsa <- convert(meanO2_Hsa, "E0", pH = 7, logaH2O = meanH2O_Hsa) * 1000
  Eh_Dme <- convert(meanO2_Dme, "E0", pH = 7, logaH2O = meanH2O_Dme) * 1000
  Eh_Bsu <- convert(meanO2_Bsu, "E0", pH = 7, logaH2O = meanH2O_Bsu) * 1000

  # Make Eh plot
  plot(range(PS), c(-300, -100), xlab = NA, ylab = NA, type = "n", xaxt = "n", xaxs = "i", yaxs = "i")
  mtext("PS", 1, 2.1, font = 2, cex = par("cex"))
  mtext("Eh (mV)", side = 2, line = 2, las = 0, cex = par("cex"), font = 2)
  axis(1, at = 1:16, labels = c(1,NA,NA,NA,5,NA,NA,NA,NA,10,NA,NA,NA,NA,15,NA))
  lines(PS, Eh_Hsa, lwd = 1.5, col = 4)
  lines(PS, Eh_Dme, lwd = 1.5, col = 4, lty = 2)
  lines(PS, Eh_Bsu, lwd = 1.5, col = 4, lty = 3)

  # Make legends
  par(mar = c(0, 0, 0, 0))

  plot.new()
  par(xpd = NA)
  legend <- c("high", "Kernel density estimate", "low", "Low-density points", "Median")
  legend("center", legend, text.font = c(3, 1, 3, 1, 1), title = "Background proteome",
    pch = c(19, 19, 19, 19, NA), lty = c(NA, NA, NA, NA, 2), col = c(blues9[c(9, 7, 5, 2)], "slategray4"), bty = "n")
  points(-0.016, 0.372, pch = 15, cex = 0.35)
  label.figure("a", font = 2, cex = 2, xfrac = 0.1, yfrac = 0.9)

  plot.new()
  legend("top", c("Phylostrata", "Biofilm", "Fly development"), pch = c(15, 16, 17), col = c(3, 5, 6), title = "Target proteins", bty = "n", inset = -0.2)
  par(xpd = FALSE)

  plot.new()
  legend("center", c("Human", "D. melanogaster", "B. subtilis"), lty = c(1, 2, 3), col = 4, lwd = 2, title = "Background proteome", bty = "n", text.font = c(1, 3, 3))
  label.figure("b", font = 2, cex = 2, xfrac = 0.9, yfrac = 0.9)

  if(pdf) {
    dev.off()
    addexif("evdevH2O4", "Ranges of chemical metrics and thermodynamic parameters for different background proteomes", "Dick (2021) (preprint)")
  }

}

# Chemical and thermodynamic analysis of B. subtilis biofilm transcriptome and proteome 20201221
evdevH2O5 <- function(pdf = FALSE, boot.R = 99) {

  # Setup plot
  if(pdf) pdf("evdevH2O5.pdf", width = 7, height = 4.5)
  par(mfrow = c(2, 3))
  par(mar = c(4, 4, 1, 1), las = 1, mgp = c(3, 0.8, 0))

  # Read the overall amino acid compositions calculated from Futo et al., 2020 data
  devodir <- system.file("extdata/evdevH2O/devodata", package = "JMDplots")
  aa <- read.csv(file.path(devodir, "FOK+21_mean_aa.csv"), as.is = TRUE)
  # Identify rows with transcriptome and proteome data
  isT <- aa$organism == "transcriptome"
  isP <- aa$organism == "proteome"
  # Identify stages with proteomic data
  iP <- match(aa$protein[isP], aa$protein[isT])

  # Function to plot confidence intervals and means 20210708
  plot_CI_and_mean <- function(X) {
    # Confidence interval for transcriptome
    cix <- c(1:11, 11:1)
    ciy <- c(X$low[isT], rev(X$high[isT]))
    # col = 3 with transparency
    polygon(cix, ciy, col = "#61d04f50", border = NA)
    # Confidence interval for proteome
    cix <- c(iP, rev(iP))
    ciy <- c(X$low[isP], rev(X$high[isP]))
    # col = 4 with transparency
    polygon(cix, ciy, col = "#2297e650", border = NA)
    # Points for transcriptome and proteome
    lines(1:11, X$mean[isT], type = "b", col = 3, pch = 19, cex = 1.3)
    lines(iP, X$mean[isP], type = "b", col = 4, pch = 15, cex = 1.3)
  }

  # Plot A: protein length
  X <- getFOK21("length", boot.R = boot.R)
  plot(c(1, 11), range(X$mean, X$low, X$high), type = "n", xlab = NA, ylab = "Protein length", xaxt = "n", font.lab = 2)
  # Make rotated labels (modified from https://www.r-bloggers.com/rotated-axis-labels-in-r-plots/)
  text(x = 1:11, y = par()$usr[3] - 1.5 * strheight("A"), labels = aa$protein[isT], srt = 45, adj = 1, xpd = TRUE)
  axis(1, at = 1:11, labels = NA)
  abline(v = c(5, 9), lty = 3, lwd = 1.5, col = "gray40")
  plot_CI_and_mean(X)
  legend("bottomleft", c("Proteome", "Transcriptome"), lty = c(1, 1), pch = c(15, 19), col = c(4, 3), pt.cex = 1.3, bty = "n")
  label.figure("a", cex = 1.5, xfrac = 0.025, font = 2)

  # Plot B: nH2O
  X <- getFOK21("H2O", boot.R = boot.R)
  plot(c(1, 11), range(X$mean, X$low, X$high), type = "n", xlab = NA, ylab = nH2Olab, xaxt = "n", font.lab = 2)
  text(x = 1:11, y = par()$usr[3] - 1.5 * strheight("A"), labels = aa$protein[isT], srt = 45, adj = 1, xpd = TRUE)
  axis(1, at = 1:11, labels = NA)
  abline(v = c(5, 9), lty = 3, lwd = 1.5, col = "gray40")
  plot_CI_and_mean(X)
  label.figure("b", cex = 1.5, xfrac = 0.025, font = 2)

  # Plot C: ZC
  X <- getFOK21("ZC", boot.R = boot.R)
  plot(c(1, 11), range(X$mean, X$low, X$high), type = "n", xlab = NA, ylab = ZClab, xaxt = "n", font.lab = 2)
  text(x = 1:11, y = par()$usr[3] - 1.5 * strheight("A"), labels = aa$protein[isT], srt = 45, adj = 1, xpd = TRUE)
  axis(1, at = 1:11, labels = NA)
  abline(v = c(5, 9), lty = 3, lwd = 1.5, col = "gray40")
  plot_CI_and_mean(X)
  label.figure("c", cex = 1.5, xfrac = 0.025, font = 2)

  # Plot D: logaH2O
  par(mgp = c(2.8, 0.8, 0))
  datadir <- system.file("extdata/evdevH2O/MaximAct", package = "JMDplots")
  T_H2O <- read.csv(file.path(datadir, "transcriptome_H2O_Bsu.csv"), as.is = TRUE)[, -1]
  P_H2O <- read.csv(file.path(datadir, "proteome_H2O_Bsu.csv"), as.is = TRUE)[, -1]
  # Get mean values of logaH2O
  T_meanH2O <- colMeans(T_H2O)
  P_meanH2O <- colMeans(P_H2O)
  plot(c(1, 11), range(c(T_meanH2O, P_meanH2O)), type = "n", xlab = "Biofilm stage", ylab = logaH2Olab, xaxt = "n", font.lab = 2)
  text(x = 1:11, y = par()$usr[3] - 1.5 * strheight("A"), labels = aa$protein[isT], srt = 45, adj = 1, xpd = TRUE)
  axis(1, at = 1:11, labels = NA)
  abline(v = c(5, 9), lty = 3, lwd = 1.5, col = "gray40")
  abline(h = 0, lty = 4, lwd = 1.5, col = "slategray4")
  lines(1:11, T_meanH2O, type = "b", col = 3, pch = 19, cex = 1.3)
  lines(iP, P_meanH2O, type = "b", col = 4, pch = 15, cex = 1.3)
  label.figure("d", cex = 1.5, xfrac = 0.025, font = 2)

  # Plot E: logfO2
  T_O2 <- read.csv(file.path(datadir, "transcriptome_O2_Bsu.csv"), as.is = TRUE)[, -1]
  P_O2 <- read.csv(file.path(datadir, "proteome_O2_Bsu.csv"), as.is = TRUE)[, -1]
  # Get mean values of logfO2
  T_meanO2 <- colMeans(T_O2)
  P_meanO2 <- colMeans(P_O2)
  plot(c(1, 11), range(c(T_meanO2, P_meanO2)), type = "n", xlab = "Biofilm stage", ylab = logfO2lab, xaxt = "n", font.lab = 2)
  text(x = 1:11, y = par()$usr[3] - 1.5 * strheight("A"), labels = aa$protein[isT], srt = 45, adj = 1, xpd = TRUE)
  axis(1, at = 1:11, labels = NA)
  abline(v = c(5, 9), lty = 3, lwd = 1.5, col = "gray40")
  abline(h = 0, lty = 4, lwd = 1.5, col = "slategray4")
  lines(1:11, T_meanO2, type = "b", col = 3, pch = 19, cex = 1.3)
  lines(iP, P_meanO2, type = "b", col = 4, pch = 15, cex = 1.3)
  label.figure("e", cex = 1.5, xfrac = 0.025, font = 2)

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
  lines(1:11, T_Eh, type = "b", col = 3, pch = 19, cex = 1.3)
  lines(iP, P_Eh, type = "b", col = 4, pch = 15, cex = 1.3)

  # Add lines for logaH2O = 0  20210713
  T_pe0 <- 0.25 * T_meanO2 - pH - 0.5 * logK
  P_pe0 <- 0.25 * P_meanO2 - pH - 0.5 * logK
  T_Eh0 <- convert(T_pe0, "Eh") * 1000
  P_Eh0 <- convert(P_pe0, "Eh") * 1000
  lines(1:11, T_Eh0, lty = 2, col = 3, pch = 19, cex = 1.3)
  lines(iP, P_Eh0, lty = 2, col = 4, pch = 15, cex = 1.3)

  label.figure("f", cex = 1.5, xfrac = 0.025, font = 2)

  if(pdf) {
    dev.off()
    addexif("evdevH2O5", "Chemical and thermodynamic analysis of B. subtilis biofilm transcriptome and proteome", "Dick (2021) (preprint)")
  }
}

# Organismal water content, proteomic nH2O, and optimal logaH2O for fruit fly development 20210116
evdevH2O6 <- function(pdf = FALSE, boot.R = 99) {

  # Setup plot
  if(pdf) pdf("evdevH2O6.pdf", width = 10, height = 6)
  layout(matrix(c(1,1,1,1, 2,2,2,2, 3,3,3,3, 4,4,4, 5,5,5, 6,6,6, 7,7,7), byrow = TRUE, nrow = 2))
  par(mar = c(5, 4, 3.5, 2))
  par(font.lab = 2, las = 1)

  ## Plot A:
  # Developmental time points and water content values from Fig. 4 of Church and Robertson, 1966
  # doi:10.1002/jez.1401620309
  # Hours 110, 120, 130 are for P1, P2, Adult
  Hours <- c(13.47, 17.4, 23.6, 30.5, 36.17, 48.05, 53.82, 60.02, 72.64, 84.06, 95.17, 103.95, 110, 120, 130)
  Water <- c(79.57, 75.23, 75.33, 81.69, 82.4, 82.81, 76.81, 76.77, 79.58, 76.19, 76.1, 76.01, 66.72, 66.38, 72.04)
  plot(Hours, Water, type = "b", pch = 19, xlim = c(0, 140), xlab = "Hours after hatching                       ", ylab = "Water content (%)",
       xaxs = "i", xaxt = "n", cex = 1.5)
  axis(1, at = seq(0, 100, 20))
  axis(1, at = 110, "P1")
  axis(1, at = 120, "P2")
  axis(1, at = 130, "  Adult")
  title("Organismal water content\n(Church and Robertson, 1966)", font.main = 1)
  label.figure("a", cex = 1.6, font = 2)

  ## Plot B:
  # Read mean amino acid compositions for developmental time points
  devodir <- system.file("extdata/evdevH2O/devodata", package = "JMDplots")
  aa <- read.csv(file.path(devodir, "CBS+17_mean_aa.csv"), as.is = TRUE)
#  # Plot nH2O for model proteins
#  plot(H2OAA(aa), type = "b", xaxt = "n", xlab = "Developmental stage", ylab = nH2Olab, cex = 1.5)

  # Get mean nH2O and bootstrap confidence interval for weighted mean 20210708
  H2O <- getCBS17(boot.R = boot.R)
  # Setup plot
  plot(c(1, 17), range(c(H2O$mean, H2O$low, H2O$high)), type = "n", xaxt = "n", xlab = "Developmental stage", ylab = nH2Olab, cex = 1.5)
  # Reorder points to make shaded CI area with polygon() 20210707
  cix <- c(1:17, 17:1)
  ciy <- c(H2O$low, rev(H2O$high))
  polygon(cix, ciy, col = "lightgray", border = NA)
  points(1:17, H2O$mean, type = "b", cex = 1.5)

  labels <- gsub("p", "P", gsub("e", "E", aa$protein))
  # Label "f" and "m" as sub-labels of "Ay" and "A" 20210713
  labels <- gsub("A", "", gsub("Ay", "", labels))
  text(x = 1:13, y = par()$usr[3] - 1.5 * strheight("A"), labels = labels[1:13], srt = 45, adj = 1, xpd = TRUE)
  text(x = 14:17, y = par()$usr[3] - 1.5 * strheight("A"), labels = labels[14:17], xpd = TRUE)
  axis(1, at = 1:17, labels = NA)
  axis(1, at = c(14.5, 16.5), tick = FALSE, labels = c("Ay", "A"))
  title("Developmental proteome\n(Proteomic data from Casas-Vila et al., 2017)", font.main = 1)
  label.figure("b", cex = 1.6, font = 2)

  ## Plot C:
  # Read optimal logaH2O and logfO2 for differentially expressed proteins of Fabre et al., 2019
  datadir <- system.file("extdata/evdevH2O/MaximAct", package = "JMDplots")
  H2O_embryo <- read.csv(file.path(datadir, "fly_embryo_H2O_Dme.csv"), as.is = TRUE)
  O2_embryo <- read.csv(file.path(datadir, "fly_embryo_O2_Dme.csv"), as.is = TRUE)
  H2O_adult <- read.csv(file.path(datadir, "fly_adult_H2O_Dme.csv"), as.is = TRUE)
  O2_adult <- read.csv(file.path(datadir, "fly_adult_O2_Dme.csv"), as.is = TRUE)
  # Get mean values
  H2O_embryo <- median(colMeans(H2O_embryo)[-1])
  O2_embryo <- median(colMeans(O2_embryo)[-1])
  H2O_adult <- median(colMeans(H2O_adult)[-1])
  O2_adult <- median(colMeans(O2_adult)[-1])

  # MaximAct results for developmental proteome of Casas-Vila et al., 2017  20210403
  # Make logaH2O plot
  H2O <- read.csv(file.path(datadir, "fly_development_H2O_Dme.csv"), as.is = TRUE)[, -1]
  meanH2O <- colMeans(H2O)
  plot(meanH2O, ylim = range(0, meanH2O), type = "b", xaxt = "n", xlab = "Developmental stage", ylab = logaH2Olab, cex = 1.5)
  abline(h = 0, lty = 4, lwd = 1.5, col = "slategray4")
  lines(c(1, 4), rep(H2O_embryo, 2), lty = 1, lwd = 2, col = 4)
  lines(c(14, 17), rep(H2O_adult, 2), lty = 1, lwd = 2, col = 2)
  # Make rotated labels (modified from https://www.r-bloggers.com/rotated-axis-labels-in-r-plots/)
  text(x = 1:13, y = par()$usr[3] - 1.5 * strheight("A"), labels = labels[1:13], srt = 45, adj = 1, xpd = TRUE)
  text(x = 14:17, y = par()$usr[3] - 1.5 * strheight("A"), labels = labels[14:17], xpd = TRUE)
  axis(1, at = 1:17, labels = NA)
  axis(1, at = c(14.5, 16.5), tick = FALSE, labels = c("Ay", "A"))
  text(9, 1, "DP", font = 2)
  text(2.5, 1.55, "DEP\nembryo", font = 2)
  text(15.5, 0.65, "DEP\nadult", font = 2)
  label.figure("c", cex = 1.6, font = 2)
  title("Developmental proteome (DP) and\nDifferentially expressed proteins (DEP)", font.main = 1)

  ## Plot D:
  # ZC and nH2O of differentially expressed proteins in dataset of Fabre et al., 2019  20210401

  # Get amino acid composition of differentially expressed proteins
  pd <- pdat_fly("FKL+19_protein")
  # Get ZC and nH2O for all proteins
  ZC <- ZCAA(pd$pcomp$aa)
  nH2O <- H2OAA(pd$pcomp$aa)
  # Get ZC and nH2O for proteins with higher expression in embryo and adult
  ZC_embryo <- ZC[!pd$up2]
  ZC_adult <- ZC[pd$up2]
  nH2O_embryo <- nH2O[!pd$up2]
  nH2O_adult <- nH2O[pd$up2]

  # First two plots: boxplots for ZC and nH2O 20210403
  # Make boxplot for ZC
  ZCdat <- list(embryo = ZC_embryo, adult = ZC_adult)
  boxplot(ZCdat, ylab = ZClab, outpch = 21, outbg = "#55555580", outcol = NA)
  # Add p-value
  ZC_pvalue <- wilcox.test(ZC_embryo, ZC_adult)$p.value
  legend <- bquote(italic(P) == .(signif(ZC_pvalue, 3)))
  legend("topright", legend = legend, bty = "n")
  label.figure("d", cex = 1.7, yfrac = 0.96, xfrac = 0.05, font = 2)
  # Make boxplot for nH2O
  nH2Odat <- list(embryo = nH2O_embryo, adult = nH2O_adult)
  boxplot(nH2Odat, ylab = nH2Olab, outpch = 21, outbg = "#55555580", outcol = NA)
  # Add p-value
  nH2O_pvalue <- wilcox.test(nH2O_embryo, nH2O_adult)$p.value
  legend <- bquote(italic(P) == .(signif(nH2O_pvalue, 3)))
  legend("bottomleft", legend = legend, bty = "n")

  # Add title for first two plots
  par(xpd = NA)
  text(-0.4, -0.3, "Differentially expressed proteins\n(Proteomic data from Fabre et al., 2019)", cex = 1.3)
  text(-0.4, -1.35, bquote(italic(N) == .(sum(!pd$up2))~"with higher expression in embryos"), cex = 1.3)
  text(-0.4, -1.42, bquote(italic(N) == .(sum(pd$up2))~"with higher expression in adults"), cex = 1.3)
  par(xpd = FALSE)

  ## Plot E:
  # Point symbols to use for embryo and adult
  # embryo: filled blue circle
  # adult: filled red square
  pch <- c(19, 15)
  col <- c(4, 2)
  cex <- c(0.8, 0.9)

  # Plot median ZC and nH2O of differentially expressed proteins
  plot(c(-0.145, -0.135), c(-0.76, -0.72), xlab = ZClab, ylab = nH2Olab, type = "n")
  ZCmed <- c(median(ZC_embryo), median(ZC_adult))
  nH2Omed <- c(median(nH2O_embryo), median(nH2O_adult))
  points(ZCmed, nH2Omed, pch = pch, col = col, cex = cex * 2)
  text(ZCmed, nH2Omed - 0.0025, c("embryo", "adult"))
  label.figure("e", cex = 1.7, yfrac = 0.96, xfrac = 0.05, font = 2)

  # Add title for second two plots
  par(xpd = NA)
  text(-0.131, -0.714, "Median values for differentially expressed proteins", cex = 1.3)
  par(xpd = FALSE)

  # Make logaH2O-logfO2 plot
  plot(c(-71, -70.2), c(0, 2), xlab = logfO2lab, ylab = logaH2Olab, type = "n")
  abline(h = 0, lty = 4, lwd = 1.5, col = "slategray4")
  O2 <- c(O2_embryo, O2_adult)
  H2O <- c(H2O_embryo, H2O_adult)
  points(O2, H2O, pch = pch, col = col, cex = cex * 2)
  text(O2, H2O - 0.13, c("embryo", "adult"))

  if(pdf) {
    dev.off()
    addexif("evdevH2O6", "Organismal water content, proteomic nH2O, and optimal logaH2O for fruit fly development", "Dick (2021) (preprint)")
  }

}

# Evolution of protein ZC in eukaryotic lineages 20211103
evdevH2O7 <- function(pdf = FALSE, H2O = FALSE) {

  # Setup figure
  if(pdf) pdf("evdevH2O7.pdf", width = 10, height = 6)
  mat <- matrix(c(1,1,1,1,1,1, 2,2,2,2, 3,3,3,3,3, 5,5,5,5,5, 4,4,4,4,4, 6,6,6,6,6), nrow = 10)
  layout(mat, widths = c(2, 1, 1))
  par(mar = c(4, 4, 3.5, 1), mgp = c(2.5, 1, 0))

  ## Panel A: gene ages from Liebeskind et al. (2016)

  # Read list of reference proteomes
  datadir <- system.file("extdata/evdevH2O/LMM16", package = "JMDplots")
  refprot <- read.csv(file.path(datadir, "reference_proteomes.csv"))
  # Put HUMAN last so it's more visible
  irp <- 1:nrow(refprot)
  ishuman <- refprot$OSCODE == "HUMAN"
  irp <- c(irp[!ishuman], irp[ishuman])
  refprot <- refprot[irp, ]

  # Start ZC or nH2O plot
  if(!H2O) plot(c(1, 9), c(-0.18, -0.02), xlab = "Gene Age", ylab = ZClab, type = "n", font.lab = 2)
  if(H2O) plot(c(1, 9), c(-0.9, -0.65), xlab = "Gene Age", ylab = nH2Olab, type = "n", font.lab = 2)
  # Add drop line at gene age 5 (Opisthokonta)
  abline(v = 5, lty = 2, col = "gray40")
  # Loop over proteomes
  for(i in 1:nrow(refprot)) {
    # Read modeAge and ZC/nH2O values
    dat <- read.csv(file.path(datadir, "metrics", paste0(refprot$OSCODE[i], ".csv.xz")))
    # Get mean ZC/nH2O for each modeAge
    modeAge <- 1:max(dat$modeAge)
    if(!H2O) X <- sapply(modeAge, function(Age) mean(subset(dat, modeAge == Age)$ZC))
    if(H2O) X <- sapply(modeAge, function(Age) mean(subset(dat, modeAge == Age)$nH2O))
    # Add lines to plot
    col <- "#99999980"
    lwd <- 1
    if(refprot$OSCODE[i] == "HUMAN") {
      col <- 2
      lwd <- 2
    }
    lines(modeAge, X, col = col, lwd = lwd)
  }
  # Add text to indicate divergence at Opisthokonta
  if(H2O) y <- -0.65 else y <- -0.03
  text(3.95, y, "Common ancestors", adj = c(0.5, 0.5))
  text(5.95, y, "Lineages diverge", adj = c(0.5, 0.5))
  # Add labels for divergence times (Kumar et al., 2017)
  par(xpd = NA)
  if(H2O) y <- -0.625 else y <- -0.003
  text(1, y, 4290, srt = 45) # Cellular organisms
  text(4, y, 2101, srt = 45) # Eukaryota
  text(5, y, 1105, srt = 45) # Eukaryota
  text(6, y, 948, srt = 45)  # Eumetazoa
  text(7, y, 615, srt = 45)  # Vertebrata
  text(8, y, 177, srt = 45)  # Mammalia
  par(xpd = FALSE)
  legend("left", "Human", lty = 1, col = 2, lwd = 2, bty = "n")
  title("Divergence times for human lineage (Mya)", line = 2.5, font.main = 1)
  label.figure("a", font = 2, cex = 1.6)

  # Assemble age groups for legend
  modeAges <- read.csv(file.path(datadir, "modeAges.csv"))
  legend <- sapply(1:9, function(i) {
    agetab <- table(modeAges[, paste0("X", i)])
    paste0(i, ": ", paste0(names(agetab), " (", agetab, ")", collapse = ", "))
  })
  # Start plot for legend
  plot.new()
  par(mar = c(1, 1, 1, 1))
  par(xpd = NA)
  legend("center", legend, bty = "n", y.intersp = 2)
  par(xpd = FALSE)

  ## Panel B: Protein families from James et al. (2021)  20211221

  par(mar = c(4, 4, 3.5, 1))

  for(set in c("pfam_plant_nontrans", "pfam_plant_trans", "pfam_animal_nontrans", "pfam_animal_trans")) {

    # Read amino acid composition from file
    datadir <- system.file("extdata/evdevH2O/JWN+21", package = "JMDplots")
    file <- file.path(datadir, paste0(set, "_AA.csv.xz"))
    AA <- read.csv(file)

    # Order data by Age
    AA <- AA[order(AA$protein), ]
    Age <- AA$protein
    # Calculate ZC/nH2O
    if(!H2O) X <- ZC(protein.formula(AA))
    if(H2O) X <- H2OAA(AA)

    # Make boxplots and regression lines
    # Some parts adapted from https://github.com/MaselLab/ProteinEvolution/blob/master/Figures/BoxAndWhiskerPlots_LinearModelSlopes_MetricsVsAge.py
    LinearModel <- lm(X ~ Age)
    Intercept <- coef(LinearModel)["(Intercept)"]
    Slope <- coef(LinearModel)["Age"]
    Pvalue <- summary(LinearModel)$coefficients[2, 4]
    R2 <- summary(LinearModel)$r.squared

    if(H2O) ylab <- nH2Olab else ylab <- ZClab
    plot(Age, X, type = "n", xlab = "Age (Mya)", ylab = ylab, xlim = rev(range(Age)), font.lab = 2)
    boxplot(X ~ Age, add = TRUE, at = unique(Age), boxfill = "lightblue", xaxt = "n", yaxt = "n",
            position = "dodge", varwidth = TRUE, boxwex = 200, outpch = 21, outcex = 0.4, outbg = "#66666680", outcol = NA)
    abline(Intercept, Slope, lwd = 2, col = "dodgerblue")

    # Create title
    main <- gsub("pfam_p", "P", set)
    main <- gsub("pfam_a", "A", main)
    main <- gsub("_nontrans", " non-transmembrane", main)
    main <- gsub("_trans", " transmembrane", main)
    main <- paste0(main, " (", nrow(AA), ")")
    title(main, line = 2.2, font.main = 1, cex.main = 1.1)

    # Add legend
    # NOTE: The ages are revered (x-axis), so slope is multiplied by -1
    stattext <- bquote("Slope"==.(signif(-Slope * 1e3, digits = 2))*"/Gya, "~italic(P)==.(signif(Pvalue, digits = 2)))
    title(stattext, line = 1, font.main = 1, cex.main = 1)

    if(set == "pfam_plant_nontrans") label.figure("b", font = 2, cex = 1.6)

  }

  if(pdf) dev.off()
}

# Calculate optimal logaH2O and logfO2 for various datasets 20210402
# History: Phylostrata 20201218, B. subtilis biofilm 20201221
runMaximAct <- function(dataset = "TPPG17", seed = 1:100, nbackground = 2000, res = 256, bg_organism = "Hsa", writeFiles = TRUE) {

  message(paste("runMaximAct: target", dataset, "with", length(seed), "samples of", bg_organism, "background"))

  # Process 'dataset' argument
  if(dataset %in% c("TPPG17", "LMM16")) {
    if(dataset == "TPPG17") xlab <- "Trigos phylostrata"
    if(dataset == "LMM16") xlab <- "Liebeskind gene ages"
    # Get mean amino acid compositions for phylostrata
    AA_target <- getphyloaa(dataset)$aa
    O2 <- c(-72, -65)
    H2O <- c(-2, 3)
  } else if(dataset %in% c("transcriptome", "proteome")) {
    xlab <- paste("Biofilm", dataset)
    # Read amino acid compositions of overall proteins in each biofilm stage
    devodir <- system.file("extdata/evdevH2O/devodata", package = "JMDplots")
    aa <- read.csv(file.path(devodir, "FOK+21_mean_aa.csv"), as.is = TRUE)
    AA_target <- aa[aa$organism == dataset, ]
    O2 <- c(-72, -65)
    H2O <- c(-2, 5)
  } else if(dataset %in% c("fly_embryo", "fly_adult")) {
    # Get amino acid composition of differentially expressed proteins in fly adults vs embryos 20210402
    pd <- pdat_fly("FKL+19_protein")
    aain <- pd$pcomp$aa
    if(missing(seed)) seed <- 1:20
    # To optimize activities for individual proteins instead of average compositions,
    # we need to extend the ranges of logfO2 and logaH2O
    O2 <- c(-80, -62)
    H2O <- c(-10, 17)
    xlab <- "protein"
    # Down-expressed: Higher in embryo
    if(dataset == "fly_embryo") AA_target <- aain[!pd$up2, ]
    # Up-expressed: Higher in adult
    if(dataset == "fly_adult") AA_target <- aain[pd$up2, ]
  } else if(dataset == "fly_development") {
    # Read mean amino acid compositions for developmental proteome of Casas-Vila et al., 2017 20210403
    devodir <- system.file("extdata/evdevH2O/devodata", package = "JMDplots")
    AA_target <- read.csv(file.path(devodir, "CBS+17_mean_aa.csv"), as.is = TRUE)
    xlab <- "time point"
    O2 <- c(-72, -68)
    H2O <- c(-2, 4)
  }

  # Get different background proteins if specified 20210712
  # Human (Hsa) is the default background in MaximAct()
  if(identical(bg_organism, "Hsa")) AA_background <- NULL
  else if(identical(bg_organism, "Sce")) AA_background <- yeast.aa()
  else if(identical(bg_organism, "Eco")) AA_background <- read.csv(system.file("extdata/organisms/UP000000625_83333.csv.xz", package = "JMDplots"), as.is = TRUE)
  else if(identical(bg_organism, "Dme")) AA_background <- read.csv(system.file("extdata/organisms/UP000000803_7227.csv.xz", package = "JMDplots"), as.is = TRUE)
  else if(identical(bg_organism, "Bsu")) AA_background <- read.csv(system.file("extdata/organisms/UP000001570_224308.csv.xz", package = "JMDplots"), as.is = TRUE)
  else if(identical(bg_organism, "Mja")) {
    AA_background <- read.csv(system.file("extdata/organisms/UP000000805_243232.csv.xz", package = "JMDplots"), as.is = TRUE)
    # We use all the proteins (1787) without sampling
    seed = 1
    # Offset logfO2 and logaH2O
    O2 <- O2 + 7
    H2O <- H2O - 2
  }
  else stop(paste("unrecognized background organism:", bg_organism))
  bgname <- paste0("_", bg_organism)

  # Set resolution 20210714
  O2 <- c(O2, res)
  H2O <- c(H2O, res)

  # Start plot
  if(writeFiles) {
    png(paste0(dataset, bgname, ".png"), width = 1000, height = 1000)
    par(cex = 2)
  }
  # Run MaximAct()
  MA <- MaximAct(AA_target, seed = seed, nbackground = nbackground, xlab = xlab, O2 = O2, H2O = H2O, AA_background = AA_background)
  if(writeFiles) {
    # Round and save values
    O2vals <- round(MA$O2, 3)
    H2Ovals <- round(MA$H2O, 3)
    write.csv(O2vals, paste0(dataset, "_O2", bgname, ".csv"), row.names = FALSE, quote = FALSE)
    write.csv(H2Ovals, paste0(dataset, "_H2O", bgname, ".csv"), row.names = FALSE, quote = FALSE)
    # Close plot
    dev.off()
  }
  # Return output of MaximAct()
  MA

}

# Example of protein chemical formula, formation reaction, and equilibrium constant 20210115
# Add logQ 20210719
LYSC_example <- function() {
  ip <- pinfo("LYSC_CHICK")
  pl <- protein.length(ip)

  # Calculate the per-residue formula, rounded to 3 decimal places
  residue.formula <- round(protein.formula(ip) / pl, 3)
  # Format it as a text object
  formula <- as.chemical.formula(residue.formula)
  # Calculate per-residue ΔGf° at 25 °C and 1 bar
  G <- suppressMessages(protein.OBIGT(ip)$G / pl)
  # Add the residue as a new species
  suppressMessages(mod.OBIGT("LYSC_residue", formula = formula, G = G))

  # Calculate properties of formation reaction from QEC basis species
  basis("QEC")
  sres <- suppressMessages(subcrt("LYSC_residue", 1, T = 25))

  # Print results
  message("Per-residue chemical formula of LYSC_CHICK")
  print(formula)
  message("Stoichiometry of formation reaction")
  print(sres$reaction[, 1:3])
  message("logK of formation reaction (25 \u00B0C, 1 bar)")
  print(round(sres$out$logK, 2))

  # Sanity check: we get the same equilibrium constant starting with the whole formula
  logK.protein <- suppressMessages(subcrt("LYSC_CHICK", 1, T = 25)$out$logK)
  logK.residue <- logK.protein / pl
  stopifnot(all.equal(sres$out$logK, logK.residue, tol = 0.01, scale = 1))

  # Calculate logQ 20210719
  # Use activities of amino acids and H2O set in basis("QEC"), but with logfO2 = -70
  basis("O2", -70)
  logQ.residue_round <- suppressMessages(subcrt("LYSC_residue", 1, T = 25, logact = 0)$out$logQ)

  # The reaction auto-balance is not very precise, so let's do it another way
  # For activity of residue = 1, activity of protein is 1/length
  loga.protein <- log10(1/pl)
  logQ.protein <- suppressMessages(subcrt("LYSC_CHICK", 1, T = 25, logact = loga.protein)$out$logQ)
  logQ.residue <- logQ.protein / pl
  # Sanity check 2: logQ is similar for the two methods above
  stopifnot(all.equal(logQ.residue, logQ.residue_round, tol = 0.02, scale = 1))

  message("logQ of formation reaction (logfO2 = -70, logaH2O = 0)")
  print(round(logQ.residue, 2))
  message("log(K/Q)")
  logKQ <- logK.residue - logQ.residue
  print(round(logKQ, 2))

  # Sanity check 3: same value for log(K/Q) is returned by affinity()
  species("LYSC_CHICK", loga.protein)
  a <- suppressMessages(affinity())
  A.protein <- a$values[[1]][1]
  A.residue <- A.protein / pl
  stopifnot(all.equal(logKQ, A.residue))
}

############################
### UNEXPORTED FUNCTIONS ###
############################

# Mean ZC and nH2O of phylostrata 20191122
plotphylo <- function(var = "ZC", PS_source = "TPPG17", memo = NULL, xlab = "PS", boot.R = 99) {
  if(is.null(memo)) {
    dat <- read.csv(system.file("extdata/phylostrata/TPPG17.csv.xz", package = "canprot"), as.is = TRUE)
    if(PS_source == "LMM16") {
      dat <- read.csv(system.file("extdata/phylostrata/LMM16.csv.xz", package = "canprot"), as.is = TRUE)
      colnames(dat)[c(1,3)] <- c("Entry", "Phylostrata")
      # remove entries that have ENSP instead of UniProt IDs
      dat <- dat[!grepl("^ENSP", dat$Entry), ]
    }
    # Update old UniProt IDs
    dat <- check_IDs(dat, "Entry")
    # Remove genes with no UniProt mapping 20210718
    dat <- dat[!is.na(dat$Entry), ]
    # Run protcomp and suppress warning about duplicated IDs 20210718
    pcomp <- suppressWarnings(protcomp(dat$Entry))
    AA <- pcomp$aa
  } else {
    dat <- memo$dat
    AA <- memo$AA
  }
  # Use AA functions (H2OAA, ZCAA) because protcomp no longer returns these values 20201216
  X <- switch(var, ZC = ZCAA(AA), nH2O = H2OAA(AA), nAA = protein.length(AA), n = NA, Cost = CostAA(AA))
  # Get mean chemical metrics for each phylostratum
  PS <- sort(unique(dat$Phylostrata))
  high.X <- low.X <- cum.X <- mean.X <- numeric()
  for(p in PS) {
    if(var %in% c("ZC", "nH2O", "nAA", "Cost")) {
      # Point mean
      this.X <- X[dat$Phylostrata == p]
      mean.X <- c(mean.X, mean(this.X))
      # Confidence interval from bootstrap 20210707
      set.seed(1234)
      samplemean <- function(x, i) mean(x[i])
      boot.X <- boot::boot(this.X, samplemean, R = boot.R)
      ci.X <- boot::boot.ci(boot.X, conf = 0.95, type = "perc")
      low.X <- c(low.X, ci.X$percent[4])
      high.X <- c(high.X, ci.X$percent[5])
      # Cumulative mean
      cum.X <- c(cum.X, mean(X[dat$Phylostrata <= p]))
    }
    if(var == "n") {
      # Number of proteins in this phylostratum
      mean.X <- c(mean.X, sum(dat$Phylostrata == p))
      cum.X <- c(cum.X, sum(dat$Phylostrata <= p))
    }
  }
  ylab <- switch(var, ZC = ZClab, nH2O = nH2Olab, nAA = quote("Protein length or"~italic(n)/10), Cost = "Cost")
  if(var %in% c("ZC", "nH2O", "nAA", "Cost")) {
    ylim <- range(c(mean.X, low.X, high.X))
    # Extend the ylim to zero for protein length and number plot 20210713
    if(var == "nAA") ylim <- range(0, ylim)
    plot(range(PS), ylim, type = "n", xlab = xlab, ylab = ylab, font.lab = 2)
    # Reorder points to make shaded CI area with polygon() 20210707
    cix <- c(PS, rev(PS))
    ciy <- c(low.X, rev(high.X))
    polygon(cix, ciy, col = "lightgray", border = NA)
    # Plot the point
    points(PS, mean.X, pch = 19, type = "b", cex = 0.7)
    # Don't plot the cumulative mean 20210824
    #lines(PS, cum.X, col = 2, lty = 2)
  }
  if(var == "n") {
    # Add points for (number of proteins) / 10  20210713
    points(PS, mean.X / 10, type = "b", cex = 0.8)
  }
  # return the dat and AA for memoization 20191211
  invisible(list(dat = dat, AA = AA))
}

# Calculate metabolic cost from amino acid compositions of proteins 20211220
CostAA <- function(AAcomp) {
  # Amino acid cost from Akashi and Gojobori (2002)
  # doi:10.1073/pnas.062526999
  Cost <- c(Ala = 11.7, Cys = 24.7, Asp = 12.7, Glu = 15.3, Phe = 52.0,
            Gly = 11.7, His = 38.3, Ile = 32.3, Lys = 30.3, Leu = 27.3,
            Met = 34.3, Asn = 14.7, Pro = 20.3, Gln = 16.3, Arg = 27.3,
            Ser = 11.7, Thr = 18.7, Val = 23.3, Trp = 74.3, Tyr = 50.0)
  # find columns with names for the amino acids
  isAA <- colnames(AAcomp) %in% names(Cost)
  iAA <- match(colnames(AAcomp)[isAA], names(Cost))
  # calculate total of cost values for each protein
  sumCost <- rowSums(t(t(AAcomp[, isAA]) * Cost[iAA]))
  # divide by length of proteins to get average
  sumCost / rowSums(AAcomp[, isAA])
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
          frame.plot = axes, add2 = NULL, plot.key = TRUE, ...)
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

    if(plot.key) {
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
    }

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
