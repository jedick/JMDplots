# JMDplots/gradH2O.R
# 20190711 R file started for Goldschmidt poster
# 20190930 started move to JMDplots for salinity gradients paper
# 20200402 paper now on bioRxiv (doi:10.1101/2020.04.01.020008)

# number of reactions for each amino acid in E. coli metabolic recconstruction 20191101
# (data from Feist et al., 2007; doi:10.1038/msb4100155)
gradH2O0 <- function() {
  # define compound names used for amino acids
  # (from "compounds" table of FHR+07 spreadsheet)
  AAnames <- list(
    Ala = " ala-L ",
    Cys = " cys-L ", # but not hcys-l
    Asp = " asp-L ",
    Glu = " glu-L ",
    Phe = " phe-L ",
    Gly = " gly ",   # but not glyald, glyb, glyc, ...
    His = " his-L ",
    Ile = " ile-L ",
    Lys = " lys-L ",
    Leu = " leu-L ",
    Met = " met-L ",
    Asn = " asn-L ",
    Pro = " pro-L ",
    Gln = " gln-L ",
    Arg = " arg-L ",
    Ser = " ser-L ", # but not pser-L
    Thr = " thr-L ", # but not athr-L
    Val = " val-L ",
    Trp = " trp-L ",
    Tyr = " tyr-L "
  )
  # read reaction list
  file <- system.file("extdata/gradH2O/reaction_equations.csv", package = "JMDplots")
  dat <- read.csv(file, as.is = TRUE)
  rxn <- dat$equation
  # add spaces before and after each reaction
  rxn <- paste0(" ", rxn, " ")
  # replace [ with space (e.g. [c] after compound name)
  rxn <- gsub("[", " ", rxn, fixed = TRUE)
  # count the number of reactions with each amino acid
  sort(sapply(lapply(AAnames, grepl, rxn, fixed = TRUE), sum), decreasing = TRUE)
}

# basis species comparison 20190713 / stoichiometric hydration state 20191005
gradH2O1 <- function(pdf = FALSE) {

  # set up figure
  if(pdf) pdf("gradH2O1.pdf", width = 8, height = 4)
  par(mfrow = c(2, 4))
  par(mar = c(3.2, 3.4, 2.5, 1))
  par(mgp = c(2.2, 0.7, 0))
  par(las = 1)
  par(cex.lab = 1.2)

  # define axis labels
  nH2Olab <- expression(italic(n)[H[2] * O])
  nO2lab <- expression(italic(n)[O[2]])
  ZClab <- expression(italic(Z)[C])

  # function to plot linear model
  lmfun <- function(ZC, y, xlim = par("usr")[1:2], legend.x = NULL, ...) {
    mylm <- lm(y ~ ZC)
    lines(xlim, predict(mylm, data.frame(ZC = xlim)), ...)
    # add R-squared text
    if(!is.null(legend.x)) {
      R2 <- format(round(summary(mylm)$r.squared, 3), nsmall = 3)
      R2txt <- substitute(italic(R)^2 == R2, list(R2 = R2))
      legend(legend.x, legend = R2txt, bty = "n")
    }
    invisible(round(residuals(mylm), 3))
  }

  # function to add plot title and panel label
  mainlab <- function(main, lab) {
    title(main, font.main = 1, cex.main = 1)
    label.figure(lab, xfrac = 0.15, yfrac = 0.92, cex = 1.6)
  }

  # function to plot nH2O or nO2 vs ZC of amino acids
  aaplot <- function(ZC, y, ylab, legend.x, main, lab) {
    plot(ZC, y, type = "p", pch = aminoacids(1), xlab = ZClab, ylab = ylab)
    mainlab(main, lab)
    lmfun(ZC, y, c(-1, 1), legend.x)
  }

  scatterH2O <- function(ZC, nH2O, legend.x, main, lab) {
    smoothScatter(ZC, nH2O, xlab = ZClab, ylab = nH2Olab, colramp = colorRampPalette(c("transparent", blues9)))
    lmfun(ZC, nH2O, par("usr")[1:2], legend.x, lty = 2, lwd = 2, col = "grey40")
    mainlab(main, lab)
  }

  scatterO2 <- function(ZC, nO2, legend.x, main, lab) {
    smoothScatter(ZC, nO2, xlab = ZClab, ylab = nO2lab, colramp = colorRampPalette(c("transparent", blues9)))
    lmfun(ZC, nO2, par("usr")[1:2], legend.x, lty = 2, lwd = 2, col = "grey40")
    mainlab(main, lab)
  }

  # get names of amino acids
  aa <- aminoacids("")
  # calculate ZC of the amino acids
  ZC.aa <- ZC(info(aa, "aq"))

  # get amino acid compositions of E. coli proteins (UniProt)
  ecoli <- read.csv(system.file("/extdata/organisms/ecoli.csv.xz", package = "JMDplots"), as.is = TRUE)
  pf.ecoli <- protein.formula(ecoli)
  ZC.ecoli <- ZC(pf.ecoli)
  # get nH2O and nO2 with QCa basis species
  basis(c("glutamine", "cysteine", "acetic acid", "H2O", "O2"))
  pb.ecoli <- protein.basis(ecoli)
  pl.ecoli <- protein.length(ecoli)
  nH2O.ecoli <- pb.ecoli[, "H2O"] / pl.ecoli
  nO2.ecoli <- pb.ecoli[, "O2"] / pl.ecoli
  # check that nH2O equals that calculated with canprot::H2OAA
  nH2O.ref <- H2OAA(ecoli, basis = "QCa")
  stopifnot(all(nH2O.ecoli == nH2O.ref))

  # Calculate nH2O and nO2 of amino acids with CHNOS basis
  basis("CHNOS")
  species(aa)
  # Plot 1: nH2O-ZC of amino acids (CHNOS)
  aaplot(ZC.aa, species()$H2O, nH2Olab, "bottomleft", "Amino acids (CHNOS)", "(a)")
  # Plot 2: nO2-ZC of amino acids (CHNOS)
  aaplot(ZC.aa, species()$O2, nO2lab, "bottomright", "Amino acids (CHNOS)", "(b)")

  # Calculate nH2O and nO2 with QEC basis
  basis("QEC")
  species(aa)
  # Plot 3: nH2O-ZC of amino acids (QEC)
  aaplot(ZC.aa, species()$H2O, nH2Olab, "bottomright", "Amino acids (QEC)", "(c)")
  # Plot 4: nO2-ZC of amino acids (QEC)
  aaplot(ZC.aa, species()$O2, nO2lab, "bottomright", "Amino acids (QEC)", "(d)")

  ## Calculate nH2O and nO2 of amino acids with QCa basis
  basis(c("glutamine", "cysteine", "acetic acid", "H2O", "O2"))
  species(aa)
  # Plot 5: nH2O-ZC of amino acids (QCa)
  aaplot(ZC.aa, species()$H2O, nH2Olab, "bottomright", "Amino acids (QCa)", "(e)")
  # Plot 6: nO2-ZC of amino acids (QCa)
  aaplot(ZC.aa, species()$O2, nO2lab, "bottomright", "Amino acids (QCa)", "(f)")

  # Plot 7: nH2O-ZC of E. coli proteins (QCa)
  scatterH2O(ZC.ecoli, nH2O.ecoli, "bottomright", quote(italic(E.~coli)*" proteins (QCa)"), "(g)")
  # Plot 8: nO2-ZC of E. coli proteins (QCa)
  scatterO2(ZC.ecoli, nO2.ecoli, "bottomright", quote(italic(E.~coli)*" proteins (QCa)"), "(h)")

  if(pdf) {
    dev.off()
    addexif("gradH2O1", "Derivation of stoichiometric hydration state", "Dick et al. (2020) (preprint)")
  }
  
  # Return linear model of nH2O-ZC for amino acids
  invisible(lm(species()$H2O ~ ZC.aa))
}

# nH2O-ZC scatterplots for redox gradients and the Baltic Sea 20190713
gradH2O2 <- function(pdf = FALSE, vars = "H2O-ZC", lm = NULL) {
  if(pdf) pdf("gradH2O2.pdf", width = 12, height = 5.6)
  par(mfrow = c(1, 2))
  par(mar = c(4, 4.5, 2, 1), las = 1, cex = 1.2)
  par(cex.lab = 1.3)

  # plot 1: compare ZC and nH2O of proteins in datasets from gradox paper
  mgradox <- ppage("gradoxGS", plot.it = FALSE)
  pgradox <- ppage("gradoxGS", H2O = TRUE, plot.it = FALSE)
  pcomp(mgradox, pgradox, reorder = FALSE, yline = 3.2, vars = vars, cex.ylab = 1.5, font = 2, labdy = 0.003)
  # add proteomes from Nif-encoding genomes 20191014
  np <- NifProteomes()
  if(vars == "H2O-ZC") {
    points(np$ZC, np$nH2O, pch = 15)
    lines(np$ZC, np$nH2O, col = "dimgray", lwd = 0.8, lty = 2)
  }
  if(vars == "pIG") {
    points(np$pI, np$GRAVY, pch = 15)
    lines(np$pI, np$GRAVY, col = "dimgray", lwd = 0.8, lty = 2)
  }
  # add text labels
  if(vars == "H2O-ZC") {
    text(-0.21, -0.007, "hot spring", cex = 0.8)
    text(-0.21, -0.0095, "source", cex = 0.8)
    text(-0.157, -0.006, "photo-", adj = 0, cex = 0.8)
    text(-0.157, -0.0085, "trophic", adj = 0, cex = 0.8)
    text(-0.157, -0.011, "zone", adj = 0, cex = 0.8)
    text(-0.17, -0.015, "photosynthetic", cex = 0.8)
    text(-0.17, -0.018, "fringe", cex = 0.8)
    text(-0.127, -0.002, "> 3 mm", adj = 0, cex = 0.8)
    text(-0.138, 0.001, "3 mm", adj = 0, cex = 0.8)
    text(-0.134, 0.0035, "2 mm", adj = 0, cex = 0.8)
    text(-0.131, 0.012, "1 mm", adj = 0, cex = 0.8)
    text(-0.205, 0.011, "vent", cex = 0.8)
    text(-0.205, 0.0085, "fluids", cex = 0.8)
    text(-0.157, 0.013, "plume", cex = 0.8)
    text(-0.1615, 0.005, "sea-", adj = 0, cex = 0.8)
    text(-0.1615, 0.0025, "water", adj = 0, cex = 0.8)
    text(c(-0.188, -0.212, -0.186, -0.161), c(0.0455, 0.0344, 0.0257, 0.0248), c("Nif-D", "Nif-C", "Nif-B", "Nif-A"), adj = 0, cex = 0.8)
    text(-0.218, 0.0365, "NF", cex=0.7, font = 2)
    # add H2O-ZC fit for E. coli 20200819
    if(!is.null(lm)) {
      x <- par("usr")[1:2]
      y <- predict(lm, data.frame(ZC.aa = x))
      for(dy in seq(-0.50, -0.60, -0.01)) lines(x, y + dy, col = "gray80")
    }
  }
  if(vars == "pIG") {
    text(c(5.00, 7.49, 5.65, 5.68), c(-0.123, -0.130, -0.109, -0.066), c("Nif-D", "Nif-C", "Nif-B", "Nif-A"), adj = 0)
    text(c(5.37, 5.42, 5.57), c(-0.158, -0.176, -0.216), c("1 mm", "2 mm", "3 mm"), adj = 1)
  }
  title("Redox gradients", font.main = 1)
  label.figure("(a)", cex = 2, xfrac = 0.035)

  # plot 2: compare ZC and nH2O of proteins in Baltic Sea surface
  mbaltics <- ppage("balticsurface", plot.it = FALSE)
  pbaltics <- ppage("balticsurface", H2O = TRUE, plot.it = FALSE)
  pcomp(mbaltics, pbaltics, reorder = FALSE, yline = 3.2, labels.at = NA, vars = vars, pch = list(c(2, 2, 2, 2, 17, 17, 17, 17, 17, 17)), cex.ylab = 1.8)
  # compare ZC and nH2O of proteins in Baltic Sea 10-20 m
  mbalticd <- ppage("balticdeep", plot.it = FALSE)
  pbalticd <- ppage("balticdeep", H2O = TRUE, plot.it = FALSE)
  pcomp(mbalticd, pbalticd, reorder = FALSE, add = TRUE, lty = 3, labels.at = NA, vars = vars, pch = list(c(6, 6, 6, 6, 25, 25, 25, 25, 25)))
  # add text labels
  if(vars == "H2O-ZC") {
    text(-0.14, 0.037, "surface\n< 6 PSU")
    text(-0.135, 0.026, "surface\n> 6 PSU")
    text(-0.185, 0.031, "chl a max\n< 6 PSU")
    text(-0.18, 0.015, "chl a max\n> 6 PSU")
    # add H2O-ZC fit for E. coli 20200819
    if(!is.null(lm)) {
      x <- par("usr")[1:2]
      y <- predict(lm, data.frame(ZC.aa = x))
      for(dy in seq(-0.50, -0.60, -0.01)) lines(x, y + dy, col = "gray80")
    }
  }
  title("Baltic Sea", font.main = 1)
  label.figure("(b)", cex = 2, xfrac = 0.035)

  if(pdf) {
    dev.off()
    addexif("gradH2O2", "nH2O-ZC scatterplots for redox gradients and the Baltic Sea", "Dick et al. (2020) (preprint)")
  }
}

# nH2O for Baltic Sea metagenome and metatranscriptome in different size fractions 20190715
gradH2O3 <- function(pdf = FALSE, var = NULL) {
  if(pdf) pdf("gradH2O3.pdf", width = 6, height = 3.5)
  nvar <- ifelse(is.null(var), 1, length(var))
  if(nvar==2) par(mfrow = c(2, 3))
  else par(mfrow = c(1, 3))
  par(mar = c(3, 4, 2, 1), mgp = c(3.5, 0.7, 0), las = 1)
  par(cex.lab = 1.5)

  for(i in 1:nvar) {
    ylim <- NULL
    if(identical(var[i], "pI")) ylim <- c(5, 7.5)
    figlab <- c("(a)", "(b)", "(c)")
    if(i==2) figlab <- c("(d)", "(e)", "(f)")

    mplot("Baltic_Sea-0.1s", "iMicrobe_MGP", H2O = TRUE, plottype = "#FF000030",
          col = "red", add.title = FALSE, ylim = ylim, yline = 2.7, var = var[i], srt = NA, ilabel = c(1, 12))
    mplot("Baltic_Sea-0.1s", "SRA_MTP", H2O = TRUE, plottype = "#0000FF30",
          col = "blue", add.title = FALSE, add = TRUE, pch = 15, var = var[i])
    title(quote("0.1-0.8"~mu*m))
    legend("bottomleft", c("metagenomes", "metatranscriptomes"), pch = c(19, 15), col = c("red", "blue"), lty = c(2, 2), bty = "n")
    label.figure(figlab[1], cex = 1.7)

    mplot("Baltic_Sea-0.8s", "iMicrobe_MGP", H2O = TRUE, plottype = "#FF000030",
          col = "red", add.title = FALSE, ylim = ylim, yline = 2.7, var = var[i], srt = NA, ilabel = c(1, 12))
    mplot("Baltic_Sea-0.8s", "SRA_MTP", H2O = TRUE, plottype = "#0000FF30",
          col = "blue", add.title = FALSE, add = TRUE, pch = 15, var = var[i])
    title(quote("0.8-3.0"~mu*m))
    label.figure(figlab[2], cex = 1.7)

    mplot("Baltic_Sea-3.0s", "iMicrobe_MGP", H2O = TRUE, plottype = "#FF000030",
          col = "red", add.title = FALSE, ylim = ylim, yline = 2.7, var = var[i], srt = NA, ilabel = c(1, 12))
    mplot("Baltic_Sea-3.0s", "SRA_MTP", H2O = TRUE, plottype = "#0000FF30",
          col = "blue", add.title = FALSE, add = TRUE, pch = 15, var = var[i])
    title(quote("3.0-200"~mu*m))
    label.figure(figlab[3], cex = 1.7)
  }

  if(pdf) {
    dev.off()
    addexif("gradH2O3", "nH2O for Baltic Sea metagenome and metatranscriptome in different size fractions", "Dick et al. (2020) (preprint)")
  }
}

# nH2O vs ZC for freshwater, marine, and hypersaline environments 20191004
# make separate plot for sediments 20191006
# add Amazon River 20191007
# remove sediments and add GRAVY - pI plots 20191027
gradH2O4 <- function(pdf = FALSE) {

  if(pdf) pdf("gradH2O4.pdf", width = 8, height = 5)
  layout(matrix(1:6, nrow = 2))
  par(las = 1, mar = c(4, 4.2, 2, 1), mgp = c(2.5, 1, 0))
  par(cex.lab = 1.5)
  xlim <- c(-0.2, -0.08)
  if(getOption("basis") == "QCa") ylim <- c(-1.15, -1)
  if(getOption("basis") == "QEC") ylim <- c(-0.8, -0.72)

  # plots 1-2: Amazon river metagenome
  mout <- ppage("amazon", plot.it = FALSE)
  pout <- ppage("amazon", H2O = TRUE, plot.it = FALSE)
  pcomp(mout, pout, xlim = xlim, ylim = ylim, lty = 0, yline = 2.8, labels.at = NA)
  hullfun(mout, pout, c(1, 3), "green3", c("riverFL", "riverPA"))
  hullfun(mout, pout, c(1, 3), "blue", c("plumeFL", "plumePA"))
  text(c(-0.13, -0.13), c(0.010, -0.005), c("river", "plume"))
  legend("topleft", c("free-living", "particle-associated"), pch = c(20, 15), bty = "n", inset = c(0.05, 0))
  legend("topleft", legend = c(NA, NA), pch = c(21, 22), pt.cex = c(0.7, 1), bty = "n")
  title("Amazon River metagenomes", font.main = 1)
  label.figure("(a)", xfrac = 0.1, cex = 1.8)
  # plot 2: GRAVY - pI
  pcomp(mout, pout, lty = 0, yline = 3, vars = "pIG", labels.at = NA, cex.ylab = 0.9)
  hullfun(mout, pout, c(1, 3), "green3", c("riverFL", "riverPA"), vars = "pIG")
  hullfun(mout, pout, c(1, 3), "blue", c("plumeFL", "plumePA"), vars = "pIG")
  text(c(7.2, 8.1), c(-0.14, -0.175), c("river", "plume"))
  title("Amazon River metagenomes", font.main = 1)
  label.figure("(d)", xfrac = 0.1, cex = 1.8)

  # plots 3-4: Amazon river metatranscriptome
  pcomp(mout, pout, "MT", xlim = xlim, ylim = ylim, lty = 0, yline = 2.8, labels.at = NA)
  hullfun(mout, pout, c(2, 4), "green3", c("riverFL", "riverPA"))
  hullfun(mout, pout, c(2, 4), "blue", c("plumeFL", "plumePA"))
  text(c(-0.125, -0.12), c(0.025, 0.005), c("river", "plume"))
  title("Amazon River metatranscriptomes", font.main = 1)
  label.figure("(b)", xfrac = 0.1, cex = 1.8)
  # plot 4: GRAVY - pI
  pcomp(mout, pout, "MT", lty = 0, yline = 3, vars = "pIG", labels.at = NA, cex.ylab = 0.9)
  hullfun(mout, pout, c(2, 4), "green3", c("riverFL", "riverPA"), vars = "pIG")
  hullfun(mout, pout, c(2, 4), "blue", c("plumeFL", "plumePA"), vars = "pIG")
  text(c(8.2, 6.5), c(-0.105, -0.11), c("river", "plume"))
  title("Amazon River metatranscriptomes", font.main = 1)
  label.figure("(e)", xfrac = 0.1, cex = 1.8)

  # start plot 5: Eiler et al. (freshwater vs marine)
  moutE <- ppage("eiler", plot.it = FALSE)
  poutE <- ppage("eiler", H2O = TRUE, plot.it = FALSE)
  pcomp(moutE, poutE, xlim = xlim, ylim = ylim, lty = 0, yline = 2.8, labels.at = NA)
  hullfun(moutE, poutE, 1, "green3")
  hullfun(moutE, poutE, 2, "blue")

  # add hypersaline water data
  moutH <- ppage("hypersaline", plot.it = FALSE)
  poutH <- ppage("hypersaline", H2O = TRUE, plot.it = FALSE)
  pcomp(moutH, poutH, reorder = FALSE, add = TRUE, labdx = 0.006, labdy = 0.004)
  hullfun(moutH, poutH, 1:3, "turquoise3")
  text(c(-0.16, -0.16, -0.12), c(0.037, -0.007, -0.002), c("freshwater", "marine", "hypersaline"))
  legend("bottomright", c("lower salinity", "higher salinity"), pch = c(0, 15), col = "turquoise3", bty = "n")
  legend("bottomright", c("hypersaline datasets", "", ""), bty = "n")
  title("Freshwater - marine - hypersaline", font.main = 1)
  label.figure("(c)", xfrac = 0.1, cex = 1.8)

  # plot 6: GRAVY - pI
  pcomp(moutE, poutE, lty = 0, yline = 3, vars = "pIG", labels.at = NA, cex.ylab = 0.9)
  hullfun(moutE, poutE, 1, "green3", vars = "pIG")
  hullfun(moutE, poutE, 2, "blue", vars = "pIG")
  pcomp(moutH, poutH, reorder = FALSE, add = TRUE, vars = "pIG", labels.at = "min", labdx = -0.25)
  hullfun(moutH, poutH, 1:3, "turquoise3", vars = "pIG")
  text(c(7.5, 7.4, 6.3), c(-0.14, -0.20, -0.27), c("freshwater", "marine", "hypersaline"))
  title("Freshwater - marine - hypersaline", font.main = 1)
  label.figure("(f)", xfrac = 0.1, cex = 1.8)

  if(pdf) {
    dev.off()
    addexif("gradH2O4", "nH2O vs ZC for freshwater, marine, and hypersaline environments", "Dick et al. (2020) (preprint)")
  }
}

# scatterplots of ZC and nH2O vs GRAVY 20191117
gradH2O5 <- function(pdf = FALSE) {
  if(pdf) pdf("gradH2O5.pdf", width = 7, height = 33)
  par(las = 1, mar = c(4, 4.2, 1, 1), mgp = c(2.5, 1, 0))
  par(cex.lab = 1.3)
  scatterfun <- function(xvar, yvar, AAcomp) {
    if(yvar=="ZC") {
      y <- ZCAA(AAcomp)
      ylab <- expression(italic(Z)[C])
    }
    if(yvar=="nH2O") {
      y <- H2OAA(AAcomp)
      ylab <- expression(italic(n)[H[2] * O])
    }
    if(xvar=="GRAVY") {
      x <- GRAVY(AAcomp)
      xlab <- "GRAVY"
    }
    if(xvar=="pI") {
      x <- pI(AAcomp)
      xlab <- "pI"
    }
    # make scatterplot
    smoothScatter(x, y, xlab = xlab, ylab = ylab, colramp = colorRampPalette(c("transparent", blues9)))
    # add linear fit
    mylm <- lm(y ~ x)
    xs <- par("usr")[1:2]
    ys <- predict(mylm, data.frame(x = xs))
    lines(xs, ys, lty = 2, lwd = 2, col = "grey40")
    # add legend with R-squared value
    R2 <- format(round(summary(mylm)$r.squared, 2), nsmall = 2)
    R2txt <- substitute(italic(R)^2 == R2, list(R2 = R2))
    legend("bottomleft", legend = R2txt, bty = "n")
  }
  par(mfrow = c(1, 2))
  # get amino acid compositions of E.coli proteins (UniProt)
  ecoli <- read.csv(system.file("/extdata/organisms/ecoli.csv.xz", package = "JMDplots"), as.is = TRUE)
  scatterfun("GRAVY", "ZC", ecoli)
  label.figure("(a)", cex = 1.8)
  scatterfun("GRAVY", "nH2O", ecoli)
  label.figure("(b)", cex = 1.8)
  if(pdf) {
    dev.off()
    addexif("gradH2O5", "Scatterplots of ZC and nH2O vs GRAVY", "Dick et al. (2020) (preprint)")
  }
}

# nH2O-ZC and GRAVY-pI plots for Baltic Sea and Rodriguez-Brito et al. data 20200421
gradH2O6 <- function(pdf = FALSE) {
  if(pdf) pdf("gradH2O6.pdf", width = 8, height = 8)
  layout(matrix(1:4, nrow = 2))
  par(mar = c(4, 4.5, 2, 1), las = 1, cex = 1.2)
  par(cex.lab = 1.3)

  # Baltic Sea nH2O - ZC
  mout <- ppage("balticsurface", plot.it = FALSE)
  pout <- ppage("balticsurface", H2O = TRUE, plot.it = FALSE)
  pcomp(mout, pout, reorder = FALSE, yline = 3.2, cex.ylab = 1.5, font = 2, labdy = 0.003, labels.at = NA, xlim = c(-0.2, -0.08))
  text(c(-0.126, -0.124), c(0.04, 0.02), c("< 6 PSU", "> 6 PSU"))
  title("Baltic Sea")
  label.figure("(a)", cex = 1.7)
  # Baltic Sea GRAVY - pI
  pcomp(mout, pout, reorder = FALSE, vars = "pIG", yline = 3.2, cex.ylab = 1.5, font = 2, labdy = 0.003, labels.at = NA)
  text(c(7.5, 7.2), c(-0.11, -0.23), c("< 6 PSU", "> 6 PSU"))
  title("Baltic Sea")
  label.figure("(c)", cex = 1.7)

  # Rodriguez-Brito et al. nH2O - ZC
  mout <- ppage("socal", plot.it = FALSE)
  pout <- ppage("socal", H2O = TRUE, plot.it = FALSE)
  pcomp(mout, pout, reorder = FALSE, yline = 3.2, cex.ylab = 1.5, font = 2, labdy = 0.003, labels.at = NA, xlim = c(-0.2, -0.08))
  text(c(-0.155, -0.125, -0.115), c(-0.005, 0.015, 0.027), c("FW", "LS", "MS-HS"))
  title("Rodriguez-Brito et al.")
  label.figure("(b)", cex = 1.7)
  # Rodriguez-Brito et al. GRAVY - pI
  pcomp(mout, pout, reorder = FALSE, vars = "pIG", yline = 3.2, cex.ylab = 1.5, font = 2, labdy = 0.003, labels.at = NA)
  text(c(8.3, 6, 4.5), c(-0.18, -0.19, -0.115), c("FW", "LS", "MS-HS"))
  title("Rodriguez-Brito et al.")
  label.figure("(d)", cex = 1.7)
  if(pdf) {
    dev.off()
    addexif("gradH2O6", "nH2O-ZC and GRAVY-pI plots for Baltic Sea and Rodriguez-Brito et al. data", "Dick et al. (2020) (preprint)")
  }
}

# differential gene and protein expression, time-course and NaCl vs organic solutes 20200420
gradH2O7 <- function(pdf = FALSE) {
  if(pdf) pdf("gradH2O7.pdf", width = 8, height = 4.5)
  layout(matrix(0:11, nrow = 3, byrow = TRUE), widths = c(0.2, 1, 1, 1), heights = c(0.2, 1, 1))
  # add titles
  par(mar = c(0, 0, 0, 0))
  plot.new()
  text(0.57, 0.5, "All compiled datasets for bacteria", font = 2)
  plot.new()
  text(0.57, 0.5, "Time-course experiments", font = 2)
  plot.new()
  text(0.57, 0.5, "NaCl or organic solutes", font = 2)
  plot.new()
  text(0.5, 0.6, "Proteins coded by\ndifferentially expressed genes", srt = 90, font = 2)
  par(mar = c(4, 4, 0.2, 1), mgp = c(2.5, 1, 0))

  # function to plot an arrow partway along a line
  mkarrow <- function(row1, ct, frac = 0.5) {
    # row1 is the starting point
    x1 <- ct$ZC.diff[row1]
    y1 <- ct$nH2O.diff[row1]
    # the next row is the end of the full line (not the arrow)
    x2 <- ct$ZC.diff[row1 + 1]
    y2 <- ct$nH2O.diff[row1 + 1]
    # calculate slope
    m <- (y2 - y1) / (x2 - x1)
    # calculate value of x and y on the line (arrow tip)
    x <- x1 + (x2 - x1) * frac
    y <- y1 + m * (x - x1)
    # calculate value of x0 and y0 (start of arrow - just a short line)
    x0 <- x1 + (x2 - x1) * frac * 0.99
    y0 <- y1 + m * (x0 - x1)
    # draw the arrow
    suppressWarnings(arrows(x0, y0, x, y, length = 0.1, angle = 20, lwd = 2))
  }

  # function to make diff plot with arrows and points
  mkdiff <- function(ct, ndat, pch = 21, ...) {
    # make an empty plot, add lines, then add points
    diffplot(ct, pch = NA, pt.text = NA, contour = FALSE)
    n <- 0
    for(i in 1:length(ndat)) {
      idat <- n + 1:ndat[i]
      lines(ct$ZC.diff[idat], ct$nH2O.diff[idat])
      lapply(head(idat, -1), mkarrow, ct = ct)
      n <- n + ndat[i]
    }
    par(bg = "white")
    diffplot(ct, pch = pch, add = TRUE, contour = FALSE, cex.text = 0.8, ...)
  }

  # plot A: transcriptomes compilation
  osmotic_gene <- read.csv(system.file(paste0("vignettes/osmotic_gene_", getOption("basis"), ".csv"), package = "JMDplots"))
  diffplot(osmotic_gene, pt.text = NA, contour = FALSE, cex = 1.5)
  points(mean(osmotic_gene$ZC.diff), mean(osmotic_gene$nH2O.diff), pch = 19, cex = 1.7)
  legend("topleft", c("dataset", "mean"), pch = c(1, 19), bty = "n")
  label.figure("(a)", cex = 1.8, xfrac = 0.12, yfrac = 1.05)
  # print the p-values
  p.ZC <- round(t.test(osmotic_gene$ZC.down, osmotic_gene$ZC.up, paired = TRUE)$p.value, 3)
  p.nH2O <- round(t.test(osmotic_gene$nH2O.down, osmotic_gene$nH2O.up, paired = TRUE)$p.value, 3)
  print(paste("number of transcriptomics datasets:", nrow(osmotic_gene)))
  print(paste("p-values for transcriptomics:", p.ZC, "(ZC),", p.nH2O, "(nH2O)"))

  # plot C: transcriptomics: increasing time
  Ttime <- list(
    KKG = c("KKG+14_Gene_30min", "KKG+14_Gene_80min", "KKG+14_Gene_310min"),
    SLM = c("SLM+14_5", "SLM+14_30", "SLM+14_60"),
    FRH = c("FRH+15_NaCl_1h", "FRH+15_NaCl_6h", "FRH+15_NaCl_24h"),
    HLL = c("HLL17_45min", "HLL17_14h"),
    QHT = c("QHT+13_Gene.24.h", "QHT+13_Gene.48.h", "QHT+13_Gene.72.h")
  )
  comptab <- osmotic_gene[match(unlist(Ttime), osmotic_gene$dataset), ]
  ndat <- sapply(Ttime, length)
  mkdiff(comptab, ndat)
  label.figure("(c)", cex = 1.8, xfrac = 0.12, yfrac = 1.05)

  # plot E: transcriptomics: NaCl or organic solutes
  Tsolute <- list(
    KSA = c("KSA+02_NaCl", "KSA+02_sorbitol"),
    HZP = c("HZP+05_HSS", "HZP+05_HOS"),
    KLB = c("KLB+15_trans-NaCl", "KLB+15_trans-suc"),
    FRH_1 = c("FRH+15_NaCl_1h", "FRH+15_glycerol_1h"),
    FRH_6 = c("FRH+15_NaCl_6h", "FRH+15_glycerol_6h"),
    SBB = c("SBB+09_NaCl", "SBB+09_Sucrose"),
    WGB = c("WGB+13_N", "WGB+13_U")
  )
  comptab <- osmotic_gene[match(unlist(Tsolute), osmotic_gene$dataset), ]
  ndat <- sapply(Tsolute, length)
  mkdiff(comptab, ndat, pt.text = LETTERS[1:sum(ndat)], pch = c(21, 22))
  legend("topleft", c("NaCl", "organic"), pch = c(21, 22), bty = "n")
  label.figure("(e)", cex = 1.8, xfrac = 0.12, yfrac = 1.05)

  par(mar = c(0, 0, 0, 0))
  par(xpd = NA)
  plot.new()
  text(0.5, 0.6, "Differentially expressed proteins", srt = 90, font = 2)
  par(mar = c(4, 4, 0.2, 1), mgp = c(2.5, 1, 0))
  par(xpd = FALSE)

  # plot B: proteomes compilation
  osmotic_bact <- read.csv(system.file(paste0("vignettes/osmotic_bact_", getOption("basis"), ".csv"), package = "canprot"))
  diffplot(osmotic_bact, pt.text = NA, contour = FALSE, cex = 1.5)
  points(mean(osmotic_bact$ZC.diff), mean(osmotic_bact$nH2O.diff), pch = 19, cex = 1.7)
  legend("topleft", c("dataset", "mean"), pch = c(1, 19), bty = "n")
  label.figure("(b)", cex = 1.8, xfrac = 0.12, yfrac = 1.05)
  # print the p-values
  p.ZC <- round(t.test(osmotic_bact$ZC.down, osmotic_bact$ZC.up, paired = TRUE)$p.value, 3)
  p.nH2O <- round(t.test(osmotic_bact$nH2O.down, osmotic_bact$nH2O.up, paired = TRUE)$p.value, 3)
  print(paste("number of proteomics datasets:", nrow(osmotic_bact)))
  print(paste("p-values for proteomics:", p.ZC, "(ZC),", p.nH2O, "(nH2O)"))
  # print total number of studies
  osmotic_gene_studies <- sapply(strsplit(osmotic_gene$dataset, "_"), "[", 1)
  osmotic_bact_studies <- sapply(strsplit(osmotic_bact$dataset, "_"), "[", 1)
  nstudies <- length(unique(c(osmotic_gene_studies, osmotic_bact_studies)))
  print(paste("total number of studies (transcriptomics and proteomics):", nstudies))

  # plot D: proteomics: increasing time
  Ptime <- list(
    KKG = c("KKG+14_Protein_30min", "KKG+14_Protein_80min", "KKG+14_Protein_310min"),
    QHT = c("QHT+13_Protein.24.h", "QHT+13_Protein.48.h")
  )
  comptab <- osmotic_bact[match(unlist(Ptime), osmotic_bact$dataset), ]
  ndat <- sapply(Ptime, length)
  mkdiff(comptab, ndat, pt.text = c("a", "b", "c", "l", "m"))
  label.figure("(d)", cex = 1.8, xfrac = 0.12, yfrac = 1.05)

  # plot F: proteomics: NaCl or organic solutes
  Psolute <- list(
    KLB = c("KLB+15_prot-NaCl", "KLB+15_prot-suc"),
    SKV = c("SKV+16_Osmotic.stress.glucose_LB", "SKV+16_Glucose_LB")
  )
  comptab <- osmotic_bact[match(unlist(Psolute), osmotic_bact$dataset), ]
  ndat <- sapply(Psolute, length)
  mkdiff(comptab, ndat, pt.text = c("O", "P", "Q", "R"), pch = c(21, 22))
  label.figure("(f)", cex = 1.8, xfrac = 0.12, yfrac = 1.05)

  if(pdf) {
    dev.off()
    addexif("gradH2O7", "differential gene and protein expression, time-course and NaCl vs organic solutes", "Dick et al. (2020) (preprint)")
  }
}

# differentially expressed proteins in halophiles under hypo- and hyperosmotic stress
# adapted from canprot/hyperosmotic.Rmd 20190717-20191007
# add GRAVY and pI plot 20191028
# use different symbols for eukaryotes and add halophilic bacteria and archaea 20191102-20191103
# separate plots for halophiles and non-halophiles (no eukaryotes) 20200424
gradH2O8 <- function(pdf = FALSE) {
  if(pdf) pdf("gradH2O8.pdf", width = 5, height = 4)
  layout(matrix(c(1, 1, 2, 2, 0, 3, 3, 0), nrow = 2, byrow = 2))
  par(mar = c(3.2, 3.2, 2, 1))
  par(mgp = c(1.9, 0.6, 0))
  par(cex.lab = 1.2)

  # get proteomic data for halophiles
  osmotic_halo <- read.csv(system.file(paste0("vignettes/osmotic_halo_", getOption("basis"), ".csv"), package = "canprot"))
  # use triangles for hyperosmotic, squares for hypoosmotic
  pch <- ifelse(grepl("hypoosmotic", osmotic_halo$tags), 0, 2)
  # make nH2O-ZC plot
  diffplot(osmotic_halo, pch = pch, contour = FALSE, cex.text = 0.75, labtext = NA)
  legend("topleft", c("hypoosmotic", "hyperosmotic"), pch = c(0, 2), bty = "n")
  title("Halophile proteomics")
  label.figure("(a)", cex = 1.7, yfrac = 0.93)
  # make GRAVY-pI plot
  diffplot(osmotic_halo, vars = c("pI", "GRAVY"), pch = pch, contour = FALSE, cex.text = 0.75, labtext = NA)
  legend("topleft", c("hypoosmotic", "hyperosmotic"), pch = c(0, 2), bty = "n")
  title("Halophile proteomics")
  label.figure("(b)", cex = 1.7, yfrac = 0.93)

  # get proteomic data for non-halophiles in hyperosmotic stress
  osmotic_bact <- read.csv(system.file(paste0("vignettes/osmotic_bact_", getOption("basis"), ".csv"), package = "canprot"))
  # combine halophile and non-halophile data for hyperosmotic stress
  osmotic_halo <- osmotic_halo[osmotic_halo$tags!="hypoosmotic", ]
  alldat <- rbind(osmotic_bact, osmotic_halo)
  pch <- c(rep(1, nrow(osmotic_bact)), rep(2, nrow(osmotic_halo)))
  # make GRAVY-pI plot
  diffplot(alldat, vars = c("pI", "GRAVY"), pt.text = NA, contour = FALSE, pch = pch, labtext = NA, cex = 1.5)
  points(mean(alldat$pI.diff), mean(alldat$GRAVY.diff), pch = 19, cex = 1.7)
  legend("topleft", c("non-halophiles", "halophiles", "mean"), pch = c(1, 2, 19), bty = "n")
  title("All hyperosmotic proteomic data")
  label.figure("(c)", cex = 1.7, yfrac = 0.93)
  # print the p-values
  p.pI <- round(t.test(alldat$pI.down, alldat$pI.up, paired = TRUE)$p.value, 3)
  p.GRAVY <- round(t.test(alldat$GRAVY.down, alldat$GRAVY.up, paired = TRUE)$p.value, 3)
  print(paste("p-values for all hyperosmotic proteomic data:", p.pI, "(pI),", p.GRAVY, "(GRAVY)"))

  if(pdf) {
    dev.off()
    addexif("gradH2O8", "Differentially expressed proteins in halophiles under hypo- and hyperosmotic stress", "Dick et al. (2020) (preprint)")
  }
}

# calculate ZC and nH2O of proteomes encoding different Nif homologs (Poudel et al., 2018)
# 20191014
NifProteomes <- function() {
  # read file with Nif genome classifications and taxids
  Niffile <- system.file("extdata/gradH2O/Nif_homolog_genomes.csv", package = "JMDplots")
  Nif <- read.csv(Niffile, as.is = TRUE)
  # drop NA taxids
  Nif <- Nif[!is.na(Nif$taxid), ]
  # read refseq data
  RSfile <- system.file("extdata/refseq/protein_refseq.csv.xz", package = "JMDplots")
  refseq <- read.csv(RSfile, as.is = TRUE)
  # the Nif types, arranged from anaerobic to aerobic
  types <- c("Nif-D", "Nif-C", "Nif-B", "Nif-A")
  # assemble the compositional metrics
  ZC <- ZC.SD <- nH2O <- nH2O.SD <- GRAVY <- GRAVY.SD <- pI <- pI.SD <- numeric()
  for(type in types) {
    # get the taxids for genomes with this type of Nif
    iNif <- Nif$Type == type
    taxid <- Nif$taxid[iNif]
    # remove duplicated taxids 20191018
    taxid <- taxid[!duplicated(taxid)]
    # get the row number in the refseq data frame
    irefseq <- match(taxid, refseq$organism)
    # include only organisms with at least 1000 protein sequences
    i1000 <- refseq$chains[irefseq] >= 1000
    irefseq <- irefseq[i1000]
    #print(paste(type, "represented by", length(irefseq), "nonredundant genomes with at least 1000 protein sequences"))
    # get the amino acid composition from refseq
    AAcomp <- refseq[irefseq, ]
    # calculate ZC and nH2O
    ZC <- c(ZC, mean(ZCAA(AAcomp)))
    ZC.SD <- c(ZC.SD, sd(ZCAA(AAcomp)))
    nH2O <- c(nH2O, mean(H2OAA(AAcomp)))
    nH2O.SD <- c(nH2O.SD, sd(H2OAA(AAcomp)))
    GRAVY <- c(GRAVY, mean(GRAVY(AAcomp)))
    GRAVY.SD <- c(GRAVY.SD, sd(GRAVY(AAcomp)))
    pI <- c(pI, mean(pI(AAcomp)))
    pI.SD <- c(pI.SD, sd(pI(AAcomp)))
  }
  # return values
  list(types = types, ZC = ZC, ZC.SD = ZC.SD, nH2O = nH2O, nH2O.SD = nH2O.SD, GRAVY = GRAVY, GRAVY.SD = GRAVY.SD, pI = pI, pI.SD = pI.SD)
}

############################
### UNEXPORTED FUNCTIONS ###
############################

# function to add convex hulls 20191007
hullfun <- function(mout, pout, istudy, basecol, group = NULL, vars = "ZC") {
  x <- y <- numeric()
  for(ist in istudy) {
    thisx <- mout[[ist]]$AA
    thisy <- pout[[ist]]$AA
    if(vars == "GRAVY") thisx <- mout[[ist]]$GRAVY
    if(vars == "pI") thisx <- mout[[ist]]$pI
    if(vars == "pIG") {
      thisx <- mout[[ist]]$pI
      thisy <- pout[[ist]]$GRAVY
    }
    if(!is.null(group)) {
      thisx <- thisx[mout[[ist]]$group %in% group]
      thisy <- thisy[mout[[ist]]$group %in% group]
    }
    x <- c(x, thisx)
    y <- c(y, thisy)
  }
  i <- chull(x, y)
  r <- as.numeric(col2rgb(basecol))
  col <- rgb(r[1], r[2], r[3], 80, maxColorValue=255)
  polygon(x[i], y[i], col = col, border = NA)
}
