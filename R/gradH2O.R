# JMDplots/gradH2O.R
# make plots for Goldschmidt poster 20190711
# added to JMDplots starting 20190930

#source("mplot.R")

## load required packages
#library(canprot)
#library(CHNOSZ)

# basis species comparison 20190713 / relative water content 20191005
gradH2O1 <- function(pdf = FALSE) {

  # check if E. coli data is available
  if(!exists("ecoli")) stop("please load amino acid compositions of E. coli proteins with data(ecoli)")

  # set up figure
  if(pdf) pdf("gradH2O1.pdf", width = 6, height = 6)
  par(mfrow = c(3, 3))
  par(mar = c(3.2, 3.2, 2.5, 1))
  par(mgp = c(2.2, 0.7, 0))
  par(las = 1)

  # define axis labels
  nH2Olab <- expression(italic(n)[H[2] * O])
  nO2lab <- expression(italic(n)[O[2]])
  ZClab <- expression(italic(Z)[C])

  # function to plot linear model
  lmfun <- function(ZC, nH2O, legend.x = NULL, ...) {
    mylm <- lm(nH2O ~ ZC)
    lines(c(-1, 1), predict(mylm, data.frame(ZC = c(-1, 1))), ...)
    # add R-squared text
    if(!is.null(legend.x)) {
      R2 <- round(summary(mylm)$r.squared, 2)
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
    lmfun(ZC, y, legend.x)
  }

  scatterfun <- function(ZC, nH2O, main, lab) {
    smoothScatter(ZC, nH2O, xlab = ZClab, ylab = nH2Olab, colramp = colorRampPalette(c("transparent", blues9)))
    lmfun(ZC, nH2O, lty = 2, lwd = 2, col = "grey40")
    mainlab(main, lab)
  }

  # get names of amino acids
  aa <- aminoacids("")
  # calculate ZC of the amino acids
  ZC.aa <- ZC(info(aa, "aq"))

  # get amino acid composition of human proteins (UniProt)
  human <- get("human_base", canprot)
  pf.human <- protein.formula(human)
  ZC.human <- ZC(pf.human)

  # get amino acid compositions of E. coli proteins (UniProt)
  pf.ecoli <- protein.formula(ecoli)
  ZC.ecoli <- ZC(pf.ecoli)

  # calculate nH2O and nO2 of amino acids with CHNOS basis
  basis("CHNOS")
  species(aa)

  # plot 1: nH2O-ZC of amino acids (CHNOS)
  aaplot(ZC.aa, species()$H2O, nH2Olab, "bottomleft", "Amino acids (CHNOS)", "A")

  # plot 2: nO2-ZC of amino acids (CHNOS)
  aaplot(ZC.aa, species()$O2, nO2lab, "bottomright", "Amino acids (CHNOS)", "B")

  # calculate nH2O and nO2 with QEC basis
  basis("QEC")
  species(aa)

  # plot 3: nO2-ZC of amino acids (QEC)
  aaplot(ZC.aa, species()$O2, nO2lab, "bottomright", "Amino acids (QEC)", "C")

  # plot 4: nH2O-ZC of amino acids (QEC)
  # save the residuals here
  rQEC <- aaplot(ZC.aa, species()$H2O, nH2Olab, "bottomright", "Amino acids (QEC)", "D")
  names(rQEC) <- aminoacids(3)

  # plot 5: nH2O-ZC of human proteins (QEC)
  nH2O.human <- H2OAA(human, basis = "QEC")
  scatterfun(ZC.human, nH2O.human, "Human proteins (QEC)", "E")

  # plot 6: nH2O-ZC of E. coli proteins (QEC)
  nH2O.ecoli <- H2OAA(ecoli, basis = "QEC")
  scatterfun(ZC.ecoli, nH2O.ecoli, "E. coli proteins (QEC)", "F")

  # plot 7: nH2O-ZC of amino acids (rQEC)
  aaplot(ZC.aa, rQEC, nH2Olab, NULL, "Amino acids (rQEC)", "G")

  # plot 8: nH2O-ZC of human proteins (rQEC)
  nH2O.human <- H2OAA(human, basis = "rQEC")
  scatterfun(ZC.human, nH2O.human, "Human proteins (rQEC)", "H")

  # plot 9: nH2O-ZC of E. coli proteins (rQEC)
  nH2O.ecoli <- H2OAA(ecoli, basis = "rQEC")
  scatterfun(ZC.ecoli, nH2O.ecoli, "E. coli proteins (rQEC)", "I")

  if(pdf) invisible(dev.off())
  ## output value of rQEC for checking code of H2OAA()
  #rQEC
}

# ZC for selected redox gradients 20190715
gradH2O2 <- function(pdf = FALSE) {
  if(pdf) pdf("gradH2O2.pdf", width = 6, height = 2.3)
  par(mfrow = c(1, 3))
  par(mar = c(4, 4, 4, 1), mgp = c(2.5, 0.7, 0))
  mplot("Bison_Pool", "IMG_MGP", add.label = FALSE, plottype = "#FF000030", col = "red")
  label.figure("A", cex = 2, yfrac = 0.9)
  mplot("Diffuse_Vents", "SRA_MGP", add.label = FALSE, plottype = "#FF000030", col = "red")
  label.figure("B", cex = 2, yfrac = 0.9)
  mplot("Guerrero_Negro", "IMG_MGP", add.label = FALSE, plottype = "#FF000030", col = "red")
  label.figure("C", cex = 2, yfrac = 0.9)
  if(pdf) invisible(dev.off())
}

# nH2O-ZC scatterplots for redox gradients and Baltic Sea 20190713
gradH2O3 <- function(pdf = FALSE) {
  if(pdf) pdf("gradH2O3.pdf", width = 13, height = 5.6)
  par(mfrow = c(1, 3))
  par(mar = c(4, 4.5, 1, 1), las = 1, cex = 1.2)
  # compare ZC and nH2O of proteins in datasets from gradox paper
  mgradox <- ppage("gradoxGS", plot.it = FALSE)
  pgradox <- ppage("gradoxGS", H2O = TRUE, plot.it = FALSE)
  pcomp(mgradox, pgradox, type = "both", reorder = FALSE, yline = 3.5)
  # overlay general trend from human proteins (printed in basis_comparison())
  lines(c(-1, 1), c(-0.61, -0.98), lty = 2, lwd = 3, col = "grey40")
  legend("topleft", c("redox", "gradients"), bty = "n", text.font = 2)
  # add legend for environment type
  legend("topright", c("vent fluids", "plume", "seawater", "hot spring", "phototrophic", "mat > 3 mm", "mat 1-3 mm"),
         pch = c(19, 19, 19, 15, 15, 15, 15), col = c("red", "purple1", "purple1", "orange", "orange", "green3", "green3"), bty = "n")
  # overlay symbols for seawater, phototrophic and mat surface
  legend(-0.1627, -0.7176, lty=0, lwd=0, bty="n", pt.cex=1.6, pt.lwd=1,
           pch = c(1, 1, 1, 0, 0, 0, 0),
           legend=c("", "", "", "", "", "", ""),
           col=c(NA, NA, "purple1", NA, "green3", NA, "green3"))
  label.figure("A", cex = 2, xfrac = 0.035)
  # compare ZC and nH2O of proteins in Baltic Sea surface
  mbaltics <- ppage("balticsurface", plot.it = FALSE)
  pbaltics <- ppage("balticsurface", H2O = TRUE, plot.it = FALSE)
  pcomp(mbaltics, pbaltics, type = "both", reorder = FALSE, yline = 3.5)
  lines(c(-1, 1), c(-0.61, -0.98), lty = 2, lwd = 3, col = "grey40")
  legend("topleft", c("Baltic Sea", "surface"), bty = "n", text.font = 2)
  label.figure("B", cex = 2, xfrac = 0.035)
  # compare ZC and nH2O of proteins in Baltic Sea deep
  mbalticd <- ppage("balticdeep", plot.it = FALSE)
  pbalticd <- ppage("balticdeep", H2O = TRUE, plot.it = FALSE)
  pcomp(mbalticd, pbalticd, type = "both", reorder = FALSE, yline = 3.5)
  lines(c(-1, 1), c(-0.61, -0.98), lty = 2, lwd = 3, col = "grey40")
  legend("topleft", c("Baltic Sea", "10-20 m"), bty = "n", text.font = 2)
  # add legend for particle size
  legend("topright", legend = as.expression(c(quote("0.1-0.8"~mu*m))),
         pch = c(17), col = c("black"), bty = "n")
  label.figure("C", cex = 2, xfrac = 0.035)
  if(pdf) invisible(dev.off())
}

# nH2O for Baltic Sea size fractions 20190715
gradH2O4 <- function(pdf = FALSE) {
  if(pdf) pdf("gradH2O4.pdf", width = 6, height = 2.5)
  par(mfrow = c(1, 3))
  par(mar = c(5, 4, 1, 1), mgp = c(3, 0.7, 0), las = 1)

  mplot("Baltic_Sea-0.1s", "iMicrobe_MGP", H2O = TRUE, plottype = "#FF000030", col = "red", add.title = FALSE, ylim = c(-0.8, -0.7), yline = 2.7)
  mplot("Baltic_Sea-0.1s", "SRA_MTP", H2O = TRUE, plottype = "#0000FF30", col = "blue", add.title = FALSE, add = TRUE, pch = 1)
  legend("topright", legend = quote("0.1-0.8"~mu*m), bty = "n")
  label.figure("A", cex = 2, xfrac = 0.04)

  mplot("Baltic_Sea-0.8s", "iMicrobe_MGP", H2O = TRUE, plottype = "#FF000030", col = "red", add.title = FALSE, ylim = c(-0.8, -0.7), yline = 2.7)
  mplot("Baltic_Sea-0.8s", "SRA_MTP", H2O = TRUE, plottype = "#0000FF30", col = "blue", add.title = FALSE, add = TRUE, pch = 1)
  legend("topright", legend = quote("0.8-3.0"~mu*m), bty = "n")
  label.figure("B", cex = 2, xfrac = 0.035)

  mplot("Baltic_Sea-3.0s", "iMicrobe_MGP", H2O = TRUE, plottype = "#FF000030", col = "red", add.title = FALSE, ylim = c(-0.8, -0.7), yline = 2.7)
  mplot("Baltic_Sea-3.0s", "SRA_MTP", H2O = TRUE, plottype = "#0000FF30", col = "blue", add.title = FALSE, add = TRUE, pch = 1)
  legend("topright", legend = quote("3.0-200"~mu*m), bty = "n")
  label.figure("C", cex = 2, xfrac = 0.035)

  if(pdf) invisible(dev.off())
}

