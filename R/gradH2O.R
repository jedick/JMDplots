# JMDplots/gradH2O.R
# make plots for Goldschmidt poster 20190711
# added to JMDplots starting 20190930

#source("mplot.R")

## load required packages
#library(canprot)
#library(CHNOSZ)

# basis species comparison, from canprot/vignettes/basis_comparison.Rmd 20190713
gradH2O1 <- function(pdf = FALSE) {
  if(pdf) pdf("gradH2O1.pdf", width = 10, height = 10)
  # define labels used in the plot
  nH2Olab <- expression(italic(n)[H[2] * O])
  nO2lab <- expression(italic(n)[O[2]])
  ZClab <- expression(italic(Z)[C])
  QEClab <- CHNOSZ::syslab(c("glutamine", "glutamic acid", "cysteine", "H2O", "O2"))
  CHNOSlab <- CHNOSZ::syslab(c("CO2", "NH3", "H2S", "H2O", "O2"))
  # get amino acid composition of human proteins (UniProt)
  aa <- get("human_base", canprot)
  protein.formula <- CHNOSZ::protein.formula(aa)
  ZC <- CHNOSZ::ZC(protein.formula)
  # set up figure
  par(mfrow = c(2, 2))
  par(mar = c(3.2, 3.2, 2.5, 1))
  par(cex = 1.7)
  par(mgp = c(2.2, 0.7, 0))
  par(las = 1)
  # plot the per-residue compositions of the proteins projected into different sets of basis species
  for(basis in c("QEC", "CHNOS")) {
    CHNOSZ::basis(basis)
    protein.basis <- CHNOSZ::protein.basis(aa)
    protein.length <- CHNOSZ::protein.length(aa)
    residue.basis <- protein.basis / protein.length
    # nO2 vs ZC
    smoothScatter(ZC, residue.basis[, "O2"], xlab = ZClab, ylab = nO2lab, colramp = colorRampPalette(c("transparent", blues9)))
    if(basis=="QEC") figlab <- "A" else figlab <- "C"
    label.figure(figlab, yfrac = 0.88, cex = 1.3)
    # nH2O vs ZC
    smoothScatter(ZC, residue.basis[, "H2O"], xlab = ZClab, ylab = nH2Olab, colramp = colorRampPalette(c("transparent", blues9)))
    if(basis=="QEC") figlab <- "B" else figlab <- "D"
    label.figure(figlab, yfrac = 0.88, cex = 1.3)
    # add titles
    if(basis=="QEC") mtext(QEClab, outer = TRUE, cex = 1.7, line = -1.2)
    if(basis=="CHNOS") mtext(CHNOSlab, outer = TRUE, cex = 1.7, line = -15.7)
    # add linear fit for QEC basis 20190713
    if(basis=="QEC") {
      lm.QEC <- lm(residue.basis[, "H2O"] ~ ZC)
      nH2O.pred <- predict.lm(lm.QEC, data.frame(ZC = c(-1, 1)))
      lines(c(-1, 1), nH2O.pred, lty = 2, lwd = 3, col = "grey40")
    }
  }
  if(pdf) invisible(dev.off())
}

# plot ZC for selected redox gradients 20190715
gradH2O2 <- function(pdf = FALSE) {
  if(pdf) pdf("gradH2O2.pdf", width = 6, height = 2.3)
  par(mfrow = c(1, 3))
  par(mar = c(4, 4, 4, 1), mgp = c(2.5, 0.7, 0))
  mplot("Bison_Pool", "IMG_MGP", add.label = FALSE)
  label.figure("A", cex = 2, yfrac = 0.9)
  mplot("Diffuse_Vents", "SRA_MGP", add.label = FALSE)
  label.figure("B", cex = 2, yfrac = 0.9)
  mplot("Guerrero_Negro", "IMG_MGP", add.label = FALSE)
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

