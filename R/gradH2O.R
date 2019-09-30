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
      print(nH2O.pred)
    }
  }
  if(pdf) dev.off()
}

