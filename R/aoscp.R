# source code to make figures in the paper
# Average oxidation state of carbon in proteins
# by Jeffrey M. Dick, 2014-07-30
# added to JMDplots package 20191018

# package dependencies:
# png (used in yeast())
# plotrix (using in potential())

# ZC of amino acids vs ZC of RNA codons and hydropathy index of amino acids
aoscp1 <- function(pdf = FALSE) {
  # start plot
  if(pdf) pdf("aoscp1.pdf", width=10, height=5, family="Times")
  par(mfrow=c(1, 2))

  ## the codon stuff
  file <- system.file("/extdata/aoscp/codons.csv", package = "JMDplots")
  dat <- read.csv(file, as.is=TRUE)
  # take only first two nucleobases, then remove duplicates
  dat$codon <- sapply(dat$codon, substr, 1, 2)
  dat <- unique(dat)
  # tabulate counts of first two nucleobases and get them into matrix form
  codon <- sapply(sapply(dat$codon, strsplit, ""), table)
  codon <- sapply(codon, as.matrix)
  codon <- sapply(codon, t)
  # rbind them, filling in missing columns
  codon <- lapply(codon, function(x) {
    this <- t(matrix(numeric(4), dimnames=list(c("A", "C", "G", "U"))))
    icol <- match(colnames(x), colnames(this))
    this[, icol] <- x
    this
  })
  codon <- do.call("rbind", codon)
  # calculate ZC of the codons
  ZC.codon <- ZC(nucleic.formula(codon))
  # calculate ZC of the amino acids
  aa.formula <- info(info(dat$aminoacid))$formula
  # the "stop" are NA - give them an undefined ZC
  aa.formula[is.na(aa.formula)] <- "C0"
  ZC.aminoacid <- ZC(aa.formula)
  # make plot
  plot(ZC.codon, ZC.aminoacid, cex=2.5,
    xlab=expression(italic(Z)[C]~"in first two bases of RNA codon"),
    ylab=expression(italic(Z)[C]~"in amino acid"))
  # add labels for amino acids
  aa <- aminoacids(1, which=dat$aminoacid)
  iF <- which(aa=="F")
  iP <- which(aa=="P")
  iE <- which(aa=="E")
  iFPE <- c(iF, iP, iE)
  text(ZC.codon[-iFPE], ZC.aminoacid[-iFPE], aa[-iFPE])
  # plot F, P and E to not overlap with other symbols
  text(ZC.codon[iF]+0.04, ZC.aminoacid[iF]-0.03, aa[iF])
  text(ZC.codon[iP]+0.04, ZC.aminoacid[iP]+0.03, aa[iP])
  text(ZC.codon[iE]+0.04, ZC.aminoacid[iE], aa[iE])
  # identify codons for amino acids that appear multiple times
  dy <- -0.12
  cex <- 0.8
  # the serine and arginine codons that have different first two bases
  iUC <- match("UC", dat$codon)
  text(ZC.codon[iUC], ZC.aminoacid[iUC]+dy, "UC", cex=cex)
  iGG <- match("GG", dat$codon)
  text(ZC.codon[iGG], ZC.aminoacid[iGG]+dy, "GG", cex=cex)
  # the serine and arginine codons that share the first two bases
  iAG_Ser <- which(dat$codon=="AG" & dat$aminoacid=="serine")
  text(ZC.codon[iAG_Ser], ZC.aminoacid[iAG_Ser]+dy, "AG", cex=cex)
  iAG_Arg <- which(dat$codon=="AG" & dat$aminoacid=="arginine")
  text(ZC.codon[iAG_Arg], ZC.aminoacid[iAG_Arg]+dy, "AG", cex=cex)
  # label plot
  label.figure("a", paren = TRUE, italic = TRUE)

  ## the hydropathy stuff
  aa <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")
  hydropathy <- c(1.8, 2.5, -3.5, -3.5, 2.8, -0.4, -3.2, 4.5, -3.9, 3.8, 1.9, -3.5, -1.6, -3.5, -4.5, -0.8, -0.7, 4.2, -0.9, -1.3)
  # the properties of the amino acids in the thermodynamic database
  aa.name <- aminoacids("", aa)
  aa.dat <- info(info(aa.name))
  # get amino acid formula, ZC
  aa.form <- aa.dat$formula
  aa.ZC <- ZC(aa.form)
  # plot ZC vs hydropathy
  plot(hydropathy, aa.ZC, cex=2.5,
    xlab="hydropathy index of amino acid",
    ylab=expression(italic(Z)[C]~"in amino acid")
  )
  # add labels - drop N and Q for clarity
  iNQ <- aa %in% c("N", "Q")
  text(hydropathy[!iNQ], aa.ZC[!iNQ], aa[!iNQ])
  # label plot
  label.figure("b", paren = TRUE, italic = TRUE)

  # done!
  if(pdf) invisible(dev.off())
  ## we can print some statistics
  #print(summary(lm(ZC.aminoacid ~ ZC.codon)))
  #print(summary(lm(aa.ZC ~ hydropathy)))
}

