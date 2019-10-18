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

# histograms of ZC of all human proteins and human membrane proteins
aoscp2 <- function(pdf = FALSE) {
  # read the HUMAN data file
  file <- system.file("/extdata/aoscp/ZC_HUMAN.csv.xz", package = "JMDplots")
  HUMAN <- read.csv(file)
  # only take sequences containing at least 50 amino acids
  i50.HUMAN <- HUMAN$length >= 50
  ZC.HUMAN <- HUMAN$ZC[i50.HUMAN]
  # read the membrane data file
  file <- system.file("/extdata/aoscp/ZC_membrane.csv", package = "JMDplots")
  membrane <- read.csv(file)
  # only take sequences containing at least 50 amino acids
  i50.membrane <- membrane$length >= 50
  ZC.membrane <- membrane$ZC[i50.membrane]
  # setup plot
  if(pdf) pdf("aoscp2.pdf", width=6.6, height=4.4, family="Times")
  layout(matrix(1:4, 2, 2, byrow=TRUE), widths=c(3, 1.5))
  par(mgp=c(2.3, 1, 0))
  par(mar=c(4, 4, 1, 1))
  xlim <- range(c(ZC.HUMAN, ZC.membrane))
  # compute breakpoints: start with a range beyond all ZC values
  breaks <- seq(-1, 1, by=0.02)
  # which of these breaks are within the range for membrane ZC
  i.membrane <- which(breaks > min(ZC.membrane) & breaks < max(ZC.membrane))
  # grow it by one break to include the extreme values
  i.membrane <- c(i.membrane[1]-1, i.membrane, rev(i.membrane)[1]+1)
  breaks.membrane <- breaks[i.membrane]
  # repeat for HUMAN
  i.HUMAN <- which(breaks > min(ZC.HUMAN) & breaks < max(ZC.HUMAN))
  i.HUMAN <- c(i.HUMAN[1]-1, i.HUMAN, rev(i.HUMAN)[1]+1)
  breaks.HUMAN <- breaks[i.HUMAN]
  # HUMAN plot: histogram
  ZC.lab <- expression(italic(Z)[C])
  hist(ZC.HUMAN, breaks=breaks.HUMAN, xlim=xlim, xlab=ZC.lab, main="")
  legend("topright", legend="all human proteins        ", bty="n")
  label.figure("a", paren = TRUE, italic = TRUE)
  # HUMAN QQ plot
  # to reduce file size, plot only ~10% of the points between the 5% and 95% quantiles
  len <- length(ZC.HUMAN)
  # indices for 10% of all the points
  isamp <- sample(sample(len), round(len/10))
  # ensure all points beyond the 5% and 95% quantiles are included
  isamp <- unique(c(isamp, 1:round(len/20), round(19*len/20):len))
  # calculate the quantiles (full data set), but don't plot them
  qqn <- qqnorm(ZC.HUMAN, plot.it=FALSE)
  # now make the plot with selection of points
  plot(sort(qqn$x)[isamp], sort(qqn$y)[isamp], xlab="Theoretical Quantiles", ylab=ZC.lab, type="p", cex=0.15)
  # draw line through 1st and 3rd quartiles
  probs <- c(25, 75)/100
  qqline(ZC.HUMAN, probs=probs)
  points(quantile(qqn$x, probs), quantile(qqn$y, probs), pch=3, cex=1.5)
  label.figure("b", paren = TRUE, italic = TRUE)
  # membrane plot: histogram
  hist(ZC.membrane, breaks=breaks.membrane, xlim=xlim, xlab=ZC.lab, main="")
  legend("topright", legend="human membrane proteins", bty="n")
  label.figure("c", paren = TRUE, italic = TRUE)
  # membrane QQ plot
  len <- length(ZC.membrane)
  isamp <- sample(sample(len), round(len/10))
  isamp <- unique(c(isamp, 1:round(len/20), round(19*len/20):len))
  qqn <- qqnorm(ZC.membrane, plot.it=FALSE)
  plot(sort(qqn$x)[isamp], sort(qqn$y)[isamp], xlab="Theoretical Quantiles", ylab=ZC.lab, type="p", cex=0.15)
  qqline(ZC.membrane, probs=probs)
  points(quantile(qqn$x, probs), quantile(qqn$y, probs), pch=3, cex=1.5)
  label.figure("d", paren = TRUE, italic = TRUE)
  if(pdf) invisible(dev.off())
  ## run a t-test
  #print(t.test(ZC.HUMAN, ZC.membrane))
}

