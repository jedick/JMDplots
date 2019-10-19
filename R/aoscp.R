# source code to make figures in the paper
# Average oxidation state of carbon in proteins
# by Jeffrey M. Dick, 2014-07-30
# added to JMDplots package 20191018

# package dependencies:
# png (used in aoscp3())
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

# draw yeast cell color-coded with median ZC of proteins in different locations
aoscp3 <- function(png=FALSE, outline=FALSE) {
  # set 'png' to TRUE to make the base plot (no labels)
  # set 'outline' to TRUE to skip plotting the cell components
  # (in order to quickly test outline, labels and color bar)
  # ZC values from reference proteome based on SGD
  ZC <- list(
    bud=yeast.ZC("B"),
    cytoplasm=yeast.ZC("C"),
    bud.neck=yeast.ZC("K"),
    mitochondrion=yeast.ZC("M"),
    vacuole=yeast.ZC("V"),
    nucleus=yeast.ZC("N"),
    nucleolus=yeast.ZC("U"),
    ER=yeast.ZC("E"),
    Golgi=yeast.ZC("G"),
    extracellular=yeast.ZC("X"),
    vacuolar.membrane=yeast.ZC("O"),
    ER.membrane=yeast.ZC("D"),
    plasma.membrane=yeast.ZC("P"),
    nuclear.membrane=yeast.ZC("L"),
    inner.membrane=yeast.ZC("R"),
    outer.membrane=yeast.ZC("T"),
    Golgi.membrane=yeast.ZC("A")
  )
  # get the median for plotting; also print mean and n for table in paper
  ZC.mean <- sapply(ZC, mean)
  ZC.median <- sapply(ZC, median)
  #for(i in 1:length(ZC)) print(paste(names(ZC)[i], "n", length(ZC[[i]]),
  #  "median", round(ZC.median[i], 3), "mean", round(ZC.mean[i], 3)))
  # convert ZC median in a given min/max range to the interval [0, 1]
  ZCmin <- -0.21
  ZCmax <- -0.095
  ZC01 <- (ZC.median-ZCmin)/(ZCmax-ZCmin)
  # adjust bias so cytoplasm is close to white
  cR <- colorRamp(c("black", "red", "white", "blue"), bias=1.085)
  col <- cR(ZC01)
  # round the values so cytoplasm is 255 255 255
  col <- round(col)
  # to make base plot (PNG)
  if(png) {
    png("aoscp3.png", width=1384, height=785)
    plot.window(c(0, 1), c(0, 1))
    par(mar=c(0, 0, 0, 0))
    plot.new()
    # get dimensions from outline PNG
    file <- system.file("/extdata/aoscp/cell/outline.png", package = "JMDplots")
    img <- png::readPNG(file)
    if(outline) ZC <- list(outline=NULL, extracellular=NULL)
    # loop over locations
    for(j in 1:length(ZC)) {
      file <- system.file(paste0("/extdata/aoscp/cell/", names(ZC)[j], ".png"), package = "JMDplots")
      img <- png::readPNG(file)
      if(!outline) {
        # points are where the alpha is not zero
        isthere <- img[, , 4]!=0
        # set the colors
        for(i in 1:3) {
          imgnew <- img[, , i]
          imgnew[isthere] <- col[j, i]/255
          img[, , i] <- imgnew
        }
      }
      pu <- par("usr")
      rasterImage(img, pu[1], pu[3], pu[2], pu[4])
    }
    invisible(dev.off())
  } else {
    # to make labeled plot (PDF)
    pdf("aoscp3.pdf", width=12, height=8, family="Times")
    layout(t(matrix(c(1, 1, 1, 1, 1, 2))))
    # plot base
    img <- readPNG("aoscp3.png")
    # scale dimensions 5/4
    d <- dim(img[, , 1])
    r <- d[1] / d[2] * 5/4
    y1 <- (1 - r)/2
    y2 <- y1 + r
    plot.new()
    rasterImage(img, 0, y1, 1, y2)
    # draw labels
    # cex <- 4   # for PNG
    cex <- 3
    text(0.02, 0.85, "endoplasmic\nreticulum", cex=cex, adj=0)
    lines(c(0.095, 0.21), c(0.795, 0.6), lwd=2)
    text(0.35, 0.9, "nucleus", cex=cex)
    lines(c(0.34, 0.37), c(0.87, 0.64), lwd=2)
    text(0.53, 0.88, "nucleolus", cex=cex)
    lines(c(0.5, 0.4), c(0.85, 0.58), lwd=2)
    text(0.66, 0.77, "bud\nneck", cex=cex, adj=0)
    lines(c(0.69, 0.7), c(0.71, 0.57), lwd=2)
    text(0.83, 0.76, "bud tip", cex=cex)
    lines(c(0.83, 0.8), c(0.73, 0.55), lwd=2)
    text(0.29, 0.08, "cytoplasm", cex=cex)
    lines(c(0.29, 0.32), c(0.11, 0.3), lwd=2)
    text(0.47, 0.1, "vacuole", cex=cex)
    lines(c(0.45, 0.44), c(0.12, 0.31), lwd=2)
    text(0.67, 0.11, "mitochondrion", cex=cex)
    lines(c(0.64, 0.61), c(0.14, 0.38), lwd=2)
    text(0.81, 0.16, "extracellular", cex=cex)
    lines(c(0.81, 0.83), c(0.18, 0.3), lwd=2)
    text(0.02, 0.15, "Golgi\napparatus", cex=cex, adj=0)
    lines(c(0.09, 0.18), c(0.21, 0.38), lwd=2)
    # draw color legend
    plot.new()
    plot.window(xlim=c(0, 1), ylim=c(ZCmin, ZCmax))
    par(las=1, mar=c(6, 0, 4, 10), cex=cex/2)
    axis(4, at=seq(-0.21, -0.095, 0.023))
    # make color bar
    rst <- array(1, c(101, 1, 4))
    col <- cR(seq(1, 0, -0.01))/255
    rst[, , 1] <- col[, 1]
    rst[, , 2] <- col[, 2]
    rst[, , 3] <- col[, 3]
    # use angle=180 to get the right orientation
    rasterImage(rst, -1, ZCmin, 1, ZCmax)
    title(main=expression(italic(Z)[C]))
    invisible(dev.off())
  }
}

# ZC and Eh ranges in yeast and ER-cytoplasm electron-transfer scheme
aoscp4 <- function(pdf = FALSE) {
  ## start plot
  if(pdf) pdf("aoscp4.pdf", width=8, height=4, family="Times")
  par(mfrow=c(1, 2))
  par(mgp=c(2.8, 1, 0))
  par(mar=c(4, 4, 1, 1))
  par(las = 1)
  ## make ZC-Eh plot for cell compartments
  dat <- data.frame(
    location = c("C", "M", "E", "V", "X"),
    Eh.min = c(-320, -360, -208, -160, -150),
    Eh.max = c(-240, -255, -133, 100, 160)
  )
  plot(c(-350, 0), c(-0.25, -0.05), type="n", 
    xlab="Eh, mV", ylab=expression(italic(Z)[C]))
  for(i in 1:nrow(dat)) {
    ZC <- yeast.ZC(dat$location[i])
    rect(dat$Eh.min[i], fivenum(ZC)[2], dat$Eh.max[i], fivenum(ZC)[4])
  }
  # it's easier just to place the labels manually
  text(-320, -0.086, "cytoplasm", adj=0)
  text(-360, -0.196, "mitochondrion", adj=0)
  text(-178, -0.245, "ER", adj=1)
  text(-120, -0.228, "vacuole", adj=0)
  text(-150, -0.048, "extracellular", adj=0)
  # add plot label (a)
  par(xpd=NA)
  label.figure("a", paren = TRUE, italic = TRUE)
  ## to plot protein electron donating scheme ER-cytoplasm
  par(mar=c(0, 1, 0, 0))
  plot.new()
  plot.window(c(-2, 2), c(-1.7, 1.7))
  draw.arc(-0.5, -1.2, deg1=10, deg2=170)
  draw.arc(-0.5, 1.2, deg1=190, deg2=350)
  # separate and label compartments
  lines(c(-0.5, -0.5), c(-1.5, 1.5), lty=2)
  text(-1.5, 0, "ER")
  text(0.5, 0, "cytoplasm")
  # label Eh
  text(-1.5, 1.2, "-208 mV")
  text(0.5, 1.2, "-320 mV")
  text(1, 1.2, "Eh (GSH)", adj=0)
  # label ZC
  text(-1.5, -1.2, "-0.19")
  text(0.5, -1.2, "-0.14")
  text(1, -1.2, expression(italic(Z)[C]~"(protein)"), adj=0)
  # add electron arrow
  arrows(-0.9, -0.4, -0.1, 0.4, length=0.2)
  text(-0.3, 0, expression(italic(e)^-phantom()))
  # add plot label (b)
  par(xpd=NA)
  label.figure("b", paren = TRUE, italic = TRUE)
  ## done!
  if(pdf) invisible(dev.off())
}

# average oxidation state of carbon in proteins from different organisms
# adapted from ?protein.formula
aoscp5 <- function(pdf = FALSE) {
  if(pdf) pdf("aoscp5.pdf", width=10, height=6, family="Times")
  par(mar=c(7, 4.5, 2, 2))
  par(las = 1)
  # get amino acid compositions of microbial proteins 
  # generated from the RefSeq database 
  file <- system.file("/extdata/aoscp/protein_refseq61.csv.xz", package = "JMDplots")
  aa <- read.csv(file, as.is = TRUE)
  # calculate ZC
  ip <- add.protein(aa)
  pf <- protein.formula(ip)
  zc <- ZC(pf)
  ## save the calculated values in out/ZC_refseq.csv
  #name <- gsub("]$","",unlist(lapply(lapply(strsplit(aa$ref, ")[", fixed=TRUE), rev), "[", 1)))
  #out <- data.frame(taxid=aa$organism, name=name, length=aa$abbrv, ZC=round(zc, 4))
  #write.csv(out, "out/ZC_refseq.csv", row.names=FALSE, quote=2)
  # only use those organisms with a minimum
  # number of sequenced bases
  ibig <- as.numeric(aa$abbrv) >= 50000
  ip <- ip[ibig]
  aa <- aa[ibig, ]
  zc <- zc[ibig]
  # the organism names we search for
  # "" matches all organisms
  terms <- c("Natr", "Halo", "Rhodo", "Acido", "Methylo",
    "Chloro", "Nitro", "Desulfo", "Geo", "Methano",
    "Thermo", "Pyro", "Sulfo",
    "Streptomyces", "Mycobacterium", "Pseudomonas", "Escherichia", "Salmonella", 
    "Vibrio", "Bacteroides", "Lactobacillus", "Staphylococcus", "Streptococcus",
    "Listeria", "Bacillus", "Clostridium", "Mycoplasma", "Buchnera", "")
  nterms <- length(terms)
  tps <- thermo()$protein$ref[ip]
  plot(0, 0, xlim=c(1, nterms), ylim=c(-0.35, -0.05), pch="",
    ylab=expression(italic(Z)[C]),
    xlab="", xaxt="n", mar=c(6, 3, 1, 1))
  # set seed for reproducible jitter
  set.seed(101)
  for(i in 1:nterms) {
    it <- grep(terms[i], tps)
    zct <- zc[it]
    #print(paste(terms[i], round(mean(zct), 3)))
    if(i < 14) factor <- 1 else factor <- 0.5
    points(jitter(rep(i, length(zct)), factor), zct, pch=20, cex=0.7)
  }
  terms[terms==""] <- paste("all", length(ip))
  axis(1, 1:nterms, terms, las=2)
  if(pdf) invisible(dev.off())
}

