# JMDplots/gradH2O.R
# R file started for Goldschmidt poster 20190711
# moved to JMDplots for salinity gradients paper starting 20190930

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

# ZC for selected redox gradients 20190715
gradH2O1 <- function(pdf = FALSE) {
  if(pdf) pdf("gradH2O1.pdf", width = 7, height = 5)
  par(mfrow = c(2, 2))
  par(mar = c(4, 4, 2, 1), mgp = c(2.5, 0.7, 0))
  par(cex.lab = 1)
  mplot("Bison_Pool", "IMG_MGP", add.label = FALSE, plottype = "#FF000030", col = "red")
  label.figure("A", cex = 1.7, yfrac = 0.93)
  mplot("Diffuse_Vents", "SRA_MGP", add.label = FALSE, plottype = "#FF000030", col = "red")
  label.figure("B", cex = 1.7, yfrac = 0.93)
  mplot("Guerrero_Negro", "IMG_MGP", add.label = FALSE, plottype = "#FF000030", col = "red", ylim = c(-0.150, -0.125))
  label.figure("C", cex = 1.7, yfrac = 0.93)
  # add proteomes from Nif-encoding genomes 20191014
  # get mean and SD of ZC values
  np <- NifProteomes()
  ZCmean <- np$ZC
  ZChi <- ZCmean + np$ZC.SD
  ZClo <- ZCmean - np$ZC.SD
  nn <- length(ZCmean)
  # set up plot
  plot(0, 0, xlim = c(1, nn), ylim = c(-0.24, -0.13), xlab = "Nif type", ylab = NA, xaxt = "n")
  mtext(quote(italic(Z)[C]), side = 2, line = 2, las = 0)
  axis(1, 1:nn, np$type)
  # draw polygon (filled area)
  polygon(c(1:nn, nn:1), c(ZChi, rev(ZClo)), col = "#FF000030", border = NA)
  # draw lines and points
  lines(1:nn, ZCmean, col = "red", lty = 2)
  points(1:nn, ZCmean, pch = 19, col = "red")
  # add title
  title("Nif-bearing genomes (NF)", font.main = 1)
  label.figure("D", cex = 1.7, yfrac = 0.93)
  # done!
  if(pdf) {
    dev.off()
    addexif("gradH2O1", "Carbon oxidation state for selected redox gradients", "Dick et al. (2020) (preprint)")
  }
}

# nH2O-ZC scatterplots for redox gradients and the Baltic Sea 20190713
gradH2O2 <- function(pdf = FALSE, vars = "H2O-ZC") {
  if(pdf) pdf("gradH2O2.pdf", width = 12, height = 5.6)
  par(mfrow = c(1, 2))
  par(mar = c(4, 4.5, 2, 1), las = 1, cex = 1.2)
  par(cex.lab = 1.5)

  # plot 1: compare ZC and nH2O of proteins in datasets from gradox paper
  mgradox <- ppage("gradoxGS", plot.it = FALSE)
  pgradox <- ppage("gradoxGS", H2O = TRUE, plot.it = FALSE)
  pcomp(mgradox, pgradox, reorder = FALSE, yline = 3.2, vars = vars, cex.ylab = 1.8, font = 2, labdy = 0.003)
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
  }
  if(vars == "pIG") {
    text(c(5.00, 7.49, 5.65, 5.68), c(-0.123, -0.130, -0.109, -0.066), c("Nif-D", "Nif-C", "Nif-B", "Nif-A"), adj = 0)
    text(c(5.37, 5.42, 5.57), c(-0.158, -0.176, -0.216), c("1 mm", "2 mm", "3 mm"), adj = 1)
  }
  title("Redox gradients", font.main = 1)
  label.figure("A", cex = 2, xfrac = 0.035)

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
  }
  title("Baltic Sea", font.main = 1)
  label.figure("B", cex = 2, xfrac = 0.035)

  if(pdf) {
    dev.off()
    addexif("gradH2O2", "nH2O-ZC scatterplots for redox gradients and the Baltic Sea", "Dick et al. (2020) (preprint)")
  }
}

# nH2O for Baltic Sea metagenome and metatranscriptome in different size fractions 20190715
gradH2O3 <- function(pdf = FALSE, var = NULL) {
  if(pdf) pdf("gradH2O3.pdf", width = 7, height = 3)
  nvar <- ifelse(is.null(var), 1, length(var))
  if(nvar==2) par(mfrow = c(2, 3))
  else par(mfrow = c(1, 3))
  par(mar = c(5, 4, 1, 1), mgp = c(3, 0.7, 0), las = 1)
  par(cex.lab = 1.5)

  for(i in 1:nvar) {
    ylim <- NULL
    if(identical(var[i], "pI")) ylim <- c(5, 7.5)
    figlab <- c("A", "B", "C")
    if(i==2) figlab <- c("D", "E", "F")

    mplot("Baltic_Sea-0.1s", "iMicrobe_MGP", H2O = TRUE, plottype = "#FF000030", col = "red", add.title = FALSE, ylim = ylim, yline = 2.7, var = var[i])
    mplot("Baltic_Sea-0.1s", "SRA_MTP", H2O = TRUE, plottype = "#0000FF30", col = "blue", add.title = FALSE, add = TRUE, pch = 15, var = var[i])
    legend("topright", legend = quote("0.1-0.8"~mu*m), bty = "n")
    legend("bottomleft", c("metagenomes", "metatranscriptomes"), pch = c(19, 15), col = c("red", "blue"), lty = c(2, 2), bty = "n")
    label.figure(figlab[1], cex = 1.8, xfrac = 0.04)

    mplot("Baltic_Sea-0.8s", "iMicrobe_MGP", H2O = TRUE, plottype = "#FF000030", col = "red", add.title = FALSE, ylim = ylim, yline = 2.7, var = var[i])
    mplot("Baltic_Sea-0.8s", "SRA_MTP", H2O = TRUE, plottype = "#0000FF30", col = "blue", add.title = FALSE, add = TRUE, pch = 15, var = var[i])
    legend("topright", legend = quote("0.8-3.0"~mu*m), bty = "n")
    label.figure(figlab[2], cex = 1.8, xfrac = 0.035)

    mplot("Baltic_Sea-3.0s", "iMicrobe_MGP", H2O = TRUE, plottype = "#FF000030", col = "red", add.title = FALSE, ylim = ylim, yline = 2.7, var = var[i])
    mplot("Baltic_Sea-3.0s", "SRA_MTP", H2O = TRUE, plottype = "#0000FF30", col = "blue", add.title = FALSE, add = TRUE, pch = 15, var = var[i])
    legend("topright", legend = quote("3.0-200"~mu*m), bty = "n")
    label.figure(figlab[3], cex = 1.8, xfrac = 0.035)
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
  ylim <- c(-0.04, 0.04)

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
  label.figure("A", xfrac = 0.1, cex = 1.8)
  # plot 2: GRAVY - pI
  pcomp(mout, pout, lty = 0, yline = 3, vars = "pIG", labels.at = NA, cex.ylab = 0.9)
  hullfun(mout, pout, c(1, 3), "green3", c("riverFL", "riverPA"), vars = "pIG")
  hullfun(mout, pout, c(1, 3), "blue", c("plumeFL", "plumePA"), vars = "pIG")
  text(c(7.2, 8.1), c(-0.14, -0.175), c("river", "plume"))
  title("Amazon River metagenomes", font.main = 1)
  label.figure("D", xfrac = 0.1, cex = 1.8)

  # plots 3-4: Amazon river metatranscriptome
  pcomp(mout, pout, "MT", xlim = xlim, ylim = ylim, lty = 0, yline = 2.8, labels.at = NA)
  hullfun(mout, pout, c(2, 4), "green3", c("riverFL", "riverPA"))
  hullfun(mout, pout, c(2, 4), "blue", c("plumeFL", "plumePA"))
  text(c(-0.125, -0.12), c(0.025, 0.005), c("river", "plume"))
  title("Amazon River metatranscriptomes", font.main = 1)
  label.figure("B", xfrac = 0.1, cex = 1.8)
  # plot 4: GRAVY - pI
  pcomp(mout, pout, "MT", lty = 0, yline = 3, vars = "pIG", labels.at = NA, cex.ylab = 0.9)
  hullfun(mout, pout, c(2, 4), "green3", c("riverFL", "riverPA"), vars = "pIG")
  hullfun(mout, pout, c(2, 4), "blue", c("plumeFL", "plumePA"), vars = "pIG")
  text(c(8.2, 6.5), c(-0.105, -0.11), c("river", "plume"))
  title("Amazon River metatranscriptomes", font.main = 1)
  label.figure("E", xfrac = 0.1, cex = 1.8)

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
  label.figure("C", xfrac = 0.1, cex = 1.8)

  # plot 6: GRAVY - pI
  pcomp(moutE, poutE, lty = 0, yline = 3, vars = "pIG", labels.at = NA, cex.ylab = 0.9)
  hullfun(moutE, poutE, 1, "green3", vars = "pIG")
  hullfun(moutE, poutE, 2, "blue", vars = "pIG")
  pcomp(moutH, poutH, reorder = FALSE, add = TRUE, vars = "pIG", labels.at = "min", labdx = -0.25)
  hullfun(moutH, poutH, 1:3, "turquoise3", vars = "pIG")
  text(c(7.5, 7.4, 6.3), c(-0.14, -0.20, -0.27), c("freshwater", "marine", "hypersaline"))
  title("Freshwater - marine - hypersaline", font.main = 1)
  label.figure("F", xfrac = 0.1, cex = 1.8)

  if(pdf) {
    dev.off()
    addexif("gradH2O4", "nH2O vs ZC for freshwater, marine, and hypersaline environments", "Dick et al. (2020) (preprint)")
  }
}

# mean differences of compositional metrics for differentially expressed proteins in osmotic stress
# adapted from canprot/hyperosmotic.Rmd 20190717-20191007
# add GRAVY and pI plot 20191028
# use different symbols for eukaryotes and add more bacteria and archaea (osmotic2) 20191102-20191103
gradH2O5 <- function(pdf = FALSE) {
  if(pdf) pdf("gradH2O5.pdf", width = 8, height = 5.2)
  layout(matrix(c(0,0,1,1,1,1,1,1,0,0, 2,2,2,2,2,3,3,3,3,3), nrow = 2, byrow = TRUE), heights = c(0.3, 1))
  par(cex = 1)

  # make legend
  par(mar = c(0.5, 1, 0, 1))
  plot.new()
  plot.window(c(1, 12), c(1, 5.3))
  par(xpd = NA)
  points(rep(0.5, 3), 1:3, pch = c(15, 2, 1), cex = c(0.5, 1.2, 1.5), col = "slategray3")
  points(rep(2.5, 3), 2:4, pch = c(2, 1, 0), cex = 2)
  text(rep(2.5, 3), 2:4, c("E", "C", "A"), cex = 0.75)
  text(c(0.3, 2.7), c(5.3, 5.3), c("Dick (2017)", "This study"), adj = c(0.5, 1))
  text(4, 5.3, "Domain, method, conditions", adj = c(0, 1))
  text(4, 4, "Bacteria and Archaea, proteomics, optimal vs hypoosmotic", adj = 0)
  text(4, 3, "Bacteria and Archaea, proteomics, hyperosmotic", adj = 0)
  text(4, 2, "Bacteria and Archaea, transcriptomics, hyperosmotic", adj = 0)
  text(4, 1, "Eukaryotes, proteomics, hyperosmotic", adj = 0)
  z <- MASS::kde2d(rep(c(0, 0.5, 1, 1.5, 2, 2.5, 3), 5), rep(c(2.6, 2.8, 3, 3.2, 3.4), each = 7))
  contour(z, drawlabels = FALSE, levels = 0.22, lty = 2, add = TRUE, col = "maroon3")
  # adjust plot parameters for main plots
  par(mar = c(4, 4, 1, 1), mgp = c(2.5, 1, 0))
  par(xpd = FALSE)

  # generate compositional table
  datasets1 <- pdat_osmotic()
  pdat1 <- lapply_canprot(datasets1, function(dataset) {
    pdat_osmotic(dataset, basis = "rQEC")
  })
  # use filled square for eukaryotes, open circle for bacteria, open triangle for gene expression 20191102
  pch1 <- rep(15, length(datasets1))
  igene1 <- grepl("trans", datasets1)
  pch1[igene1] <- 2
  ibact1 <- (grepl("KKG\\+12", datasets1) | grepl("KLB\\+15", datasets1)) & !igene1
  pch1[ibact1] <- 1
  pt1 <- rep(NA, length(datasets1))
  col1 <- rep("slategray3", length(datasets1))
  cex1 <- rep(0.5, length(datasets1))
  cex1[ibact1] <- 1.5
  cex1[igene1] <- 1.2

  # get new data (added for this paper)
  datasets2 <- pdat_osmotic2()
  pdat2 <- lapply_canprot(datasets2, function(dataset) {
    pdat_osmotic2(dataset, basis = "rQEC")
  })
  # use open circle except open triangle for gene expression and open square for hypo-osmotic conditions 20191103
  pch2 <- rep(1, length(datasets2))
  igene2 <- grepl("trans", datasets2)
  pch2[igene2] <- 2
  ihypo2 <- datasets2 %in% c("LRB+09_2.6", "ZLZ+16_10", "LLYL17_0", "JSP+19_LoS")
  pch2[ihypo2] <- 0
  pt2 <- LETTERS[1:length(datasets2)]
  col2 <- rep("black", length(datasets1))
  cex2 <- rep(2, length(datasets2))
  # include non-eukaryotes and non-transcriptomes in kde2d contours
  contour <- c(ibact1, !(igene2 | ihypo2))

  # put the data together
  pdat <- c(pdat1, pdat2)
  datasets <- c(datasets1, datasets2)
  pch <- c(pch1, pch2)
  pt.text <- c(pt1, pt2)
  col <- c(col1, col2)
  cex <- c(cex1, cex2)

  # get the nH2O - ZC values
  comptab1 <- lapply_canprot(pdat, get_comptab, plot.it = FALSE, mfun = "mean")
  # get the GRAVY - pI values
  comptab2 <- lapply_canprot(pdat, get_comptab, var1 = "pI", var2 = "GRAVY", plot.it = FALSE, mfun = "mean")
  # make the plots
  diffplot(comptab1, pt.text = pt.text, pch = pch, col = col,
           contour = contour, cex = cex, col.contour = "maroon3", cex.text = 0.75)
  label.figure("A", cex = 1.7)
  diffplot(comptab2, vars = c("pI", "GRAVY"), pt.text = pt.text, pch = pch, col = col,
           contour = contour, cex = cex, col.contour = "maroon3", cex.text = 0.75)
  label.figure("B", cex = 1.7)
  if(pdf) {
    dev.off()
    addexif("gradH2O5", "Mean differences of compositional metrics for differentially expressed proteins in osmotic stress", "Dick et al. (2020) (preprint)")
  }
}

# basis species comparison 20190713 / relative water content 20191005
gradH2O6 <- function(pdf = FALSE) {

  # set up figure
  if(pdf) pdf("gradH2O6.pdf", width = 6, height = 6)
  par(mfrow = c(3, 3))
  par(mar = c(3.2, 3.2, 2.5, 1))
  par(mgp = c(2, 0.7, 0))
  par(las = 1)
  par(cex.lab = 1.2)

  # define axis labels
  nH2Olab <- expression(italic(n)[H[2] * O])
  nO2lab <- expression(italic(n)[O[2]])
  ZClab <- expression(italic(Z)[C])

  # function to plot linear model
  lmfun <- function(ZC, nH2O, legend.x = NULL, ...) {
    mylm <- lm(nH2O ~ ZC)
    xlim <- par("usr")[1:2]
    # for amino acid plots, draw line from -1 to +1
    if(!is.null(legend.x)) xlim <- c(-1, 1)
    lines(xlim, predict(mylm, data.frame(ZC = xlim)), ...)
    # add R-squared text
    if(!is.null(legend.x)) {
      R2 <- format(round(summary(mylm)$r.squared, 2), nsmall = 2)
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
  ecoli <- read.csv(system.file("/extdata/organisms/ecoli.csv.xz", package = "JMDplots"), as.is = TRUE)
  pf.ecoli <- protein.formula(ecoli)
  ZC.ecoli <- ZC(pf.ecoli)

  # calculate nH2O and nO2 of amino acids with CHNOS basis
  basis("CHNOS")
  species(aa)

  # plot 1: nH2O-ZC of amino acids (CHNOS)
  aaplot(ZC.aa, species()$H2O, nH2Olab, "bottomleft", "Amino acids (CHNOS)", "A")

  # plot 2: nO2-ZC of amino acids (CHNOS)
  par(mgp = c(2.2, 0.7, 0), xpd = NA)
  aaplot(ZC.aa, species()$O2, nO2lab, "bottomright", "Amino acids (CHNOS)", "B")

  # calculate nH2O and nO2 with QEC basis
  basis("QEC")
  species(aa)

  # plot 3: nO2-ZC of amino acids (QEC)
  aaplot(ZC.aa, species()$O2, nO2lab, "bottomright", "Amino acids (QEC)", "C")

  # plot 4: nH2O-ZC of amino acids (QEC)
  par(mgp = c(2, 0.7, 0))
  # save the residuals here
  # subtract 0.355 so mean of human proteins = 0 20191114
  # (see canprot::H2OAA)
  rQEC <- aaplot(ZC.aa, species()$H2O, nH2Olab, "bottomright", "Amino acids (QEC)", "D") - 0.355
  names(rQEC) <- aminoacids(3)

  # plot 5: nH2O-ZC of human proteins (QEC)
  par(mgp = c(2.2, 0.7, 0))
  nH2O.human <- H2OAA(human, basis = "QEC")
  scatterfun(ZC.human, nH2O.human, "Human proteins (QEC)", "E")

  # plot 6: nH2O-ZC of E. coli proteins (QEC)
  nH2O.ecoli <- H2OAA(ecoli, basis = "QEC")
  scatterfun(ZC.ecoli, nH2O.ecoli, quote(italic(E.~coli)*" proteins (QEC)"), "F")

  # plot 7: nH2O-ZC of amino acids (rQEC)
  par(mgp = c(2, 0.7, 0))
  aaplot(ZC.aa, rQEC, nH2Olab, "bottomright", "Amino acids (rQEC)", "G")

  # plot 8: nH2O-ZC of human proteins (rQEC)
  par(mgp = c(2.2, 0.7, 0))
  nH2O.human <- H2OAA(human, basis = "rQEC")
  scatterfun(ZC.human, nH2O.human, "Human proteins (rQEC)", "H")

  # plot 9: nH2O-ZC of E. coli proteins (rQEC)
  nH2O.ecoli <- H2OAA(ecoli, basis = "rQEC")
  scatterfun(ZC.ecoli, nH2O.ecoli, quote(italic(E.~coli)*" proteins (rQEC)"), "I")

  if(pdf) {
    dev.off()
    addexif("gradH2O6", "Basis species comparison and stoichiometric hydration state", "Dick et al. (2020) (preprint)")
  }
  ## output value of rQEC for checking code of H2OAA()
  #rQEC
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

# Supplementary Figure S1: scatterplots of GRAVY vs ZC and nH2O 20191117
gradH2OS1 <- function(pdf = FALSE) {
  if(pdf) pdf("gradH2OS1.pdf", width = 7, height = 7)
  # get amino acid compositions of human and E.coli proteins (UniProt)
  human <- get("human_base", canprot)
  ecoli <- read.csv(system.file("/extdata/organisms/ecoli.csv.xz", package = "JMDplots"), as.is = TRUE)
  ZC <- ZCAA(human)
  nH2O <- H2OAA(human)
  GRAVY <- GRAVY(human)
  scatterfun <- function(xvar, yvar, AAcomp) {
    if(xvar=="ZC") {
      x <- ZCAA(AAcomp)
      xlab <- expression(italic(Z)[C])
    }
    if(xvar=="nH2O") {
      x <- H2OAA(AAcomp)
      xlab <- expression(italic(n)[H[2] * O])
    }
    if(yvar=="GRAVY") {
      y <- GRAVY(AAcomp)
      ylab <- "GRAVY"
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
  par(mfrow = c(2, 2))
  scatterfun("ZC", "GRAVY", human)
  title("Human proteins", font.main = 1)
  scatterfun("nH2O", "GRAVY", human)
  title("Human proteins (rQEC)", font.main = 1)
  scatterfun("ZC", "GRAVY", ecoli)
  title(quote(italic(E.~coli)*" proteins"), font.main = 1)
  scatterfun("nH2O", "GRAVY", ecoli)
  title(quote(italic(E.~coli)*" proteins (rQEC)"), font.main = 1)
  if(pdf) {
    dev.off()
    addexif("gradH2OS1", "Scatterplots of GRAVY vs ZC and nH2O", "Dick et al. (2020) (preprint)")
  }
}

# Supplementary Figure S2 20191028 (provisional, not in final paper)
gradH2OS2 <- function(pdf = FALSE) {
  if(pdf) pdf("gradH2OS2.pdf", width = 12, height = 5.6)
  gradH2O2(vars = "pIG")
  if(pdf) {
    dev.off()
    addexif("gradH2OS2", "GRAVY-pI scatterplots for redox gradients and the Baltic Sea", "Dick et al. (2020) (preprint)")
  }
}

# Supplementary Figure S3 20191028 (provisional, not in final paper)
gradH2OS3 <- function(pdf = FALSE) {
  if(pdf) pdf("gradH2OS3.pdf", width = 6, height = 2.5)
  gradH2O3(var = c("GRAVY", "pI"))
  if(pdf) {
    dev.off()
    addexif("gradH2OS3", "GRAVY and pI for Baltic Sea metagenome and metatranscriptome in different size fractions", "Dick et al. (2020) (preprint)")
  }
}

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

