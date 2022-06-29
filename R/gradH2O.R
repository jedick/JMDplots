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
  if(pdf) pdf("gradH2O1.pdf", width = 8, height = 4.3)
  layout(matrix(c(1, 3, 2, 4, 5, 5), nrow = 2), widths = c(1, 1, 2))
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
    R2 <- summary(mylm)$r.squared
    if(!is.null(legend.x)) {
      rR2 <- format(round(R2, 3), nsmall = 3)
      R2txt <- substitute(italic(R)^2 == R2, list(R2 = rR2))
      legend(legend.x, legend = R2txt, bty = "n")
    }
    invisible(R2)
  }

  # function to add plot title and panel label
  mainlab <- function(main, lab, xfrac) {
    title(main, font.main = 1, cex.main = 1)
    label.figure(lab, xfrac = xfrac, yfrac = 0.8, cex = 1.6)
  }

  # function to plot nH2O or nO2 vs ZC of amino acids
  aaplot <- function(ZC, y, ylab, legend.x, main, lab, xfrac = 0.06) {
    plot(ZC, y, type = "p", pch = aminoacids(1), xlab = ZClab, ylab = ylab)
    mainlab(main, lab, xfrac)
    lmfun(ZC, y, c(-1, 1), legend.x)
  }

  # get names of amino acids
  aa <- aminoacids("")
  # calculate ZC of the amino acids
  ZC.aa <- ZC(info(aa, "aq"))

  # Calculate nH2O and nO2 of amino acids with CHNOS basis
  basis("CHNOS")
  species(aa)
  # Plot 1: nH2O-ZC of amino acids (CHNOS)
  CHNOS.H2O.R2 <- aaplot(ZC.aa, species()$H2O, nH2Olab, "bottomleft", "", "(a)")
  # Plot 2: nO2-ZC of amino acids (CHNOS)
  CHNOS.O2.R2 <- aaplot(ZC.aa, species()$O2, nO2lab, "bottomright", "", "(b)", 0.04)
  mtext(quote(list(CO[2], NH[3], H[2]*S, H[2]*O, O[2])~"("*bold(CHNOS)*")                 "), line = 0.5, adj = 1, cex = 0.8)

  # Calculate nH2O and nO2 with QEC basis
  basis("QEC")
  species(aa)
  # Plot 3: nH2O-ZC of amino acids (QEC)
  aaplot(ZC.aa, species()$H2O, nH2Olab, "bottomright", "", "(c)")
  # Plot 4: nO2-ZC of amino acids (QEC)
  aaplot(ZC.aa, species()$O2, nO2lab, "bottomright", "", "(d)", 0.04)
  mtext(quote(list(Glutamine, "glutamic acid", cysteine, H[2]*O, O[2])~"("*bold(QEC)*")     "), line = 0.5, adj = 1, cex = 0.8)

  # Plot 5: R2 of nH2O-ZC and nO2-ZC fits
  par(mar = c(4, 4.7, 2.5, 1))
  par(mgp = c(3, 0.7, 0))
  plotbasisfun()
  file <- system.file("extdata/gradH2O/AAbasis.csv", package = "JMDplots")
  AAbasis <- read.csv(file, as.is = TRUE)
  nbasis <- sum(!is.na(AAbasis$O2.R2))
  # add 1 for CHNOS 20200821
  nbasis <- nbasis + 1
  # add point for CHNOS basis 20200821
  points(CHNOS.O2.R2, CHNOS.H2O.R2, pch = 19, cex = 0.6, col = 2)
  text(CHNOS.O2.R2, CHNOS.H2O.R2, "CHNOS", adj = -0.15, font = 2)
  # label QEC 20201011
  QEC <- AAbasis[AAbasis$abbrv == "CEQ", ]
  text(0.98, 0.06, "QEC", font = 2)
  lines(c(QEC$O2.R2 + 0.004, 0.94), c(QEC$H2O.R2 + 0.004, 0.055))
  # label MWY 20201015
  MWY <- AAbasis[AAbasis$abbrv == "MWY", ]
  points(MWY$O2.R2, MWY$H2O.R2, pch = 19, cex = 0.6, col = 2)
  text(0.97, 0.002, "MWY")
  title(paste(nbasis, "combinations of basis species"), font.main = 1)
  label.figure("(e)", cex = 1.6, yfrac = 0.94)

  if(pdf) {
    dev.off()
    addexif("gradH2O1", "Comparison of different sets of basis species", "Dick et al. (2020)")
  }
  
  # Return linear model of nH2O-ZC for amino acids
  invisible(lm(species()$H2O ~ ZC.aa))
}

# Schematic of nH2O and ZC calculations 20200826
gradH2O2 <- function(pdf = FALSE) {
  # Setup plot
  if(pdf) pdf("gradH2O2.pdf", width = 5, height = 2.5)
  plot.new()
  par(mar = c(0, 0, 0, 0))
  par(usr = c(0, 1, 0.4, 1))
  par(xpd = NA)

  # Make headings
  text(0.30, 0.95, "Elemental composition")
  text(0.80, 0.95, "Amino acid composition")
  text(0.05, 0.80, "Basis\nspecies", srt = 90, cex = 0.9)
  text(0.05, 0.65, quote(italic(n)[H[2]*O]))
  text(0.05, 0.52, quote(italic(Z)[C]))

  # Use LYSC_CHICK
  protein <- "LYSC_CHICK"
  f <- protein.formula("LYSC_CHICK")
  fexpr <- expr.species(as.chemical.formula(f))

  # Write basis species
  basis(c("glutamine", "glutamic acid", "cysteine", "H2O", "O2"))
  reaction <- subcrt(protein, -1)$reaction
  text(0.30, 0.90, bquote(.(fexpr) == phantom(.)), cex = 0.8)
  for(i in 2:6) {
    dy <- (i - 2) * -0.03
    coeff <- formatC(reaction$coeff[i], 1, format = "f")
    if(i > 4) fexpr <- expr.species(reaction$formula[i])
    else fexpr <- bquote(.(expr.species(reaction$formula[i]))~"("*.(reaction$name[i])*")")
    text(0.22, 0.86 + dy, coeff, cex = 0.7, adj = 1)
    text(0.24, 0.86 + dy, fexpr, cex = 0.7, adj = 0)
  }
  #text(0.12, 0.86 + dy, "+", cex = 0.7)

  # Write nH2O equation
  # Divide nH2O by protein length
  pl <- protein.length(protein)
  coeff <- formatC(reaction$coeff[5], 1, format = "f")
  nH2O <- formatC(reaction$coeff[5] / pl, 3, format = "f")
  nH2Oeqn <- bquote(frac(.(coeff), .(pl)~"(protein length)") == bold(.(nH2O)))
  text(0.3, 0.65, nH2Oeqn, cex = 0.75)

  # Write ZC calculation
  df <- as.data.frame(f)
  ZC <- formatC(ZC(f), 3, format = "f")
  ZCeqn <- bquote(frac(- .(df$H) + 3*"("*.(df$N)*")" + 2*"("*.(df$O)*")" + 2*"("*.(df$S)*")", .(df$C)) == bold(.(ZC)))
  text(0.35, 0.52, ZCeqn, cex = 0.75)
  text(0.31, 0.45, "(Equation 1)")

  # Write amino acid composition
  xAA <- seq(0.65, 0.95, length.out = 10)
  AA <- pinfo(pinfo("LYSC_CHICK"))[, 6:25]
  text(xAA, 0.9, aminoacids()[1:10], cex = 0.75)
  text(xAA, 0.87, AA[1:10], cex = 0.75)
  text(xAA, 0.82, aminoacids()[11:20], cex = 0.75)
  text(xAA, 0.79, AA[11:20], cex = 0.75)

  # Write Equations using AA composition
  text(0.75, 0.65, bquote(phantom(.) %<-% "Equation 2"))
  text(0.75, 0.52, bquote(phantom(.) %<-% "Equation 3"))
  mid <- 0.603
  rect(0.001, 0.401, mid, 0.999, border = "gray60")
  rect(mid, 0.401, 0.999, 0.999, border = "gray60")

  # done!
  if(pdf) {
    dev.off()
    addexif("gradH2O2", "Schematic of nH2O and ZC calculations", "Dick et al. (2020)")
  }
  
}

# nH2O-ZC scatterplots for redox gradients and the Baltic Sea 20190713
gradH2O3 <- function(pdf = FALSE, vars = "H2O-ZC") {
  if(pdf) pdf("gradH2O3.pdf", width = 12, height = 5.6)
  layout(matrix(c(1, 2, 3), nrow = 1), widths = c(3, 2, 2))
  par(mar = c(3.5, 4.5, 2, 1), las = 1, cex = 1.2)
  par(mgp = c(2.5, 1, 0))
  par(cex.lab = 1.3)

  # define axis labels
  nH2Olab <- expression(italic(n)[H[2] * O])
  ZClab <- expression(italic(Z)[C])
  # set y-axis limit depending on basis species 20201011
  if(getOption("basis") == "QCa") ylim <- c(-1.10, -1.02)
  if(getOption("basis") == "QEC") ylim <- c(-0.78, -0.70)

  # plot 1: compare ZC and nH2O of proteins in datasets from gradox paper
  plot(c(-0.22, -0.12), ylim, xlab = ZClab, ylab = NA, type = "n")
  mtext(nH2Olab, side = 2, las = 0, line = 3.2, cex = 1.6)
  lmlines()
  mgradox <- ppage("gradoxGS", plot.it = FALSE)
  pgradox <- ppage("gradoxGS", H2O = TRUE, plot.it = FALSE)
  if(getOption("basis") == "QCa") {labdx <- NULL; labdy <- 0.003}
  if(getOption("basis") == "QEC") {labdx <- c(0, 0, -0.0125); labdy <- c(0.003, 0.003, 0.002)}
  pcomp(mgradox, pgradox, reorder = FALSE, add = TRUE, vars = vars, font = 2, labdx = labdx, labdy = labdy)
  # add proteomes from Nif-encoding genomes 20191014
  np <- NifProteomes()
  if(vars == "H2O-ZC") {
    points(np$ZC.mean, np$nH2O.mean, pch = 15)
    lines(np$ZC.mean, np$nH2O.mean, col = "dimgray", lwd = 0.8, lty = 2)
    # add study label and outline most oxidized sample 20200823
    points(np$ZC.mean[4], np$nH2O.mean[4], pch = 0, cex = 1.6)
    text(np$ZC.mean[4], np$nH2O.mean[4], "NF", adj = c(0.5, -1), cex = 0.7, font = 2)
  }
  if(vars == "pIG") {
    points(np$pI.mean, np$GRAVY.mean, pch = 15)
    lines(np$pI.mean, np$GRAVY.mean, col = "dimgray", lwd = 0.8, lty = 2)
  }
  # add text labels
  if(vars == "H2O-ZC") {
    if(getOption("basis") == "QCa") {
      text(-0.213, -1.0885, "hot spring", cex = 0.8, srt = 12)
      text(-0.213, -1.091, "source", cex = 0.8, srt = 12)
      text(-0.172, -1.062, "phototrophic", cex = 0.8, srt = 12)
      text(-0.172, -1.0645, "zone", cex = 0.8, srt = 12)
      text(-0.127, -0.002, "> 3 mm", adj = 0, cex = 0.8)
      text(-0.1395, -1.078, "3 mm", adj = 0, cex = 0.8)
      text(-0.135, -1.0665, "2 mm", adj = 0, cex = 0.8)
      text(-0.132, -1.0555, "1 mm", adj = 0, cex = 0.8)
      text(-0.182, -1.094, "vent", cex = 0.8, srt = 12)
      text(-0.182, -1.097, "  fluids", cex = 0.8, srt = 12)
      text(-0.151, -1.079, "plume", cex = 0.8, srt = 12)
      text(-0.151, -1.0815, "and", cex = 0.8, srt = 12)
      text(-0.151, -1.0840, "  seawater", cex = 0.8, srt = 12)
      text(c(-0.196, -0.222, -0.191, -0.155), c(-1.0455, -1.0645, -1.055, -1.039), c("Nif-D", "Nif-C", "Nif-B", "Nif-A"), adj = 0, cex = 0.8)
      text(-0.218, 0.0365, "NF", cex=0.7, font = 2)
    }
    if(getOption("basis") == "QEC") {
      text(-0.210, -0.7535, "hot spring", cex = 0.8)
      text(-0.210, -0.7565, "source", cex = 0.8)
      text(-0.154, -0.773, "phototrophic", cex = 0.8)
      text(-0.154, -0.776, "zone", cex = 0.8)
      text(-0.191, -0.767, "photosynthetic", cex = 0.8)
      text(-0.191, -0.770, "fringe", cex = 0.8)
      lines(c(-0.156, -0.156), c(-0.771, -0.765))
      text(-0.131, -0.776, "> 3 mm", cex = 0.8)
      text(-0.129, -0.765, "2 mm", cex = 0.8)
      text(-0.126, -0.757, "1 mm", cex = 0.8)
      text(-0.180, -0.746, "vent fluids", cex = 0.8, srt = -23)
      text(-0.156, -0.751, "plume", cex = 0.8)
      text(-0.146, -0.746, "seawater", cex = 0.8)
      lines(c(-0.1475, -0.1475), c(-0.748, -0.7625))
      text(c(-0.197, -0.222, -0.177, -0.151), c(-0.706, -0.708, -0.728, -0.7365), c("Nif-D", "Nif-C", "Nif-B", "Nif-A"), adj = 0, cex = 0.8)
    }
  }
  if(vars == "pIG") {
    text(c(5.00, 7.49, 5.65, 5.68), c(-0.123, -0.130, -0.109, -0.066), c("Nif-D", "Nif-C", "Nif-B", "Nif-A"), adj = 0)
    text(c(5.37, 5.42, 5.57), c(-0.158, -0.176, -0.216), c("1 mm", "2 mm", "3 mm"), adj = 1)
  }
  title("Redox gradients", font.main = 1)
  label.figure("(a)", cex = 1.7, xfrac = 0.04)

  # plot 2: compare ZC and nH2O of proteins in Baltic Sea surface
  plot(c(-0.18, -0.12), ylim, xlab = ZClab, ylab = NA, type = "n")
  mtext(nH2Olab, side = 2, las = 0, line = 3.2, cex = 1.6)
  lmlines()
  mbaltics <- ppage("balticsurface", plot.it = FALSE)
  pbaltics <- ppage("balticsurface", H2O = TRUE, plot.it = FALSE)
  pcomp(mbaltics, pbaltics, reorder = FALSE, add = TRUE, labels.at = NA, vars = vars, pch = list(c(17, 17, 17, 17, 19, 19, 19, 19, 19, 19)))
  # add text labels
  if(vars == "H2O-ZC") {
    if(getOption("basis") == "QCa") {
      text(-0.138, -1.031, "< 6 PSU")
      text(-0.139, -1.060, "> 6 PSU")
    }
    if(getOption("basis") == "QEC") {
      text(-0.140, -0.724, "< 6 PSU")
      text(-0.137, -0.747, "> 6 PSU")
    }
  }
  title("Baltic Sea", font.main = 1, line = 1.2)
  title("surface", font.main = 1, line = 0.2)
  label.figure("(b)", cex = 1.7, xfrac = 0.035)

  # plot 3: compare ZC and nH2O of proteins in Baltic Sea 10-20 m
  plot(c(-0.18, -0.12), ylim, xlab = ZClab, ylab = NA, type = "n")
  mtext(nH2Olab, side = 2, las = 0, line = 3.2, cex = 1.6)
  lmlines()
  mbalticd <- ppage("balticdeep", plot.it = FALSE)
  pbalticd <- ppage("balticdeep", H2O = TRUE, plot.it = FALSE)
  pcomp(mbalticd, pbalticd, reorder = FALSE, add = TRUE, labels.at = NA, vars = vars, pch = list(c(17, 17, 17, 17, 19, 19, 19, 19, 19)))
  # add text labels
  if(vars == "H2O-ZC") {
    if(getOption("basis") == "QCa") {
      text(-0.151, -1.052, "< 6 PSU")
      text(-0.136, -1.080, "> 6 PSU")
    }
    if(getOption("basis") == "QEC") {
      text(-0.148, -0.73, "< 6 PSU")
      text(-0.140, -0.75, "> 6 PSU")
    }
  }
  title("Baltic Sea", font.main = 1, line = 1.2)
  title("chl a max", font.main = 1, line = 0.2)
  label.figure("(c)", cex = 1.7, xfrac = 0.035)

  if(pdf) {
    dev.off()
    addexif("gradH2O3", "nH2O-ZC scatterplots for redox gradients and the Baltic Sea", "Dick et al. (2020)")
  }
}

# nH2O for Baltic Sea metagenome and metatranscriptome in different size fractions 20190715
gradH2O4 <- function(pdf = FALSE, var = NULL) {
  if(pdf) pdf("gradH2O4.pdf", width = 6, height = 3.5)
  nvar <- ifelse(is.null(var), 1, length(var))
  if(nvar==2) par(mfrow = c(2, 3))
  else par(mfrow = c(1, 3))
  par(mar = c(3, 4, 2, 1), mgp = c(3.5, 0.7, 0), las = 1)
  par(cex.lab = 1.5)

  for(i in 1:nvar) {
    if(identical(var[i], "pI")) ylim <- c(5, 7.5) else {
      if(getOption("basis") == "QCa") ylim <- NULL
      if(getOption("basis") == "QEC") ylim <- c(-0.785, -0.715)
    }
    figlab <- c("(a)", "(b)", "(c)")
    if(i==2) figlab <- c("(d)", "(e)", "(f)")

    mplot("Baltic_Sea-0.1s", "iMicrobe_MGP", H2O = TRUE, plottype = "#FF000030",
          col = "red", add.title = FALSE, ylim = ylim, yline = 2.7, var = var[i], srt = NA, ilabel = c(1, 12))
    mplot("Baltic_Sea-0.1s", "SRA_MTP", H2O = TRUE, plottype = "#0000FF30",
          col = "blue", add.title = FALSE, add = TRUE, pch = 15, var = var[i])
    title(quote("0.1-0.8"~mu*m))
    legend("bottomleft", c("metagenomes", "metatranscriptomes"), pch = c(19, 15), col = c("red", "blue"), lty = c(2, 2), bty = "n")
    mtext(quote(phantom(.) %->% phantom("salinity") %->% phantom(.)), 1, 1.2, cex = 0.7)
    mtext("higher\nsalinity", 1, 1.7, cex = 0.7)
    label.figure(figlab[1], cex = 1.7, xfrac = 0.17, yfrac = 0.965)

    mplot("Baltic_Sea-0.8s", "iMicrobe_MGP", H2O = TRUE, plottype = "#FF000030",
          col = "red", add.title = FALSE, ylim = ylim, yline = 2.7, var = var[i], srt = NA, ilabel = c(1, 12))
    mplot("Baltic_Sea-0.8s", "SRA_MTP", H2O = TRUE, plottype = "#0000FF30",
          col = "blue", add.title = FALSE, add = TRUE, pch = 15, var = var[i])
    mtext(quote(phantom(.) %->% phantom("salinity") %->% phantom(.)), 1, 1.2, cex = 0.7)
    mtext("higher\nsalinity", 1, 1.7, cex = 0.7)
    title(quote("0.8-3.0"~mu*m))
    label.figure(figlab[2], cex = 1.7, xfrac = 0.17, yfrac = 0.965)

    mplot("Baltic_Sea-3.0s", "iMicrobe_MGP", H2O = TRUE, plottype = "#FF000030",
          col = "red", add.title = FALSE, ylim = ylim, yline = 2.7, var = var[i], srt = NA, ilabel = c(1, 12))
    mplot("Baltic_Sea-3.0s", "SRA_MTP", H2O = TRUE, plottype = "#0000FF30",
          col = "blue", add.title = FALSE, add = TRUE, pch = 15, var = var[i])
    mtext(quote(phantom(.) %->% phantom("salinity") %->% phantom(.)), 1, 1.2, cex = 0.7)
    mtext("higher\nsalinity", 1, 1.7, cex = 0.7)
    title(quote("3.0-200"~mu*m))
    label.figure(figlab[3], cex = 1.7, xfrac = 0.17, yfrac = 0.965)
  }

  if(pdf) {
    dev.off()
    addexif("gradH2O4", "nH2O for Baltic Sea metagenome and metatranscriptome in different size fractions", "Dick et al. (2020)")
  }
}

# nH2O vs ZC for freshwater, marine, and hypersaline environments 20191004
# make separate plot for sediments 20191006
# add Amazon River 20191007
# remove sediments and add GRAVY - pI plots 20191027
gradH2O5 <- function(pdf = FALSE) {

  if(pdf) pdf("gradH2O5.pdf", width = 8, height = 5)
  layout(matrix(1:6, nrow = 2))
  par(las = 1, mar = c(4, 4.2, 2, 1), mgp = c(2.5, 0.8, 0))
  par(cex.lab = 1.5)
  xlim <- c(-0.2, -0.08)
  if(getOption("basis") == "QCa") ylim <- c(-1.15, -1)
  if(getOption("basis") == "QEC") ylim <- c(-0.8, -0.72)
  nH2Olab <- expression(italic(n)[H[2] * O])
  ZClab <- expression(italic(Z)[C])

  # plots 1-2: Amazon river metagenome
  mout <- ppage("amazon", plot.it = FALSE)
  pout <- ppage("amazon", H2O = TRUE, plot.it = FALSE)
  plot(xlim, ylim, xlab = ZClab, ylab = NA, type = "n")
  mtext(nH2Olab, side = 2, las = 0, line = 2.9)
  lmlines(0.02)
  pcomp(mout, pout, add = TRUE, lty = 0, labels.at = NA)
  hullfun(mout, pout, c(1, 3), "green3", c("riverFL", "riverPA"))
  hullfun(mout, pout, c(1, 3), "blue", c("plumeFL", "plumePA"))
  if(getOption("basis") == "QCa") {
    text(c(-0.129, -0.132), c(-1.055, -1.095), c("river", "plume"))
    rect(-0.20, -1.05, -0.18, -0.995, col = "white", border = NA)
    rect(-0.20, -1.038, -0.10, -0.995, col = "white", border = NA)
    legend("topleft", c("free-living", "particle-associated"), pch = c(20, 15), col = "blue", bty = "n", inset = c(0.2, 0))
    legend("topleft", legend = c(NA, NA), pch = c(21, 22), col = "green3", pt.cex = c(0.7, 1), bty = "n", inset = c(0.14, 0))
    text(c(-0.18, -0.17), c(-1.027, -1.027), c("river", "plume"), srt = 40, adj = 1)
  }
  if(getOption("basis") == "QEC") {
    text(c(-0.128, -0.130), c(-0.755, -0.775), c("river", "plume"))
    rect(-0.181, -0.733, -0.09, par("usr")[4] - 0.0003, col = "white", border = NA)
    rect(-0.11, -0.737, -0.094, par("usr")[4] - 0.0003, col = "white", border = NA)
    legend("topright", legend = c(NA, NA), pch = c(21, 22), col = "green3", pt.cex = c(0.7, 1), bty = "n", inset = c(0.1, 0))
    legend("topright", legend = c(NA, NA), pch = c(20, 15), col = "blue", bty = "n", inset = c(0.17, 0))
    legend("topright", c("free-living", "particle-associated"), adj = 1, bty = "n", inset = c(-0.2, 0))
    text(c(-0.107, -0.097), c(-0.735, -0.735), c("plume", "river"), srt = 40, adj = 1)
  }
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
  plot(xlim, ylim, xlab = ZClab, ylab = NA, type = "n")
  mtext(nH2Olab, side = 2, las = 0, line = 2.9)
  lmlines(0.02)
  pcomp(mout, pout, "MT", add = TRUE, lty = 0, labels.at = NA)
  hullfun(mout, pout, c(2, 4), "green3", c("riverFL", "riverPA"))
  hullfun(mout, pout, c(2, 4), "blue", c("plumeFL", "plumePA"))
  if(getOption("basis") == "QCa") text(c(-0.125, -0.12), c(-1.035, -1.073), c("river", "plume"))
  if(getOption("basis") == "QEC") text(c(-0.166, -0.122), c(-0.745, -0.774), c("river", "plume"))
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
  plot(xlim, ylim, xlab = ZClab, ylab = NA, type = "n")
  mtext(nH2Olab, side = 2, las = 0, line = 2.9)
  lmlines(0.02)
  pcomp(moutE, poutE, add = TRUE, lty = 0, labels.at = NA)
  hullfun(moutE, poutE, 1, "green3")
  hullfun(moutE, poutE, 2, "blue")
  # add hypersaline water data
  moutH <- ppage("hypersaline", plot.it = FALSE)
  poutH <- ppage("hypersaline", H2O = TRUE, plot.it = FALSE)
  pcomp(moutH, poutH, reorder = FALSE, add = TRUE, labdx = c(-0.002, 0.030, 0.056), labdy = c(-0.005, -0.011, 0.004))
  hullfun(moutH, poutH, 1:3, "turquoise3")
  if(getOption("basis") == "QCa") {
    text(c(-0.18, -0.148, -0.10), c(-1.042, -1.095, -1.088), c("freshwater", "marine", "hypersaline"))
    rect(-0.148, -1.15, -0.076, -1.11, col = "white", border = NA)
  }
  if(getOption("basis") == "QEC") {
    text(c(-0.158, -0.183), c(-0.726, -0.758), c("freshwater", "marine"))
    rect(-0.148, par("usr")[3] + 0.0003, par("usr")[2] - 0.0004, -0.779, col = "white", border = NA)
  }
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
    addexif("gradH2O5", "nH2O vs ZC for freshwater, marine, and hypersaline environments", "Dick et al. (2020)")
  }
}

# nH2O-ZC and GRAVY-pI plots for Baltic Sea and Rodriguez-Brito et al. data 20200421
gradH2O6 <- function(pdf = FALSE) {
  if(pdf) pdf("gradH2O6.pdf", width = 8, height = 8)
  layout(matrix(1:4, nrow = 2))
  par(mar = c(4, 4.5, 2, 1), las = 1, cex = 1.2)
  par(cex.lab = 1.3)
  nH2Olab <- expression(italic(n)[H[2] * O])
  ZClab <- expression(italic(Z)[C])
  if(getOption("basis") == "QCa") ylim <- c(-1.12, -1.02)
  if(getOption("basis") == "QEC") ylim <- c(-0.78, -0.7)

  # Baltic Sea nH2O - ZC
  mout <- ppage("balticsurface", plot.it = FALSE)
  pout <- ppage("balticsurface", H2O = TRUE, plot.it = FALSE)
  plot(c(-0.2, -0.08), ylim, xlab = ZClab, ylab = NA, type = "n")
  mtext(nH2Olab, side = 2, las = 0, line = 3.2, cex = 1.5)
  lmlines(0.02)
  pcomp(mout, pout, add = TRUE, reorder = FALSE, font = 2, labdy = 0.003, labels.at = NA)
  if(getOption("basis") == "QCa") text(c(-0.128, -0.128), c(-1.035, -1.075), c("< 6 PSU", "> 6 PSU"))
  if(getOption("basis") == "QEC") text(c(-0.128, -0.128), c(-0.725, -0.753), c("< 6 PSU", "> 6 PSU"))
  title("Baltic Sea", font.main = 1)
  label.figure("(a)", cex = 1.7)
  # Baltic Sea GRAVY - pI
  pcomp(mout, pout, reorder = FALSE, vars = "pIG", yline = 3.2, cex.ylab = 1.5, font = 2, labdy = 0.003, labels.at = NA)
  text(c(7.5, 7.2), c(-0.11, -0.23), c("< 6 PSU", "> 6 PSU"))
  title("Baltic Sea", font.main = 1)
  label.figure("(c)", cex = 1.7)

  # Rodriguez-Brito et al. nH2O - ZC
  mout <- ppage("socal", plot.it = FALSE)
  pout <- ppage("socal", H2O = TRUE, plot.it = FALSE)
  plot(c(-0.2, -0.08), ylim, xlab = ZClab, ylab = NA, type = "n")
  mtext(nH2Olab, side = 2, las = 0, line = 3.2, cex = 1.5)
  lmlines(0.02)
  pcomp(mout, pout, add = TRUE, reorder = FALSE, font = 2, labdy = 0.003, labels.at = NA)
  if(getOption("basis") == "QCa") text(c(-0.185, -0.185, -0.155), c(-1.08, -1.06, -1.035), c("FW", "LS", "MS-HS"))
  if(getOption("basis") == "QEC") text(c(-0.185, -0.128, -0.125), c(-0.76, -0.755, -0.725), c("FW", "LS", "MS-HS"))
  title("Fish ponds and salterns", font.main = 1)
  label.figure("(b)", cex = 1.7)
  # Rodriguez-Brito et al. GRAVY - pI
  pcomp(mout, pout, reorder = FALSE, vars = "pIG", yline = 3.2, cex.ylab = 1.5, font = 2, labdy = 0.003, labels.at = NA)
  text(c(8.3, 6, 4.5), c(-0.18, -0.19, -0.115), c("FW", "LS", "MS-HS"))
  title("Fish ponds and salterns", font.main = 1)
  label.figure("(d)", cex = 1.7)
  if(pdf) {
    dev.off()
    addexif("gradH2O6", "nH2O-ZC and GRAVY-pI plots for Baltic Sea and Rodriguez-Brito et al. data", "Dick et al. (2020)")
  }
}

# differential gene and protein expression; time-course experiments and NaCl or organic solutes 20200420
gradH2O7 <- function(pdf = FALSE) {
  if(pdf) pdf("gradH2O7.pdf", width = 8, height = 6)
  mat <- matrix(c(1,1,1,1,1,1,1,1, 2,2,2,2, 3,3,3, 4,4,4, 5,5,5, 6,6,6), nrow = 2, byrow = TRUE)
  layout(mat, heights = c(3, 2))
  par(mar = c(4, 4, 3, 1), mgp = c(2.5, 1, 0))
  par(cex.lab = 1.3)

  # function to plot an arrow partway along a line
  mkarrow <- function(row1, comptab, icol, col = 1, frac = 0.5) {
    # row1 is the starting point
    x1 <- comptab[row1, icol]
    y1 <- comptab$nH2O.diff[row1]
    # the next row is the end of the full line (not the arrow)
    x2 <- comptab[row1 + 1, icol]
    y2 <- comptab$nH2O.diff[row1 + 1]
    # calculate slope
    m <- (y2 - y1) / (x2 - x1)
    # calculate value of x and y on the line (arrow tip)
    x <- x1 + (x2 - x1) * frac
    y <- y1 + m * (x - x1)
    # calculate value of x0 and y0 (start of arrow - just a short line)
    x0 <- x1 + (x2 - x1) * frac * 0.99
    y0 <- y1 + m * (x0 - x1)
    # draw the arrow
    suppressWarnings(arrows(x0, y0, x, y, length = 0.1, angle = 20, lwd = 2, col = col))
  }

  # function to add arrows and points for Delta nH2O vs logtime or solute
  DnH2O <- function(comptab, ndat, pch = 21, pt.text = NULL, column = "logtime", arrows = TRUE, col = 1) {
    n <- 0
    icol <- grep(column, colnames(comptab))
    for(i in 1:length(ndat)) {
      idat <- n + 1:ndat[i]
      lines(comptab[idat, icol], comptab$nH2O.diff[idat], col = col)
      points(comptab[idat, icol], comptab$nH2O.diff[idat], pch = pch, bg = "white", cex = 2.5, col = col)
      # add labels
      if(is.null(pt.text)) labels <- letters[idat] else labels <- pt.text[idat]
      text(comptab[idat, icol], comptab$nH2O.diff[idat], labels, cex = 0.85, col = col)
      if(arrows) lapply(head(idat, -1), mkarrow, comptab = comptab, icol = icol, col = col)
      n <- n + ndat[i]
    }
  }

  # Read CSV files with results of chemical analysis for differential expression
  osmotic_gene <- read.csv(system.file(paste0("vignettes/osmotic_gene_", getOption("basis"), ".csv"), package = "JMDplots"))
  osmotic_bact <- read.csv(system.file(paste0("vignettes/osmotic_bact_", getOption("basis"), ".csv"), package = "canprot"))
  DnH2Olab <- quote(Delta*italic(n)[H[2]*O])
  log10timelab <- quote(log[10]*("time, minutes"))

  # plot A: time-course experiments
  if(getOption("basis") == "QCa") ylim <- c(-0.12, 0.12)
  if(getOption("basis") == "QEC") ylim <- c(-0.07, 0.10)
  plot(c(0.5, 4), ylim, xlab = log10timelab, ylab = DnH2Olab, type = "n")
  abline(h = 0, lty = 2, col = "gray40")
  # transcriptomic data
  Ttime <- list(
    KKG = c("KKG+14_Gene_30min", "KKG+14_Gene_80min", "KKG+14_Gene_310min"),
    SLM = c("SLM+14_5", "SLM+14_30", "SLM+14_60"),
    FRH = c("FRH+15_NaCl_1h", "FRH+15_NaCl_6h", "FRH+15_NaCl_24h"),
    HLL = c("HLL17_45min", "HLL17_14h"),
    QHT = c("QHT+13_Gene.24.h", "QHT+13_Gene.48.h", "QHT+13_Gene.72.h")
  )
  comptab <- osmotic_gene[match(unlist(Ttime), osmotic_gene$dataset), ]
  ndat <- sapply(Ttime, length)
  # make nH2O-time plot 20200824
  Ttime <- c(30, 80, 310, 5, 30, 60, 1*60, 6*60, 24*60, 45, 14*60, 24*60, 48*60, 72*60)
  comptab <- cbind(comptab, logtime = log10(Ttime))
  DnH2O(comptab, ndat)
  # proteomic data
  Ptime <- list(
    KKG = c("KKG+14_Protein_30min", "KKG+14_Protein_80min", "KKG+14_Protein_310min"),
    QHT = c("QHT+13_Protein.24.h", "QHT+13_Protein.48.h")
  )
  comptab <- osmotic_bact[match(unlist(Ptime), osmotic_bact$dataset), ]
  ndat <- sapply(Ptime, length)
  Ptime <- c(30, 80, 310, 24*60, 48*60)
  comptab <- cbind(comptab, logtime = log10(Ptime))
  DnH2O(comptab, ndat, pch = 22, pt.text = c("a", "b", "c", "l", "m"), col = 4)
  legend("topleft", c("transcriptomic", "proteomic"), pch = c(21, 22), col = c(1, 4), pt.cex = 1.5, bty = "n")
  label.figure("(a)", cex = 1.8)
  mtext("Time-course experiments", line = 1)

  # plot B: NaCl or organic solutes
  if(getOption("basis") == "QCa") ylim <- c(-0.07, 0.11)
  plot(c(-0.2, 1.2), ylim, xlab = "Solute", ylab = DnH2Olab, type = "n", xaxt = "n")
  axis(1, at = c(0, 1), labels = c("NaCl", "Organic"))
  abline(h = 0, lty = 2, col = "gray40")
  # transcriptomic data
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
  Tsolute <- as.numeric(!grepl("NaCl", comptab$description)) - 0.05
  Tsolute[1] <- Tsolute[1] + 0.08
  Tsolute[4] <- Tsolute[4] - 0.1
  Tsolute[8] <- Tsolute[8] + 0.1
  Tsolute[11] <- Tsolute[11] + 0.1
  Tsolute[13] <- Tsolute[13] - 0.1
  comptab <- cbind(comptab, solute = Tsolute)
  DnH2O(comptab, ndat, pt.text = LETTERS[1:sum(ndat)], column = "solute", arrows = FALSE)
  # proteomic data
  Psolute <- list(
    KLB = c("KLB+15_prot-NaCl", "KLB+15_prot-suc"),
    SKV = c("SKV+16_Osmotic.stress.glucose_LB", "SKV+16_Glucose_LB")
  )
  comptab <- osmotic_bact[match(unlist(Psolute), osmotic_bact$dataset), ]
  ndat <- sapply(Psolute, length)
  Psolute <- c(0, 1, 0, 1) + 0.05
  comptab <- cbind(comptab, solute = Psolute)
  DnH2O(comptab, ndat, pch = 22, pt.text = c("E", "F", "O", "P"), column = "solute", arrows = FALSE, col = 4)
  label.figure("(b)", cex = 1.8)
  mtext("NaCl or organic solutes", line = 1)

  # Plot C: transcriptomes nH2O-ZC
  # Add the p-values to the axis labels
  p.ZC <- formatC(t.test(osmotic_gene$ZC.down, osmotic_gene$ZC.up, paired = TRUE)$p.value, 3, format = "f")
  p.nH2O <- formatC(t.test(osmotic_gene$nH2O.down, osmotic_gene$nH2O.up, paired = TRUE)$p.value, 3, format = "f")
  labtext <- c(bquote(italic(p)==.(p.ZC)), bquote(italic(p)==.(p.nH2O)))
  diffplot(osmotic_gene, pt.text = NA, contour = FALSE, labtext = labtext, cex = 1.2)
  points(mean(osmotic_gene$ZC.diff), mean(osmotic_gene$nH2O.diff), pch = 19, cex = 1.7)
  legend("topleft", c("non-halophiles", "mean"), pch = c(1, 19), bty = "n")
  label.figure("(c)", cex = 1.8, yfrac = 0.85)

  # Plot D: transcriptomes GRAVY-pI
  p.pI <- formatC(t.test(osmotic_gene$pI.down, osmotic_gene$pI.up, paired = TRUE)$p.value, 3, format = "f")
  p.GRAVY <- formatC(t.test(osmotic_gene$GRAVY.down, osmotic_gene$GRAVY.up, paired = TRUE)$p.value, 3, format = "f")
  labtext <- c(bquote(italic(p)==.(p.pI)), bquote(italic(p)==.(p.GRAVY)))
  diffplot(osmotic_gene, vars = c("pI", "GRAVY"), pt.text = NA, contour = FALSE, labtext = labtext, cex = 1.2)
  points(mean(osmotic_gene$pI.diff), mean(osmotic_gene$GRAVY.diff), pch = 19, cex = 1.7)
  label.figure("(d)", cex = 1.8, yfrac = 0.85)
  mtext("Proteins coded by differentially expressed genes", adj = 1, line = 1.5)

  # Get proteomic data for halophiles
  osmotic_halo <- read.csv(system.file(paste0("vignettes/osmotic_halo_", getOption("basis"), ".csv"), package = "canprot"))
  # Combine halophile and non-halophile data for hyperosmotic stress proteomics
  osmotic_halo <- osmotic_halo[osmotic_halo$tags!="hypoosmotic", ]
  alldat <- rbind(osmotic_bact, osmotic_halo)
  pch <- c(rep(0, nrow(osmotic_bact)), rep(2, nrow(osmotic_halo)))
  col <- c(rep(4, nrow(osmotic_bact)), rep(2, nrow(osmotic_halo)))

  # Plot E: proteomes nH2O-ZC
  p.ZC <- formatC(t.test(alldat$ZC.down, alldat$ZC.up, paired = TRUE)$p.value, 3, format = "f")
  p.nH2O <- formatC(t.test(alldat$nH2O.down, alldat$nH2O.up, paired = TRUE)$p.value, 3, format = "f")
  labtext <- c(bquote(italic(p)==.(p.ZC)), bquote(italic(p)==bold(.(p.nH2O))))
  diffplot(alldat, pt.text = NA, contour = FALSE, labtext = labtext, cex = 1.2, pch = pch, col = col)
  points(mean(alldat$ZC.diff), mean(alldat$nH2O.diff), pch = 19, cex = 1.7)
  label.figure("(e)", cex = 1.8, yfrac = 0.85)

  # Plot F: proteomes GRAVY-pI
  p.pI <- formatC(t.test(alldat$pI.down, alldat$pI.up, paired = TRUE)$p.value, 3, format = "f")
  p.GRAVY <- formatC(t.test(alldat$GRAVY.down, alldat$GRAVY.up, paired = TRUE)$p.value, 3, format = "f")
  labtext <- c(bquote(italic(p)==.(p.pI)), bquote(italic(p)==bold(.(p.GRAVY))))
  diffplot(alldat, vars = c("pI", "GRAVY"), pt.text = NA, contour = FALSE, labtext = labtext, cex = 1.2, pch = pch, col = col)
  points(mean(alldat$pI.diff), mean(alldat$GRAVY.diff), pch = 19, cex = 1.7)
  legend("topleft", c("non-halophiles", "halophiles", "mean"), pch = c(0, 2, 19), col = c(4, 2, 1), bty = "n")
  label.figure("(f)", cex = 1.8, yfrac = 0.85)
  mtext("Differentially expressed proteins             ", adj = 1, line = 1.5)

  if(pdf) {
    dev.off()
    addexif("gradH2O7", "differential gene and protein expression; time-course experiments and NaCl or organic solutes", "Dick et al. (2020)")
  }
}

# Calculate ZC and nH2O of proteomes encoding different Nif homologs (Poudel et al., 2018)  20191014
# Also return individual ZC values of all proteomes 20220531
# Also return amino acid composition of all proteomes 20220603
NifProteomes <- function() {
  # Read file with amino acid compositions
  AAfile <- system.file("extdata/utegig/Nif_homolog_AA.csv", package = "JMDplots")
  AA <- read.csv(AAfile, as.is = TRUE)
  # The Nif types, arranged from anaerobic to aerobic
  types <- c("Nif-D", "Nif-C", "Nif-B", "Nif-A")
  # Assemble the chemical metrics
  ZC.mean <- ZC.sd <- nH2O.mean <- nH2O.sd <- GRAVY.mean <- GRAVY.sd <- pI.mean <- pI.sd <- numeric()
  ZClist <- list()
  for(i in 1:length(types)) {
    type <- types[i]
    # Get the taxids for genomes with this Nif type
    iAA <- AA$protein == type
    taxid <- AA$organism[iAA]
    # Get the amino acid compositions of all genomes with this Nif type
    AAcomp <- AA[iAA, ]
    # Calculate ZC and nH2O
    ZC <- ZCAA(AAcomp)
    ZC.mean <- c(ZC.mean, mean(ZC))
    ZC.sd <- c(ZC.sd, sd(ZC))
    H2O <- H2OAA(AAcomp)
    nH2O.mean <- c(nH2O.mean, mean(H2O))
    nH2O.sd <- c(nH2O.sd, sd(H2O))
    GRAVY <- GRAVY(AAcomp)
    GRAVY.mean <- c(GRAVY.mean, mean(GRAVY))
    GRAVY.sd <- c(GRAVY.sd, sd(GRAVY))
    pI <- pI(AAcomp)
    pI.mean <- c(pI.mean, mean(pI))
    pI.sd <- c(pI.sd, sd(pI))
    # Store ZC values 20220531
    ZClist[[i]] <- ZC
  }
  # Return values
  names(ZClist) <- types
  list(types = types, ZC.mean = ZC.mean, ZC.sd = ZC.sd, nH2O.mean = nH2O.mean, nH2O.sd = nH2O.sd,
       GRAVY.mean = GRAVY.mean, GRAVY.sd = GRAVY.sd, pI.mean = pI.mean, pI.sd = pI.sd,
       ZClist = ZClist, AA = AA)
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

# Plot R-squared of nH2O-ZC and nO2-ZC fits 20200813
plotbasisfun <- function(zoom = FALSE) {
  # read data
  file <- system.file("extdata/gradH2O/AAbasis.csv", package = "JMDplots")
  AAbasis <- read.csv(file, as.is = TRUE)

  # set up plot
  xlab <- quote(italic(R)^2~"of"~italic(n)[O[2]] - italic(Z)[C]~"fits")
  ylab <- quote(italic(R)^2~"of"~italic(n)[H[2]*O] - italic(Z)[C]~"fits")
  plot(c(0, 1), c(0, 0.8), xlab = xlab, ylab = ylab, type = "n")
  x <- extendrange(c(0.9, 0.93))
  y <- extendrange(c(0, 0.03))
  abbrv <- AAbasis$abbrv
  O2 <- AAbasis$O2.R2
  H2O <- AAbasis$H2O.R2
  # use larger blue point for QEC 20201011
  cex <- rep(0.3, length(O2))
  col <- rep(1, length(O2))
  cex[AAbasis$abbrv == "CEQ"] <- 0.6
  col[AAbasis$abbrv == "CEQ"] <- 4
  points(O2, H2O, cex = cex, col = col, pch = 19)
}

# Add nH2O-ZC guidelines parallel to regression for amino acids 20200819
lmlines <- function(step = 0.01) {
  if(FALSE) {
    # Calculate ZC of the amino acids
    aa <- aminoacids("")
    ZC.aa <- ZC(info(aa, "aq"))
    # Load amino acids with QCa or QEC basis 20200914
    if(getOption("basis") == "QCa") basis(c("glutamine", "cysteine", "acetic acid", "H2O", "O2"))
    if(getOption("basis") == "QEC") basis(c("glutamine", "glutamic acid", "cysteine", "H2O", "O2"))
    species(aa)
    # Make linear regression
    lm <- lm(species()$H2O ~ ZC.aa)
    coef <- coef(lm)
    # Clear species!
    reset()
  } else {
    # Use previously computed intercept and slope 20200920
    if(getOption("basis") == "QCa") coef <- c(-0.4830396, 0.1579203)
    if(getOption("basis") == "QEC") coef <- c(-0.1242780, -0.3088251)
  }
  x <- par("usr")[1:2]
  y <- coef[1] + coef[2] * x
  for(dy in seq(-0.48, -0.74, -step)) lines(x, y + dy, col = "gray80")
  # Add box so ends of lines don't cover plot edges 20201007
  box()
}

