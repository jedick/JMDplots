# JMDplots/bison.R
# Plots from hot spring (Bison Pool) papers (2011, 2013)
# Code moved from CHNOSZ/vignettes/hotspring.Rnw 20200712

# Measured temperature and pH
bison1 <- function() {

  par(mfrow = c(1, 3), mar = c(3.5, 3.5, 1, 1), mgp = c(2.4, 1, 0), las = 1)

  # T plot
  plot(extendrange(distance), extendrange(bison.T), xlab = "Distance (m)", ylab = axis.label("T"), type = "n")
  lines(xpoints, Tfun(xpoints))
  points(distance, bison.T, pch = 21, bg = "white", cex = 2)
  points(distance, bison.T, pch = 21, bg = site.cols.alpha, cex = 2)
  text(distance, bison.T, 1:5, cex = 0.8)

  # pH plot
  plot(extendrange(distance), extendrange(bison.pH), xlab = "Distance (m)", ylab = "pH", type = "n")
  lines(xpoints, pHfun(xpoints))
  points(distance, bison.pH, pch = 21, bg = "white", cex = 2)
  points(distance, bison.pH, pch = 21, bg = site.cols.alpha, cex = 2)
  text(distance, bison.pH, 1:5, cex = 0.8)

  # O2 plot
  plot(extendrange(distance), extendrange(bison.O2), xlab = "Distance (m)", ylab = quote(O[2]~"(mg/L)"), type = "n")
  lines(xpoints, O2fun(xpoints))
  points(distance, bison.O2, pch = 21, bg = "white", cex = 2)
  points(distance, bison.O2, pch = 21, bg = site.cols.alpha, cex = 2)
  text(distance, bison.O2, 1:5, cex = 0.8)

}

# Carbon oxidation state of proteins
bison2 <- function(plots = 1:2, add.titles = TRUE) {
  if(length(plots) == 2) par(mfrow = c(1, 2), las = 1)

  if(1 %in% plots) {
    # Proteins groups by functional annotations (2011 plot)
    plot(0, 0, xlim = c(-0.5, 5), ylim = c(-0.35, -0.10), xlab = "Site", xaxt = "n", ylab = axis.label("ZC"))
    axis(1, at = 1:5)
    col <- c(4, rep(1, 20))
    lwd <- c(4, rep(1, 20))
    lty <- c(1, rep(2, 20))
    clab <- c("hydrolase", "overall", "protease", "oxidoreductase", "transport", "membrane", "permease")
    pf.annot <- protein.formula(aa.annot)
    ZC.annot <- ZC(pf.annot)
    for(i in 1:length(classes)) {
      lines(1:5, ZC.annot[(1:5)+5*(i-1)], col = col[i], lwd = lwd[i], lty = lty[i])
      if(classes[i]=="overall") text(0.95, ZC.annot[1+5*(i-1)], "all proteins", adj = 1, font = 2)
      else if(classes[i] %in% clab) text(0.95, ZC.annot[1+5*(i-1)], classes[i], adj = 1)
    }
    if(add.titles) title(main = "Annotations")
  }

  if(2 %in% plots) {
    # Proteins grouped by phyla (2013 plot)
    # Colorful revision made on 20171217 (moved here from CHNOSZ/demo/bison.R)
    pf.phyla <- protein.formula(aa.phyla)
    ZC.phyla <- ZC(pf.phyla)
    plot(0, 0, xlim = c(1, 5), ylim = c(-0.23, -0.14), xlab = "Site", ylab = NA)
    mtext(axis.label("ZC"), side = 2, line = 3, las = 0)
    for(i in 1:length(phyla.abc)) {
      # skip Euryarchaeota because it occurs at one location, on top of Dein.-Thermus and Firmicutes
      if(phyla.abc[i]=="Euryarchaeota") next
      # which of the model proteins correspond to this phylum
      iphy <- which(aa.phyla$organism==phyla.abc[i])
      # the locations (of 1, 2, 3, 4, 5) where this phylum is found
      ilocs <- match(aa.phyla$protein[iphy], sitenames)
      # the plotting symbol: determined by alphabetical position of the phylum
      points(ilocs, ZC.phyla[iphy], pch = i-1, cex = 1.2)
      # lines to connect the phyla
      lines(ilocs, ZC.phyla[iphy], type = "c", col = phyla.cols[i], lwd = 2)
    }
    text(c(4.75, 2.0, 4.0, 4.0, 4.0, 2.0, 3.0, NA, 2.9, 1.3, 3.0),
         c(-0.146, -0.224, -0.161, -0.184, -0.145, -0.201, -0.144, NA, -0.176, -0.158, -0.192),
         phyla.abbrv, cex = 0.9)
    if(add.titles) title(main = "Major phyla")
  }

}

# Chemical affinities
bison3 <- function() {
  setup.basis()
  ip.annot <- add.protein(aa.annot)
  species("overall", sitenames)
  species(1:5, 0)
  a <- affinity(T = bison.T, pH = bison.pH, H2 = get.logaH2(bison.T))
  # divide by protein length to get per-residue affinities
  pl <- protein.length(ip.annot[1:5])
  a.res <- t(as.data.frame(a$values))/pl
  colnames(a.res) <- paste0("site", 1:5)
  rownames(a.res) <- paste0("reaction", 1:5)
  a.res
}

# Relative stabilities along a temperature and chemical gradient
bison4 <- function() {
  setup.basis()
  add.protein(aa.annot)
  species("overall", sitenames)
  species(1:5, 0)
  Tlim <- c(50, 100)
  par(mfrow = c(1, 2))
  # first plot
  a <- affinity(T = Tlim, H2 = c(-7, -4))
  diagram(a, fill = NULL, names = as.character(1:5), normalize = TRUE)
  lines(Tlim, get.logaH2(Tlim), lty = 3, col = 4, lwd = 2)
  text(68, -6.8, "Equation 2")
  # second plot
  species(1:5, -3)
  xT <- Tfun(xpoints)
  xpH <- pHfun(xpoints)
  xH2 <- get.logaH2(xT)
  a <- affinity(T = xT, pH = xpH, H2 = xH2)
  a$vars[1] <- "Distance, m"
  a$vals[[1]] <- xpoints
  e <- equilibrate(a, normalize = TRUE)
  diagram(e, legend.x = NULL, lty = 1, lwd = 2)
  diagram(e, legend.x = NULL, col = site.cols, lwd = 2, add = TRUE)
  # labels for the lines
  text(c(4, 8, 14, 19, 18), c(-2.81, -2.89, -3.06, -2.85, -2.94), 1:5)
}

# Comparing old and new group additivity parameters
bison5 <- function() {
  par(mfrow = c(2, 3))
  for(j in 1:2) {
    # use old parameters for first row and current ones for second row
    reset()
    if(j==1) {
      OldAAfile <- "extdata/OBIGT/OldAA.csv"
      add.OBIGT(system.file(OldAAfile, package = "JMDplots"))
    }
    # setup basis species and proteins
    setup.basis()
    ip.annot <- add.protein(aa.annot)
    # make the plots
    for(annot in c("overall", "transferase", "synthase")) {
      ip <- ip.annot[aa.annot$protein==annot]
      a <- affinity(T = c(50, 100), H2 = c(-7, -4), iprotein = ip)
      diagram(a, fill = NULL, names = as.character(1:5), normalize = TRUE)
      # add logaH2-T line
      lines(par("usr")[1:2], get.logaH2(par("usr")[1:2]), lty=3, col = 4, lwd = 2)
      # add a title
      title(main=annot)
    }
  }
}

# Metastable equilibrium model for relative abundances
bison6 <- function(plot.it = TRUE) {
  ip.phyla <- add.protein(aa.phyla)
  if(plot.it) layout(matrix(1:6, ncol=3), heights=c(2, 1))
  equil.results <- list()
  for(i in 1:5) {
    # get the equilibrium degrees of formation and the optimal logaH2
    ae <- alpha.equil(i, ip.phyla)
    equil.results[[i]] <- ae
    if(i %in% c(1, 3, 5) & plot.it) {
      iphy <- match(colnames(ae$alpha), phyla.abc)
      # top row: equilibrium degrees of formation
      thermo.plot.new(xlim = range(ae$H2vals), ylim = c(0, 0.5), xlab = axis.label("H2"),
        ylab=expression(alpha[equil]), yline = 2, cex.axis = 1, mgp = c(1.8, 0.3, 0))
      these.cols <- phyla.cols[match(colnames(ae$alpha), phyla.abc)]
      for(j in 1:ncol(ae$alpha)) {
        lines(ae$H2vals, ae$alpha[, j], lwd = 1.5)
        lines(ae$H2vals, ae$alpha[, j], lty = phyla.lty[iphy[j]], col = these.cols[j])
        ix <- seq(1, length(ae$H2vals), length.out = 11)
        ix <- head(tail(ix, -1), -1)
        points(ae$H2vals[ix], ae$alpha[, j][ix], pch = iphy[j]-1, col = these.cols[j])
      }
      title(main=paste("site", i))
      legend("topleft", legend = phyla.abbrv[iphy], pch = ".", bg = "white", lty = 1, lwd = 1.5)
      legend("topleft", legend = rep("", length(iphy)), pch = iphy-1, lty = phyla.lty[iphy], col = these.cols, bty = "n")
      # bottom row: Gibbs energy of transformation and position of minimum
      thermo.plot.new(xlim = range(ae$H2vals), ylim = c(0, 1/log(10)), xlab = axis.label("H2"),
        ylab = expr.property("DGtr/2.303RT"), yline = 2, cex.axis = 1, mgp = c(1.8, 0.3, 0))
      lines(ae$H2vals, ae$DGtr)
      abline(v = ae$logaH2.opt, lty = 2, lwd = 2, col = 3)
      abline(v=get.logaH2(bison.T[i]), lty = 3, lwd = 2, col = 4)
      if(i==1) legend("bottomleft", lty = c(3, 2), lwd = c(2, 2), col = c(4, 3),
        bg = "white", legend = c("Equation 2", "Optimized"))
    }
  }
  invisible(equil.results)
}

# Activity of hydrogen comparison
bison7 <- function(equil.results) {
  par(mfrow = c(1, 2), mar = c(3.5, 3.5, 1, 1), mgp = c(2.5, 1, 0), las = 1)

  # Potential of silver-silver chloride electrode in saturated KCl
  # (Bard et al., 1985; http://www.worldcat.org/oclc/12106344)
  E.AgAgCl <- function(T) {
    0.23737 - 5.3783e-4 * T - 2.3728e-6 * T^2 - 2.2671e-9 * (T+273)
  }
  # Meter readings at Bison Pool and Mound Spring
  ORP <- c(-258, -227, -55, -58, -98, -41)
  T.ORP <- c(93.9, 87.7, 75.7, 70.1, 66.4, 66.2)
  pH.ORP <- c(8.28, 8.31, 7.82, 7.96, 8.76, 8.06)
  # Convert ORP to Eh, then pe, then logaH2
  Eh <- ORP/1000 + E.AgAgCl(T.ORP)
  pe <- convert(Eh, "pe", T = convert(T.ORP, "K"))
  logK.ORP <- subcrt(c("e-", "H+", "H2"), c(-2, -2, 1), T = T.ORP)$out$logK
  logaH2.ORP <- logK.ORP - 2*pe - 2*pH.ORP

  # Sulfide and sulfate concentrations from 2005
  loga.HS <- log10(c(4.77e-6, 2.03e-6, 3.12e-7, 4.68e-7, 2.18e-7))
  loga.SO4 <- log10(c(2.10e-4, 2.03e-4, 1.98e-4, 2.01e-4, 1.89e-4))
  # Convert sulfide/sulfate ratio to logaH2
  logK.S <- subcrt(c("HS-", "H2O", "SO4-2", "H+", "H2"), c(-1, -4, 1, 1, 4), T=bison.T)$out$logK
  logaH2.S <- (logK.S + bison.pH - loga.SO4 + loga.HS) / 4

  # Convert to log molarity (log activity) then to logaH2
  logaO2 <- log10(bison.O2/1000/32)
  logK <- subcrt(c("O2", "H2", "H2O"), c(-0.5, -1, 1), T=bison.T)$out$logK
  logaH2.O <- 0 - 0.5*logaO2 - logK

  # 2011 plot
  xlab <- axis.label("T")
  ylab <- axis.label("H2")
  Tlim <- c(50, 100)
  plot(Tlim, get.logaH2(Tlim), xlim=Tlim, ylim=c(-45,0),
    xlab=xlab, ylab=ylab, type="l", lty=3, col = 4, lwd = 2)
  points(T.ORP, logaH2.ORP, pch=15)
  lines(T.ORP, logaH2.ORP, lty=2)
  points(bison.T, logaH2.O, pch=16)
  lines(bison.T, logaH2.O, lty=2)
  points(bison.T, logaH2.S, pch=17)
  lines(bison.T, logaH2.S, lty=2)
  llab <- c("Equation 2", "ORP", "dissolved oxygen", "sulfate/sulfide")
  text(c(65, 80, 80, 74), c(-4, -25, -40, -11), llab)

  # 2013 plot
  plot(Tlim, get.logaH2(Tlim), xlim=Tlim, ylim=c(-11,-2),
    xlab=xlab, ylab=ylab, type="l", lty=3, col = 4, lwd = 2)
  lines(bison.T, sapply(equil.results, "[", "logaH2.opt"), lty=2, col = 3, lwd = 2)
  points(bison.T, sapply(equil.results, "[", "logaH2.opt"), pch=21, bg="white", col = 3, lwd = 2)
  text(90, -5.3, "Equation 2")
  text(64, -9, "Optimized metastable\nequilibrium model", adj=0)
}

# Comparison of model and observed abundances
bison8 <- function(equil.results) {
  layout(matrix(c(1, 2, 3, 4, 5, 6), nrow = 2, byrow = TRUE), widths = c(2, 2, 2))
  par(mar = c(2.5, 0, 2.5, 0))
  plot.new()
  legend("topright", pch = 19, legend = phyla.abbrv, bty = "n", cex = 1.5, col = phyla.cols, text.col = 0)
  legend("topright", pch = 0:11, legend = phyla.abbrv, bty = "n", cex = 1.5)
  lim <- c(-6, -0.5)
  equil.opt <- a.blast <- alpha.blast()
  for(iloc in 1:5) {
    a.equil <- equil.results[[iloc]]
    iopt <- match(a.equil$logaH2.opt, a.equil$H2vals)
    ae.opt <- a.equil$alpha[iopt, ]
    these.cols <- phyla.cols[match(colnames(a.equil$alpha), phyla.abc)]
    # which are these phyla in the alphabetical list of phyla
    iphy <- match(names(ae.opt), phyla.abc)
    equil.opt[iloc, iphy] <- ae.opt
    mar <- c(2.5, 4.0, 2.5, 1)
    thermo.plot.new(xlab = expression(log[2]*alpha[obs]), ylab = expression(log[2]*alpha[equil]),
      xlim = lim, ylim = lim, mar = mar, cex = 1, yline = 1.5)
    # add points and 1:1 line
    points(log2(a.blast[iloc, iphy]), log2(ae.opt), pch = 19, col = these.cols)
    points(log2(a.blast[iloc, iphy]), log2(ae.opt), pch = iphy-1)
    lines(lim, lim, lty = 2)
    title(main = paste("site", iloc))
    # within-plot legend: DGtr
    DGexpr <- as.expression(quote(Delta*italic(G[tr])/italic(RT) == phantom()))
    DGval <- format(round(2.303*a.equil$DGtr[iopt], 3), nsmall = 3)
    legend("bottomright", bty = "n", legend = c(DGexpr, DGval))
  }
  invisible(equil.opt)
}

### UNEXPORTED OBJECTS ###

# names for sites 1-5
sites <- c("N", "S", "R", "Q", "P")
sitenames <- paste("bison", sites, sep = "")

# measured T and pH values
bison.T <- c(93.3, 79.4, 67.5, 65.3, 57.1)
bison.pH <- c(7.350, 7.678, 7.933, 7.995, 8.257)

# Dissolved oxygen measurements (mg/L)
bison.O2 <- c(0.173, 0.776, 0.9, 1.6, 2.8)

# distance and fitted T and pH values
distance <- c(0, 6, 11, 14, 22)
Tfun <- splinefun(distance, bison.T, method = "mono")
pHfun <- splinefun(distance, bison.pH, method = "mono")
O2fun <- splinefun(distance, bison.O2, method = "mono")
xpoints <- seq(0, 22, length.out = 128)

# read the amino acid compositions
aa.annot <- read.csv(system.file("extdata/bison/DS11.csv", package = "JMDplots"), as.is = TRUE)
aa.phyla <- read.csv(system.file("extdata/bison/DS13.csv", package = "JMDplots"), as.is = TRUE)
# functional annotations
classes <- unique(aa.annot$protein)
# the names of the phyla in alphabetical order (except Deinococcus-Thermus at end)
phyla.abc <- sort(unique(aa.phyla$organism))[c(1:7,9:11,8)]
# an abbreviation for Dein.-Thermus
phyla.abbrv <- phyla.abc
phyla.abbrv[[11]] <- "Dein.-Thermus"
# colors modified from Wu and Eisen, 2008 (doi:10.1186/gb-2008-9-10-r151)
phyla.cols <- c("#f48ba5", "#f2692f", "#cfdd2a",
  "#962272", "#87c540", "#66c3a2", "#12a64a", "#f58656",
  "#ee3237", "#25b7d5", "#3953a4")
phyla.lty <- c(1:6, 1:5)

# colors for sites (added 20200712)
# rev(hcl.colors(5, "Temps"))
site.cols <- c("#CF597E", "#E99F69", "#EAE29C", "#6CC382", "#089392")
# rev(hcl.colors(5, "Temps", 0.5))
site.cols.alpha <- c("#CF597E80", "#E99F6980", "#EAE29C80", "#6CC38280", "#08939280")

# function to load basis species
setup.basis <- function() {
  basis(c("HCO3-", "H2O", "NH3", "HS-", "H2", "H+"))
  basis(c("HCO3-", "NH3", "HS-", "H+"), c(-3, -4, -7, -7.933))
}

# function for model logaH2 (linear function of T)
get.logaH2 <- function(T) -11 + T * 3/40

# Function to return the fractional abundances based on BLAST counts
alpha.blast <- function() {
  out <- xtabs(ref ~ protein + organism, aa.phyla)
  # put it in correct order, then turn counts into fractions
  out <- out[c(1,5:2), c(1:7,9:11,8)]
  out <- out/rowSums(out)
  return(out)
}

# Function to calculate metastable equilibrium degrees of formation of proteins
# (normalized to residues) as a function of logaH2 for a specified location
alpha.equil <- function(i = 1, ip.phyla) {
  # order the names and counts to go with the alphabetical phylum list
  iloc <- which(aa.phyla$protein==sitenames[i])
  iloc <- iloc[order(match(aa.phyla$organism[iloc], phyla.abc))]
  # set up basis species, with pH specific for this location
  setup.basis()
  basis("pH", bison.pH[i])
  # calculate metastable equilibrium activities of the residues
  a <- affinity(H2 = c(-11, -1, 101), T = bison.T[i], iprotein = ip.phyla[iloc])
  e <- equilibrate(a, loga.balance = 0, as.residue = TRUE)
  # remove the logarithms to get relative abundances
  a.residue <- 10^sapply(e$loga.equil, c)
  colnames(a.residue) <- aa.phyla$organism[iloc]
  # the BLAST profile
  a.blast <- alpha.blast()
  # calculate Gibbs energy of transformation (DGtr) and find optimal logaH2
  iblast <- match(colnames(a.residue), colnames(a.blast))
  # Vectorize lists of loga1 and Astar (into matrices) for DGtr()
  loga1 <- sapply(e$loga.equil, c)
  loga2 <- log10(a.blast[i, iblast])
  Astar <- sapply(e$Astar, c)
  DGtr <- DGtr(loga1, loga2, Astar)
  logaH2.opt <- e$vals$H2[which.min(DGtr)]
  # return the calculated activities, logaH2 range, DGtr values, and optimal logaH2
  return(list(alpha = a.residue, H2vals = a$vals[[1]], DGtr = DGtr, logaH2.opt = logaH2.opt))
}

# Step-by-step example of protein metastable equilibrium calculation, after Dick and Shock, 2011
protein.equil <- function(protein, T = 25, loga.protein = 0, digits = 4) {
  out <- character()
  mymessage <- function(...) {
    message(...)
    text <- paste(list(...), collapse = " ")
    out <<- c(out, text)
  }
  # Show the individual steps in calculating metastable equilibrium among proteins
  mymessage("protein.equil: temperature from argument is ", T, " degrees C")
  # Display units
  E_units <- E.units()
  mymessage("protein.equil: energy units is ", E_units)
  TK <- convert(T, "K")
  # Get the amino acid compositions of the proteins
  aa <- pinfo(pinfo(protein))
  # Get some general information about the proteins
  pname <- paste(aa$protein, aa$organism, sep = "_")
  plength <- protein.length(aa)
  # Use thermo()$basis to decide whether to ionize the proteins
  thermo <- thermo()
  ionize.it <- FALSE
  iword <- "nonionized"
  bmat <- basis.elements()
  if("H+" %in% rownames(bmat)) {
    ionize.it <- TRUE
    iword <- "ionized"
    pH <- -thermo$basis$logact[match("H+", rownames(bmat))]
    mymessage("protein.equil: pH from thermo$basis is ", pH)
  }
  # Tell the user whose [Met] is in thermo$OBIGT
  info.Met <- info(info('[Met]', "aq"))
  mymessage("protein.equil: [Met] is from reference ", info.Met$ref1)
  ## First set of output: show results of calculations for a single protein
  mymessage("protein.equil [1]: first protein is ", pname[1], " with length ", plength[1])
  # Standard Gibbs energies of basis species
  G0basis <- unlist(suppressMessages(subcrt(thermo$basis$ispecies, T = T, property = "G")$out))
  # Coefficients of basis species in formation reactions of proteins
  protbasis <- suppressMessages(protein.basis(aa, T = T))
  # Sum of standard Gibbs energies of basis species in each reaction
  G0basissum <- colSums(t(protbasis) * G0basis)
  # Standard Gibbs energies of nonionized proteins
  G0prot <- unlist(suppressMessages(subcrt(pname, T = T, property = "G")$out))
  # Standard Gibbs energy of formation reaction of nonionized protein, E_units/mol
  G0protform <- G0prot - G0basissum
  mymessage("protein.equil [1]: reaction to form nonionized protein from basis species has G0(", E_units, "/mol) of ", signif(G0protform[1], digits))
  if(ionize.it) {
    # Standard Gibbs energy of ionization of protein, J/mol
    G0ionization <- suppressMessages(ionize.aa(aa, property = "G", T = T, pH = pH))[1, ]
    # Standard Gibbs energy of ionization of protein, E_units/mol
    if(E_units == "cal") G0ionization <- convert(G0ionization, "cal")
    mymessage("protein.equil [1]: ionization reaction of protein has G0(", E_units, "/mol) of ", signif(G0ionization[1], digits))
    # Standard Gibbs energy of formation reaction of ionized protein, E_units/mol
    G0protform <- G0protform + G0ionization
  }
  # Standard Gibbs energy of formation reaction of non/ionized residue equivalents, dimensionless
  # Gas constant
  #if(E_units == "cal") R <- 1.9872  # gas constant, cal K^-1 mol^-1
  #if(E_units == "J") R <- 8.314445  # = 1.9872 * 4.184 J K^-1 mol^-1  20220325
  if(E_units == "J") R <- 8.314463  # https://physics.nist.gov/cgi-bin/cuu/Value?r 20230630
  if(E_units == "cal") R <- 8.314463 / 4.184
  G0res.RT <- G0protform/R/TK/plength
  mymessage("protein.equil [1]: per residue, reaction to form ", iword, " protein from basis species has G0/RT of ", signif(G0res.RT[1], digits))
  # Coefficients of basis species in formation reactions of residues
  resbasis <- suppressMessages(protein.basis(aa, T = T, normalize = TRUE))
  # logQstar and Astar/RT
  logQstar <- colSums(t(resbasis) * - thermo$basis$logact)
  mymessage("protein.equil [1]: per residue, logQstar is ", signif(logQstar[1], digits))
  Astar.RT <- -G0res.RT - log(10)*logQstar
  mymessage("protein.equil [1]: per residue, Astar/RT = -G0/RT - 2.303logQstar is ", signif(Astar.RT[1], digits))
  if(!is.numeric(protein)) mymessage("protein.equil [1]: not comparing calculations with affinity() because 'protein' is not numeric")
  else {
    # For **Astar** we have to set the activities of the proteins to zero, not loga.protein!
    a <- suppressMessages(affinity(iprotein = protein, T = T, loga.protein = 0))
    aAstar.RT <- log(10) * as.numeric(a$values) / plength
    mymessage("check it!       per residue, Astar/RT calculated using affinity() is ", signif(aAstar.RT[1], digits))
    if(!isTRUE(all.equal(Astar.RT, aAstar.RT, check.attributes = FALSE)))
      stop("Bug alert! The same value for Astar/RT cannot be calculated manually as by using affinity()")
  }
  if(length(pname) == 1) mymessage("protein.equil [all]: all done... give me more than one protein for equilibrium calculations")
  else {
    ## Next set of output: equilibrium calculations
    mymessage("protein.equil [all]: lengths of all proteins are ", paste(plength, collapse = " "))
    mymessage("protein.equil [all]: Astar/RT of all residue equivalents are ", paste(signif(Astar.RT, digits), collapse = " "))
    expAstar.RT <- exp(Astar.RT)
    sumexpAstar.RT <- sum(expAstar.RT)
    mymessage("protein.equil [all]: sum of exp(Astar/RT) of all residue equivalents is ", signif(sumexpAstar.RT, digits))
    # Boltzmann distribution
    alpha <- expAstar.RT / sumexpAstar.RT    
    mymessage("protein.equil [all]: equilibrium degrees of formation (alphas) of residue equivalents are ", paste(signif(alpha, digits), collapse = " "))
    # Check with equilibrate()
    if(is.numeric(protein)) {
      loga.equil.protein <- unlist(suppressMessages(equilibrate(a, normalize = TRUE))$loga.equil)
      # Here we do have to convert from logarithms of activities of proteins to degrees of formation of residue equivalents
      a.equil.residue <- plength*10^loga.equil.protein
      ealpha <- a.equil.residue/sum(a.equil.residue)
      mymessage("check it!     alphas of residue equivalents from equilibrate() are ", paste(signif(ealpha, digits), collapse = " "))
      if(!isTRUE(all.equal(alpha, ealpha, check.attributes = FALSE)))
        stop("Bug alert! The same value for alpha cannot be calculated manually as by using equilibrate()")
    }
    # Total activity of residues
    loga.residue <- log10(sum(plength * 10^loga.protein))
    mymessage("protein.equil [all]: for activity of proteins equal to 10^", signif(loga.protein, digits), ", total activity of residues is 10^", signif(loga.residue, digits))
    # Equilibrium activities of residues
    loga.residue.equil <- log10(alpha*10^loga.residue)
    mymessage("protein.equil [all]: log10 equilibrium activities of residue equivalents are ", paste(signif(loga.residue.equil, digits), collapse = " "))
    # Equilibrium activities of proteins
    loga.protein.equil <- log10(10^loga.residue.equil/plength)
    mymessage("protein.equil [all]: log10 equilibrium activities of proteins are ", paste(signif(loga.protein.equil, digits), collapse = " "))
    # Check with equilibrate()
    if(is.numeric(protein)) {
      eloga.protein.equil <- unlist(suppressMessages(equilibrate(a, loga.balance = loga.residue, normalize = TRUE))$loga.equil)
      mymessage("check it!    log10 eq'm activities of proteins from equilibrate() are ", paste(signif(eloga.protein.equil, digits), collapse = " "))
      if(!isTRUE(all.equal(loga.protein.equil, eloga.protein.equil, check.attributes = FALSE)))
        stop("Bug alert! The same value for log10 equilibrium activities of proteins cannot be calculated manually as by using equilibrate()")
    }
  }
  return(out)
}

