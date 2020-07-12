# JMDplots/bison.R
# Plots from hot spring (Bison Pool) papers (2011, 2013)
# Code moved from CHNOSZ/hotspring.Rnw 20200712

# This is used in the vignette to reproduce the 2011 calculations
#add.obigt("OldAA")

# Measured temperature and pH
bison1 <- function() {
  distance <- c(0, 6, 11, 14, 22)
  par(mfrow = c(1, 2), mar = c(4, 4, 3, 2))
  xpoints <- seq(0, 22, length.out = 128)
  # T plot
  plot(distance, bison.T, xlab = "distance, m", ylab = axis.label("T"))
  Tfun <- splinefun(distance, bison.T, method = "mono")
  lines(xpoints, Tfun(xpoints))

  # pH plot
  plot(distance, bison.pH, xlab = "distance, m", ylab = "pH")
  pHfun <- splinefun(distance, bison.pH, method = "mono")
  lines(xpoints, pHfun(xpoints))
}

# Carbon oxidation state of proteins
bison2 <- function() {
  par(mfrow = c(1, 2))

  # 2011 plot
  ylab <- expression(bar(italic(Z))[C])
  plot(0, 0, xlim = c(-0.5, 5), ylim = c(-0.27, -0.11), xlab = "location", xaxt = "n", ylab = ylab)
  axis(1, at = 1:5)
  col <- c("green", rep("black", 20))
  lwd <- c(3, rep(1, 20))
  clab <- c("hydrolase", "overall", "protease", "oxidoreductase", "transport", "membrane", "permease")
  pf.annot <- protein.formula(aa.annot)
  ZC.annot <- ZC(pf.annot)
  for(i in 1:length(classes)) {
    lines(1:5, ZC.annot[(1:5)+5*(i-1)], col = col[i], lwd = lwd[i])
    if(classes[i] %in% clab) text(0.8, ZC.annot[1+5*(i-1)], classes[i], adj = 1)
  }
  title(main = "Annotations")

  # 2013 plot
  pf.phyla <- protein.formula(aa.phyla)
  ZC.phyla <- ZC(pf.phyla)
  # set up plot
  plot(0, 0, xlim = c(1, 5), ylim = c(-0.27, -0.11), xlab = "location", ylab = ylab)
  for(i in 1:length(phyla.abc)) {
    # which of the model proteins correspond to this phylum
    iphy <- which(aa.phyla$organism==phyla.abc[i])
    # the locations (of 1, 2, 3, 4, 5) where this phylum is found
    ilocs <- match(aa.phyla$protein[iphy], sitenames)
    # the plotting symbol: determined by alphabetical position of the phylum
    points(ilocs, ZC.phyla[iphy], pch = i-1, cex = 1.2)
    # a line to connect same phyla occurring at adjacent sites
    inlocs <- rep(NA, 5)
    inlocs[ilocs] <- ilocs
    lines(inlocs, ZC.phyla[iphy][match(1:5, ilocs)])
  }
  legend("bottomright", pch = 0:10, legend = phyla.abbrv, bg = "white", cex = 0.9)
  title(main = "Major phyla")
}

# Chemical affinities
bison3 <- function() {
  setup.basis()
  ip.annot <- add.protein(aa.annot)
  pl <- protein.length(ip.annot[1:5])
  species("overall", sitenames)
  species(1:5, 0)
  a <- affinity(T = bison.T, pH = bison.pH, H2 = get.logaH2(bison.T))
  a.res <- t(as.data.frame(a$values))/pl
  colnames(a.res) <- paste0("site", 1:5)
  rownames(a.res) <- paste0("reaction", 1:5)
  a.res
}

### UNEXPORTED OBJECTS ###

# names for sites 1-5
sites <- c("N", "S", "R", "Q", "P")
sitenames <- paste("bison", sites, sep = "")

# measured T and pH values
bison.T <- c(93.3, 79.4, 67.5, 65.3, 57.1)
bison.pH <- c(7.350, 7.678, 7.933, 7.995, 8.257)

# read the amino acid compositions
aa.annot <- read.csv(system.file("extdata/protein/DS11.csv", package = "CHNOSZ"), as.is = TRUE)
aa.phyla <- read.csv(system.file("extdata/protein/DS13.csv", package = "CHNOSZ"), as.is = TRUE)
# functional annotations
classes <- unique(aa.annot$protein)
# the names of the phyla in alphabetical order (except Deinococcus-Thermus at end)
phyla.abc <- sort(unique(aa.phyla$organism))[c(1:7,9:11,8)]
# an abbreviation for Dein.-Thermus
phyla.abbrv <- phyla.abc
phyla.abbrv[[11]] <- "Dein.-Thermus"
# colors modified from Wu and Eisen, 2008
phyla.cols <- c("#f48ba5", "#f2692f", "#cfdd2a",
  "#962272", "#87c540", "#66c3a2", "#12a64a", "#f58656",
  "#ee3237", "#25b7d5", "#3953a4")
phyla.lty <- c(1:6, 1:5)

# function to load basis species
setup.basis <- function() {
  basis(c("HCO3-", "H2O", "NH3", "HS-", "H2", "H+"))
  basis(c("HCO3-", "NH3", "HS-", "H+"), c(-3, -4, -7, -7.933))
}

# function for model logaH2 (linear function of T)
get.logaH2 <- function(T) -11 + T * 3/40

