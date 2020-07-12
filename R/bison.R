# JMDplots/bison.R
# Plots from hot spring (Bison Pool) papers (2011, 2013)
# Code moved from CHNOSZ/hotspring.Rnw 20200712

# This is used in the vignette to reproduce the 2011 calculations
#add.obigt("OldAA")

# Measured temperature and pH
bison1 <- function() {
  par(mfrow = c(1, 2), mar = c(3.5, 3.5, 1, 1), mgp = c(2.5, 1, 0), las = 1)
  # T plot
  plot(extendrange(distance), extendrange(bison.T), xlab = "Distance, m", ylab = axis.label("T"), type = "n")
  lines(xpoints, Tfun(xpoints))
  points(distance, bison.T, pch = 21, bg = "white", cex = 2)
  # col <- rev(hcl.colors(5, "Temps", 0.5))
  col <- c("#CF597E80", "#E99F6980", "#EAE29C80", "#6CC38280", "#08939280")
  points(distance, bison.T, pch = 21, bg = col, cex = 2)
  text(distance, bison.T, 1:5, cex = 0.8)

  # pH plot
  plot(extendrange(distance), extendrange(bison.pH), xlab = "Distance, m", ylab = "pH", type = "n")
  lines(xpoints, pHfun(xpoints))
  points(distance, bison.pH, pch = 21, bg = "white", cex = 2)
  points(distance, bison.pH, pch = 21, bg = col, cex = 2)
  text(distance, bison.pH, 1:5, cex = 0.8)
}

# Carbon oxidation state of proteins
bison2 <- function() {
  par(mfrow = c(1, 2), las = 1)

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
    if(classes[i]=="overall") text(0.8, ZC.annot[1+5*(i-1)], "ALL PROTEINS", adj = 1, font = 2)
    else if(classes[i] %in% clab) text(0.8, ZC.annot[1+5*(i-1)], classes[i], adj = 1)
  }
  title(main = "Annotations")

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
  title(main = "Major phyla")

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
  lines(Tlim, get.logaH2(Tlim), lty = 3)
  # second plot
  species(1:5, -3)
  xT <- Tfun(xpoints)
  xpH <- pHfun(xpoints)
  xH2 <- get.logaH2(xT)
  a <- affinity(T = xT, pH = xpH, H2 = xH2)
  a$vars[1] <- "Distance, m"
  a$vals[[1]] <- xpoints
  e <- equilibrate(a, normalize = TRUE)
  diagram(e, legend.x = NULL)
  legend("bottom", lty = 1:5, legend = 1:5, bty = "n", cex = 0.6)
}

# Comparing old and new group additivity parameters
bison5 <- function() {
  par(mfrow = c(2, 3))
  for(j in 1:2) {
    # use old parameters for first row and current ones for second row
    reset()
    if(j==1) add.OBIGT("OldAA")
    # setup basis species and proteins
    setup.basis()
    ip.annot <- add.protein(aa.annot)
    # make the plots
    for(annot in c("overall", "transferase", "synthase")) {
      ip <- ip.annot[aa.annot$protein==annot]
      a <- affinity(T = c(50, 100), H2 = c(-7, -4), iprotein = ip)
      diagram(a, fill = NULL, names = as.character(1:5), normalize = TRUE)
      # add logaH2-T line
      lines(par("usr")[1:2], get.logaH2(par("usr")[1:2]), lty=3)
      # add a title
      title(main=annot)
    }
  }
}

### UNEXPORTED OBJECTS ###

# names for sites 1-5
sites <- c("N", "S", "R", "Q", "P")
sitenames <- paste("bison", sites, sep = "")

# measured T and pH values
bison.T <- c(93.3, 79.4, 67.5, 65.3, 57.1)
bison.pH <- c(7.350, 7.678, 7.933, 7.995, 8.257)

# distance and fitted T and pH values
distance <- c(0, 6, 11, 14, 22)
Tfun <- splinefun(distance, bison.T, method = "mono")
pHfun <- splinefun(distance, bison.pH, method = "mono")
xpoints <- seq(0, 22, length.out = 128)

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

