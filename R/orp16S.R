# JMDplots/orp16S.R
# Plots for paper on Zc-ORP correlations 20210827

# Group studies by environment types 20210828
envirotype <- list(
  "River & Seawater" = c("MLL+18", "HXZ+20", "GSBT20_Prefilter", "GSBT20_Postfilter", "WHL+21", "ZLH+22", "ZZL+21", "LWJ+21", "GZL21", "RARG22"),
  "Lake & Pond" = c("SAR+13", "LLC+19", "BCA+21", "HLZ+18", "BWD+19", "IBK+22", "NLE+21", "MTC21", "SPA+21"),
  "Groundwater" = c("KLM+16", "WLJ+16", "ZDW+19", "DJK+18", "SRM+19", "APV+20", "YHK+20", "ZCZ+21", "MGW+22", "MCR+22"),
  # NOTE: Keep Geothermal in 4th location to get red color 20210904
  "Geothermal" = c("PCL+18_Acidic", "PCL+18_Alkaline", "GWS+20", "PBU+20", "MWY+21"),
  "Hyperalkaline" = c("SBP+20", "RMB+17", "CTS+17", "KSR+21", "PSB+21", "NTB+21"),
  "Sediment" = c("ZML+17", "BSPD17", "RKN+17", "HDZ+19", "OHL+18_DNA", "WHLH21a", "RSS+18", "CLS+19", "HSF+19", "ZHZ+19",
                 "LMBA21_2017", "HSF+22", "ZZLL21", "WFB+21", "HCW+22", "WKG+22"),
  "Soil" = c("MLL+19", "BMOB18", "WHLH21", "CWC+20", "PSG+20", "LJC+20", "DTJ+20", "RKSK22", "DLS21_Bulk", "WKP+22",
             "CKB+22", "CLZ+22")
)
# Turn the list into a data frame for easier lookup 20210904
envirodat <- do.call(rbind, lapply(seq_along(envirotype), function(i) data.frame(study = envirotype[[i]], groupnum = i)))
envirodat <- cbind(envirodat, group = names(envirotype)[envirodat$groupnum])

# Select color palette 20210914
orp16Scol <- palette.colors(n = length(envirotype), palette = "Classic Tableau", alpha = 0.75)

# Read Eh7 - Zc data
EZdat <- read.csv(system.file("extdata/orp16S/EZdat.csv", package = "JMDplots"))
# Read linear fit coefficients
EZlm <- read.csv(system.file("extdata/orp16S/EZlm.csv", package = "JMDplots"))

##################################
### Functions for main figures ###
##################################

# Figure 1: Thermodynamic model for the relationship between carbon oxidation state of reference proteomes and redox potential 20221213
orp16S_1 <- function(pdf = FALSE) {
  # Setup plot
  if(pdf) pdf("Figure_1.pdf", width = 7, height = 6, bg = "white")
  mat <- matrix(c(1,1,1,2,2,2, 1,1,1,3,3,3, 4,4,5,5,6,6), nrow = 3, byrow = TRUE)
  layout(mat)

  # Species with optimal growth temperatures from 30-40 °C
  species <- c(
    # Class I
    "Methanococcus vannielii", "Methanococcus maripaludis", "Methanococcus voltae", "Methanobrevibacter smithii",
    # Class II
    #"Methanosarcina barkeri", "Methanoplanus limicola", "Methanofollis liminatans"
    "Methanofollis liminatans"
  )
  # Get reference proteomes
  refseq <- read.csv(system.file("extdata/RefDB/RefSeq/genome_AA.csv.xz", package = "JMDplots"))
  irefseq <- match(species, refseq$ref)
  aa <- refseq[irefseq, ]
  ip <- add.protein(aa, as.residue = TRUE)

  # Plot theoretical (equilibrium) Eh-pH diagram
  basis(c("glutamine", "glutamic acid", "cysteine", "H2O", "H+", "e-"), c(-3.2, -4.5, -3.6, 0, 0, 0))
  res <- 385
  a <- affinity(pH = c(3, 11, res), Eh = c(-0.5, 0.1, res), iprotein = ip)
  # Use abbreviations for species
  names <- paste0("M", substr(sapply(strsplit(species, " "), "[", 2), 1, 2))
  srt <- c(0, -40, -40, -40, 0)
  d <- diagram(a, names = names, ylab = "Eh (V)", srt = srt, mar = c(3.1, 3.1, 3, 1))
  abline(v = 5, lty = 2, col = 2)
  abline(v = 7, lty = 2, col = 1)
  abline(v = 9, lty = 2, col = 4)
  title("Metastable equilibrium (Met. Equil.)\namong selected reference proteomes", font.main = 1)
  label.figure("a", cex = 1.7, font = 2)

  # Plot theoretical Zc as a function of Eh and Eh7
  par(mar = c(3.1, 4.1, 2, 1))
  par(yaxs = "r")
  # Calculate Zc
  Zc <- Zc(aa)
  # Get Zc values at discrete pHs
  pHs <- c(5, 7, 9)
  adjs <- c(0.6, 0.8, 1.1)
  Zcvals <- lapply(pHs, function(pH) {
    ipH <- d$vals$pH == pH
    ispecies_pH <- d$predominant[ipH, ]
    Zc[ispecies_pH]
  })
  names(Zcvals) <- pHs
  col <- c(2, 1, 4)
  # Loop over Eh and Eh7
  for(xvar in c("Eh", "Eh7")) {
    plot(c(-0.5, 0.1), extendrange(unlist(Zcvals)), xlab = paste(xvar, "(V)"), ylab = "", type = "n")
    mtext(axis.label("Zc"), side = 2, line = 2.5, cex = par("cex"), las = 0)
    # Apply offset to distinguish same Zc values at different pH
    dy <- c(-0.0005, 0, 0.0005)
    for(i in 1:length(pHs)) {
      xvals <- d$vals$Eh
      # Calculate Eh7
      if(xvar == "Eh7") xvals <- xvals + -0.05916 * (7 - pHs[i])
      # Dashed lines between Zc
      lines(xvals, Zcvals[[i]] + dy[i], col = col[i], lty = 2)
      # Solid lines at single Zc
      uZc <- unique(Zcvals[[i]])
      lapply(unique(Zcvals[[i]]), function(uZc) {
        iuZc <- Zcvals[[i]] == uZc
        lines(xvals[iuZc], Zcvals[[i]][iuZc] + dy[i], col = col[i], lwd = 2)
      })
      if(xvar == "Eh") {
        # Plot pH labels
        iEh <- which(Zcvals[[i]] == Zc[3])[1]
        text(xvals[iEh], -0.19, paste("pH =", pHs[i]), adj = adjs[i])
      }
    }
    if(xvar == "Eh") title("Met. Equil. at three pH values", font.main = 1)
    if(xvar == "Eh7") title("Met. Equil. with Eh7 correction", font.main = 1)
    if(xvar == "Eh") label.figure("b", cex = 1.7, font = 2)
  }

  # Plot theoretical Zc as a function of Eh7 with simulated noise added 20221220
  par(xaxs = "r")
  # Use second pH value defined above (i.e. pH = 7)
  i <- 2
  pH <- pHs[i]
  # Calculate Eh7 (xvals)
  xvals <- round(d$vals$Eh, 8)
  xvals <- xvals + -0.05916 * (7 - pH)
  # Use coarser resolution for scatterplot
  ipoint <- seq(1, 385, length.out = 97)
  # Use values in the range [-0.3, -0.1]
  iEh7 <- which(xvals >= -0.26 & xvals <= -0.14)
  iuse <- intersect(ipoint, iEh7)
  xvals <- xvals[iuse]
  # Get Zc values
  Zc <- Zcvals[[i]][iuse]

  # Plot 1: Theoretical Zc vs Eh7
  plot(xvals, Zc, xlab = "Eh7 (V)", ylab = "", pch = 19, cex = 0.5, col = "#00000080")
  mtext(axis.label("Zc"), side = 2, line = 2.5, cex = par("cex"), las = 0)
  add.linear.local(xvals, Zc, legend = "topleft", cex = 1, inset = c(-0.08, -0.05), with.N = FALSE, scale = 1)
  title("Met. Equil. sampled at pH = 7", font.main = 1)
  label.figure("c", cex = 1.7, font = 2)

  # Use RefSeq amino acid compositions and taxon names
  species_aa <- refseq
  taxa <- read.csv(system.file("extdata/RefDB/RefSeq/taxonomy.csv.xz", package = "JMDplots"), as.is = TRUE)
  # Make sure the data tables have consistent taxids
  stopifnot(all(species_aa$organism == taxa$taxid))
  # Keep only Bacteria and Archaea classified at species level 20220104
  isspecies <- !is.na(taxa$species)
  ivirus <- taxa$superkingdom == "Viruses"
  ivirus[is.na(ivirus)] <- FALSE
  species_aa <- species_aa[isspecies & !ivirus, ]
  taxa <- taxa[isspecies & !ivirus, ]
  # Keep only species with at least 500 sequences
  ilow <- species_aa$chains < 500
  species_aa <- species_aa[!ilow, ]
  taxa <- taxa[!ilow, ]

  # Plot 2: Random Zc vs Eh7
  set.seed(42)
  # Randomly sample rows from RefSeq data frame
  ispecies <- sample(1:nrow(species_aa), length(Zc))
  randvals <- Zc(species_aa[ispecies, ])
  ylim <- range(c(randvals, -0.08))
  plot(xvals, randvals, xlab = "Eh7 (V)", ylab = "", ylim = ylim, pch = 19, cex = 0.5, col = "#00000080")
  mtext(axis.label("Zc"), side = 2, line = 2.5, cex = par("cex"), las = 0)
  add.linear.local(xvals, randvals, legend = "topleft", cex = 1, inset = c(-0.08, -0.05), with.N = FALSE, scale = 1)
  title("Random species from RefSeq", font.main = 1)

  # Plot 3: Theoretical + random Zc vs Eh7
  Zc_shaped <- (Zc + 4*randvals) / 5
  ylim <- range(c(Zc_shaped, -0.08))
  plot(xvals, Zc_shaped, xlab = "Eh7 (V)", ylab = "", ylim = ylim, pch = 19, cex = 0.5, col = "#00000080")
  mtext(axis.label("Zc"), side = 2, line = 2.5, cex = par("cex"), las = 0)
  add.linear.local(xvals, Zc_shaped, legend = "topleft", cex = 1, inset = c(-0.08, -0.05), with.N = FALSE, scale = 1)
  title("20% Met. Equil. + 80% Random  ", font.main = 1)

  if(pdf) dev.off()

}


# Zc of reference proteomes compared with oxygen tolerance and with metaproteomes 20221228
orp16S_2 <- function(pdf = FALSE) {

  # Setup plot
  if(pdf) pdf("Figure_2.pdf", width = 8, height = 3.5, bg = "white")
  mat <- matrix(1:3, nrow = 1)
  layout(mat, widths = c(0.9, 1.4, 0.9))
  par(mgp = c(2.5, 1, 0))

  # Panel B: Zc in obligate anaerobic and aerotolerant genera 20221017
  par(mar = c(4, 4, 2.5, 1))
  # Read table from Million and Raoult, 2018
  dat <- read.csv(system.file("extdata/orp16S/MR18_Table_S1.csv", package = "JMDplots"))
  # Clean up names
  dat$Genus.name <- gsub(" ", "", dat$Genus.name)
  # Get amino acid compositions for genera
  aa <- taxon_AA$RefSeq
  # Calculate Zc
  values <- Zc(aa)
  ylim <- c(-0.24, -0.095)
  # Match genus names to RefSeq
  iref <- match(dat$Genus.name, aa$organism)
  # Print coverage information
  nna <- sum(is.na(iref))
  print(paste(nna, "genera not matched to RefSeq"))
  # Report archaeal genera 20230107
  genera <- aa[iref, ]$organism
  taxa <- read.csv(system.file("extdata/RefDB/RefSeq/taxonomy.csv.xz", package = "JMDplots"), as.is = TRUE)
  itax <- match(genera, taxa$genus)
  taxa <- taxa[itax, ]
  archaeal_genera <- taxa$genus[taxa$superkingdom == "Archaea"]
  print(paste("Archaeal genera:", paste(archaeal_genera, collapse = " ")))
  # Get values for obligate anaerobes and aerotolerant genera
  values <- values[iref]
  values <- list(
    Anaerobe = values[dat$Obligate.anerobic.prokaryote == 1],
    Aerotolerant = values[dat$Obligate.anerobic.prokaryote == 0]
  )
  boxplot(values, ylab = quote(italic(Z)[C]~"of genus reference proteome"), ylim = ylim, names = NA)
  axis(1, 1, labels = "Strictly\nAnaerobic", tick = FALSE, line = 0.7)
  axis(1, 2, labels = "Aerotolerant", tick = FALSE, line = 0.3)
  # Add p-values
  ptext <- bquote(italic(P) == .(signif(t.test(values$Anaerobe, values$Aerotolerant)$p.value, 3)))
  legend("bottomright", legend = ptext, bty = "n")
  # Show number of genera below names
  n <- c(
    sum(dat$Obligate.anerobic.prokaryote == 1 & !is.na(iref)),
    sum(dat$Obligate.anerobic.prokaryote == 0 & !is.na(iref))
  )
  l1 <- bquote(italic(N)==.(n[1]))
  l2 <- bquote(italic(N)==.(n[2]))
  axis(1, 1, labels = l1, tick = FALSE, line = 1.8)
  axis(1, 2, labels = l2, tick = FALSE, line = 1.8)
  # Add mean difference to title
  md <- round(diff(sapply(values, mean, na.rm = TRUE)), 4)
  main <- paste("Mean difference =", md)
  title(main, font.main = 1, cex.main = 1, line = -0.9)
  title("Oxygen tolerance", font.main = 1)
  label.figure("a", font = 2, cex = 1.5, xfrac = 0.04)

  # Panel A: Zc of community reference proteomes vs metaproteomes 20221222
  par(mar = c(4, 4, 2.5, 1))
  # Calculate chemical metrics
  dat <- getMP_orp16S()
  Zclim <- c(-0.20, -0.10)
  plot(Zclim, Zclim, xlab = quote(italic(Z)[C]~"of metaproteome"), ylab = quote(italic(Z)[C]~"of community reference proteome"), type = "n")
  lines(Zclim, Zclim, lty = 2, col = 8)
  points(dat$Zc_MP, dat$Zc_16S, pch = dat$pch, bg = dat$bg, col = dat$col)
  add.linear.local(dat$Zc_MP, dat$Zc_16S, legend = "topleft", with.N = FALSE, scale = 1, nounits = TRUE)
  title(hyphen.in.pdf("Metaproteome - Reference proteome comparison"), font.main = 1)
  label.figure("b", font = 2, cex = 1.5, xfrac = 0.03)
  # Make legend for datasets
  idup <- duplicated(dat$Name)
  # Get dataset names and number of samples
  names <- dat$Name[!idup]
  nsamp <- sapply(lapply(names, `==`, dat$Name), sum)
  names <- paste0(names, " (", nsamp, ")")
  dat <- dat[!idup, ]
  opar <- par(mar = c(4, 0, 4, 0))
  plot.new()
  legend("left", sapply(names, hyphen.in.pdf), pch = dat$pch, col = dat$col, pt.bg = dat$bg, bty = "n")
  par(opar)

  if(pdf) dev.off()

}


# Figure 3: Methods overview and Winogradsky columns 20221110
orp16S_3 <- function(pdf = FALSE) {

  if(pdf) pdf("Figure_3.pdf", width = 6/(8/15), height = 4, bg = "white")

  # Figure 1a: Schematic of data and methods 20210830
  layout(matrix(1:3, nrow = 1), widths = c(8, 3.5, 3.5))
  par(cex = 1)

  # Define colors; adapted from
  # https://cdn.elifesciences.org/author-guide/tables-colour.pdf
  Purple <- "#9E86C9"
  Purple80 <- "#9E86C980"
  PurpleText <- "#9C27B0"
  Blue80 <- "#90CAF980"
  BlueText <- "#366BFB"
  Red80 <- "#E5737380"
  RedText <- "#D50000"
  Orange80 <- "#FFB74D80"
  OrangeText <- "#f7831c"
  Gray <- "#E6E6E6"

  # Setup plot
  par(mar = c(0, 0, 0, 0))
  plot(c(5, 95), c(1, 80), type = "n", axes = FALSE, xlab = "", ylab = "")
  # Uncomment this as a guide for making the 'grid' one 20210927
  #box(lwd = 2)
  grid.roundrect(0.5*(8/15), 0.5, 0.49, 0.99, gp = gpar(fill = "azure"))

  dy <- -2.5

  # Plot lines to go behind shapes
  lines(c(20, 50), c(20, 60)+dy)
  lines(c(20, 50), c(60, 60)+dy)
  lines(c(50, 80), c(60, 60)+dy)
  # Add arrows along lines 20210927
  arrows(20, 20+dy, 20*0.35 + 50*0.65, (20*0.35 + 60*0.65)+dy, length = 0.1)
  arrows(20, 60+dy, 37, 60+dy, length = 0.1)
  arrows(50, 60+dy, 66, 60+dy, length = 0.1)

  # Plot shapes and text for biological methods
  text(20, 79+dy, "Biological Methods", col = RedText, font = 2)
  for(bg in c("white", Red80)) points(20, 60+dy, pch = 21, cex = 17, bg = bg)
  text(20, 66+dy, "RefSeq", font = 2, col = RedText)
  text(20, 57+dy, "Reference\nproteomes\nof taxa")
  for(bg in c("white", Red80)) points(20, 20+dy, pch = 21, cex = 17, bg = bg)
  text(20, 25+dy, "16S + RDP", font = 2, col = RedText)
  text(20, 18+dy, "Taxonomic\nabundances")

  # Plot shapes and text for chemical methods
  text(80, 79+dy, "Chemical Methods", col = BlueText, font = 2)
  for(bg in c("white", Blue80)) points(80, 60+dy, pch = 22, cex = 16, bg = bg)
  text(80, 67+dy, quote(bolditalic(Z)[bold(C)]), col = BlueText, cex = 1.2)
  text(80, 57.5+dy, "Carbon\noxidation\nstate")
  # Show multiple physicochemical variables 20210927
  # Function to draw rectangle at x,y with width and height w,h
  myrect <- function(x, y, w, h, ...) rect(x - w/2, y - h/2, x + w/2, y + h/2, ...)
  # T, Eh, pH, O2
  for(col in c("white", Blue80)) myrect(73, 16+dy, 5, 6, col = col)
  text(73, 17.3+dy, "pH", col = BlueText, adj = c(0.5, 1))
  for(col in c("white", Blue80)) myrect(80, 16+dy, 7, 7, col = col)
  text(80, 17.3+dy, "Eh", col = BlueText, adj = c(0.5, 1))
  for(col in c("white", Blue80)) myrect(87.5, 16+dy, 6, 6, col = col)
  text(87.5, 17.3+dy, "T", font = 3, col = BlueText, adj = c(0.5, 1))
  for(col in c("white", Blue80)) myrect(80, 8+dy, 5, 6, col = col)
  text(80, 8+dy, quote(O[2]), col = BlueText)
  # Add Eh7 and lines to source data 20221008
  for(col in c("white", Blue80)) myrect(80, 25.5+dy, 7, 7, col = col)
  text(80, 25.5+dy, "Eh7", font = 2, col = BlueText)
  lines(c(73, 87.5), c(19.5, 19.5))
  lines(c(73, 73), c(16.5, 19.5))
  lines(c(80, 80), c(17, 19.5))
  lines(c(87.5, 87.5), c(16.5, 19.5))

  # Plot inner rectangle and text
  third <- 100/3
  # Uncomment this to make the original rectangle (as a guide for the 'grid' one)
  #for(bg in c("white", Orange80)) rect(1.2*third, 50+dy, 1.8*third, 70+dy, col = bg)
  # Use this to get rounded corners 20210927
  for(fill in c("white", Orange80)) grid.roundrect(0.5*(8/15), 0.7, 0.11, 0.23, gp = gpar(fill = fill))
  text(50, 67+dy, quote(bold(C[bolditalic(c)]*H[bolditalic(h)]*N[bolditalic(n)]*O[bolditalic(o)]*S[bolditalic(s)])), col = OrangeText)
  text(50, 58+dy, "Community\nreference\nproteomes")

  # Plot arrows and text labels
  arrows(1*third - 1, 20+dy, 1*third + 6, 20+dy, code = 1, lty = 1, length = 0.1)
  arrows(2*third - 4, 16+dy, 2*third + 2.5, 16+dy, code = 2, lty = 1, length = 0.1)
  text(50, 18+dy, "Community and\nenvironmental data")
  arrows(80, 31+dy, 80, 46.5+dy, code = 3, lwd = 1.5, length = 0.1, col = BlueText)
  text(78, 46+dy, "Thermodynamic prediction", font = 2, cex = 0.9, adj = c(1, 1))
  text(78, 46+dy, "\nPositive correlation between\ncarbon oxidation state\nand redox potential", font = 3, cex = 0.9, adj = c(1, 1))
  label.figure("a", font = 2, cex = 1.5, xfrac = 0.015, yfrac = 0.97)

  ## Figure 1b: Chemical depth profiles in Winogradsky columns 20210829
  par(mar = c(4, 4, 1, 1), mgp = c(2.5, 1, 0))

  # ORP measurements from Diez-Ercilla et al. (2019)
  depth <- c(25, 22, 18, 14, 12, 9, 5)
  ORP <- c(-229, -241, -229, -201, 302, 501, 641) / 1000 # Units: V
  pH <- c(5.5, 5.5, 5.5, 4.7, 3.0, 2.7, 2.6)
  SWIdepth <- 11.3
  # Shift depth values so SWI is at zero
  MODdepth <- depth - SWIdepth
  plot(ORP, MODdepth, ylim = c(14.7, -6.3), type = "b", lty = 2, xlab = "Eh or Eh7 (V)", ylab = "Depth (cm)", xlim = c(-0.4, 0.8), las = 1)
  abline(h = 0, lty = 2, col = "darkgray", lwd = 1.5)
  # Calculate and plot Eh7 20210926
  Eh7 <- ORP + -0.05916 * (7 - pH)
  lines(Eh7, MODdepth, type = "b", pch = 19)
  text(0.08, -3.5, "Eh7")
  text(0.72, -3.5, "Eh")
  label.figure("b", font = 2, cex = 1.5, yfrac = 0.97)

  # Get 16S metadata and chemical metrics for Rundell et al. (2014) experiments
  metrics.in <- getmetrics_orp16S("RBW+14")
  mdat <- getmdat_orp16S("RBW+14", metrics.in)
  metrics <- mdat$metrics
  # Get Zc values for each layer
  layers <- c("12 cm", "8 cm", "4 cm", "SWI", "Top")
  Zc <- lapply(layers, function(layer) metrics$Zc[mdat$metadata$layer == layer])
  # Make boxplots
  boxplot(Zc, horizontal = TRUE, show.names = FALSE, xlab = axis.label("Zc"), ylim = c(-0.18, -0.145), yaxs = "i", col = "azure3")
  axis(2, 1:5, labels = layers, las = 1)
  # Add sample sizes
  par(xpd = NA)
  for(i in 1:5) {
    label <- bquote(italic(N) == .(length(Zc[[i]])))
    text(-0.1835, i - 0.3, label, adj = 1, cex = 0.9)
  }
  par(xpd = FALSE)
  label.figure("c", font = 2, cex = 1.5, yfrac = 0.97)

  # Reset layout to make orp16S_2 in the examples run nicely 20211011
  if(pdf) dev.off() else layout(1)

}


# Figure 4: Sample locations and Eh-pH diagram 20221110
orp16S_4 <- function(pdf = FALSE) {

  ## Panel A: Sample locations on world map

  if(pdf) pdf("Figure_4.pdf", width = 26, height = 15+9-1, bg = "white")
  mat <- matrix(c(1,1,1,1, 1,2,3,1, 0,2,3,0), nrow = 3, byrow = TRUE)
  layout(mat, heights = c(14, 1, 8), widths = c(7, 9, 3, 7))
  par(cex = 1)

  # Coordinates for orp16S datasets
  file <- tempfile()
  # Write spaces here (but don't save them in the file) to make this easier to read
  writeLines(con = file, text = gsub(" ", "", c(
    # Column names
    "study, latitude, longitude",
    ## River & Seawater (comments indicate source of coordinates from paper or SRA metadata)
    "MLL+18, 22.20, 113.09", # Fig. 1
    "HXZ+20, 16.52, 111.77", # Table 1
    # GSBT20 see below
    "WHL+21, 30.76, 115.36", # SAMN13242327
    "ZLH+22, 24.179, 98.902", # SAMN16122993
    "ZZL+21, 22.77, 113.79", # SAMN16964962
    "LWJ+21, 36.569, 120.976", # SAMN18558469
    # GZL21 see below
    "RARG22, 18.5, -99.5", # Materials and methods 
    ## Lake & Pond
    "SAR+13, 1.96, -157.33", # Materials and methods
    "LLC+19, 24.82, 118.15", # SAMN04549101
    "BCA+21, 46.3615, 25.0509", # SAMN07409474
    "HLZ+18, 24.795, 118.138", # SAMN07638080
    "BWD+19, 47.120571, -88.545425", # SAMN09980099
    "IBK+22, 53.1516, 13.0262", # SAMN15366194
    "NLE+21, 32.833, 35.583", # SAMEA7280991
    "MTC21, 41.396, -0.058", # SAMN10716683
    "SPA+21, 45.8126, 8.7401", # SAMN17524543
    ## Geothermal (Hot Spring)
    "PCL+18_Acidic, -38.5, 176.0", # Fig. 1
    "GWS+20, 30.12, 101.94", # SAMN13430433
    "PBU+20, 54.4395, 160.144194", # SAMN14538724
    # MWY+21 see below
    ## Hyperalkaline
    "SBP+20, 38.862, -122.414", # SAMN03850954
    "RMB+17, 22.9052, 58.6606", # SAMN05981641
    "CTS+17, 10.94323, -85.63485", # SAMN06226041
    "KSR+21, 44.264340, 8.46442", # SAMN17101425
    "PSB+21, 38.8621, -122.4304", # SAMN17252996
    "NTB+21, 22.881, 58.701", # SAMN19998441
    ## Soil - put this group before Groundwater and Sediment for clearer visualization in GBA 20210927
    "MLL+19, 26.1, 112.5", # Materials and methods
    "BMOB18, 40.60842, -74.19258", # SAMN07828017  ### Laboratory
    "WHLH21, 37.53, 105.03", # Materials and methods
    "CWC+20, 28.226, 116.898", # Materials and methods  ### Laboratory
    "PSG+20, 36.61, -119.53", # Web search for Parlier, CA   ### Mesocosm
    # LJC+20 see below
    "DTJ+20, 26.45, 111.52", # SAMN14332759   ### Laboratory
    "RKSK22, 31.97283, -81.0304", # SAMN16678415
    "DLS21_Bulk, 39.39, -75.44", # SAMN17245435  ### Mesocosm
    "WKP+22, 53.29, 17.79", # SAMN23457780
    "CKB+22, 37.8635, 138.9426", # Wikipedia Niigata University   ### Laboratory
    "CLZ+22, 25.1, 113.63", # Materials and Methods   ### Laboratory
    ## Groundwater
    "KLM+16, 42.99, -82.30", # SAMN04423023
    "WLJ+16, 41, 107", # SAMN05938707
    "ZDW+19, 30.18, 113.61", # SAMN07528610
    "DJK+18, 39.44, -82.22", # Materials and methods
    "SRM+19, 12.67417, 101.3889", # Materials and methods
    "APV+20, 20.12, -99.23", # Materials and methods
    "YHK+20, 51.209467, 10.791968", # SAMEA5714424
    "ZCZ+21, 45.21, 9.57", # Table 1
    "MGW+22, -36.849, 174.769", # Wikipedia University of Auckland
    "MCR+22, -38.386099, 144.865068", # SAMN29926991
    ## Sediment
    "ZML+17, 22.494, 114.029", # Table 1
    "BSPD17, 57.89297, 16.5855", # SAMN05163191   ### Laboratory
    "RKN+17, 62.03, 29.98", # SAMN05933570   ### Laboratory
    "HDZ+19, 29.901, 113.52435", # SAMN05990289
    "OHL+18_DNA, 56.02, 10.26", # Experimental procedures
    "WHLH21a, 23.52, 113.495", # Materials and methods
    "RSS+18, 82, -71", # Introduction
    "CLS+19, 32.22, 118.83", # SAMN08683376
    "HSF+19, 47.803, 16.709", # methods
    "ZHZ+19, 23.130, 113.671", # Materials and methods   ### Laboratory
    "LMBA21_2017, 43.42, -2.7", # Fig. 1
    "HSF+22, -9.42979, 46.49524", # SAMN14343437
    "ZZLL21, 22.68, 113.97", # SAMN16964887
    "WFB+21, 56.440893, -2.863194", # methods   ### Mesocosm
    "HCW+22, 31.3, 119.98", # Materials and methods   ### Laboratory
    "WKG+22, 53.164, 17.704" # Materials and methods
  )))

  # This reads the data
  coords <- read.csv(file, as.is = TRUE)
  # Make sure the study keys are correct 20220522
  stopifnot(all(coords$study %in% unlist(envirotype)))

  # Start map
  # https://cran.r-project.org/web/packages/oce/vignettes/map_projections.html
  par(mar = c(2, 0.5, 0, 0.5))
  # We don't need data(coastlineWorld) ... it's the default map 20211003
  # A color between azure2 and azure3 20220516
  azure23 <- rgb(t((col2rgb("azure2") + col2rgb("azure3")) / 2), maxColorValue = 255)
  mapPlot(col = azure23, projection = "+proj=wintri", border = "white", drawBox = FALSE)

  # Add Great Lakes
  # https://www.sciencebase.gov/catalog/item/530f8a0ee4b0e7e46bd300dd
  # https://gis.stackexchange.com/questions/296170/r-shapefile-transform-longitude-latitude-coordinates
  hydro_pdir <- system.file("extdata/orp16S/hydro_p", package = "JMDplots")
  for(lake in c("Erie", "Huron", "Michigan", "Ontario", "Superior", "StClair")) {
    file <- file.path(hydro_pdir, paste0("hydro_p_Lake", lake, ".shp"))
    tx <- st_read(file, quiet = TRUE)
    # https://mhweber.github.io/AWRA_2020_R_Spatial/coordinate-reference-systems.html

    ## library(sf) (version 1.0.8) gave:
    ## Linking to GEOS 3.10.2, GDAL 3.4.2, PROJ 9.0.0; sf_use_s2() is TRUE
    ## This was working with those versions
    #tx_ll <- st_transform(tx, "+proj=longlat +ellps=GRS80 +datum=NAD83")

    # library(sf) (version 1.0.9) gave:
    # Linking to GEOS 3.11.1, GDAL 3.6.1, PROJ 9.1.1; sf_use_s2() is TRUE
    # Drop datum setting to workaround "GDAL Error 6: Cannot find coordinate operations ..."
    tx_ll <- st_transform(tx, "+proj=longlat +ellps=GRS80")

    gl.coords <- st_coordinates(tx_ll)
    mapPolygon(gl.coords[, "X"], gl.coords[, "Y"], col = "white", border = NA)
  }

  # Get colors for studies
  icol <- envirodat$groupnum[match(coords$study, envirodat$study)]
  # Identify studies that use samples from laboratory or mesocosm experiments
  lab <- c(
    # Sediment
    "BSPD17", "RKN+17", "ZHZ+19", "WFB+21", "HCW+22",
    # Soil
    "BMOB18", "CWC+20", "PSG+20", "DTJ+20",
    "DLS21_Bulk", "CKB+22", "CLZ+22"
  )
  pch <- ifelse(coords$study %in% lab, 15, 19)
  # Use smaller points for high-density regions 20210915
  cex <- ifelse(coords$study %in% c(
    "ZZL+21", "MLL+18", "ZML+17", "ZZLL21", "ZHZ+19", "WHLH21a", # GD-HK-MO GBA
    "HDZ+19", # Hubei
    "WKG+22" # Poland
  ), 1.5, 2.5)
  # Plot sample locations
  mapPoints(coords$longitude, coords$latitude, pch = pch, col = orp16Scol[icol], cex = cex)

  # Plot transects 20210929
  # Coordinates for East Asia Paddy Soil dataset are from
  # Sourcedata.xlsx from https://doi.org/10.6084/m9.figshare.12622829
  dat <- getmdat_orp16S("LJC+20")
  mapPoints(dat$longitude, dat$latitude, col = orp16Scol[7], lwd = 2)
  # Coordinates for Southern Tibetan Plateau dataset are from Table 1 of MWY+21
  dat <- getmdat_orp16S("MWY+21")
  latlon <- paste(dat$latitude, dat$longitude)
  isuniq <- !duplicated(latlon)
  dat <- dat[isuniq, ]
  mapPoints(dat$longitude, dat$latitude, col = orp16Scol[4], lwd = 2)
  # Coordinates for Three Gorges Reservoir are from Table S1 of GZL21
  dat <- getmdat_orp16S("GZL21")
  dat <- dat[!is.na(dat$"ORP (mV)"), ]
  latlon <- paste(dat$Latitude, dat$Longitude)
  isuniq <- !duplicated(latlon)
  dat <- dat[isuniq, ]
  mapPoints(dat$Longitude, dat$Latitude, col = orp16Scol[1], lwd = 1)
  # Coordinates for Port Microbes are from https://github.com/rghannam/portmicrobes/data/metadata/pm_metadata.csv
  dat <- getmdat_orp16S("GSBT20")
  # Use first sample for each port
  dat <- dat[!duplicated(dat$Port), ]
  mapPoints(dat$Longitude, dat$Latitude, col = orp16Scol[1], pch = 17)

  # Add legend
  ienv = c(1, 2, 4, 5, 3, 6, 7)
  par(xpd = NA)
  legend("bottomleft", names(envirotype)[ienv], pch = 19, col = orp16Scol[ienv], bty = "n", cex = 2, pt.cex = 3, inset = c(0, -0.03))
  ltext <- c("Field sites", "Laboratory or mesocosm", "Smaller symbols for", "densely sampled areas", "Port sites", "Open symbols for transects")
  legend("bottomright", ltext, pch = c(19, 15, 20, NA, 17, 1), col = c(1, 1, 1, NA, orp16Scol[1], 1),
         bty = "n", cex = 2, pt.cex = c(3, 3, 3, 3, 2, 2), inset = c(0, -0.03))
  par(xpd = FALSE)
  label.figure("a", font = 2, cex = 3.8)

  ## Panel B: Eh-pH diagram for all environment types 20220516
  par(cex = 2)
  # Get data for unique samples
  ssr <- paste(EZdat$study, EZdat$sample, EZdat$Run, sep = "_")
  idup <- duplicated(ssr)
  thisdat <- EZdat[!idup, ]
  # Remove NA pH
  thisdat <- thisdat[!is.na(thisdat$pH), ]
  # Get colors
  ienv <- match(thisdat$envirotype, names(envirotype))
  col <- orp16Scol[ienv]
  # Make plot
  par(mar = c(4, 4, 1, 1))
  plot(thisdat$pH, thisdat$Eh / 1000, col = col, pch = 19, cex = 0.5, xlab = "pH", ylab = "Eh (V)",
       xlim = c(-2, 14), ylim = c(-0.6, 1), xaxs = "i", yaxs = "i", axes = FALSE)
  box()
  axis(1, at = seq(-2, 14, 2))
  axis(2, at = round(seq(-0.6, 1, 0.2), 1), gap.axis = 0)
  # Add outline from Bass Becking et al. (1960)
  BKM60 <- read.csv(system.file("extdata/orp16S/BKM60.csv", package = "JMDplots"))
  BKM60$Eh <- BKM60$Eh / 1000 # Units: V
  segment <- NULL
  top <- subset(BKM60, segment == "top")
  lines(top$pH, top$Eh)
  bottom <- subset(BKM60, segment == "bottom")
  lines(bottom$pH, bottom$Eh)
  left <- subset(BKM60, segment == "left")
  for(i in 1:(nrow(left)/2)) {
    ii <- c(i*2-1, i*2)
    lines(left$pH[ii], left$Eh[ii])
  }
  right <- subset(BKM60, segment == "right")
  for(i in 1:(nrow(right)/2)) {
    ii <- c(i*2-1, i*2)
    lines(right$pH[ii], right$Eh[ii])
  }
  O2 <- subset(BKM60, segment == "O2")
  lines(O2$pH, O2$Eh, lwd = 2)
  H2 <- subset(BKM60, segment == "H2")
  lines(H2$pH, H2$Eh, lwd = 2)
  text(0, 0.055, quote(H^"+"), cex = 1.2, srt = -30)
  text(-0.5, -0.035, quote(H[2]), cex = 1.2, srt = -30)
  text(5.5, 0.84, quote(H[2]*O), cex = 1.2, srt = -30)
  text(6, 0.93, quote(O[2]), cex = 1.2, srt = -30)
  label.figure("b", font = 2, cex = 1.9)

  # Add legend
  par(mar = c(4, 1, 1, 1))
  plot.new()
  ienv = c(1, 2, 4, 5, 3, 6, 7)
  ltext <- names(envirotype)[ienv]
  legend("left", ltext, pch = 19, col = orp16Scol[ienv], bty = "n", xpd = NA)

  if(pdf) dev.off()

}


# Figure 5: Associations between Eh7 and Zc at local scales 20220517
orp16S_5 <- function(pdf = FALSE) {
  if(pdf) pdf("Figure_5.pdf", width = 8, height = 6, bg = "white")
  mat <- matrix(c(1,1,1,1, 2,2,2,2, 3,3,3,3,
                  0,0,0,0,0,0,0,0,0,0,0,0,
                  4,4,4,4,4, 0, 6,6,6, 7,7,7,
                  5,5,5,5,5, 0, 6,6,6, 7,7,7),
                nrow = 4, byrow = TRUE)
  layout(mat, heights = c(3, 0.75, 3, 3))
  par(mar = c(3, 4, 3, 1))
  par(mgp = c(2.5, 1, 0))

  ## Panel A: Analysis of selected datasets 20211003
  # Daya Bay (Sediment)
  plotEZ("WHLH21a", "Bacteria", groupby = "Position", groups = c("Surface", "Middle", "Bottom"),
    legend.x = "bottomleft", title.line = NULL, dxlim = c(-20, 0), slope.legend = "right")
  title("Daya Bay\n(Sediment bacteria)", font.main = 1)
  label.figure("a", font = 2, cex = 2.3, yfrac = 0.9)
  # Bay of Biscay (Sediment)
  plotEZ("LMBA21_2017", "Bacteria", groupby = "Season", groups = c("Summer", "Winter"),
    legend.x = "bottomright", title.line = NULL, dxlim = c(0, 170), dylim = c(-0.01, 0), slope.legend = "bottom")
  title("Bay of Biscay\n(Sediment bacteria)", font.main = 1)
  # Hunan Soil (Soil)
  plotEZ("MLL+19", "Bacteria", groupby = "Type", groups = c("Upland", "Paddy", "Sediment"),
    title.line = NULL, dylim = c(0, 0.01), slope.legend = "top")
  title("Hunan Province\n(Soil and sediment bacteria)", font.main = 1)

  ## Panel B: Slope vs log10(number of samples) for all datasets
  par(mar = c(3.6, 4, 1, 1))
  # Use Bacteria only 20210913
  lmbac <- EZlm[EZlm$lineage == "Bacteria", ]
  # Get color according to environment group
  env <- envirodat[match(lmbac$study, envirodat$study), ]
  # Get range for included samples 20210905
  i1 <- c(1, 2, 4, 5)
  j1 <- env$groupnum %in% i1
  i2 <- c(3, 6, 7)
  j2 <- env$groupnum %in% i2
  xlim <- range(log10(lmbac$nsamp[j1 | j2]))
  # Get range of slopes (absolute value)
  ymaxabs <- max(abs(lmbac$slope[j1 | j2]))
  ylim <- c(-ymaxabs*1.2, ymaxabs)
  # River & seawater, lake & pond, geothermal, hyperalkaline
  plot(xlim, ylim, type = "n", xlab = quote(log[10]~"(Number of samples)"), ylab = quote("Slope of linear fit"~(1/V)))
  abline(h = 0, lty = 2, lwd = 1.5, col = "gray50")
  abline(h = c(-0.1, 0.1), lty = 3, lwd = 1.5, col = "gray50")
  x1 <- log10(lmbac$nsamp[j1])
  y1 <- lmbac$slope[j1]
  col1 <- orp16Scol[env$groupnum[j1]]
  # Make bigger points for larger slopes 20221230
  cex1 <- ifelse(abs(y1) > 0.1, 1.5, 1)
  points(x1, y1, pch = 19, col = col1, cex = cex1)
  # Add legend
  ltext <- names(envirotype)[i1[1:2]]
  legend("bottomleft", ltext, pch = 19, col = orp16Scol[i1[1:2]])
  ltext <- names(envirotype)[i1[3:4]]
  legend("bottomright", ltext, pch = 19, col = orp16Scol[i1[3:4]])
  title("Linear regressions for bacterial communities\nin each dataset", font.main = 1, xpd = NA, line = 0.7)
  label.figure("b", font = 2, cex = 2.3, yfrac = 1.05)
  # Groundwater, sediment, soil
  plot(xlim, ylim, type = "n", xlab = quote(log[10]~"(Number of samples)"), ylab = quote("Slope of linear fit"~(1/V)))
  abline(h = 0, lty = 2, lwd = 1.5, col = "gray50")
  abline(h = c(-0.1, 0.1), lty = 3, lwd = 1.5, col = "gray50")
  x2 <- log10(lmbac$nsamp[j2])
  y2 <- lmbac$slope[j2]
  col2 <- orp16Scol[env$groupnum[j2]]
  cex2 <- ifelse(abs(y2) > 0.1, 1.5, 1)
  points(x2, y2, pch = 19, col = col2, cex = cex2)
  ltext <- names(envirotype)[i2]
  legend("bottomright", ltext, pch = 19, col = orp16Scol[i2], bg = "white")

  ## Panel C: Distinctions in carbon oxidation state estimated for different geothermal areas 20210930
  # Use Geothermal datasets
  i <- 4
  studies <- envirotype[[i]]
  # Assign colors
  # Use blue for neutral/alkaline
  col <- rep(orp16Scol[1], length(studies))
  # Use red for acidic
  iacid <- studies %in% c("PCL+18_Acidic", "LMG+20")
  col[iacid] <- orp16Scol[4]
  # Use gray for sediments
  ised <- studies %in% c("PBU+20", "MWY+21")
  col[ised] <- "gray"

  # Offset for labels 20211012
  dx <- list(
    c(0.04, -0.19, -0.22, 0.04, -0.18),
    c(-0.38, -0.13, -0.12, NA, -0.20)
  )
  dy <- list(
    c(0, -0.015, -0.03, 0, -0.01),
    c(-0.002, 0.03, -0.007, NA, 0.007)
  )

  # Loop over Bacteria and Archaea
  lineages <- c("Bacteria", "Archaea")
  for(k in 1:2) {
    dat <- EZlm[EZlm$lineage == lineages[k], ]
    plot(c(-0.45, 0.40), c(-0.24, -0.12), type = "n", xlab = "Eh7 (V)", ylab = cplab$Zc)
    for(j in seq_along(studies)) {
      with(dat[dat$study == studies[j], ], {
        # Divide by 1000 to convert slope from V-1 to mV-1 (used by lm fit) 20220520
        y <- function(x) intercept + slope * x / 1000
        Eh7 <- c(Eh7min, Eh7max)
        Zc <- y(Eh7)
        if(length(Eh7) > 0) {
          lines(Eh7 / 1000, Zc, col = col[j], lwd = 2)
          # Add number to identify dataset
          text(tail(Eh7 / 1000, 1) + dx[[k]][j], tail(Zc, 1) + dy[[k]][j], j)
        }
      })
    }
    title(paste0("Geothermal\n", tolower(lineages[k])), font.main = 1, xpd = NA, line = 0.7)
    if(k==1) {
      label.figure("c", font = 2, cex = 2.3, yfrac = 1.025)
      # Add legend
      lacid <- paste0("Acidic (", paste(which(iacid), collapse = ", "), ")")
      lneut <- "Circumneutral to"
      lalk <- paste0("alkaline (", paste(which(col == orp16Scol[1]), collapse = ", "), ")")
      lsed <- paste0("Sediment (", paste(which(ised), collapse = ", "), ")")
      legend("topleft", c(lacid, lneut, lalk, lsed), col = c(orp16Scol[4], orp16Scol[1], NA, "gray"), lwd = 2, bty = "n")
    }
  }

  if(pdf) dev.off()

}


# Figure 6: Associations between Eh7 and Zc at a global scale 20210828
orp16S_6 <- function(pdf = FALSE, EMP_primers = FALSE) {

  if(pdf) {
    if(EMP_primers) pdf("Figure_S3.pdf", width = 10, height = 7, bg = "white")
    else pdf("Figure_6.pdf", width = 10, height = 7, bg = "white")
  }
  mat <- matrix(c(
    16, 16, rep(1:7, each = 6), 17, 17,
    16, 16, rep(8:14, each = 6), 18, 18,
    rep(15, 46),
    rep(0, 46),
    rep(0, 5), rep(19, 16), rep(20, 16), rep(21, 9)
  ), nrow = 5, byrow = TRUE)
  layout(mat, heights = c(2, 2, 0.5, 0.5, 4))
  par(mar = c(4, 4, 1, 1))
  par(mgp = c(2.5, 1, 0))

  ## Take only datasets that use EMP primers 20221004
  if(EMP_primers) {
    uses_EMP_primers <- c(
      "MLL+18", "LWJ+21", "RARG22",                           # River & Seawater
      "MTC21",                                                # Lake & Pond
      "PCL+18_Acidic", "PCL+18_Alkaline", "GWS+20", "MWY+21", # Geothermal
      "SBP+20", "RMB+17", "PSB+21",                           # Hyperalkaline
      "WLJ+16", "ZDW+19", "DJK+18", "APV+20", "MGW+22",       # Groundwater
      "OHL+18_DNA", "ZHZ+19",                                 # Sediment
      "MLL+19", "PSG+20", "RKSK22", "CKB+22"                  # Soil
    )
    EZdat <- EZdat[EZdat$study %in% uses_EMP_primers, ]
    # Make sure we didn't miss any studies (check for typos ...)
    stopifnot(all(uses_EMP_primers %in% EZdat$study))
  }

  ## Panel A: Scatterplots and fits for Bacteria and Archaea in each environment
  par(mar = c(0, 0, 1, 0))
  xlim <- range(EZdat$Eh7)
  ylim <- range(EZdat$Zc, -0.1)
  eedat <- EZdat[EZdat$lineage == "Bacteria", ]
  global.slopes <- list()
  global.slopes$Bacteria <- eachenv(eedat, xlim = xlim, ylim = ylim, lineage = "Bacteria")
  par(mar = c(1, 0, 0, 0))
  eedat <- EZdat[EZdat$lineage == "Archaea", ]
  global.slopes$Archaea <- eachenv(eedat, xlim = xlim, ylim = ylim, lineage = "Archaea")
  # Add labels
  plot.new()
  text(0.5, -0.5, "Eh7 (V)", cex = 1.2, xpd = NA)
  plot.new()
  text(0.2, 0.5, cplab$Zc, cex = 1.2, srt = 90, xpd = NA)
  label.figure("a", font = 2, cex = 2, xfrac = 0.2, yfrac = 0.97)
  plot.new()
  text(0.3, 0.5, "Bacteria", srt = 90, xpd = NA)
  plot.new()
  text(0.3, 0.5, "Archaea", srt = 90, xpd = NA)
 
  ## Panel B: Scatterplots and fits for Bacteria and Archaea in all environments 20210914
  # Start plot for Bacteria
  par(mar = c(4, 4, 1, 1))
  plot(c(-0.50, 0.65), range(EZdat$Zc), type = "n", xlab = "Eh7 (V)", ylab = cplab$Zc)
  # Use Bacteria only
  thisdat <- EZdat[EZdat$lineage == "Bacteria", ]
  # Add linear fit; include number of studies in legend 20210925
  nstudy <- length(unique(thisdat$study))
  add.linear.global(thisdat$Eh7, thisdat$Zc, nstudy)
  # Add points
  eachenv(thisdat, add = TRUE, do.linear = FALSE)
  title("Bacteria", font.main = 1, line = 0.5, xpd = NA)
  label.figure("b", font = 2, cex = 2, xfrac = 0.05, yfrac = 1)

  # Now do Archaea
  plot(c(-0.50, 0.65), range(EZdat$Zc), type = "n", xlab = "Eh7 (V)", ylab = cplab$Zc)
  thisdat <- EZdat[EZdat$lineage == "Archaea", ]
  nstudy <- length(unique(thisdat$study))
  add.linear.global(thisdat$Eh7, thisdat$Zc, nstudy, inset = c(0, 0.05))
  eachenv(thisdat, add = TRUE, do.linear = FALSE)
  title("Archaea", font.main = 1, line = 0.5, xpd = NA)

  # Add legend
  par(mar = c(4, 1, 1, 1))
  plot.new()
  ienv = c(1, 2, 4, 5, 3, 6, 7)
  ltext <- names(envirotype)[ienv]
  legend("left", ltext, pch = 19, col = orp16Scol[ienv], bty = "n")

  if(pdf) dev.off()
  # Return slopes 20220611
  invisible(global.slopes)

}


###########################################
### Functions for supplementary figures ###
###########################################

# Figure S1: Zc-Eh scatterplots for all studies 20210827
# This also creates files EZdat (Eh and Zc values) and
# EZlm (linear fits) for use by other plotting functions
orp16S_S1 <- function(pdf = FALSE) {

  # Setup figure
  if(pdf) pdf("Figure_S1.pdf", width = 9, height = 12, bg = "white")
  par(mfrow = c(4, 3))

  results <- c(

    message("\nRiver & Seawater"),
    plotEZ("MLL+18", "Bacteria", groupby = "Season", groups = c("Summer", "Winter"), legend.x = "bottomright"),
    plotEZ("HXZ+20", "Bacteria", groupby = "Station", groups = c("SYBL", "C4")),
    plotEZ("GSBT20_Prefilter", "two", groupby = "Region", groups = c("West Coast U.S.", "Great Lakes", "East Coast U.S.", "Europe", "Asia"),
           legend.x = "bottomright", dylim = c(-0.015, 0)),
    plotEZ("GSBT20_Postfilter", "two", groupby = "Region", groups = c("West Coast U.S.", "Great Lakes", "East Coast U.S.", "Europe", "Asia"),
           legend.x = "bottomright", dylim = c(-0.007, 0), dxlim = c(0, 75)),
    plotEZ("WHL+21", "Bacteria", groupby = "Season", groups = c("Spring", "Summer", "Autumn", "Winter"), legend.x = "bottomleft"),
    plotEZ("ZLH+22", "Bacteria"),
    plotEZ("ZZL+21", "Bacteria", groupby = "Location", groups = c("Main Stream", "Animal Farm", "Hospital", "WWTP", "Tributary"),
           legend.x = "bottomright", dxlim = c(0, 100)),
    plotEZ("LWJ+21", "two", groupby = "Type", groups = c("Freshwater", "Freshwater Plastic", "Seawater", "Seawater Plastic"), dylim = c(0, 0.01)),
    plotEZ("GZL21", "Bacteria", groupby = "Type", groups = c("Surface water", "Middle water", "Bottom water"), legend.x = "bottomleft"),
    plotEZ("RARG22", "two", groupby = "Group", groups = c("AMD", "Abiotic Treatment", "Biotic Treatment", "Stream"), dylim = c(0, 0.007)),

    message("\nLake & Pond"),
    plotEZ("SAR+13", "two", groupby = "Zone", groups = c("Photic-oxic", "Transition", "Anoxic")),
    plotEZ("LLC+19", "Bacteria", groupby = "Size", groups = c("Free-living", "Particle-associated")),
    plotEZ("BCA+21", "Bacteria", groupby = "Month", groups = c("Jul", "Nov", "Feb", "Apr")),
    plotEZ("HLZ+18", "Bacteria", groupby = "Type", groups = c("Reservoir", "Pond"), legend.x = "bottomright"),
    plotEZ("BWD+19", "Bacteria", groupby = "Cover", groups = c("Ice", "Ice Free"), legend.x = "bottomright"),
    plotEZ("IBK+22", "two", groupby = "Land Use", groups = c("Arable", "Forest", "Grassland")),
    plotEZ("NLE+21", "Bacteria", groupby = "Year", groups = c("2017", "2018"), legend.x = "bottomleft"),
    plotEZ("MTC21", "two", groupby = "Year", groups = c(2012, 2013, 2014)),
    plotEZ("SPA+21", "Bacteria", groupby = "Depth", groups = c("Epi", "Secchi", "Meso"), legend.x = "bottomleft"),

    message("\nGeothermal"),
    plotEZ("PCL+18_Acidic", "two", legend.x = "bottomright"),
    plotEZ("PCL+18_Alkaline", "two"),
    plotEZ("GWS+20", "two", groupby = "Hydrothermal Field", groups = c("Batang", "Litang", "Kangding")),
    plotEZ("PBU+20", "Bacteria", groupby = "Type", groups = c("Cauldron", "Sampling Pit", "Spring", "Geyser Valley (Control)"), legend.x = "bottomright"),
    plotEZ("MWY+21", "two", groupby = "Location", groups = c("Quseyongba", "Moluojiang", "Daggyai", "Quzhuomu"), legend.x = "bottomright"),

    message("\nHyperalkaline"),
    plotEZ("SBP+20", "Bacteria", groupby = "pH Group", groups = c("< 10", "> 10")),
    plotEZ("RMB+17", "two", groupby = "pH Group", groups = c("< 10", "> 10")),
    plotEZ("CTS+17", "two", groupby = "Type", groups = c("River", "Well", "Spring"), legend.x = "bottomleft"),
    plotEZ("KSR+21", "Bacteria", groupby = "Location", groups = c("Lerone", "Branega", "Branega Creek Water")),
    plotEZ("PSB+21", "Bacteria", groupby = "O2 range", groups = c("> 0.5 mg/L", "0.2-0.5 mg/L", "< 0.2 mg/L"), dxlim = c(-100, 0), dylim = c(0, 0.01)),
    plotEZ("NTB+21", "two", groupby = "Well", groups = c("BA1A", "BA1D")),

    message("\nGroundwater"),
    plotEZ("KLM+16", "Bacteria", groupby = "Day", groups = c(-1, 246, 448, 671)),
    plotEZ("WLJ+16", "two", groupby = "Arsenic", groups = c("As < 100 \u03BCg/L", "As > 100 \u03BCg/L")),
    plotEZ("ZDW+19", "two", groupby = "Season", groups = c("Non-monsoon", "Spring", "Monsoon", "Autumn"), legend.x = "bottomright"),
    plotEZ("DJK+18", "two", groupby = "Aquifer", groups = c("Athens", "Greene", "Licking"), legend.x = "bottomleft"),
    plotEZ("SRM+19", "Bacteria", groupby = "Land Use", groups = c("Agriculture", "Community", "Landfill", "Mine")),
    plotEZ("APV+20", "two", groupby = "Type", groups = c("Canal", "Piezometer", "Well", "Spring")),
    plotEZ("YHK+20", "Bacteria", groupby = "Location", groups = c("Upper Hillslope", "Middle Slope", "Lower Footslope"), dylim = c(0, 0.005)),
    plotEZ("ZCZ+21", "Bacteria", groupby = "Location", groups = c("LO", "CR1", "MN", "VA", "BS", "CR2"), legend.x = "topright"),
    plotEZ("MGW+22", "two", groupby = "Region", groups = c("Auckland", "Canterbury", "Taupo", "Wellington"), legend.x = "bottomleft", dylim = c(-0.006, 0)),
    plotEZ("MCR+22", "two", groupby = "Location", groups = c("Upstream", "WWTP", "Farm", "Swamp"), legend.x = "bottomleft"),

    message("\nSediment"),
    plotEZ("ZML+17", "two", groupby = "Depth", groups = c("\u2264 5 cm", "\u2265 10 cm", "\u2265 20 cm"), legend.x = "bottomleft", dylim = c(-0.003, 0)),
    plotEZ("BSPD17", "Bacteria", groupby = "Type", groups = c("Water", "Sediment")),
    plotEZ("RKN+17", "two", groupby = "Treatment", groups = c("CH4", "CH4+Fe", "CH4+Mn"), legend.x = "bottomleft"),
    plotEZ("HDZ+19", "Bacteria", groupby = "Type", groups = c("Water", "Sediment")),
    plotEZ("OHL+18_DNA", "two", groupby = "Location", groups = c("Kal\u00f8 Vig", "Lake Constance", "Norsminde Fjord")),
    plotEZ("WHLH21a", "Bacteria", groupby = "Position", groups = c("Surface", "Middle", "Bottom"), legend.x = "bottomleft"),
    plotEZ("RSS+18", "Bacteria", groupby = "Site", groups = c("Deep Hole", "Snowgoose Bay", "John's Island", "Skeleton Lake"), dylim = c(0, 0.005)),
    plotEZ("CLS+19", "two", groupby = "Type", groups = c("Water", "Sediment"), legend.x = "bottomright", dylim = c(-0.002, 0)),
    plotEZ("HSF+19", "Bacteria", groupby = "Type", groups = c("Water", "Sediment-Water Interface", "Sediment"), legend.x = "topright"),
    plotEZ("ZHZ+19", "two", groupby = "Treatment", groups = c("Original", "Nitrate-reducing", "Ferric-reducing", "Sulfate-reducing", "Methanogenic"),
           dylim = c(0, 0.005)),
    plotEZ("LMBA21_2017", "two", groupby = "Season", groups = c("Summer", "Winter"), legend.x = "bottomright", dylim = c(-0.005, 0)),
    plotEZ("HSF+22", "Bacteria", groupby = "Location", groups = c("West Lagoon", "North Lagoon", "South Lagoon", "Cinq Cases"),
           legend.x = "bottomright", dxlim = c(0, 100)),
    plotEZ("ZZLL21", "Bacteria", groupby = "Type", groups = c("Main stream", "Animal farm", "Hospital", "WWTP", "Tributary"),
           legend.x = "bottomleft", dylim = c(-0.005, 0)),
    plotEZ("WFB+21", "Bacteria", groupby = "Treatment", groups = c("C. volutator", "H. diversicolor", "Cv & Hd", "MPB", "Manual turbation"),
           dylim = c(0, 0.002)),
    plotEZ("HCW+22", "Bacteria", groupby = "Condition", groups = c("Static", "Weak", "Strong"), dylim = c(0, 0.007)),
    plotEZ("WKG+22", "Bacteria", groupby = "Season", groups = c("Spring", "Summer", "Autumn", "Winter"), dylim = c(0, 0.002)),

    message("\nSoil"),
    plotEZ("MLL+19", "two", groupby = "Type", groups = c("Upland", "Paddy", "Sediment")),
    plotEZ("BMOB18", "two", groupby = "Treatment", groups = c("Acetate", "No amendment", "Pre-incubation"), dxlim = c(-50, 0)),
    plotEZ("WHLH21", "Bacteria", groupby = "Stage", groups = c("Cyanobacteria", "Cyanolichen", "Chlorolichen", "Moss"), legend.x = "bottomright"),
    plotEZ("CWC+20", "Bacteria", groupby = "Management", groups = c("Flooding", "Draining"), legend.x = "bottomright"),
    plotEZ("PSG+20", "two", groupby = "Treatment", groups = c("Initial", "NCC", "RB", "RGP", "TP"), legend.x = "bottomright"),
    plotEZ("LJC+20", "Bacteria", groupby = "MAT", groups = c(">= 21.5 \u00b0C", "< 21.5 \u00b0C"), dylim = c(0, 0.005)),
    plotEZ("DTJ+20", "Bacteria", groupby = "Zone", groups = c("Bulk Soil", "Mature", "Elongation", "Tip")),
    plotEZ("RKSK22", "two", groupby = "Compartment", groups = c("Bulk sediment", "Rhizosphere", "Root"), legend.x = "bottomright"),
    plotEZ("DLS21_Bulk", "Bacteria", groupby = "Treatment", groups = c("Control", "Char", "Silicate", "Husk")),
    plotEZ("WKP+22", "Bacteria", groupby = "Type", groups = c("Intercropping", "Monoculture"), legend.x = "bottomleft"),
    plotEZ("CKB+22", "two", groupby = "Treatment", groups = c("FM 5 Mg/ha", "FM 10 Mg/ha", "RS 5 Mg/ha", "RS 10 Mg/ha"), legend.x = "bottomright"),
    plotEZ("CLZ+22", "Bacteria")

  )
  # Done plotting!
  if(pdf) dev.off()

  # Assemble all data
  EZdat <- do.call(rbind, results[names(results) == "EZdat"])
  # Assemble all linear model results
  coefficients <- sapply(results[names(results) == "EZlm"], "[", "coefficients")
  intercept <- sapply(coefficients, "[[", "(Intercept)")
  slope <- sapply(coefficients, "[[", "Eh7")
  model <- lapply(results[names(results) == "EZlm"], "[[", "model")
  Eh7lim <- sapply(lapply(model, "[", "Eh7"), "range")
  Eh7min <- Eh7lim[1, ]
  Eh7max <- Eh7lim[2, ]
  # Get MOE95 (95% CI) 20221008
  CI95 <- lapply(results[names(results) == "EZlm"], "confint", parm = "Eh7", level = 0.95)
  MOE95 <- sapply(lapply(CI95, range), diff) / 2
  # Get study name etc.
  study <- sapply(results[names(results) == "study"], "[")
  name <- sapply(results[names(results) == "name"], "[")
  envirotype <- sapply(results[names(results) == "envirotype"], "[")
  lineage <- sapply(results[names(results) == "lineage"], "[")
  nsamp <- sapply(model, nrow)
  pearson.r <- unlist(lapply(results[names(results) == "pearson"], "[[", "estimate"))
  # Note: slope and MOE95 are mutiplied by 1000 to convert from mV-1 to V-1
  EZlm <- data.frame(study, name, envirotype, lineage, nsamp, Eh7min, Eh7max,
    slope = signif(slope * 1000, 6), MOE95 = signif(MOE95 * 1000, 6), intercept = signif(intercept, 6), pearson.r = signif(pearson.r, 6))
  # Save data and results to files
  write.csv(EZdat, "EZdat.csv", row.names = FALSE, quote = FALSE)
  write.csv(EZlm, "EZlm.csv", row.names = FALSE, quote = 2)

}


# Figure S2: Compare regressions with Eh7, Eh, and O2 20220517
orp16S_S2 <- function(pdf = FALSE) {

  if(pdf) pdf("Figure_S2.pdf", width = 6, height = 8, bg = "white")
  mat <- matrix(1:8, ncol = 2)
  layout(mat, heights = c(2, 2, 2, 1))

  # Use only samples with non-NA O2
  hasO2dat <- EZdat[!is.na(EZdat$O2_umol_L), ] 
  ylim <- range(hasO2dat$Zc, -0.09)

  # Loop over domains
  for(lineage in c("Bacteria", "Archaea")) {
    par(mar = c(4, 4, 1.5, 1))
    thisdat <- hasO2dat[hasO2dat$lineage == lineage, ]
    # Zc-Eh7 plot
    plot(c(-0.50, 0.70), ylim, type = "n", xlab = "Eh7 (V)", ylab = cplab$Zc)
    # Add linear fit; include number of studies in legend 20210925
    nstudy <- length(unique(thisdat$study))
    add.linear.global(thisdat$Eh7, thisdat$Zc, nstudy)
    # Add points
    eachenv(thisdat, add = TRUE, do.linear = FALSE, xvar = "Eh7")
    title(lineage, font.main = 1)
    # Zc-Eh plot
    plot(c(-0.50, 0.70), ylim, type = "n", xlab = "Eh (V)", ylab = cplab$Zc)
    add.linear.global(thisdat$Eh, thisdat$Zc, nstudy, xvar = "Eh")
    eachenv(thisdat, add = TRUE, do.linear = FALSE, xvar = "Eh")
    # Zc-O2 plot
    plot(c(0, 750), ylim, type = "n", xlab = quote(O[2]~"("*mu*"mol"/"liter"*")"), ylab = cplab$Zc)
    add.linear.global(thisdat$O2_umol_L, thisdat$Zc, nstudy, xvar = "O2", scale = 1)
    eachenv(thisdat, add = TRUE, do.linear = FALSE, xvar = "O2")
    # Add legend
    par(mar = c(0, 0, 0, 0))
    plot.new()
    ienv = c(1, 2, 4, 5, 3, 6, 7)
    ltext <- names(envirotype)[ienv]
    # Add number of samples in each environment 20220518
    nsamp <- table(thisdat$envirotype)[ltext]
    nsamp[is.na(nsamp)] <- 0
    ltext <- paste0(ltext, " (", nsamp, ")")
    legend("bottom", ltext, pch = 19, col = orp16Scol[ienv], bty = "n")
  }

  if(pdf) dev.off()

}

### NOTE: Figure S3 is made with orp16S_4(EMP_primers = FALSE)

##########################################
### Functions for tables and datasets ####
##########################################

# Function to gather info (study name, number of samples with Bacteria and Archaea,
# T, pH, Eh, and Eh7 ranges) for filling in Table S1 20220513
# Rewritten to use precomputed values from EZdat 20220611
orp16S_info <- function(study) {
  # Get data and linear regression for this dataset
  thisdat <- EZdat[EZdat$study == study, ]
  thislm <- EZlm[EZlm$study == study, ]
  # Print dataset name
  print(name <- thislm$name[1])
  # Number of samples with Bacteria and Archaea
  nBac <- sum(thisdat$lineage == "Bacteria")
  nArc <- sum(thisdat$lineage == "Archaea")
  print(paste0("nBac: ", nBac, "; nArc: ", nArc))
  # Format the T, pH, Eh, Eh7 ranges
  if(all(is.na(thisdat$T))) T <- NA else T <- round(unique(range(thisdat$T, na.rm = TRUE)), 1)
  if(!all(is.na(T))) T <- paste(T, collapse = " to ")
  if(all(is.na(thisdat$pH))) pH <- NA else pH <- round(unique(range(thisdat$pH, na.rm = TRUE)), 1)
  if(!all(is.na(pH))) pH <- paste(pH, collapse = " to ")
  if(all(is.na(thisdat$Eh))) Eh <- NA else Eh <- round(unique(range(thisdat$Eh, na.rm = TRUE)))
  if(!all(is.na(Eh))) Eh <- paste(Eh, collapse = " to ")
  if(all(is.na(thisdat$Eh7))) Eh7 <- NA else Eh7 <- round(unique(range(thisdat$Eh7, na.rm = TRUE)))
  if(!all(is.na(Eh7))) Eh7 <- paste(Eh7, collapse = " to ")
  infotxt <- paste0(study, ": T ", T, ", pH ", pH, ", Eh ", Eh, ", Eh7 ", Eh7)
  print(infotxt)
  # Read linear fit coefficients
  bacslope <- NA
  for(lineage in c("Bacteria", "Archaea")) {
    idat <- thislm$lineage == lineage
    if(any(idat)) {
      slope <- thislm$slope[idat]
      if(slope > 0) slopetxt <- "positive"
      if(slope < 0) slopetxt <- "negative"
      if(lineage == "Bacteria") bacslope <- slopetxt
      slopetxt <- paste("Slope of Zc-Eh7 correlation for", lineage, "is", slopetxt)
      print(slopetxt)
    }
  }
  # Return values for Table S1 invisibly 20220520
  if(!is.na(bacslope)) bacslope <- strsplit(bacslope, " ")[[1]][1]
  out <- data.frame(study = study, name = name, nBac = nBac, nArc = nArc, T = T, pH = pH, Eh = Eh, Eh7 = Eh7, bacslope = bacslope)
  invisible(out)
}

# Table of regression slopes (positive/negative) and p-values 20220516
orp16S_T1 <- function(samesign = FALSE) {
  # All environment types
  envirotype <- unique(EZlm$envirotype)
  # Loop over Bacteria and Archaea
  out <- lapply(c("Bacteria", "Archaea"), function(lineage) {
    # Initialize table
    out <- matrix(nrow = length(envirotype), ncol = 4)
    # Loop over environment type
    for(ienv in 1:length(envirotype)) {
      # Get regression results
      thisdat <- EZlm[EZlm$lineage == lineage & EZlm$envirotype == envirotype[ienv], ]
      # Filter to keep only datasets with same sign of minimum and maximum slope within 95% CI 20221022
      if(samesign) thisdat <- thisdat[sign(thisdat$slope + thisdat$MOE95) == sign(thisdat$slope - thisdat$MOE95), ]
      # Number of datasets and numbers with positive and negative slopes
      ndat <- nrow(thisdat)
      npos <- sum(thisdat$slope > 0)
      nneg <- sum(thisdat$slope < 0)
      # Enter counts into table
      out[ienv, 1] <- ndat
      out[ienv, 2] <- npos
      out[ienv, 3] <- nneg
    }
    # Add a row for totals
    out <- rbind(out, colSums(out))
    # Calculate binomial probabilities
    for(i in 1:nrow(out)) if(out[i, 1] != 0) out[i, 4] <- binom.test(out[i, 2], out[i, 1], p = 0.5, alternative = "greater")$p.value
    out[1:length(envirotype), 4] <- round(out[1:length(envirotype), 4], 3)
    out[nrow(out), 4] <- round(out[nrow(out), 4], 5)
    out
  })
  out <- do.call(cbind, out)
  # Add names
  rownames(out) <- c(envirotype, "Total")
  colnames(out) <- paste(rep(c("N", "Pos", "Neg", "p"), 2), rep(c("Bac", "Arc"), each = 4), sep = "_")
  out
}

### NOTE: Dataset S1 is the EZdat.csv file made with orp16S_S1()
### NOTE: Dataset S2 is the EZlm.csv file made with orp16S_S1()

# Dataset S3: Most abundant genera at high and low Eh7 20221006
orp16S_D3 <- function(mincount = 100) {

  # Get amino acid compositions of taxa compiled from RefSeq sequences
  taxon_AA <- taxon_AA$RefSeq
  # Calculate Zc for all taxa
  Zc <- Zc(taxon_AA)

  # Use EZlm as a template for our output
  out <- lapply(1:nrow(EZlm), function(i) {

    study <- EZlm$study[i]
    lineage <- EZlm$lineage[i]

    Q1.genus <- NA
    Q1.Zc <- NA
    Q4.genus <- NA
    Q4.Zc <- NA
    mapperc <- 0

    # Get RDP counts, mapping to NCBI taxonomy, and chemical metrics
    studyfile <- gsub("_.*", "", study)
    datadir <- system.file("extdata/orp16S/RDP", package = "JMDplots")
    RDPfile <- file.path(datadir, paste0(studyfile, ".tab.xz"))
    # If there is no .xz file, look for a .tab file 20210607
    if(!file.exists(RDPfile)) RDPfile <- file.path(datadir, paste0(studyfile, ".tab"))
    # Use try() to catch error for no mapped sequences for Archaea
    RDP <- try(read_RDP(RDPfile, lineage = lineage, mincount = mincount), silent = TRUE)
    if(!inherits(RDP, "try-error")) {

      # Calculate metrics to make sure we get the same samples used in the analysis for the paper
      map <- map_taxa(RDP, refdb = "RefSeq")
      metrics <- getmetrics_orp16S(study, lineage = lineage, mincount = mincount)
      mdat <- getmdat_orp16S(study, metrics)
      metadata <- mdat$metadata
      # Drop rows with NA Eh7
      ina <- is.na(metadata$Eh7)
      metadata <- metadata[!ina, ]

      # Keep only genus-level classifications
      igenus <- RDP$rank == "genus"
      RDP <- RDP[igenus, , drop = FALSE]
      map <- map[igenus]
      # Calculate percentage of mapped genera
      is.mapped <- !is.na(map)
      totcount <- sum(RDP[, -(1:4), drop = FALSE])
      mapcount <- sum(RDP[is.mapped, -(1:4), drop = FALSE])
      mapperc <- round(mapcount / totcount * 100, 2)
      # Keep only the mapped taxa
      RDP <- RDP[is.mapped, , drop = FALSE]
      map <- map[is.mapped]

      # Get the lower and upper quartiles for Eh7
      Eh7 <- metadata$Eh7
      Qval <- quantile(Eh7, probs = c(0, 0.25, 0.5, 0.75, 1))
      # Use try() to catch error if 'breaks' are not unique (e.g. only 3 values for Eh7)
      Qind <- try(cut(Eh7, Qval, labels = 1:4, include.lowest = TRUE), silent = TRUE)
      if(inherits(Qind, "try-error")) {
        Q1 <- Eh7 <= Qval[2]
        Q4 <- Eh7 >= Qval[4]
      } else {
        Q1 <- which(Qind == 1)
        Q4 <- which(Qind == 4)
      }

      # Find the most abundant genus for Q1 and Q4 (low and high Eh7)
      Q1run <- metadata$Run[Q1]
      Q1count <- rowSums(RDP[, Q1run, drop = FALSE])
      Q1max <- which.max(Q1count)
      Q1map <- map[Q1max]
      Q4run <- metadata$Run[Q4]
      Q4count <- rowSums(RDP[, Q4run, drop = FALSE])
      Q4max <- which.max(Q4count)
      Q4map <- map[Q4max]

      # Get genus names and Zc for Q1 and Q4
      Q1.genus <- taxon_AA$organism[Q1map]
      Q1.Zc <- round(Zc[Q1map], 6)
      Q4.genus <- taxon_AA$organism[Q4map]
      Q4.Zc <- round(Zc[Q4map], 6)

    }

    out <- data.frame(mapperc, Q1.genus, Q1.Zc, Q4.genus, Q4.Zc)
    out

  })

  out <- do.call(rbind, out)
  # Prepend study, name, envirotype, lineage columns from EZlm
  out <- cbind(EZlm[, 1:4], out)
  # Save to file
  write.csv(out, "Dataset_S3.csv", row.names = FALSE, quote = 2)

}

#################################
### Data processing functions ###
#################################

# Get metadata for a study, appending columns for pch and col 20200914
getmdat_orp16S <- function(study, metrics = NULL, dropNA = TRUE, size = NULL, quiet = TRUE) {
  # Read metadata file
  # Remove suffix after underscore 20200929
  studyfile <- gsub("_.*", "", study)
  datadir <- system.file("extdata/orp16S", package = "JMDplots")
  file <- file.path(datadir, "metadata", paste0(studyfile, ".csv"))
  metadata <- read.csv(file, as.is = TRUE, check.names = FALSE)

  if(dropNA) {
    # Exclude samples with NA name 20200916
    iname <- match("name", tolower(colnames(metadata)))
    noname <- is.na(metadata[, iname])
    if(any(noname)) {
      if(!quiet) print(paste0("getmdat_orp16S [", study, "]: dropping ", sum(noname), " samples with NA name"))
      metadata <- metadata[!is.na(metadata[, iname]), ]
    }
  }
  # Use NULL pch as flag for unavailable dataset 20210820
  pch <- NULL
  infotext <- NULL

  ## Datasets for orp16S paper 20211003
  shortstudy <- study
  if(study == "RBW+14") {
    type <- rep("reducing", nrow(metadata))
    type[metadata$layer == "Top"] <- "oxidizing"
    type[metadata$layer == "SWI"] <- "transition"
    pch <- sapply(type, switch, oxidizing = 24, transition = 20, reducing = 25)
    col <- sapply(type, switch, oxidizing = 4, transition = 1, reducing = 2)
  }
  if(grepl("PCL\\+18", study)) {
    # PCL+18, PCL+18_Acidic, PCL+18_Alkaline
    Type <- sapply(strsplit(study, "_"), "[", 2)
    if(!is.na(Type)) metadata <- metadata[metadata$Type == Type, ]
    shortstudy <- "PCL+18"
  }
  if(grepl("LMBA21", study)) {
    # LMBA21, LMBA21_2013, LMBA21_2015, LMBA21_2016, LMBA21_2017, LMBA21_2018, LMBA21_2019
    Year <- sapply(strsplit(study, "_"), "[", 2)
    if(!is.na(Year)) metadata <- metadata[metadata$Year == Year, ]
    shortstudy <- "LMBA21"
  }
  if(grepl("DLS21", study)) {
    # DLS21, DLS21_Bulk, DLS21_Rhizosphere, DLS21_Plaque
    source <- sapply(strsplit(study, "_"), "[", 2)
    if(!is.na(source)) {
      if(source == "Bulk") metadata <- metadata[metadata$Source == "bulk soil", ]
      if(source == "Rhizosphere") metadata <- metadata[metadata$Source == "rhizosphere soil", ]
      if(source == "Plaque") metadata <- metadata[metadata$Source == "iron plaque", ]
    }
    shortstudy <- "DLS21"
  }
  if(grepl("GSBT20", study)) {
    # GSBT20, GSBT20_Prefilter, GSBT20_Postfilter
    filter <- sapply(strsplit(study, "_"), "[", 2)
    if(!is.na(filter)) {
      if(filter == "Postfilter") metadata <- metadata[grepl("post", metadata$Sample), ]
      if(filter == "Prefilter") metadata <- metadata[!grepl("post", metadata$Sample), ]
    }
    shortstudy <- "GSBT20"
  }
  if(grepl("OHL\\+18", study)) {
    # OHL+18, OHL+18_DNA, OHL+18_cDNA
    Molecule <- sapply(strsplit(study, "_"), "[", 2)
    if(!is.na(Molecule)) metadata <- metadata[metadata$Molecule == Molecule, ]
    shortstudy <- "OHL+18"
  }
  if(grepl("WHLH21", study)) {
    # WHLH21, WHLH21_Jun, WHLH21_Sep
    Month <- sapply(strsplit(study, "_"), "[", 2)
    if(!is.na(Month)) metadata <- metadata[grepl(Month, metadata$LibraryName), ]
    shortstudy <- "WHLH21"
  }

  ## Additional datasets for comparison with metaproteomes 20221222
  if(study %in% c(
    # 20220829 Manus Basin Inactive Chimney, M.B. Active Chimneys
    "MPB+19", "RYP+14",
    # 20221028 Soda Lake Biomats, Mock Communities, Saanich Inlet
    "KTS+17", "KTS+17.mock", "HTZ+17"
  )) {
    pch <- 20
    col <- 1
  }

  # Further processing of metadata for orp16S datasets 20210820
  if(shortstudy %in% c(
    "MLL+19", "HXZ+20", "BCA+21", "RMB+17", "SBP+20", "NTB+21", "MWY+21", "SAR+13", "CTS+17", "HDZ+19",
    "ZHZ+19", "YHK+20", "SRM+19", "HLZ+18", "PSG+20", "KSR+21", "ZCZ+21", "ZZL+21", "PBU+20", "MLL+18",
    "BWD+19", "WHL+21", "HSF+19", "ZML+17", "DTJ+20", "WFB+21", "KLM+16", "LMBA21", "BSPD17", "CWC+20",
    "BMOB18", "LJC+20", "DLS21",  "ZZLL21", "GWS+20", "CLS+19", "GZL21",  "LLC+19", "NLE+21", "APV+20",
    "WHLH21a", "PCL+18", "GSBT20", "SPA+21", "IBK+22", "HSF+22", "HCW+22", "WKP+22", "CKB+22", "WHLH21",
    "RSS+18", "PSB+21", "RKSK22", "WLJ+16", "DJK+18", "OHL+18", "ZLH+22", "LWJ+21", "MGW+22", "RKN+17",
    "ZDW+19", "WKG+22", "CLZ+22", "RARG22", "MCR+22", "MTC21"
  )) {
    # Get Eh or ORP values (uses partial name matching, can match a column named "Eh (mV)")
    Eh <- metadata$Eh
    Ehname <- "Eh"
    if(is.null(Eh)) {
      Eh <- metadata$ORP
      Ehname <- "ORP"
    }
    # Also match a column starting with "redox" (case-insensitive)
    iEh <- grep("^redox", colnames(metadata), ignore.case = TRUE)
    if(length(iEh) > 0) {
      Eh <- metadata[, iEh[1]]
      Ehname <- strsplit(colnames(metadata)[iEh[1]], " ")[[1]][1]
    }
    if(is.null(Eh)) stop("can't find Eh or ORP column")

    # Append Ehorig column (normalized column name used for exporting data) 20210821
    metadata <- cbind(metadata, Ehorig = Eh)
    
    # Get temperature for dEh/dpH calculation 20210829
    iT <- match("T", colnames(metadata))  # matches "T"
    if(is.na(iT)) iT <- grep("^T\\ ", colnames(metadata))[1]  # matches "T (°C)" but not e.g. "Treatment"
    if(is.na(iT)) iT <- grep("^Temp", colnames(metadata))[1]  # matches "Temperature (°C)"
    if(!is.na(iT)) {
      T <- metadata[, iT]
      Ttext <- paste(round(range(na.omit(T)), 1), collapse = " to ")
    } else {
      T <- rep(25, nrow(metadata))
      Ttext <- "assumed 25"
    }

    # Adjust Eh to pH 7 20210828
    # Find pH column
    ipH <- match("pH", colnames(metadata))
    if(!is.na(ipH)) {
      # pH values
      pH <- metadata[, ipH]
      pHtext <- paste(round(range(na.omit(pH)), 1), collapse = " to ")
      ## Eh(mV)-pH slope at 25 °C
      #dEhdpH <- -59.2
      # Find dEh/dpH as a function of T 20210829
      TK <- T + 273.15
      R <- 0.0083147
      F <- 96.4935
      dEhdpH <- - (log(10) * R * TK) / F * 1000
      # Difference to pH 7
      pHdiff <- 7 - pH
      # Adjust to pH 7
      Eh7 <- round(Eh + pHdiff * dEhdpH, 2)
    } else {
      Eh7 <- Eh
      pHtext <- "assumed 7"
    }
    metadata <- cbind(metadata, Eh7 = Eh7)
    # Print message about T, pH and Eh ranges
    Ehtext <- paste(round(range(na.omit(Eh))), collapse = " to ")
    Eh7text <- paste(round(range(na.omit(Eh7))), collapse = " to ")
    infotext <- paste0(study, ": T ", Ttext, ", pH ", pHtext, ", Eh ", Ehtext, ", Eh7 ", Eh7text)
    if(!quiet) print(infotext)

    # Get O2 values 20220517
    O2_umol_L <- NA
    iO2 <- grep("^O2", colnames(metadata))[1]
    if(is.na(iO2)) iO2 <- grep("^DO", colnames(metadata))[1]
    if(is.na(iO2)) iO2 <- grep("^oxygen", tolower(colnames(metadata)))[1]
    if(!is.na(iO2)) {
      # Get the units
      O2name <- colnames(metadata)[iO2]
      # Remove text up to and including left parenthesis, and right parenthesis
      units <- gsub(".*\\(|\\)", "", O2name)
      if(units %in% c("\u03BCmol/L", "\u03BCmol L-1", "umol/L", "umol L-1", "umol kg-1", "\u03BCM")) O2_umol_L <- metadata[, iO2]
      else if(units %in% c("mM", "mmol L-1")) O2_umol_L <- metadata[, iO2] * 1000
      else if(units %in% c("mg/L", "mg L-1")) O2_umol_L <- metadata[, iO2] * 1000 / 31.9988
      # This isn't needed now but leave it here for future reference 20220610
      # Source: USGS Office of Water Quality Technical Memorandum, 2011.03,
      #         Change to solubility equations for oxygen in water
      #         https://water.usgs.gov/water-resources/memos/
      #else if(units %in% c("mL/L", "mL L-1")) O2_umol_L <- metadata[, iO2] * 1.42905 * 1000 / 31.9988
      else if(units == "%") {
        # Calculate logK for O2(gas) = O2(aq)
        logK <- rep(NA, nrow(metadata))
        not_NA <- !is.na(T)
        logK[not_NA] <- subcrt(c("O2", "O2"), c("gas", "aq"), c(-1, 1), T = T[not_NA])$out$logK
        # Calculate saturated oxygen concentration from logK = logaO2(aq) - logfO2(gas)
        # Assume fugacity is partial pressure of oxygen in sea-level atmosphere
        logaO2 <- logK + log10(0.21 * 1.01325)
        # Remove logarithm and convert mol to umol
        O2_umol_L_sat <- 10^logaO2 * 1e6
        # Multiply by percent to get concentration
        O2_umol_L <- metadata[, iO2] / 100 * O2_umol_L_sat
      }
      else if(identical(units, O2name)) warning(paste0(study, ": no units given for ", O2name))
      else stop(paste0(study, ": unrecognized units for oxygen concentration: ", units))
      # Round value 20221011
      O2_umol_L <- round(O2_umol_L, 2)
    }
    metadata <- cbind(metadata, O2_umol_L = O2_umol_L)

    # Cluster samples by Eh7 value
    if(nrow(metadata) > 3) {
      # Keep NA values out of clusters
      ina <- is.na(Eh7)
      # Divide data into 3 clusters
      cl <- kmeans(Eh7[!ina], 3, nstart = 100)
      # Find the clusters with the lowest and highest means
      imin <- which.min(cl$centers[, 1])
      imax <- which.max(cl$centers[, 1])
      imid <- (1:3)[-c(imin, imax)]
      # Name the clusters
      cluster <- rep(NA, nrow(metadata))
      cluster[!ina][cl$cluster == imin] <- "reducing"
      cluster[!ina][cl$cluster == imid] <- "transition"
      cluster[!ina][cl$cluster == imax] <- "oxidizing"
      # Get pch and col for plot
      # Use is.null to let previous assignments take precedence (for SVH+19 in geo16S) 20211008
      if(is.null(pch)) {
        pch <- sapply(cluster, switch, oxidizing = 24, transition = 21, reducing = 25, NA)
        col <- sapply(cluster, switch, oxidizing = 4, transition = 1, reducing = 2, NA)
      }
      if(!quiet) {
        # Print message about Eh7 or ORP ranges of clusters
        ro <- range(na.omit(Eh7[cluster=="oxidizing"]))
        no <- sum(na.omit(cluster=="oxidizing"))
        message(paste0("getmdat_orp16S: oxidizing Eh7 ", paste(ro, collapse = " to "), " mV (n = ", no, ")"))
        rt <- range(na.omit(Eh7[cluster=="transition"]))
        nt <- sum(na.omit(cluster=="transition"))
        message(paste0("getmdat_orp16S: transition Eh7 ", paste(rt, collapse = " to "), " mV (n = ", nt, ")"))
        rr <- range(na.omit(Eh7[cluster=="reducing"]))
        nr <- sum(na.omit(cluster=="reducing"))
        message(paste0("getmdat_orp16S: reducing Eh7 ", paste(rr, collapse = " to "), " mV (n = ", nr, ")"))
      }
      # Append cluster names to data frame
      metadata <- cbind(metadata, cluster)
    } else {
      # For very small datasets, don't create clusters but still assign pch and col
      pch <- rep(21, length(Eh7))
      pch[which.max(Eh7)] <- 24
      pch[which.min(Eh7)] <- 25
      col <- rep(1, length(Eh7))
      col[which.max(Eh7)] <- 4
      col[which.min(Eh7)] <- 2
    }
  }

  if(is.null(pch)) stop(paste(study, "metadata file exists, but not set up for processing"))

  metadata <- cbind(metadata, pch, col)
  # Use the infotext as an attribute (was previously used by orp16S_info) 20220513
  attr(metadata, "infotext") <- infotext
  # Return both metadata and metrics, if provided 20220506
  if(is.null(metrics)) metadata else {
    # Keep metadata only for samples with metrics 20201006
    metadata <- metadata[metadata$Run %in% metrics$Run, ]
    # Keep specified number of random samples 20220509
    if(!is.null(size)) if(nrow(metadata) > size) {
      isample <- sample(1:nrow(metadata), size = size)
      metadata <- metadata[isample, ]
    }
    # Put metrics in same order as metadata 20220505
    imet <- match(metadata$Run, metrics$Run)
    metrics <- metrics[imet, ]
    # Insert sample column in metrics
    # Use first column name starting with "sample" or "Sample" 20210818
    sampcol <- grep("^sample", colnames(metadata), ignore.case = TRUE)[1]
    metrics <- data.frame(Run = metrics$Run, sample = metadata[, sampcol], nH2O = metrics$nH2O, Zc = metrics$Zc)
    list(metadata = metadata, metrics = metrics)
  }
}

# Function to calculate metrics for a given study 20220506
getmetrics_orp16S <- function(study, mincount = 100, quiet = TRUE, ...) {
  # Remove suffix after underscore 20200929
  studyfile <- gsub("_.*", "", study)
  datadir <- system.file("extdata/orp16S/RDP", package = "JMDplots")
  RDPfile <- file.path(datadir, paste0(studyfile, ".tab.xz"))
  # If there is no .xz file, look for a .tab file 20210607
  if(!file.exists(RDPfile)) RDPfile <- file.path(datadir, paste0(studyfile, ".tab"))
  RDP <- read_RDP(RDPfile, mincount = mincount, quiet = quiet, ...)
  map <- map_taxa(RDP, refdb = "RefSeq", quiet = quiet)
  get_metrics(RDP, refdb = "RefSeq", map)
}

############################
### UNEXPORTED FUNCTIONS ###
############################

# Calculate and plot linear regression and return coefficients 20210920
add.linear.global <- function(xvals, Zc, nstudy = NA, xvar = "Eh7", legend.x = "topright", inset = 0, scale = 1000) {
  text.col <- "black"
  line.col <- "gray62"
  # Use lighter colors for environments with few datasets 20210925
  if(!is.na(nstudy)) if(nstudy < 5) {
    text.col <- "gray60"
    line.col <- "gray70"
  }
  # Show number of studies and exit if there are zero 20210926
  if(!is.na(nstudy)) {
    legend("topleft", legend = nstudy, bty = "n", text.font = 2, text.col = text.col, inset = c(-0.05, 0))
    if(nstudy == 0) return()
  }
  # Fit data with linear model
  thislm <- lm(Zc ~ xvals)
  # Draw linear regression
  xlim <- range(xvals)
  plx <- predict(thislm, newdata = data.frame(xvals = xlim), se = T)
  lines(xlim / scale, plx$fit, col = line.col)

  # Add legend with number of samples
  Ntext <- bquote(italic(N) == .(length(Zc)))
  legend(legend.x, legend = Ntext, bty = "n", text.col = text.col, inset = inset)

  # Get the slope
  slope <- thislm$coefficients[2]
  # Get the MOE for 95% confidence interval 20221008
  CI95 <- confint(thislm, "xvals", level = 0.95)
  stopifnot(all.equal(mean(CI95), as.numeric(slope)))
  MOE95 <- diff(range(CI95)) / 2
  # Multiply by 1e3 to convert from mV-1 to V-1 or umol-1 to mmol-1
  slope <- slope * 1e3
  MOE95 <- MOE95 * 1e3
  # Round to fixed number of decimal places
  slopetxt <- formatC(slope, digits = 3, format = "f")
  MOE95txt <- formatC(MOE95, digits = 3, format = "f")
  # Units for Eh7/Eh/O2
  if(xvar %in% c("Eh7", "Eh")) slopeleg <- bquote(italic(m) == .(slopetxt) %+-% .(MOE95txt) / V)
  if(xvar == "O2") slopeleg <- bquote(italic(m) == .(slopetxt) %+-% .(MOE95txt) ~ liter / mmol)
  # Add text to plot
  legend <- as.expression(c("", slopeleg))
  legend(legend.x, legend = legend, bty = "n", text.col = text.col, inset = inset)

  # Calculate Pearson correlation 20211009
  pearson <- cor.test(xvals, Zc, method = "pearson")
  # Format correlation coefficient
  rtext <- formatC(pearson$estimate, digits = 2, format = "f")
  rtext <- bquote(italic(r) == .(rtext))
  legend("bottomright", legend = rtext, bty = "n", text.col = text.col)
  # Return slope 20220611
  slope
}

add.linear.local <- function(Eh7, Zc, col = "gray62", lwd = 1, legend = NULL, cex = 0.85, inset = 0, with.N = TRUE, scale = 1000, nounits = FALSE) {

  # Get slope and MOE95 (95% CI) for linear regression
  EZdat <- data.frame(Eh7, Zc)
  EZlm <- lm(Zc ~ Eh7, EZdat)
  slope <- EZlm$coefficients[2]
  CI95 <- confint(EZlm, "Eh7", level = 0.95)
  stopifnot(all.equal(mean(CI95), as.numeric(slope)))
  MOE95 <- diff(range(CI95)) / 2
  # Convert from mV-1 to V-1
  slope <- slope * scale
  MOE95 <- MOE95 * scale
  # Use solid or dashed line to indicate large or small slope 20210926
  if(is.na(slope)) lty <- 3 else if(abs(slope) < 0.01) lty <- 2 else lty <- 1
  # Plot regression line
  Eh7lim <- range(EZlm$model$Eh7)
  Zcpred <- predict.lm(EZlm, data.frame(Eh7 = Eh7lim))
  lines(Eh7lim / scale, Zcpred, col = col, lwd = lwd, lty = lty)
  # Calculate Pearson correlation 20220520
  pearson <- cor.test(EZdat$Eh7, EZdat$Zc, method = "pearson")

  if(!is.null(legend)) {

    # Format correlation coefficient 20221001
    rtxt <- formatC(pearson$estimate, digits = 2, format = "f")
    # Format slope and MOE95
    slopetxt <- formatC(slope, digits = 3, format = "f")
    MOE95txt <- formatC(MOE95, digits = 3, format = "f")

    if(legend == "title") {
      # Put statistics in subtitle 20221004
      main <- bquote(list(
        italic(N) == .(nrow(EZdat)),
        italic(r) == .(rtxt),
        italic(m) == .(slopetxt) %+-% .(MOE95txt) / V
      ))
      title(main = main, line = 0.75, cex.main = 1)
    } else {
      # Put statistics in plot 20221001
      ntext <- bquote(italic(N) == .(nrow(EZdat)))
      rtext <- bquote(italic(r) == .(rtxt))
      stext <- bquote(italic(m) == .(slopetxt))
      MOEtext <- bquote(phantom(xx) %+-% .(MOE95txt) / V)
      if(nounits) MOEtext <- bquote(phantom(xx) %+-% .(MOE95txt))
      if(with.N) ltext <- c(ntext, rtext, stext, MOEtext)
      if(!with.N) ltext <- c(rtext, stext, MOEtext)
      legend(legend, legend = ltext, bty = "n", cex = cex, inset = inset)
    }
  }
  # Return results for plotEZ() 20221004
  invisible(list(EZlm = EZlm, Eh7lim = Eh7lim, Zcpred = Zcpred, pearson = pearson, MOE95 = MOE95))
}

# Scatterplots for all samples in each environment type 20210913
eachenv <- function(eedat, add = FALSE, do.linear = TRUE, ienv = c(1, 2, 4, 5, 3, 6, 7), cols = orp16Scol,
  lineage = NULL, xlim = NULL, ylim = NULL, xvar = "Eh7") {
  # Get x values and scale (for mV to V conversion)
  scale <- 1000
  if(xvar == "Eh7") xvals <- eedat$Eh7
  if(xvar == "Eh") xvals <- eedat$Eh
  if(xvar == "O2") {
    xvals <- eedat$O2_umol_L
    scale <- 1
  }
  eedat <- cbind(eedat, xvals)
  # Get overall x and y limits
  if(is.null(xlim)) xlim <- range(eedat$xvals)
  if(is.null(ylim)) ylim <- range(eedat$Zc)
  # Get names of environment types
  envirotypes <- names(envirotype)
  slopes <- rep(NA, length(ienv))
  # Loop over environment types
  for(i in ienv) {
    # Start plot
    if(!add) plot(xlim / scale, ylim, type = "n", xlab = "", ylab = "", axes = FALSE)
    # Get Eh7/O2 and Zc values
    thisdat <- eedat[eedat$envirotype == envirotypes[i], ]
    xvals <- thisdat$xvals
    Zc <- thisdat$Zc
    if(do.linear) {
      # Include number of studies in legend 20210925
      nstudy <- length(unique(thisdat$study))
      # Adjust legend position for Archaea in Lake & Pond (high Zc in hypersaline lakes) 20221021
      inset <- c(0, 0)
      if(lineage == "Archaea" & envirotypes[i] == "Lake & Pond") inset <- c(0, 0.05)
      if(nrow(thisdat) > 0) slopes[[i]] <- add.linear.global(xvals, Zc, nstudy, inset = inset)
    }
    if(!isTRUE(add)) {
      # Add plot axes and title
      if(lineage == "Archaea") {
        axis(1, labels = NA)
        # Make rotated labels (modified from https://www.r-bloggers.com/rotated-axis-labels-in-r-plots/)
        x <- c(-0.40, 0, 0.60)
        text(x = x, y = par()$usr[3] - 1.5 * strheight("A"), labels = x, srt = 45, adj = 1, xpd = NA)
      }
      if(i == 1) axis(2)
      box()
      if(lineage == "Bacteria") title(envirotypes[i], font.main = 1, cex.main = 1)
    }
    # Plot points
    points(xvals / scale, Zc, pch = 19, cex = 0.2, col = cols[i])
  }
  # Return slopes 20220611
  if(do.linear) {
    names(slopes) <- envirotypes
    slopes
  }
}

# Gather values of Zc and nH2O for metaproteomes and 16S-based estimates 20220828
getMP_orp16S <- function() {

  # List sources of metaproteomic and 16S data
  studies_MP <- c("MPB+19", "PMM+18", "KTS+17", "KTS+17.mock", "HTZ+17")
  studies_16S <- c("MPB+19", "RYP+14", "KTS+17", "KTS+17.mock", "HTZ+17")
  n <- length(studies_16S)
  pchs <- c(24, 25, 12, 18, 20)
  bgs <- c(adjustcolor(c(5, 6), alpha.f = 0.815), NA, NA, NA)
  cols <- c(1, 1, 2, 3, "#00000080")

  out <- lapply(seq_along(studies_16S), function(i) {
    # Get 16S values
    metrics <- getmetrics_orp16S(studies_16S[i])
    mdat <- getmdat_orp16S(studies_16S[i], metrics)
    # Get metaproteome values
    datadir <- system.file("extdata/orp16S/metaproteome", package = "JMDplots")
    aa <- read.csv(paste0(datadir, "/", studies_MP[i], "/", studies_MP[i], "_aa.csv"))
    # Match 16S samples to metaproteome
    iaa <- match(mdat$metadata$Metaproteome, aa$organism)
    # Drop non-matching samples
    ina <- is.na(iaa)
    aa <- aa[iaa[!ina], ]
    metadata <- mdat$metadata[!ina, , drop = FALSE]
    metrics <- mdat$metrics[!ina, , drop = FALSE]
    stopifnot(all(metadata$Metaproteome == aa$organism))
    # Make data frame
    data.frame(Name = metadata$Name,
      Study_16S = studies_16S[i], Run_16S = metadata$Run, Sample_16S = metadata$Sample,
      Study_MP = studies_MP[i], Sample_MP = aa$organism,
      Zc_16S = round(metrics$Zc, 6), Zc_MP = round(Zc(aa), 6),
      nH2O_16S = round(metrics$nH2O, 6), nH2O_MP = round(nH2O(aa), 6),
      #pch = metadata$pch, col = metadata$col
      pch = pchs[i], bg = bgs[i], col = cols[i]
    )
  })

  # Return value instead of saving to file 20221018
  do.call(rbind, out)

}
