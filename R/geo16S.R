# JMDplots/geo16S.R
# Make plots for the paper:
# Geobiochemistry of microbial community proteomes inferred from 16S RNA sequences
# 20210416 Initial commit to JMDplots
# 20210527 Updated plots for RefSeq release 206

# Figure 1: Chemical compositions of taxonomic groups 20200925
geo16S1 <- function(pdf = FALSE) {

  if(pdf) pdf("geo16S1.pdf", width = 11, height = 5)
  par(mfrow = c(1, 3))

  taxacomp("majorphyla", legend.x = "bottomleft", hline = c(-0.81, -0.68))
  title("Major phyla and their classes", font.main = 1, cex.main = 1.4)
  label.figure("A", font = 2, cex = 1.6)
  # Draw lines indicating zoom area in next plot
  par(xpd = NA)
  lines(c(-0.05, -0.015), c(-0.81, -0.9), lty = 2, col = "gray40")
  lines(c(-0.05, -0.015), c(-0.68, -0.65), lty = 2, col = "gray40")
  par(xpd = FALSE)
  taxacomp("majorcellular", legend.x = "bottomleft", hline = c(-0.76, -0.71))
  title("Major cellular phyla and their classes", font.main = 1, cex.main = 1.4)
  label.figure("B", font = 2, cex = 1.6)
  par(xpd = NA)
  lines(c(-0.05, -0.015), c(-0.76, -0.81), lty = 2, col = "gray40")
  lines(c(-0.05, -0.015), c(-0.71, -0.68), lty = 2, col = "gray40")
  par(xpd = FALSE)
  taxacomp("Proteobacteria", legend.x = "topright")
  title("Proteobacterial classes and their orders", font.main = 1, cex.main = 1.4)
  label.figure("C", font = 2, cex = 1.6)

  if(pdf) dev.off()

}

# Figure 2: Natural environment datasets 20200923
geo16S2 <- function(pdf = FALSE) {

  if(pdf) pdf("geo16S2.pdf", width = 9, height = 7)
  oopar <- par(no.readonly = TRUE)
  par(mar = c(4, 4, 3, 2))
  par(mgp = c(2.5, 1, 0))
  layout(matrix(c(1,2,3, 4,9,5, 6,7,8), nrow = 3, byrow = TRUE))

  # Function to add points
  pointfun <- function(pcomp) {
    ifill <- pcomp$pch > 20
    points(pcomp$ZC[ifill], pcomp$nH2O[ifill], pch = pcomp$pch[ifill], col = 1, bg = pcomp$col[ifill])
    points(pcomp$ZC[!ifill], pcomp$nH2O[!ifill], pch = pcomp$pch[!ifill], col = pcomp$col[!ifill])
  }

  p1 <- plotcomp("BGPF13", title = FALSE, points = FALSE)
#  title("Yellowstone hot springs\nBowen De Le\u00f3n et al., 2013", font.main = 1)
  title("Yellowstone hot springs", font.main = 1)
  addhull(p1$ZC, p1$nH2O, 2, TRUE)
  pointfun(p1)
  legend <- c("Archaea", "Bacteria")
  legend("bottomleft", legend, pch = c(23, 22), col = c(1, 1), pt.bg = c(6, 5), bg = "white")

  p2 <- plotcomp("SVH+19", title = FALSE, points = FALSE)
#  title("Black Sea\nSollai et al., 2019", font.main = 1)
  title("Black Sea", font.main = 1)
  addhull(p2$ZC, p2$nH2O, "blue", TRUE, lty = 2)
  pointfun(p2)
  legend <- c("Oxic", "Suboxic", "Euxinic")
  legend("bottomright", legend, pch = c(24, 20, 25), pt.bg = c(4, 1, 2), bg = "white")

  p3 <- plotcomp("HLA+16", title = FALSE, points = FALSE)
#  title("Baltic Sea\nHerlemann et al., 2016", font.main = 1)
  title("Baltic Sea", font.main = 1)
  addhull(p3$ZC, p3$nH2O, "blue", TRUE)
  pointfun(p3)
  legend <- c("< 6 PSU", "6-20 PSU", "> 20 PSU")
  legend("bottomright", legend, pch = c(24, 20, 21), col = c(1, 1, 1), pt.bg = c(3, NA, 4), bg = "white")

  p4 <- plotcomp("MPB+17", title = FALSE, points = FALSE)
#  title("Manus Basin submarine vents\nMeier et al., 2017", font.main = 1)
  title("Manus Basin submarine vents", font.main = 1)
  addhull(p4$ZC, p4$nH2O, 2, TRUE, lty = 2)
  pointfun(p4)
  legend <- as.expression(c(quote(italic(T)~"< 50 \u00B0C"), quote(italic(T)~"> 50 \u00B0C")))
  legend("bottomleft", legend, pch = c(21, 23), col = c(1, 1), pt.bg = c(4, 2), bg = "white", title = "Water")
  legend("bottomright", c("Rock", "Fauna"), pch = c(20, 8), col = c(1, "yellow4"), bg = "white")

  p5 <- plotcomp("ZLM+16", title = FALSE, points = FALSE)
#  title("Tibetan Plateau lakes\nZhong et al., 2016", font.main = 1)
  title("Tibetan Plateau lakes", font.main = 1)
  addhull(p5$ZC, p5$nH2O, "turquoise3", TRUE, lty = 2)
  pointfun(p5)
  legend <- c("< 10 g/l", "24-99 g/l", "> 300 g/l")
  legend("topright", legend, pch = c(24, 20, 21), col = c(1, 1, 1), pt.bg = c(3, NA, 4), bg = "white")

  p6 <- plotcomp("JHM+16", title = FALSE, points = FALSE)
#  title("Lake Fryxell oxygen gradient\nJungblut et al., 2016", font.main = 1)
  title("Lake Fryxell oxygen gradient", font.main = 1)
  addhull(p6$ZC, p6$nH2O, "tan1", TRUE)
  pointfun(p6)
  legend <- c("Oxic", "Transition", "Anoxic")
  legend("bottomright", legend, pch = c(24, 20, 25), pt.bg = c(4, 1, 2), bg = "white")

  p7 <- plotcomp("HCW+13", title = FALSE, points = FALSE, ylim = c(-0.765, -0.7545))
#  title("Guerrero Negro mat layers\nHarris et al., 2013", font.main = 1)
  title("Guerrero Negro mat layers", font.main = 1)
  addhull(p7$ZC, p7$nH2O, "tan1", TRUE, lty = 2)
  pointfun(p7)
  text(c(-0.1516, -0.1571, -0.1576), c(-0.7551, -0.7608, -0.7649), c("0-1 mm", "1-2 mm", "2-3 mm"))
  legend <- c("Photic/oxic", "Low sulfide", "High sulfide")
  legend("topleft", legend, pch = c(24, 20, 25), pt.bg = c(4, 1, 2), bg = "white")

  p8 <- plotcomp("XDZ+17", title = FALSE, points = FALSE)
#  title("Qarhan Salt Lake and\nnormal soils, Xie et al., 2017", font.main = 1)
  title("Qarhan Salt Lake\nand normal soils", font.main = 1)
  addhull(p8$ZC, p8$nH2O, "turquoise3", TRUE)
  pointfun(p8)
  legend <- c("Normal", "Saline")
  legend("topright", legend, pch = c(24, 21), pt.bg = c(3, 4), bg = "white")

  # Make an index plot
  opar <- par(mar = c(2.5, 2.5, 0.5, 0.5))
  xlim <- c(-0.22, -0.09)
  ylim <- c(-0.77, -0.71)
  if(options("basis")$basis == "QCa") {
    xlim <- c(-0.22, -0.09)
    ylim <- c(-1.12, -1)
  }
  plot(xlim, ylim, xlab = "", ylab = "", type = "n")
  lmlines()
#  # Add convex hull for stratified lakes and sewater (from Fig. 3) 20210503
#  fig3 <- geo16S3(plot.it = FALSE)
#  addhull(fig3$mar$ZC, fig3$mar$nH2O, "blue")
#  addhull(fig3$ursu$ZC, fig3$ursu$nH2O, "turquoise3")
  # Add convex hulls for each dataset in this figure
  addhull(p1$ZC, p1$nH2O, 2, TRUE)
  addhull(p2$ZC, p2$nH2O, "blue", TRUE, lty = 2)
  addhull(p3$ZC, p3$nH2O, "blue", TRUE)
  addhull(p4$ZC, p4$nH2O, 2, TRUE, lty = 2)
  addhull(p5$ZC, p5$nH2O, "turquoise3", TRUE, lty = 2)
  addhull(p6$ZC, p6$nH2O, "tan1", TRUE)
  addhull(p7$ZC, p7$nH2O, "tan1", TRUE, lty = 2)
  addhull(p8$ZC, p8$nH2O, "turquoise3", TRUE)
  par(opar)

  # Add environment type labels 20210427
  par(xpd = NA)
  text(-0.397, -0.703, "Hydrothermal", cex = 1.5, font = 2, srt = 90)
  text(-0.24, -0.857, "Microbial Mats", cex = 1.5, font = 2)
  text(-0.075, -0.635, "Seawater", cex = 1.5, font = 2)
  text(0.074, -0.785, "Hypersaline", cex = 1.5, font = 2, srt = 90)
  par(xpd = FALSE)

  par(oopar)
  if(pdf) dev.off()

}

# Figure 3: Stratified lakes and seawater 20210428
geo16S3 <- function(pdf = FALSE, plot.it = TRUE) {

  if(plot.it) {
    if(pdf) pdf("geo16S3.pdf", width = 7, height = 9)
    mat <- matrix(c(1,1,1,1, 2,3,4,5, 6,7,8,9, 10,11,12,13), byrow = TRUE, nrow = 4)
    layout(mat, heights = c(1, 10, 10, 10))
    # Make legend
    par(mar = c(0, 0, 0, 0))
    plot.new()
    legend <- as.expression(c(quote(italic(Z)[C]), quote(O[2]), quote(NO[3]^"-" / NO[2]^"-")))
    legend("top", legend = legend, lty = c(1, 1, 2), lwd = 1.5, col = c(1, 2, 4), pch = c(21, NA, NA), pt.bg = "white", ncol = 3, bty = "n")
    # Setup plot parameters
    par(mgp = c(1.8, 0.5, 0), mar = c(3, 3, 3, 1))
  }

  # Identify datasets to plot
  study <- c("SVH+19", "MZG+20", "MZG+20", "GBL+15", "GBL+15",
             "HXZ+20", "HXZ+20", "HXZ+20", "HXZ+20",
             "BCA+20", "BCA+20", "BCA+20", "BCA+20")
  column <- c("study", "lake", "lake", "size", "size",
              "station", "station", "station", "station",
              "month", "month", "month", "month")
  ID <- c("SVH+19", "Lake Zug", "Lake Lugano", "0.2-1.6micron", "1.6-30micron",
          "SYBL", "SYBL", "C4", "C4",
          "Jul", "Nov", "Feb", "Apr")
  title <- c("Black Sea", "Lake Zug", "Lake\nLugano", "ETNP", "ETNP",
             # Use leading or trailing space to flag NO3-/NO2- plots
             "Blue Hole\n", " Blue Hole\n", "Blue Hole\n", "Blue Hole\n ",
             "Ursu Lake\n", "Ursu Lake\n", "Ursu Lake\n", "Ursu Lake\n")
  subtitle <- c("", "", "", "", "",
                "Inside", "Inside", "Outside (C4)", "Outside (C4)",
                "July 2015", "November 2015", "February 2016", "April 2016")
  titlesub <- paste(title, subtitle)
  # Make objects to hold all ZC and nH2O values (for convex hull in Figure 2)
  marZC <- marnH2O <- numeric()
  ursuZC <- ursunH2O <- numeric()
  for(i in 1:length(study)) {
    # ZC range for plots
    if(study[i] == "BCA+20") ZClim <- c(-0.175, -0.145) else ZClim <- c(-0.17, -0.14)
    # Get the metadata and compositional metrics for this study
    # Keep all rows for higher-resolution O2 measurements
    mdat <- getmdat(study[i], dropNA = FALSE)
    metrics <- getmetrics(study[i])
    # Use depths < 500 m (excludes 500-2000 m Black Sea to better visualize shallower trends)
    mdat <- mdat[mdat$depth < 500, ]
    # Get the rows matching the ID
    iID <- mdat[, column[i]] == ID[i]
    mdat <- mdat[iID, ]
    # Sort the data by depth
    alldat <- mdat <- mdat[order(mdat$depth), ]
    # Now exclude NA samples
    mdat <- mdat[!is.na(mdat$name), ]
    depth <- mdat$depth
    # Get the ZC and nH2O values
    imet <- match(mdat$Run, metrics$Run)
    ZC <- metrics$ZC[imet]
    nH2O <- metrics$nH2O[imet]

    if(plot.it) {
      # Reverse y-axis (depth)
      if(study[i] == "BCA+20") ylim <- c(11, 0) else ylim <- rev(range(depth))
      # Determine whether the title has changed
      newplot <- TRUE
      if(i > 1) if(titlesub[i]==titlesub[i-1]) newplot <- FALSE
      if(newplot) {
        plot(ZC, depth, xlim = ZClim, ylim = ylim, xlab = axis.label("ZC"), ylab = "Depth (m)", type = "b")
      } else {
        # Add to plot if the title hasn't changed
        points(ZC, depth, type = "b", pch = 0)
      }

      if(newplot) {
        # Add title in lower right
        if(grepl("Outside", subtitle[i])) {
          text(ZClim[1], ylim[1], title[i], adj = c(0, 0), font = 2)
          text(ZClim[1], ylim[1], subtitle[i], adj = c(0, 0))
        } else {
          text(ZClim[2], ylim[1], title[i], adj = c(1, 0), font = 2)
          text(ZClim[2], ylim[1], subtitle[i], adj = c(1, 0))
        }
        # Plot O2 concentrations or NO3-/NO2- ratio
        nc <- nchar(title[i])
        if(substr(title[i], 1, 1) == " " | substr(title[i], nc, nc) == " ") {
          what <- "NO3.NO2"
          xlim <- c(0, 250)
          col <- 4
          lty <- 2
          # Calculate NO3- / NO2- ratio 20210511
          NO3.NO2 <- alldat$`NO3- (umol L-1)` / alldat$`NO2- (umol L-1)`
          alldat <- cbind(alldat, NO3.NO2)
        } else {
          what <- "O2"
          if(study[i] == "BCA+20") xlim <- c(0, 25) else xlim <- c(0, 220)
          col <- 2
          lty <- 1
        }
        icol <- grep(paste0("^", what), colnames(alldat))
        # Remove NA values
        alldat <- alldat[!is.na(alldat[, icol]), ]
        depth <- alldat$depth
        par(new = TRUE)
        plot(alldat[, icol], depth, col = col, lty = lty, type = "l", axes = FALSE, xlab = "", ylab = "", xlim = xlim, ylim = ylim)
        # Add second axis labels
        xlab <- colnames(alldat)[icol]
        if(xlab == "O2 (umol kg-1)") xlab <- quote(O[2]~"(\u00B5mol kg"^-1*")")
        if(xlab == "O2 (umol L-1)") xlab <- quote(O[2]~"(\u00B5mol L"^-1*")")
        if(xlab == "O2 (mg L-1)") xlab <- quote(O[2]~"(mg L"^-1*")")
        if(xlab == "NO3.NO2") xlab <- quote(NO[3]^"-" / NO[2]^"-"~"(mol/mol)")
        axis(3)
        mtext(xlab, side = 3, line = 1.7, cex = par("cex"))
        # Extra labels for ETNP
        if(title[i]=="ETNP") {
          text(40, 78, "0.2-\n1.6 \u00B5m")
          text(140, 135, "1.6-\n30 \u00B5m")
        }
        # Restore xlim for plotting ZC
        par(new = TRUE)
        plot(0, 0, type = "n", axes = FALSE, xlab = "", ylab = "", xlim = ZClim, ylim = ylim)
      }
    } # end if(plot.it)

    # Store all non-NA ZC and nH2O values
    if(grepl("Ursu", title[i])) {
      ursuZC <- c(ursuZC, na.omit(ZC))
      ursunH2O <- c(ursunH2O, na.omit(nH2O))
    } else {
      marZC <- c(marZC, na.omit(ZC))
      marnH2O <- c(marnH2O, na.omit(nH2O))
    }
  }

  if(plot.it) {
    if(pdf) dev.off()
  }
  # Return values for marine and freshwater datasets and Ursu Lake 20210521
  list(mar = list(ZC = marZC, nH2O = marnH2O), ursu = list(ZC = ursuZC, nH2O = ursunH2O))

}

# Figure 4: Compositional differences at different taxonomic levels 20200924
geo16S4 <- function(pdf = FALSE) {

  if(pdf) pdf("geo16S4.pdf", width = 10.5, height = 7)
  mat <- matrix(c(1,2,3,4, 5,6,7,8), nrow = 2, byrow = TRUE)
  layout(mat, widths = c(1, 2, 2, 2))
  par(mar = c(4, 4, 3, 1))
  par(mgp = c(2.5, 1, 0))
  par(cex.lab = 1.1)

  # Preload data for faster running
  mdat <- getmdat("MPB+17")
  RDP <- getRDP("MPB+17", mdat = mdat)
  map <- getmap("MPB+17", RDP = RDP)

  # Make plots for Manus Basin
  p <- groupcomp("MPB+17", "ZC", "domain", pch.up = 23, pch.down = 21, ylim = c(-0.23, -0.13), xlim = c(0, 100),
    xadj = c(Bacteria = 1), yadj = c(Bacteria = -5),
    mdat = mdat, RDP = RDP, map = map
  )
  title("Manus Basin")
  text(40, -0.157, "   < 50 \u00B0C", font = 2)
  text(10, -0.157, "T", font = 4)
  text(40, -0.202, "   > 50 \u00B0C", font = 2)
  text(10, -0.202, "T", font = 4)
  p <- groupcomp("MPB+17", "ZC", "phylum", pch.up = 23, pch.down = 21, ylim = c(-0.23, -0.13),
    xadj = c(Proteobacteria = 1, Bacteroidetes = -0.18, Campilobacterota = 0.4),
    yadj = c(Bacteroidetes = 1, Campilobacterota = -3),
    mdat = mdat, RDP = RDP, map = map
  )
  title(paste0("Phylum (", round(p), "% of total)"))
  p <- groupcomp("MPB+17", "ZC", "class", pch.up = 23, pch.down = 21, ylim = c(-0.23, -0.13),
    xadj = c(Flavobacteriia = -0.17, Gammaproteobacteria = 0.2, Campylobacteria = 0.4),
    yadj = c(Flavobacteriia = 1.2, Gammaproteobacteria = 1.7, Campylobacteria = -3),
    mdat = mdat, RDP = RDP, map = map
  )
  title(paste0("Class (", round(p), "% of total)"))
  p <- groupcomp("MPB+17", "ZC", "genus", pch.up = 23, pch.down = 21, ylim = c(-0.23, -0.13), minpercent = 1,
    xadj = c(Alteromonas = 0.1, Sulfurimonas = 1.05, Alcanivorax = -1, Halomonas = -0.65, Thiogranum = 0.17, Sulfurovum = -0.1, Pseudomonas = 0.1, Pseudoalteromonas = -0.25, Acinetobacter = 0.7),
    yadj = c(Sulfurimonas = 1.5, Pseudomonas = 1.8, Sulfurovum = -0.5, Thiogranum = 5, Marinimicrobia_genera_incertae_sedis = -0.8, Alteromonas = 1.8, Acinetobacter = 1.5),
    mdat = mdat, RDP = RDP, map = map
  )
  lines(c(0, 0), c(-0.1525, -0.1425))
  title(paste0("Genus (", round(p), "% of total)"))

  # Make plots for Baltic Sea
  mdat <- getmdat("HLA+16")
  RDP <- getRDP("HLA+16", mdat = mdat)
  map <- getmap("HLA+16", RDP = RDP)
  p <- groupcomp("HLA+16", "nH2O", "domain", pch.up = 24, pch.down = 21, ylim = c(-0.78, -0.71), xlim = c(0, 100),
    xadj = c(Bacteria = 1), yadj = c(Bacteria = 1.5),
    mdat = mdat, RDP = RDP, map = map
  )
  title("Baltic Sea")
  text(40, -0.741, "PSU < 6", font = 2)
  text(40, -0.753, "PSU > 20", font = 2)
  p <- groupcomp("HLA+16", "nH2O", "phylum", pch.up = 24, pch.down = 21, ylim = c(-0.78, -0.71),
    xadj = c(Proteobacteria = -0.2, Planctomycetes = 0.1, "Cyanobacteria/Chloroplast" = 0.35),
    yadj = c(Planctomycetes = 1.8, "Cyanobacteria/Chloroplast" = 2.5),
    mdat = mdat, RDP = RDP, map = map
  )
  title(paste0("Phylum (", round(p), "% of total)"))
  p <- groupcomp("HLA+16", "nH2O", "class", pch.up = 24, pch.down = 21, ylim = c(-0.78, -0.71),
    xadj = c(Acidimicrobiia = 0.5, Gammaproteobacteria = 0.4, Flavobacteriia = -0.1, Verrucomicrobiae = -0.25, Betaproteobacteria = 0.25),
    yadj = c(Acidimicrobiia = -1, Alphaproteobacteria = -0.6, Gammaproteobacteria = 1.6, Verrucomicrobiae = 1.5, Flavobacteriia = 1, Betaproteobacteria = 2.5),
    mdat = mdat, RDP = RDP, map = map
  )
  title(paste0("Class (", round(p), "% of total)"))
  p <- groupcomp("HLA+16", "nH2O", "genus", pch.up = 24, pch.down = 21, ylim = c(-0.78, -0.71), minpercent = 1,
    xadj = c(Spartobacteria_genera_incertae_sedis = 1.05, `Candidatus Pelagibacter` = -0.05),
    yadj = c(Spartobacteria_genera_incertae_sedis = -0.5),
    mdat = mdat, RDP = RDP, map = map
  )
  title(paste0("Genus (", round(p), "% of total)"))

#  # Between-study comparisons
#  groupcomp("MPB+17", "ZC", ylim = c(-0.02, 0.02), study2 = "HLA+16")
#  mtext("Manus Basin - Baltic Sea", line = 1)
#
#  groupcomp("BGPF13", "ZC", ylim = c(-0.12, 0.02), study2 = "XDZ+17")
#  mtext("Yellowstone - Qarhan Salt Lake", line = 1)

  if(pdf) dev.off()

}

# Figure 5: Shale gas datasets 20210414
geo16S5 <- function(pdf = FALSE) {

  if(pdf) pdf("geo16S5.pdf", width = 9, height = 6)
  par(mfrow = c(2, 2))
  par(mar = c(4, 4, 1, 1))
  par(mgp = c(2.5, 1, 0))

  ## Plot A: Pennsylvania streams affected by Marcellus Shale activity 20210324

  # Data from Ulrich et al., 2018
  xlim <- c(-0.16, -0.13)
  ylim <- c(-0.755, -0.725)
  plotcomp("UKD+18.water_2014", xlim = xlim, ylim = ylim, title = FALSE, pval = FALSE)
  legend("topleft", c("MSA+", "MSA-"), pch = c(21, 1), pt.bg = c(2, 1), bg = "white", title = "Northwestern Pennsylvania")
  label.figure("A", cex = 1.5, xfrac = 0.03, font = 2)

  ## Plot B: Comparison of different studies on Pennsylvania Streams 20210327

  studies <- c("UKD+18.water_2014", "UKD+18.sediment_2014", "CUN+18", "MMA+20_spring", "MMA+20_fall")
  # Start plot
  plot(c(-0.148, -0.138), c(-0.744, -0.732), type = "n", xlab = cplab$ZC, ylab = cplab$nH2O)
  pch <- 21:25
  # Loop over studies
  for(i in 1:5) {
    mean <- plotcomp(studies[[i]], plot.it = FALSE)$mean
    points(mean$ZC.dn, mean$nH2O.dn, pch = pch[i], cex = 1.5, lwd = 2, bg = "#ffffffa0")
    lines(c(mean$ZC.dn, mean$ZC.up), c(mean$nH2O.dn, mean$nH2O.up))
    points(mean$ZC.up, mean$nH2O.up, pch = pch[i], cex = 1.8, lwd = 2, bg = "#df536ba0")
  }
  # Add labels
  text(-0.1443, -0.736, "NW PA\nwater")
  text(-0.1466, -0.7422, "NW PA\nsediment")
  text(-0.1390, -0.7379, "NE PE\nsediment")
  text(-0.1410, -0.7417, "PASF water (spring)")
  text(-0.1420, -0.7431, "PASF water (fall)")
  # Add legend
  legend("topleft", c("Lowest disturbance", "Highest disturbance"), pch = c(21, 21), pt.bg = c("#ffffffa0", "#df536ba0"), pt.cex = c(1.4, 1.7), lwd = 2, lty = NA)
  label.figure("B", cex = 1.5, xfrac = 0.03, font = 2)

  ## Plots C-D: Comparison of different studies on produced water 20210330

  # Panel C: Cluff et al., 2014
  plotcomp("CHM+14", title = FALSE, pval = FALSE)
  legend("topright", c("Injected fluids (day 0)", "Produced water (day 49 and after)"),
         pch = c(21, 21), pt.bg = c("white", 2), bg = "white", title = "Marcellus Shale")
  box()
  label.figure("C", cex = 1.5, xfrac = 0.03, font = 2)

  # Panel D: Multiple studies
  studies <- c("CHM+14", "HRR+18", "ZLF+19")
  # Start plot
  plot(c(-0.22, -0.14), c(-0.75, -0.71), type = "n", xlab = cplab$ZC, ylab = cplab$nH2O)
  pch <- 21:25
  # Loop over studies
  for(i in 1:3) {
    mean <- plotcomp(studies[[i]], plot.it = FALSE)$mean
    points(mean$ZC.dn, mean$nH2O.dn, pch = pch[i], cex = 1.5, lwd = 2, bg = "#ffffffa0")
    lines(c(mean$ZC.dn, mean$ZC.up), c(mean$nH2O.dn, mean$nH2O.up))
    points(mean$ZC.up, mean$nH2O.up, pch = pch[i], cex = 1.8, lwd = 2, bg = "#df536ba0")
  }
  # Add labels
  text(-0.165, -0.726, "Marcellus Shale")
  text(-0.204, -0.732, "Denver-Julesburg Basin")
  text(-0.191, -0.719, "Duvernay Formation")
  # Add legend
  legend("topright", c("Source water\nor injected fluids", "Produced water"), pch = c(21, 21), pt.bg = c("#ffffffa0", "#df536ba0"), pt.cex = c(1.4, 1.7), lwd = 2, lty = NA)
  label.figure("D", cex = 1.5, xfrac = 0.03, font = 2)

  if(pdf) dev.off()

}

