# JMDplots/geo16S.R
# Make plots for the paper:
# Geobiochemistry of microbial community proteomes inferred from 16S RNA sequences
# 20200416 added to JMDplots

# Figure 1: Chemical compositions of taxonomic groups 20200925
geo16S1 <- function(pdf = FALSE) {

  if(pdf) pdf("geo16S1.pdf", width = 11, height = 5)
  par(mfrow = c(1, 3))

  taxacomp("majorphyla", legend.x = "bottomleft", hline = c(-0.81, -0.68))
  title("Major phyla and their classes", font.main = 1)
  # Draw lines indicating zoom area in next plot
  par(xpd = NA)
  lines(c(-0.05, -0.015), c(-0.81, -0.9), lty = 2, col = "gray40")
  lines(c(-0.05, -0.015), c(-0.68, -0.65), lty = 2, col = "gray40")
  par(xpd = FALSE)
  taxacomp("majorcellular", legend.x = "bottomleft", hline = c(-0.76, -0.71))
  title("Major cellular phyla and their classes", font.main = 1)
  par(xpd = NA)
  lines(c(-0.05, -0.015), c(-0.76, -0.81), lty = 2, col = "gray40")
  lines(c(-0.05, -0.015), c(-0.71, -0.68), lty = 2, col = "gray40")
  par(xpd = FALSE)
  taxacomp("Proteobacteria", legend.x = "topright")
  title("Proteobacterial classes, Bacilli, and their orders", font.main = 1)

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
  title("Yellowstone hot springs\nBowen De Le\u00f3n et al., 2013", font.main = 1)
  addhull(p1$ZC, p1$nH2O, 2, TRUE)
  pointfun(p1)
  legend <- c("Archaea", "Bacteria")
  legend("bottomleft", legend, pch = c(23, 22), col = c(1, 1), pt.bg = c(6, 5), bg = "white")

  p2 <- plotcomp("SVH+19", title = FALSE, points = FALSE)
  title("Black Sea\nSollai et al., 2019", font.main = 1)
  addhull(p2$ZC, p2$nH2O, 1, TRUE, lty = 2)
  pointfun(p2)
  legend <- c("oxic", "suboxic", "euxinic")
  legend("bottomright", legend, pch = c(24, 20, 25), pt.bg = c(4, 1, 2), bg = "white")

  p3 <- plotcomp("HLA+16", title = FALSE, points = FALSE)
  title("Baltic Sea\nHerlemann et al., 2016", font.main = 1)
  addhull(p3$ZC, p3$nH2O, 1, TRUE)
  pointfun(p3)
  legend <- c("< 6 PSU", "6-20 PSU", "> 20 PSU")
  legend("bottomright", legend, pch = c(24, 20, 21), col = c(1, 1, 1), pt.bg = c(3, NA, 4), bg = "white")

  p4 <- plotcomp("MPB+17", title = FALSE, points = FALSE)
  title("Manus Basin submarine vents\nMeier et al., 2017", font.main = 1)
  addhull(p4$ZC, p4$nH2O, 2, TRUE, lty = 2)
  pointfun(p4)
  legend <- c("T < 10 \u00B0C", "T > 10 \u00B0C", "rock", "fauna")
  legend("bottomleft", legend, pch = c(21, 23, 20, 8), col = c(1, 1, 1, "yellow4"), pt.bg = c(4, 2, NA, NA), bg = "white")

  p5 <- plotcomp("ZLM+16", title = FALSE, points = FALSE)
  title("Tibetan Plateau lakes\nZhong et al., 2016", font.main = 1)
  addhull(p5$ZC, p5$nH2O, "tan1", TRUE, lty = 2)
  pointfun(p5)
  legend <- c("< 10 g/l", "24-99 g/l", "> 300 g/l")
  legend("topright", legend, pch = c(24, 20, 21), col = c(1, 1, 1), pt.bg = c(3, NA, 4), bg = "white")

  p6 <- plotcomp("JHM+16", title = FALSE, points = FALSE)
  title("Lake Fryxell oxygen gradient\nJungblut et al., 2016", font.main = 1)
  addhull(p6$ZC, p6$nH2O, 4, TRUE)
  pointfun(p6)
  legend <- c("oxic", "transition", "anoxic")
  legend("bottomright", legend, pch = c(24, 20, 25), pt.bg = c(4, 1, 2), bg = "white")

  p7 <- plotcomp("HCW+13", title = FALSE, points = FALSE, ylim = c(-0.765, -0.7545))
  title("Guerrero Negro mat layers\nHarris et al., 2013", font.main = 1)
  addhull(p7$ZC, p7$nH2O, 4, TRUE, lty = 2)
  pointfun(p7)
  text(c(-0.1525, -0.1575, -0.1583), c(-0.7547, -0.760, -0.7647), c("0-1 mm", "1-2 mm", "2-3 mm"))
  legend <- c("photic/oxic", "low sulfide", "high sulfide")
  legend("topleft", legend, pch = c(24, 20, 25), pt.bg = c(4, 1, 2), bg = "white")

  p8 <- plotcomp("XDZ+17", title = FALSE, points = FALSE)
  title("Qarhan Salt Lake soils\nXie et al., 2017", font.main = 1)
  addhull(p8$ZC, p8$nH2O, "tan1", TRUE)
  pointfun(p8)
  legend <- c("normal", "saline")
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
  # Add convex hull for stratified lakes and sewater (from Fig. 3) 20210503
  fig3 <- geo16S3(plot.it = FALSE)
  addhull(fig3$ZC, fig3$nH2O, "slategray")
  # Add convex hulls for each dataset in this figure
  addhull(p1$ZC, p1$nH2O, 2, TRUE)
  addhull(p2$ZC, p2$nH2O, 1, TRUE, lty = 2)
  addhull(p3$ZC, p3$nH2O, 1, TRUE)
  addhull(p4$ZC, p4$nH2O, 2, TRUE, lty = 2)
  addhull(p5$ZC, p5$nH2O, "tan1", TRUE, lty = 2)
  addhull(p6$ZC, p6$nH2O, 4, TRUE)
  addhull(p7$ZC, p7$nH2O, 4, TRUE, lty = 2)
  addhull(p8$ZC, p8$nH2O, "tan1", TRUE)
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
    if(pdf) pdf("geo16S3.pdf", width = 7, height = 6)
    mat <- matrix(c(1,1,1,1, 2,3,4,5, 6,7,8,9), byrow = TRUE, nrow = 3)
    layout(mat, heights = c(1, 10, 10))
    # Make legend
    par(mar = c(0, 0, 0, 0))
    plot.new()
#    legend <- as.expression(c(quote(italic(Z)[C]), quote(O[2]), "ORP"))
    legend <- as.expression(c(quote(italic(Z)[C]), quote(O[2]), quote(NO[3]^"-" / NO[2]^"-")))
    legend("top", legend = legend, lty = 1, lwd = 1.5, col = c(1, 2, 4), pch = c(21, NA, NA), pt.bg = "white", ncol = 3, bty = "n")
    # Setup plot parameters
    par(mgp = c(1.8, 0.5, 0), mar = c(3, 3, 3, 1))
  }

  # Global ZC and O2 range for plots
  ZClim <- c(-0.17, -0.14)
  # Identify datasets to plot
  study <- c("SVH+19", "MZG+20", "MZG+20", "GBL+15", "GBL+15",
             "HXZ+20", "HXZ+20", "HXZ+20", "HXZ+20")
  column <- c("study", "lake", "lake", "size", "size",
              "station", "station", "station", "station")
  ID <- c("SVH+19", "Lake Zug", "Lake Lugano", "0.2-1.6micron", "1.6-30micron",
          "SYBL", "SYBL", "C4", "C4")
  title <- c("Black Sea", "Lake Zug", "Lake\nLugano", "ETNP", "ETNP",
#             # Use leading or trailing space to flag ORP plots
             # Use leading or trailing space to flag NO3-/NO2- plots
             "Inside\nBlue Hole", " Inside\nBlue Hole", "Outside\nBlue Hole", "Outside\nBlue Hole ")
  # Make objects to hold all ZC and nH2O values (for convex hull in Figure 2)
  allZC <- allnH2O <- numeric()
  for(i in 1:length(study)) {
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
      ylim <- rev(range(depth))
      # Determine whether the title has changed
      newplot <- TRUE
      if(i > 1) if(title[i]==title[i-1]) newplot <- FALSE
      if(newplot) {
        plot(ZC, depth, xlim = ZClim, ylim = ylim, xlab = axis.label("ZC"), ylab = "Depth (m)", type = "b")
      } else {
        # Add to plot if the title hasn't changed
        points(ZC, depth, type = "b", pch = 0)
      }

      if(newplot) {
        # Add title in lower right
        if(grepl("Outside", title[i])) text(ZClim[1], ylim[1], title[i], adj = c(0, 0))
        else text(ZClim[2], ylim[1], title[i], adj = c(1, 0))
#        # Plot O2 concentrations or ORP
        # Plot O2 concentrations or NO3-/NO2- ratio
        nc <- nchar(title[i])
        if(substr(title[i], 1, 1) == " " | substr(title[i], nc, nc) == " ") {
#          what <- "ORP"
          what <- "NO3.NO2"
#          xlim <- c(-350, 100)
          xlim <- c(0, 250)
          col <- 4
          # Calculate NO3- / NO2- ratio 20210511
          NO3.NO2 <- alldat$`NO3- (umol L-1)` / alldat$`NO2- (umol L-1)`
          alldat <- cbind(alldat, NO3.NO2)
        } else {
          what <- "O2"
          xlim <- c(0, 220)
          col <- 2
        }
        icol <- grep(paste0("^", what), colnames(alldat))
        # Remove NA values
        alldat <- alldat[!is.na(alldat[, icol]), ]
        depth <- alldat$depth
        par(new = TRUE)
        plot(alldat[, icol], depth, col = col, type = "l", axes = FALSE, xlab = "", ylab = "", xlim = xlim, ylim = ylim)
        # Add second axis labels
        xlab <- colnames(alldat)[icol]
        if(xlab == "O2 (umol kg-1)") xlab <- quote(O[2]~"(\u00B5mol kg"^-1*")")
        if(xlab == "O2 (umol L-1)") xlab <- quote(O[2]~"(\u00B5mol L"^-1*")")
        if(xlab == "NO3.NO2") xlab <- quote(NO[3]^"-" / NO[2]^"-"~"(mol/mol)")
        axis(3)
        mtext(xlab, side = 3, line = 1.7, cex = par("cex"))
        # Extra labels for ETNP
        if(title[i]=="ETNP") {
          text(40, 100, "0.2-\n1.6 \u00B5m")
          text(160, 130, "1.6-\n30 \u00B5m")
        }
        # Restore xlim for plotting ZC
        par(new = TRUE)
        plot(0, 0, type = "n", axes = FALSE, xlab = "", ylab = "", xlim = ZClim, ylim = ylim)
      }
    } # end if(plot.it)

    # Store the ZC and nH2O values
    allZC <- c(allZC, ZC)
    allnH2O <- c(allnH2O, nH2O)
  }

  if(plot.it) {
    if(pdf) dev.off()
  } else list(ZC = allZC, nH2O = allnH2O)

}


# Figure 4: Shale gas datasets 20210414
geo16S4 <- function(pdf = FALSE) {

  if(pdf) pdf("geo16S4.pdf", width = 9, height = 6)
  par(mfrow = c(2, 2))
  par(mar = c(4, 4, 1, 1))
  par(mgp = c(2.5, 1, 0))

  ## Plot A: Pennsylvania streams affected by Marcellus Shale activity 20210324

  # Data from Ulrich et al., 2018
  xlim <- c(-0.16, -0.13)
  ylim <- c(-0.755, -0.725)
  plotcomp("UKD+18.water_2014", xlim = xlim, ylim = ylim, title = FALSE)
  legend("topleft", c("MSA+", "MSA-"), pch = c(21, 1), pt.bg = c(2, 1), bg = "white")
  legend("bottomright", c("Pennsylvania Streams", "Ulrich et al., 2018"), bg = "white")
  label.figure("A", cex = 1.5, xfrac = 0.03)

  ## Plot C: Comparison of different studies on Pennsylvania Streams 20210327

  studies <- c("UKD+18.water_2014", "UKD+18.sediment_2014", "CUN+18", "MMA+20_spring", "MMA+20_fall")
  # Start plot
  plot(c(-0.15, -0.138), c(-0.744, -0.732), type = "n", xlab = cplab$ZC, ylab = cplab$nH2O)
  pch <- 21:25
  # Loop over studies
  for(i in 1:5) {
    mean <- plotcomp(studies[[i]], plot.it = FALSE)$mean
    points(mean$ZC.dn, mean$nH2O.dn, pch = pch[i], cex = 1.5, lwd = 2, bg = "white")
    lines(c(mean$ZC.dn, mean$ZC.up), c(mean$nH2O.dn, mean$nH2O.up))
    points(mean$ZC.up, mean$nH2O.up, pch = pch[i], cex = 1.8, lwd = 2, bg = 2)
  }
  # Add labels
  text(-0.1455, -0.736, "Ulrich et al., 2018\n(water)")
  text(-0.1482, -0.7404, "Ulrich et al., 2018\n(sediment)")
  text(-0.139, -0.7363, "Chen See\net al., 2018")
  text(-0.1419, -0.7404, "Mumford et al., 2020\n(spring)")
  text(-0.1432, -0.7428, "Mumford et al., 2020\n(fall)")
  # Add legend
  legend("topleft", c("Lowest disturbance", "Highest disturbance"), pch = c(21, 21), pt.bg = c("white", 2), pt.cex = c(1.4, 1.7), lwd = 2, lty = NA)
  label.figure("B", cex = 1.5, xfrac = 0.03)

  ## Plots C-D: Comparison of different studies on produced water 20210330

  # Panel C: Cluff et al., 2014
  plotcomp("CHM+14", title = FALSE)
  legend("topright", c("Injected Fluids (day = 0)", "Produced Water (day >= 49)"),
         pch = c(21, 21), pt.bg = c("white", 2), bg = "white")
  legend("bottomleft", c("Marcellus Shale", "Cluff et al., 2014"), bg = "white", inset = c(-0.03, 0))
  box()
  label.figure("C", cex = 1.5, xfrac = 0.03)

  # Panel D: Multiple studies
  studies <- c("CHM+14", "HRR+18", "ZLF+19")
  # Start plot
  plot(c(-0.22, -0.14), c(-0.75, -0.71), type = "n", xlab = cplab$ZC, ylab = cplab$nH2O)
  pch <- 21:25
  # Loop over studies
  for(i in 1:3) {
    mean <- plotcomp(studies[[i]], plot.it = FALSE)$mean
    points(mean$ZC.dn, mean$nH2O.dn, pch = pch[i], cex = 1.5, lwd = 2, bg = "white")
    lines(c(mean$ZC.dn, mean$ZC.up), c(mean$nH2O.dn, mean$nH2O.up))
    points(mean$ZC.up, mean$nH2O.up, pch = pch[i], cex = 1.8, lwd = 2, bg = 2)
  }
  # Add labels
  text(-0.165, -0.726, "Cluff et al., 2014")
  text(-0.2, -0.73, "Hull et al., 2018")
  text(-0.163, -0.746, "Zhong et al., 2019")
  # Add legend
  legend("topright", c("Source Water\nor Injected Fluids", "Produced Water"), pch = c(21, 21), pt.bg = c("white", 2), pt.cex = c(1.4, 1.7), lwd = 2, lty = NA)
  label.figure("D", cex = 1.5, xfrac = 0.03)

  if(pdf) dev.off()

}

# Figure 5: Compositional differences at different taxonomic ranks 20200924
geo16S5 <- function(pdf = FALSE) {

  if(pdf) pdf("geo16S5.pdf", width = 8, height = 6)
  mat <- matrix(c(1,1,2,2,3,3, 0,4,4,5,5,0), nrow = 2, byrow = 2)
  layout(mat)
  par(mar = c(4, 4, 3, 1))
  par(mgp = c(2.5, 1, 0))

  rankcomp("MPB+17", "ZC", pch.up = 23, pch.down = 21, ylim = c(-0.03, 0.005))
  mtext("Manus Basin", line = 1.5)
  mtext("(T > 10 \u00B0C) - (T < 10 \u00B0C)", line = 0.3, cex = 0.8)

  rankcomp("HLA+16", "nH2O", pch.up = 24, pch.down = 21, ylim = c(-0.01, 0.015))
  mtext("Baltic Sea", line = 1.5)
  mtext("(PSU < 6) - (PSU > 20)", line = 0.3, cex = 0.8)

  rankcomp("XDZ+17", "nH2O", pch.up = 24, pch.down = 21, ylim = c(-0.002, 0.03))
  mtext("Qarhan Salt Lake soils", line = 1.5)
  mtext("normal - saline", line = 0.3, cex = 0.8)

  rankcomp("MPB+17", "ZC", ylim = c(-0.02, 0.02), study2 = "HLA+16")
  mtext("Manus Basin - Baltic Sea", line = 1)

  rankcomp("BGPF13", "ZC", ylim = c(-0.12, 0.02), study2 = "XDZ+17")
  mtext("Yellowstone - Qarhan Salt Lake", line = 1)

  if(pdf) dev.off()

}