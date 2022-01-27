# JMDplots/geo16S.R
# Make plots for the paper:
# Geobiochemistry of microbial community proteomes estimated from 16S RNA sequences
# 20210416 Initial commit to JMDplots
# 20210527 Updated plots for RefSeq release 206

# Figure 1: Distinct chemical parameters of reference proteomes for taxonomic groups 20200925
geo16S1 <- function(pdf = FALSE) {

  if(pdf) pdf("geo16S1.pdf", width = 11, height = 5)
  par(mfrow = c(1, 3))

  tc1 <- taxacomp("majorphyla", legend.x = "bottomleft", hline = c(-0.81, -0.68))
  title("Major phyla and their classes", font.main = 1, cex.main = 1.4)
  label.figure("A", font = 2, cex = 1.6)
  # Draw lines indicating zoom area in next plot
  par(xpd = NA)
  lines(c(-0.05, -0.015), c(-0.81, -0.9), lty = 2, col = "gray40")
  lines(c(-0.05, -0.015), c(-0.68, -0.65), lty = 2, col = "gray40")
  par(xpd = FALSE)

  tc2 <- taxacomp("majorcellular", legend.x = "bottomleft", hline = c(-0.77, -0.71))
  title("Major cellular phyla and their classes", font.main = 1, cex.main = 1.4)
  label.figure("B", font = 2, cex = 1.6)
  par(xpd = NA)
  lines(c(-0.05, -0.015), c(-0.77, -0.81), lty = 2, col = "gray40")
  lines(c(-0.05, -0.015), c(-0.71, -0.68), lty = 2, col = "gray40")
  par(xpd = FALSE)

  tc3 <- taxacomp("proteobacteria", legend.x = "topright")
  title("Proteobacterial classes and their orders", font.main = 1, cex.main = 1.4)
  label.figure("C", font = 2, cex = 1.6)

  if(pdf) dev.off()

  # Return data for Supplementary Table 20210831
  datA <- do.call(rbind, lapply(tc1, function(x) {do.call(rbind, x)}))
  datA <- cbind(plot = "A", datA)
  datB <- do.call(rbind, lapply(tc2, function(x) {do.call(rbind, x)}))
  datB <- cbind(plot = "B", datB)
  datC <- do.call(rbind, lapply(tc3, function(x) {do.call(rbind, x)}))
  datC <- cbind(plot = "C", datC)
  out <- rbind(datA, datB, datC)
  rownames(out) <- NULL
  invisible(out)

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

  p1 <- plotmet("BGPF13", title = FALSE, points = FALSE)
#  title("Yellowstone hot springs\nBowen De Le\u00f3n et al., 2013", font.main = 1)
  title("Yellowstone hot springs", font.main = 1)
  addhull(p1$ZC, p1$nH2O, 2, TRUE)
  pointfun(p1)
  legend <- c("Archaea", "Bacteria")
  legend("bottomleft", legend, pch = c(23, 22), col = c(1, 1), pt.bg = c(6, 5), bg = "white")

  p2 <- plotmet("SVH+19", title = FALSE, points = FALSE)
#  title("Black Sea\nSollai et al., 2019", font.main = 1)
  title("Black Sea", font.main = 1)
  addhull(p2$ZC, p2$nH2O, "blue", TRUE, lty = 2)
  pointfun(p2)
  legend <- c("Oxic", "Suboxic", "Euxinic")
  legend("bottomright", legend, pch = c(24, 20, 25), pt.bg = c(4, 1, 2), bg = "white")

  p3 <- plotmet("HLA+16", title = FALSE, points = FALSE)
#  title("Baltic Sea\nHerlemann et al., 2016", font.main = 1)
  title("Baltic Sea", font.main = 1)
  addhull(p3$ZC, p3$nH2O, "blue", TRUE)
  pointfun(p3)
  legend <- c("< 6", "6-20", "> 20")
  legend("bottomright", legend, pch = c(24, 20, 21), col = c(1, 1, 1), pt.bg = c(3, NA, 4), bg = "white", title = "Salinity")

  p4 <- plotmet("MPB+17", title = FALSE, points = FALSE)
#  title("Manus Basin submarine vents\nMeier et al., 2017", font.main = 1)
  title("Manus Basin submarine vents", font.main = 1)
  addhull(p4$ZC, p4$nH2O, 2, TRUE, lty = 2)
  pointfun(p4)
  legend <- as.expression(c(quote(italic(T)~"< 50 \u00B0C"), quote(italic(T)~"> 50 \u00B0C")))
  legend("bottomleft", legend, pch = c(21, 23), col = c(1, 1), pt.bg = c(4, 2), bg = "white", title = "Water")
  legend("bottomright", c("Rock", "Fauna"), pch = c(20, 8), col = c(1, "#757500C0"), bg = "white")

  p5 <- plotmet("ZLM+16", title = FALSE, points = FALSE)
#  title("Tibetan Plateau lakes\nZhong et al., 2016", font.main = 1)
  title("Tibetan Plateau lakes", font.main = 1)
  addhull(p5$ZC, p5$nH2O, "turquoise3", TRUE, lty = 2)
  pointfun(p5)
  legend <- c("< 10 g/L", "24-99 g/L", "> 300 g/L")
  legend("topright", legend, pch = c(24, 20, 21), col = c(1, 1, 1), pt.bg = c(3, NA, 4), bg = "white", title = "Salinity")

  p6 <- plotmet("JHM+16", title = FALSE, points = FALSE)
#  title("Lake Fryxell oxygen gradient\nJungblut et al., 2016", font.main = 1)
  title("Lake Fryxell oxygen gradient", font.main = 1)
  addhull(p6$ZC, p6$nH2O, "tan1", TRUE)
  pointfun(p6)
  legend <- c("Oxic", "Transition", "Anoxic")
  legend("bottomright", legend, pch = c(24, 20, 25), pt.bg = c(4, 1, 2), bg = "white")

  p7 <- plotmet("HCW+13", title = FALSE, points = FALSE, ylim = c(-0.7685, -0.7585))
#  title("Guerrero Negro mat layers\nHarris et al., 2013", font.main = 1)
  title("Guerrero Negro mat layers", font.main = 1)
  addhull(p7$ZC, p7$nH2O, "tan1", TRUE, lty = 2)
  pointfun(p7)
  text(c(-0.1512, -0.1572, -0.1576), c(-0.7587, -0.7645, -0.7684), c("0-1 mm", "1-2 mm", "2-3 mm"))
  legend <- c("Photic/oxic", "Low sulfide", "High sulfide")
  legend("topleft", legend, pch = c(24, 20, 25), pt.bg = c(4, 1, 2), bg = "white")

  p8 <- plotmet("XDZ+17", title = FALSE, points = FALSE)
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

  # Return data for Supplementary Table 20210831
  out <- rbind(p1, p2, p3, p4, p5, p6, p7, p8)
  out$ZC <- round(out$ZC, 6)
  out$nH2O <- round(out$nH2O, 6)
  invisible(out)

}

# Figure 3: Stratified lakes and seawater 20210428
geo16S3 <- function(pdf = FALSE) {

  if(pdf) pdf("geo16S3.pdf", width = 7, height = 9)
  mat <- matrix(c(1,1,1,1, 2,3,4,5, 6,7,8,9, 10,11,12,13), byrow = TRUE, nrow = 4)
  layout(mat, heights = c(1, 10, 10, 10))
  # Make legend
  par(mar = c(0, 0, 0, 0))
  plot.new()
  legend <- as.expression(c(quote(italic(Z)[C]), quote(O[2])))
  legend("top", legend = legend, lty = c(1, 1), lwd = 1.5, col = c(1, 2), pch = c(21, NA), pt.bg = "white", ncol = 3, bty = "n")
  # Setup plot metrics
  par(mgp = c(1.8, 0.5, 0), mar = c(3, 3, 3, 1))

  # Identify datasets to plot
  # "__BSMG__" is flag for Black Sea metagenome 20220119
  study <- c("SVH+19", "MZG+20", "MZG+20", "GBL+15", "GBL+15",
             "__BSMG__", NA, "HXZ+20", "HXZ+20",
             "BCA+21", "BCA+21", "BCA+21", "BCA+21")
  column <- c("study", "lake", "lake", "size", "size",
              "", NA, "Station", "Station",
              "Month", "Month", "Month", "Month")
  ID <- c("SVH+19", "Lake Zug", "Lake Lugano", "0.2-1.6micron", "1.6-30micron",
          "", NA, "SYBL", "C4",
          "Jul", "Nov", "Feb", "Apr")
  # \n are used for vertical offset from plot bottom
  title <- c("Black Sea\n", "Lake Zug\n", "Lake\nLugano\n", "ETNP\n", "ETNP\n",
             "Black Sea\n\n", NA, "Blue Hole\n\n", "Blue Hole\n\n",
             "Ursu Lake\n\n", "Ursu Lake\n\n", "Ursu Lake\n\n", "Ursu Lake\n\n")
  subtitle <- c("", "", "", "", "",
                "Metagenome\n", NA, "Inside\n", "Outside (C4)\n",
                "July 2015\n", "November 2015\n", "February 2016\n", "April 2016\n")
  titlesub <- paste(title, subtitle)

  # Make object to hold data for Supplementary Table
  out <- list()
  # Loop over studies
  for(i in 1:length(study)) {
    if(is.na(study[i])) {
      # Leave empty space 20220119
      par(c(0, 0, 0, 0))
      plot.new()
      par(xpd = NA)
      arrows(0, 0.75, -0.4, 0.75, length = 0.1)
      text(-0.4, 0.5, paste(
        "This plot is for protein sequences",
        "inferred from shotgun metagenomes.",
        "",
        "All other plots are for estimated",
        "community proteomes from 16S rRNA",
        "and reference protein sequences.",
        sep = "\n"), cex = 0.9, adj = 0)
      par(xpd = FALSE)
      next
    }
    # ZC range for plots
    if(study[i] == "BCA+21") ZClim <- c(-0.180, -0.145) else ZClim <- c(-0.170, -0.140)
    if(study[i] == "__BSMG__") {
      ## Get data for Black Sea metagenome
      ARASTdir <- system.file("extdata/geo16S/ARAST", package = "JMDplots")
      aa <- read.csv(file.path(ARASTdir, "Black_Sea_AA.csv"))
      ZC <- ZCAA(aa)
      nH2O <- H2OAA(aa)
      depth = c(50, 70, 80, 85, 90, 95, 100, 105,
        110, 130, 170, 250, 500, 1000, 2000)
      Metagenome = c("SRR12347146", "SRR12347145", "SRR12347139", "SRR12347138", "SRR12347137", "SRR12347136", "SRR12347135", "SRR12347134",
        "SRR12347133", "SRR12347132", "SRR12347144", "SRR12347143", "SRR12347142", "SRR12347141", "SRR12347140")
      ## Get metadata for O2 concentration
      alldat <- getmdat("SVH+19")
      # Check that the samples are in the right order
      stopifnot(all(aa$protein == Metagenome))
      stopifnot(all(alldat$depth == depth))
      # Put in correct data for the metagenome
      alldat$study <- "VMW+21"
      alldat$name <- "Black Sea metagenome"
      alldat$Run <- Metagenome
    } else {
      ## Get data for 16S studies
      # Get the metadata and chemical metrics for this study
      # Keep all rows for higher-resolution O2 measurements
      mdat <- getmdat(study[i], dropNA = FALSE)
      metrics <- getmetrics(study[i])
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
    }

    # Reverse y-axis (depth)
    ylim <- rev(range(depth))
    # Visualize deeper O2 concentrations in Ursu Lake
    if(study[i] == "BCA+21") ylim <- c(11, 0)
    if(study[i] %in% c("SVH+19", "__BSMG__")) {
      # Plot 1000 and 2000 m samples closer to the others 20210608
      ylim <- c(700, 50)
      depth[match(c(1000, 2000), depth)] <- c(600, 700)
    }
    # Add space at bottom of ETNP for spearman correlations 20220114
    if(study[i] == "GBL+15") ylim <- c(320, 30)
    # Determine whether the title has changed
    newplot <- TRUE
    if(i > 1) if(titlesub[i]==titlesub[i-1]) newplot <- FALSE
    if(newplot) {
      if(study[i] %in% c("SVH+19", "__BSMG__")) {
        plot(ZC, depth, xlim = ZClim, ylim = ylim, xlab = axis.label("ZC"), ylab = "Depth (m)", type = "b", yaxt = "n")
        axis(2, at = seq(100, 700, 100), labels = c(100, 200, 300, 400, 500, 1000, 2000), gap.axis = 0)
        # Plot y-axis break 20210715
        par(xpd = NA)
        rect(-0.172, 557, -0.168, 542, col = "white", border = NA)
        text(-0.1712, 542, "/", srt = 90)
        text(-0.1712, 556, "/", srt = 90)
        par(xpd = FALSE)
      } else plot(ZC, depth, xlim = ZClim, ylim = ylim, xlab = axis.label("ZC"), ylab = "Depth (m)", type = "b")
    } else {
      # Add to plot if the title hasn't changed
      points(ZC, depth, type = "b", pch = 0)
    }
    # Save ZC and depth for Spearman correlation 20220114
    ZCdepth <- list(ZC = ZC, depth = depth)

    if(newplot) {
      # Add title in lower right
      if(grepl("Outside", subtitle[i])) {
        text(ZClim[1], ylim[1], title[i], adj = c(0, 0), font = 2)
        text(ZClim[1], ylim[1], subtitle[i], adj = c(0, 0))
      } else {
        text(ZClim[2], ylim[1], title[i], adj = c(1, 0), font = 2)
        text(ZClim[2], ylim[1], subtitle[i], adj = c(1, 0))
      }
      # Plot O2 concentrations
      nc <- nchar(title[i])
      what <- "O2"
      if(study[i] == "BCA+21") xlim <- c(0, 25) else xlim <- c(0, 220)
      col <- 2
      lty <- 1
      icol <- grep("^O2", colnames(alldat))
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
      axis(3)
      mtext(xlab, side = 3, line = 1.7, cex = par("cex"))
      # Extra labels for ETNP
      if(title[i]=="ETNP\n") {
        text(44, 76, "0.2-\n1.6 \u00B5m")
        text(135, 138, "1.6-\n30 \u00B5m")
      }
      # Restore xlim for plotting ZC
      par(new = TRUE)
      plot(0, 0, type = "n", axes = FALSE, xlab = "", ylab = "", xlim = ZClim, ylim = ylim)
    }
    # Save O2 and depth for Spearman correlation 20220114
    O2depth <- list(O2 = alldat[, icol], depth = depth)

    # Add Spearman correlation 20220114
    depths <- intersect(ZCdepth$depth, O2depth$depth)
    iZC <- match(depths, ZCdepth$depth)
    iO2 <- match(depths, O2depth$depth)
    # Calculate Spearman rank correlation; use = "complete.obs" to drop pairs with NA ZC
    spearman <- cor(ZCdepth$ZC[iZC], O2depth$O2[iO2], use = "complete.obs", method = "spearman")
    rtxt <- bquote(rho == .(formatC(spearman, digits = 2, format = "f")))
    if(grepl("ETNP", title[i])) {
      # For ETNP, put the correlations next to the lines for different size fractions
      if(newplot) text(ZClim[2], ylim[1], rtxt, adj = c(3, 0.1))
      else text(ZClim[2], ylim[1], rtxt, adj = c(1.9, 0.1))
    } else {
      if(grepl("Outside", subtitle[i])) text(ZClim[1], ylim[1], rtxt, adj = c(0, 0.1))
      else text(ZClim[2], ylim[1], rtxt, adj = c(1, 0.1))
    }

    # Assemble the data for Supplementary Table
    sd <- alldat[, c("study", "name", "Run", "sample", "depth")]
    sd <- cbind(sd, "O2 (umol kg-1)" = NA, "O2 (umol L-1)" = NA, "O2 (mg L-1)" = NA, ZC = NA, nH2O = NA)
    # Names of the columns with chemical concentrations
    cnames <- c("O2 (umol kg-1)", "O2 (umol L-1)", "O2 (mg L-1)")
    for(cname in cnames) if(cname %in% colnames(alldat)) sd[, cname] <- alldat[, cname]
    if(study[i] == "__BSMG__") {
      sd$ZC <- ZC
      sd$nH2O <- nH2O
    } else {
      # Put ZC and nH2O values in correct place
      metrics <- metrics[imet, ]
      metrics <- metrics[!is.na(metrics$Run), ]
      isd <- match(metrics$Run, sd$Run)
      sd$ZC[isd] <- metrics$ZC
      sd$nH2O[isd] <- metrics$nH2O
    }

    # Replace NA name with study name
    sd.name <- na.omit(sd$name)[1]
    sd$name[is.na(sd$name)] <- sd.name
    # Use NA instead of "" for missing Run and sample name
    sd$Run[sd$Run == ""] <- NA
    sd$sample[sd$sample == ""] <- NA
    # Place data frame into output list
    out[[i]] <- sd

  } # end loop

  if(pdf) dev.off()

  # Return data for Supplementary Table 20210831
  out <- do.call(rbind, out)
  out$ZC <- round(out$ZC, 6)
  out$nH2O <- round(out$nH2O, 6)
  invisible(out)

}

# Figure 4: Shale gas datasets and ZC difference between oxidized and reduced samples 20210414
geo16S4 <- function(pdf = FALSE) {

  if(pdf) pdf("geo16S4.pdf", width = 10, height = 6)
  layout(matrix(c(1, 3, 2, 4, 5, 5), nrow = 2), widths = c(2, 2, 1.8))
  par(mar = c(4, 4, 1, 1))
  par(mgp = c(2.5, 1, 0))

  ## Plot A: Pennsylvania streams affected by Marcellus Shale activity 20210324

  # Data from Ulrich et al., 2018
  xlim <- c(-0.16, -0.13)
  ylim <- c(-0.755, -0.725)
  plotmet("UKD+18.water_2014", xlim = xlim, ylim = ylim, title = FALSE)
  legend("topleft", c("MSA-", "MSA+"), pch = c(1, 21), pt.bg = c(1, 2), bg = "white", title = "NW PA streams (2014)")
  label.figure("A", cex = 1.5, xfrac = 0.03, font = 2)

  ## Plot B: Comparison of different studies on Pennsylvania Streams 20210327

  studies <- c("UKD+18.water_2014", "UKD+18.sediment_2014", "MMA+20_spring", "MMA+20_fall")
  # Start plot
  plot(c(-0.148, -0.140), c(-0.748, -0.735), type = "n", xlab = cplab$ZC, ylab = cplab$nH2O)
  pch <- 21:25
  xadj <- c(0.5, 0.5, 0.5, -0.5)
  yadj <- c(-0.8, -0.8, -0.8, 0.5)
  outB <- list()
  # Loop over studies
  for(i in 1:4) {
    pm <- plotmet(studies[[i]], plot.it = FALSE, extracolumn = "type")
    # Determine sample groups from values of pch returned by plotmet()  20210901
    i1 <- pm$pch == 1
    i2 <- pm$pch == 21
    means <- list(ZC1 = mean(pm$ZC[i1]), ZC2 = mean(pm$ZC[i2]), nH2O1 = mean(pm$nH2O[i1]), nH2O2 = mean(pm$nH2O[i2]))
    points(means$ZC1, means$nH2O1, pch = pch[i], cex = 1.5, lwd = 2, bg = "#ffffffa0")
    lines(c(means$ZC1, means$ZC2), c(means$nH2O1, means$nH2O2))
    points(means$ZC2, means$nH2O2, pch = pch[i], cex = 1.8, lwd = 2, bg = "#df536ba0")
    # Plot number of samples next to points 20210902
    n1 <- length(pm$ZC[i1])
    n2 <- length(pm$ZC[i2])
    if(i==4) dyadj <- 0.2 else dyadj <- 0
    text(means$ZC1, means$nH2O1, n1, adj = c(xadj[i], yadj[i] - dyadj))
    text(means$ZC2, means$nH2O2, n2, adj = c(xadj[i], yadj[i] + dyadj))
    # Save values for Supplementary Table 20210901
    outB[[i]] <- pm
  }
  # Add labels
  text(-0.143, -0.7393, "NW PA\nstreams (2014)")
  text(-0.1468, -0.7433, "NW PA\nsediment (2014)")
  text(-0.1428, -0.7455, "PASF streams (spring)")
  text(-0.1432, -0.7473, "PASF streams (fall)")
  # Add legend
  legend("topleft", c("MSA- or lowest disturbance", "MSA+ or highest disturbance"), pch = c(21, 21), pt.bg = c("#ffffffa0", "#df536ba0"), pt.cex = c(1.4, 1.7), lwd = 2, lty = NA)
  label.figure("B", cex = 1.5, xfrac = 0.03, font = 2)

  ## Plots C-D: Comparison of different studies on produced water 20210330

  # Panel C: Cluff et al., 2014
  plotmet("CHM+14_injected-49", title = FALSE)
  legend("topright", c("Injected fluids (day 0)", "Produced water (day 49 and after)"),
         pch = c(21, 21), pt.bg = c("white", 2), bg = "white", title = "Marcellus Shale")
  box()
  label.figure("C", cex = 1.5, xfrac = 0.03, font = 2)

  # Panel D: Multiple studies
  studies <- c("CHM+14_injected-49", "HRR+18_injected-22", "ZLF+19_injected-18")
  # Start plot
  plot(c(-0.22, -0.14), c(-0.75, -0.71), type = "n", xlab = cplab$ZC, ylab = cplab$nH2O)
  pch <- 21:25
  outD <- list()
  # Loop over studies
  for(i in 1:3) {
    pm <- plotmet(studies[[i]], plot.it = FALSE, extracolumn = "type")
    # Determine sample groups from values of pch returned by plotmet()  20210901
    i1 <- pm$pch == 1
    i2 <- pm$pch == 21
    means <- list(ZC1 = mean(pm$ZC[i1]), ZC2 = mean(pm$ZC[i2]), nH2O1 = mean(pm$nH2O[i1]), nH2O2 = mean(pm$nH2O[i2]))
    points(means$ZC1, means$nH2O1, pch = pch[i], cex = 1.5, lwd = 2, bg = "#ffffffa0")
    lines(c(means$ZC1, means$ZC2), c(means$nH2O1, means$nH2O2))
    points(means$ZC2, means$nH2O2, pch = pch[i], cex = 1.8, lwd = 2, bg = "#df536ba0")
    # Plot number of samples next to points 20210902
    n1 <- length(pm$ZC[i1])
    n2 <- length(pm$ZC[i2])
    if(i == 3) adj <- c(-1, 1) else adj <- c(-1, 0.5)
    text(means$ZC1, means$nH2O1, n1, adj = adj)
    text(means$ZC2, means$nH2O2, n2, adj = c(0.2, -1))
    outD[[i]] <- pm
  }
  # Add labels
  text(-0.168, -0.726, "Marcellus Shale")
  text(-0.202, -0.731, "Denver-Julesburg Basin")
  text(-0.172, -0.745, "Duvernay Formation")
  # Add legend
  legend("topright", c("Source water\nor injected fluids", "Produced water"), pch = c(21, 21), pt.bg = c("#ffffffa0", "#df536ba0"), pt.cex = c(1.4, 1.7), lwd = 2, lty = NA)
  label.figure("D", cex = 1.5, xfrac = 0.03, font = 2)

  ## Panel E: ZC differences between oxidized and reduced sample subsets in various studies
  study <- c(
    "GBL+15", "JHM+16_O2", "MPB+17", "BCA+21",
    "SVH+19_O2", "MZG+20_Zug", "MZG+20_Lugano", "HXZ+20",
    "UKD+18.water_2014", "UKD+18.sediment_2014", "MMA+20_spring", "MMA+20_fall",
    "CHM+14_injected-49", "HRR+18_injected-22", "ZLF+19_injected-18",
    "SMS+12", "EH18"
  )
  description <- c(
    "ETNP (0.2-1.6 \u00B5M)", "Lake Fryxell mat", "Manus Basin (water samples)", "Ursu Lake (all months)",
    "Black Sea", "Lake Zug", "Lake Lugano", "Blue Hole (Inside)",
    "NW PA streams (2014)", "NW PA sediment (2014)", "PASF streams (spring)", "PASF streams (fall)",
    "Marcellus Shale", "Denver-Julesburg Basin", "Duvernay Formation",
    "Bison Pool", "Mono Lake"
  )
  # Description of conditions of sample groups
  cond2 <- c("anoxic", "anoxic", "> 50 \u00B0C", "anoxic",
             "anoxic", "anoxic", "anoxic", "anoxic",
             "MSA+", "MSA+", "highest", "highest",
             "PW day 49+", "PW day 22+", "PW day 18",
             "anoxic", "anoxic"
  )
  cond1 <- c("oxic", "oxic", "< 50 \u00B0C", "oxic",
             "oxic", "oxic", "oxic", "oxic",
             "MSA-", "MSA-", "lowest", "lowest",
             "IF day 0", "SW day 0", "SW day 0",
             "oxic", "oxic"
  )
  # Values of pch (from getmdat()/plotmet()) for oxidized and reduced sample groups
  pch_ox <- c(24, 24, 21, 24,
            24, 24, 24, 24,
            1, 1, 1, 1,
            1, 1, 1,
            24, 24
  )
  pch_red <- c(25, 25, 23, 25,
            25, 25, 25, 25,
            21, 21, 21, 21,
            21, 21, 21,
            25, 25
  )

  # Loop over studies and get ZC difference
  P <- n1 <- n2 <- DZC <- numeric()
  for(i in 1:length(study)) {
    # Get metrics for samples in this study
    pm <- plotmet(study[[i]], plot.it = FALSE)
    # Determine oxidized and reduced sample groups from values of pch returned by plotmet()  20210901
    pm <- pm[!is.na(pm$pch), ]
    iox <- pm$pch == pch_ox[i]
    ired <- pm$pch == pch_red[i]
    # Keep track of the mean difference of ZC
    DZC <- c(DZC, mean(pm$ZC[ired]) - mean(pm$ZC[iox]))
    # Keep track of number of samples 20210902
    n1 <- c(n1, length(pm$ZC[iox]))
    n2 <- c(n2, length(pm$ZC[ired]))
  }

  # Include numbers of samples in condition text 20210902
  condition <- paste0(cond2, " (", n2, ") - ", cond1, " (", n1, ")")

  # Order samples by Delta ZC 20220118
  orderDZC <- order(DZC)
  DZC <- DZC[orderDZC]
  description <- description[orderDZC]
  condition <- condition[orderDZC]

  # Setup plot
  nsamp <- length(DZC)
  par(mar = c(4, 12, 3, 0.5), las = 1)
  plot(c(min(DZC), 0.01), c(nsamp + 0.5, 0.5), ylim = c(nsamp + 0.5, 0.5), yaxs = "i", yaxt = "n", ylab = "", xlab = cplab$DZC, type = "n")
  title("Mean differences between oxidized                         ", font.main = 1, line = 1.6)
  title("and reduced sample groups                         ", font.main = 1, line = 0.5)
  # Add line at DZC = 0
  abline(v = 0, lty = 2, col = "gray40")
  # Use red for shale gas and hydrothermal systems
  col <- rep(1, length(description))
  col[description %in% c("Bison Pool", "Manus Basin (water samples)", "Marcellus Shale", "Denver-Julesburg Basin", "Duvernay Formation")] <- 2
  # Plot Delta ZC 
  for(i in 1:nsamp) lines(c(DZC[i], DZC[i]), c(i - 0.5, i + 0.5), lwd = 2, col = col[i])

  # Add dataset description and conditions
  axis(2, at = 1:nsamp - 0.2, labels = description, tick = FALSE)
  axis(2, at = 1:nsamp + 0.2, labels = condition, tick = FALSE, cex.axis = 0.9)
  label.figure("E", cex = 1.5, yfrac = 0.97, font = 2)

  if(pdf) dev.off()

  # Return data for Supplementary Table 20210901
  outB <- do.call(rbind, outB)
  outD <- do.call(rbind, outD)
  out <- rbind(outB, outD)
  out <- out[, !colnames(out) %in% c("pch", "col")]
  out$nH2O <- round(out$nH2O, 6)
  out$ZC <- round(out$ZC, 6)
  invisible(out)

}

# Function to plot individual datasets for geo16S5 20220112
# Use H2O = FALSE for ZC, H2O = TRUE for nH2O
MG16S <- function(which, plot.lines = TRUE, lowest.level = NULL, lineage = NULL, rm.outliers = FALSE, H2O = FALSE, cex = 1) {

  # For MG datasets analyzed in this study (others are from gradox paper)
  ARASTdir <- system.file("extdata/geo16S/ARAST", package = "JMDplots")
  # Read data for paired metagenomes and amplicon sequences expanded from Tax4Fun paper (Aßhauer et al., 2015)
  AWDM15file <- system.file("extdata/geo16S/AWDM15.csv", package = "JMDplots")
  AWDM15 <- read.csv(AWDM15file)

  # Start a plot if there isn't one 20220122
  if(is.null(dev.list())) {
    xylim <- c(-0.22, -0.08)
    xlab <- quote(italic(Z)[C]~"from shotgun metagenome or metatranscriptome")
    ylab <- quote(italic(Z)[C]~"estimated from 16S rRNA")
    plot(xylim, xylim, type = "n", xlab = xlab, ylab = ylab)
  }

  # Use semi-transparent colors 20220122
  c1 <- addalpha(1, "80")
  c2 <- addalpha(2, "b0")
  c4 <- addalpha(4, "b0")
  c5 <- addalpha(5, "b0")
  c6 <- addalpha(6, "b0")
  c8 <- addalpha(8, "b0")

  if(which == "Guerrero_Negro") {
    ## Guerrero Negro metagenome (Kunin et al., 2008)
    dat_MG <- mplot("Guerrero_Negro", "IMG_MGP", plot.it = FALSE, H2O = H2O)
    # Reverse the order because upper mat layers are plotted on the right
    metric_MG <- rev(dat_MG$AA)
    # 16S data (Harris et al., 2013)
    mdat <- getmdat("HCW+13")
    dat_16S <- getmetrics("HCW+13", mdat = mdat, lowest.level = lowest.level, lineage = lineage)
    if(H2O) metric_16S <- dat_16S$nH2O else metric_16S <- dat_16S$ZC
    # Check that the sample names are the same
    stopifnot(all.equal(rev(rownames(dat_MG$meancomp)), gsub(".*_", "", dat_16S$sample)))
    # Add lines and points
    if(plot.lines) lines(metric_MG, metric_16S)
    points(metric_MG, metric_16S, pch = 21, bg = "white", cex = cex)
    # Fill symbol for most oxidized (surface) sample
    points(metric_MG[1], metric_16S[1], pch = 21, bg = 4, cex = cex)
    # Get sample name and ID for Supplemental Table 20220125
    Sample <- mdat$sample
    Amplicon <- mdat$GenBank
    gradox_S1 <- read.csv(system.file("extdata/gradox/Table_S1.csv", package = "JMDplots"))
    Metagenome <- gradox_S1$ID[gradox_S1$study.name == "Guerrero_Negro"]
  }

  if(which == "ETNP_MG") {
    ## ETNP OMZ metagenome (Glass et al., 2015)
    dat_MG <- mplot("ETNP_OMZ", "SRA_MGP", plot.it = FALSE, H2O = H2O)
    # Reverse the order because smaller water depths are plotted on the right
    metric_MG <- rev(dat_MG$AA)
    # 16S data (Ganesh et al., 2015)
    mdat <- getmdat("GBL+15")
    dat_16S <- getmetrics("GBL+15", mdat = mdat, lowest.level = lowest.level, lineage = lineage)
    # Use smallest size fraction
    dat_16S <- subset(dat_16S, grepl("1.6micron", dat_16S$sample))
    mdat <- mdat[mdat$Run %in% dat_16S$Run, ]
    if(H2O) metric_16S <- dat_16S$nH2O else metric_16S <- dat_16S$ZC
    # Check that the sample names are the same
    stopifnot(all.equal(rev(rownames(dat_MG$meancomp)), gsub("_.*", "", dat_16S$sample)))
    # Add lines and points
    if(plot.lines) lines(metric_MG, metric_16S)
    points(metric_MG, metric_16S, pch = 21, bg = "white", cex = cex)
    # Fill symbol for surface sample
    points(metric_MG[1], metric_16S[1], pch = 21, bg = 4, cex = cex)
    # Get sample name and ID for Supplemental Table 20220125
    Sample <- paste0(mdat$depth, "m_", mdat$size)
    Amplicon <- mdat$Run
    gradox_S1 <- read.csv(system.file("extdata/gradox/Table_S1.csv", package = "JMDplots"))
    Metagenome <- gradox_S1$ID[gradox_S1$study.name == "ETNP_OMZ" & gradox_S1$type == "MG"]
  }

  if(which == "ETNP_MT") {
    ## ETNP OMZ metatranscriptome (Ganesh et al., 2015)
    dat_MT <- mplot("ETNP_OMZ", "SRA_MTP", plot.it = FALSE, H2O = H2O)
    # Reverse the order because smaller water depths are plotted on the right
    metric_MT <- rev(dat_MT$AA)
    # 16S data (Ganesh et al., 2015)
    mdat <- getmdat("GBL+15")
    dat_16S <- getmetrics("GBL+15", lowest.level = lowest.level, lineage = lineage)
    # Use smallest size fraction
    dat_16S <- subset(dat_16S, grepl("1.6micron", dat_16S$sample))
    mdat <- mdat[mdat$Run %in% dat_16S$Run, ]
    if(H2O) metric_16S <- dat_16S$nH2O else metric_16S <- dat_16S$ZC
    # Check that the sample names are the same
    stopifnot(all.equal(rev(rownames(dat_MT$meancomp)), gsub("_.*", "", dat_16S$sample)))
    # Add lines and points
    if(plot.lines) lines(metric_MT, metric_16S, col = 8)
    points(metric_MT, metric_16S, pch = 22, bg = "white", col = 8, cex = cex)
    # Fill symbol for surface sample
    points(metric_MT[1], metric_16S[1], pch = 22, bg = 4, col = 8, cex = cex)
    # For the return value
    metric_MG <- metric_MT
    # Get sample name and ID for Supplemental Table 20220125
    Sample <- paste0(mdat$depth, "m_", mdat$size)
    Amplicon <- mdat$Run
    gradox_S1 <- read.csv(system.file("extdata/gradox/Table_S1.csv", package = "JMDplots"))
    Metagenome <- gradox_S1$ID[gradox_S1$study.name == "ETNP_OMZ" & gradox_S1$type == "MT"]
  }

  if(which == "Bison_Pool") {
    ## Bison Pool metagenome (Havig et al., 2011)
    dat_MG <- mplot("Bison_Pool", "IMG_MGP", plot.it = FALSE, H2O = H2O)
    metric_MG <- dat_MG$AA
    # 16S data (Swingley et al., 2012)
    # mincount needs to be lowered from default for when lineage = "genus"
    # (site 4 (Q) has less than 200 genus-level classifications)
    mdat <- getmdat("SMS+12")
    dat_16S <- getmetrics("SMS+12", mdat = mdat, lowest.level = lowest.level, lineage = lineage, mincount = 50)
    if(H2O) metric_16S <- dat_16S$nH2O else metric_16S <- dat_16S$ZC
    # Check that the sample names are the same
    stopifnot(all.equal(rownames(dat_MG$meancomp), mdat$"Field Code"))
    # Add lines and points
    if(plot.lines) lines(metric_MG, metric_16S)
    points(metric_MG, metric_16S, pch = 21, bg = "transparent", cex = cex)
    # Fill symbol for low-T sample
    points(metric_MG[5], metric_16S[5], pch = 21, bg = c4, cex = cex)
    # Get sample name and ID for Supplemental Table 20220125
    gradox_S1 <- read.csv(system.file("extdata/gradox/Table_S1.csv", package = "JMDplots"))
    Library <- sapply(strsplit(gradox_S1$sample.description[gradox_S1$study.name == "Bison_Pool"], " "), "tail", 1)
    Sample <- paste0("Site ", mdat$Sample, " (", mdat$`Field Code`, ") (", Library, ")")
    Amplicon <- mdat$Run
    Metagenome <- gradox_S1$ID[gradox_S1$study.name == "Bison_Pool"]
  }

  if(which == "Mono_Lake") {
    ## Mono Lake metatranscriptome (Edwardson and Hollibaugh, 2017)
    dat_MT <- mplot("Mono_Lake", "SRA_MTP", plot.it = FALSE, H2O = H2O)
    # NOTE: reverse the order because smaller water depths are plotted on the right
    metric_MT <- rev(dat_MT$AA)
    # 16S data (Edwardson and Hollibaugh, 2018)
    mdat <- getmdat("EH18")
    dat_16S <- getmetrics("EH18", lowest.level = lowest.level, lineage = lineage)
    if(H2O) metric_16S <- dat_16S$nH2O else metric_16S <- dat_16S$ZC
    # Check that the sample names are the same
    stopifnot(all.equal(rev(rownames(dat_MT$meancomp)), gsub(".*_", "", dat_16S$sample)))
    # Add lines and points
    if(plot.lines) lines(metric_MT, metric_16S, pch = 0, col = 8)
    points(metric_MT, metric_16S, pch = 22, bg = "white", col = 8, cex = cex)
    # Fill symbol for surface sample
    points(metric_MT[1], metric_16S[1], pch = 22, bg = 4, col = 8, cex = cex)
    # For the return value
    metric_MG <- metric_MT
    # Get sample name and ID for Supplemental Table 20220125
    Sample <- paste0(mdat$sample)
    Amplicon <- mdat$Run
    gradox_S1 <- read.csv(system.file("extdata/gradox/Table_S1.csv", package = "JMDplots"))
    Metagenome <- gradox_S1$ID[gradox_S1$study.name == "Mono_Lake"]
  }

  if(which == "Marcellus_Shale") {
    ## Marcellus metagenomes (Daly et al., 2016) 20211218
    aa <- read.csv(file.path(ARASTdir, "Marcellus_Shale_AA.csv"))
    if(H2O) metric_MG <- H2OAA(aa) else metric_MG <- ZCAA(aa)
    # Marcellus 16S (Cluff et al., 2014)
    mdat <- getmdat("CHM+14")
    dat_16S <- getmetrics("CHM+14", mdat = mdat, lowest.level = lowest.level, lineage = lineage)
    # Time points: input, T7, T13, T82, T328
    Sample = paste("Day", c(0, 7, 13, 82, 328))
    # List run IDs here
    Metagenome = c("SRR3111417", "SRR3111625", "SRR3111724", "SRR3111729", "SRR3111737")
    Amplicon = c("SRR1184016", "SRR1184060", "SRR1184062", "SRR1184081", "SRR1184083")
    # GC content from https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR******* 20220126
    GC_MG <- c(58.4, 37.7, 57, 39.6, 53.1)
    GC_16S <- c(53.9, 55, 54.8, 51.8, 53.2)
    ## 16S replicate 2
    #Amplicon2 = c("SRR1184049", "SRR1184061", "SRR1184063", "SRR1184082", "SRR1184084")
    # Make sure metagenomes are in correct order
    stopifnot(all(aa$protein == Metagenome))
    # Get 16S runs corresponding to metagenomes
    idat <- match(Amplicon, dat_16S$Run)
    dat_16S <- dat_16S[idat, ]
    if(H2O) metric_16S <- dat_16S$nH2O else metric_16S <- dat_16S$ZC
    # Assign colors: open circle for injected fluid, gray for flowback, red for produced
    mdat <- mdat[idat, ]
    col <- rep(4, nrow(mdat))
    col[mdat$type == "flowback fluid"] <- 8
    col[mdat$type == "produced fluid"] <- 2
    if(plot.lines) type <- "b" else type <- "p"
    points(metric_MG, metric_16S, pch = 23, bg = col, type = type, cex = cex)
  }

  if(which == "Manus_Basin") {
    ## Manus Basin metagenomes (Meier et al., 2017) 20220110
    aa <- read.csv(file.path(ARASTdir, "Manus_Basin_AA.csv"))
    if(H2O) metric_MG <- H2OAA(aa) else metric_MG <- ZCAA(aa)
    # Manus Basin 16S (Meier et al., 2017)
    mdat <- getmdat("MPB+17")
    dat_16S <- getmetrics("MPB+17", mdat = mdat, lowest.level = lowest.level, lineage = lineage)
    # Samples: NSu-F2b, NSu-F5, Fw-F1b, Fw-F3, RR-F1b
    Sample = c("MNB27-NSu-F2b", "MNB29-NSu-F5", "MNB14-Fw-F1b", "MNB17-Fw-F3", "MNB45-RR-F1b")
    # List run IDs here
    Metagenome = c("ERR1679394", "ERR1679395", "ERR1679397", "ERR1679396", "ERR1679398")
    Amplicon = c("ERR1665247", "ERR1665249", "ERR1665234", "ERR1665237", "ERR1665265")
    # Make sure metagenomes are in correct order
    stopifnot(all(aa$protein == Metagenome))
    # Get 16S runs corresponding to metagenomes
    idat <- match(Amplicon, dat_16S$Run)
    dat_16S <- dat_16S[idat, ]
    if(H2O) metric_16S <- dat_16S$nH2O else metric_16S <- dat_16S$ZC
    # Assign colors: blue for < 10 degC, red for > 50 degC
    mdat <- mdat[idat, ]
    col <- rep("white", nrow(mdat))
    col[mdat$T > 50] <- c2
    col[mdat$T < 10] <- c4
    points(metric_MG, metric_16S, pch = 22, bg = col, col = c1, cex = cex)
  }

  if(which == "Black_Sea") {
    ## Black Sea metagenome (Villanueva et al., 2021) 20220115
    aa <- read.csv(file.path(ARASTdir, "Black_Sea_AA.csv"))
    if(H2O) metric_MG <- H2OAA(aa) else metric_MG <- ZCAA(aa)
    # 16S data (Sollai et al., 2019)
    mdat <- getmdat("SVH+19")
    dat_16S <- getmetrics("SVH+19", mdat = mdat)
    if(H2O) metric_16S <- dat_16S$nH2O else metric_16S <- dat_16S$ZC
    Sample = c(50, 70, 80, 85, 90, 95, 100, 105,
      110, 130, 170, 250, 500, 1000, 2000)
    Metagenome = c("SRR12347146", "SRR12347145", "SRR12347139", "SRR12347138", "SRR12347137", "SRR12347136", "SRR12347135", "SRR12347134",
      "SRR12347133", "SRR12347132", "SRR12347144", "SRR12347143", "SRR12347142", "SRR12347141", "SRR12347140")
    # Check that the samples are in the right order
    stopifnot(all(mdat$depth == Sample))
    stopifnot(all(aa$protein == Metagenome))
    # Assign colors and symbols: blue up for < 100 m, red down for >= 100 m
    col <- rep(c4, nrow(mdat))
    pch <- rep(24, nrow(mdat))
    col[mdat$depth >= 100] <- c2
    pch[mdat$depth >= 100] <- 25
    points(metric_MG, metric_16S, pch = pch, bg = col, col = c1, cex = cex)
    Amplicon <- mdat$Run
  }

  if(which == "HMP") {
    # HMP 16S
    mdat <- getmdat("HMP12")
    met <- getmetrics("HMP12", mdat = mdat, lowest.level = lowest.level, lineage = lineage)
    # HMP metagenomes
    aa <- read.csv(file.path(ARASTdir, "HMP_AA.csv"))
    # Put data in same order
    dat <- AWDM15[AWDM15$Name == "HMP", ]
    imet <- match(dat$Amplicon, met$Run)
    met <- met[imet, ]
    iaa <- match(dat$Metagenome, aa$protein)
    aa <- aa[iaa, ]
    # Make sure the 16S and metagenomes are paired correctly
    stopifnot(all(na.omit(met$Run == dat$Amplicon)))
    stopifnot(all(aa$protein == dat$Metagenome))
    # Get ZC values
    if(H2O) metric_16S <- met$nH2O else metric_16S <- met$ZC
    if(H2O) metric_MG <- H2OAA(aa) else metric_MG <- ZCAA(aa)
    # Don't plot MG with low numbers of protein fragments 20220122
    ilow <- aa$chains < 50000
    if(any(ilow)) {
      metric_MG <- metric_MG[!ilow]
      metric_16S <- metric_16S[!ilow]
      mdat <- mdat[!ilow, ]
      aa <- aa[!ilow, ]
      dat <- dat[!ilow, ]
    }
    # Remove outliers (anomalously high ZC in metagenome) 20221215
    if(rm.outliers & !H2O) {
      iout.Nasal <- dat$Body.site=="Nasal cavity" & metric_MG > -0.12
      iout.UG <- dat$Body.site=="UG tract" & metric_MG > -0.14
      iout <- iout.Nasal | iout.UG
      metric_MG <- metric_MG[!iout]
      metric_16S <- metric_16S[!iout]
      mdat <- mdat[!iout, ]
      aa <- aa[!iout, ]
      dat <- dat[!iout, ]
    }
    # Colors: blue (Skin), green (Nasal cavity), gray (Oral cavity), red (GI tract), magenta (UG tract)
    # Symbols: up triangle (skin, GI tract), circle (Oral cavity), down triangle (Nasal cavity, UG tract)
    col <- sapply(mdat$"Body site", switch, "Skin" = c5, "Nasal cavity" = c4, "Oral cavity" = c8, "GI tract" = c2, "UG tract" = c6)
    pch <- sapply(mdat$"Body site", switch, "Skin" = 24, "Nasal cavity" = 25, "Oral cavity" = 21, "GI tract" = 24, "UG tract" = 25)
    points(metric_MG, metric_16S, pch = pch, bg = col, col = c1)
    # Get sample name and ID for Supplemental Table 20220125
    Sample <- dat$Sample.name
    Amplicon <- dat$Amplicon
    Metagenome <- dat$Metagenome
  }

  if(which == "Guts") {
    # Guts 16S
    met <- getmetrics("MKK+11", lowest.level = lowest.level, lineage = lineage, mincount = 50)
    # Guts metagenomes
    aa <- read.csv(file.path(ARASTdir, "Guts_AA.csv"))
    # Make sure the 16S and metagenomes are paired correctly
    dat <- AWDM15[AWDM15$Name == "Guts", ]
    stopifnot(all(met$Run == paste0("mgm", dat$Amplicon)))
    stopifnot(all(aa$protein == paste0("mgm", dat$Metagenome)))
    # Get ZC values
    if(H2O) metric_16S <- met$nH2O else metric_16S <- met$ZC
    if(H2O) metric_MG <- H2OAA(aa) else metric_MG <- ZCAA(aa)
    points(metric_MG, metric_16S, pch = 21, bg = c2, col = c1, cex = cex)
    # Get sample name and ID for Supplemental Table 20220125
    Sample <- dat$Sample.name
    Amplicon <- paste0("mgm", dat$Amplicon)
    Metagenome <- paste0("mgm", dat$Metagenome)
  }

  if(which == "Soils") {
    # Soils 16S
    met <- getmetrics("FLA+12", lowest.level = lowest.level, lineage = lineage)
    # Soils metagenomes
    aa <- read.csv(file.path(ARASTdir, "Soils_AA.csv"))
    # Put data in same order
    dat <- AWDM15[AWDM15$Name == "Soils", ]
    imet <- match(paste0("mgm", dat$Amplicon), met$Run)
    met <- met[imet, ]
    iaa <- match(paste0("mgm", dat$Metagenome), aa$protein)
    aa <- aa[iaa, ]
    # Make sure the 16S and metagenomes are paired correctly
    stopifnot(all(met$Run == paste0("mgm", dat$Amplicon)))
    stopifnot(all(aa$protein == paste0("mgm", dat$Metagenome)))
    # Get ZC values
    if(H2O) metric_16S <- met$nH2O else metric_16S <- met$ZC
    if(H2O) metric_MG <- H2OAA(aa) else metric_MG <- ZCAA(aa)
    points(metric_MG, metric_16S, pch = 21, bg = c4, col = c1, cex = cex)
    # Get sample name and ID for Supplemental Table 20220125
    Sample <- dat$Sample.name
    Amplicon <- paste0("mgm", dat$Amplicon)
    Metagenome <- paste0("mgm", dat$Metagenome)
  }

  list(which = rep(which, length(Sample)), metric_MG = metric_MG, metric_16S = metric_16S, Sample = Sample, Amplicon = Amplicon, Metagenome = Metagenome)

}

# Comparison of protein ZC from metagenomic or metatranscriptomic data with estimates from 16S and reference sequences 20211017
# Add Marcellus, HMP, and Soils and Guts 20211218
geo16S5 <- function(pdf = FALSE) {

  if(pdf) pdf("geo16S5.pdf", width = 10, height = 6)
  mat <- matrix(c(1,1, 2,2, 3,3, 0, 4,4, 5,5, 0), nrow = 2, byrow = TRUE)
  layout(mat)
  par(mar = c(4, 4, 2, 1), mgp = c(2.8, 1, 0))
  xylim <- c(-0.22, -0.12)
  xlimHMP <- c(-0.22, -0.08)

  ### Panel A: Comparisons with metagenomes analyzed by Dick et al. (2019)

  # Start plot A
  xlab <- quote(italic(Z)[C]~"from shotgun metagenome or metatranscriptome")
  ylab <- quote(italic(Z)[C]~"estimated from 16S rRNA")
  plot(xylim, xylim, type = "n", xlab = xlab, ylab = ylab)
  lines(xylim, xylim, lty = 2, col = "gray40")
  out1 <- MG16S("Guerrero_Negro")
  # Label points for upper 3 layers
  text(c(-0.135, -0.1382, -0.1422), c(-0.1472, -0.1598, -0.1618), c("1", "2", "3"), cex = 0.8)
  text(-0.135, -0.174, "Guerrero\nNegro\nmat")
  out2 <- MG16S("ETNP_MG")
  out3 <- MG16S("ETNP_MT")
  text(-0.172, -0.160, "ETNP water")
  out4 <- MG16S("Bison_Pool", cex = 1.4, plot.lines = FALSE)
  # Make arrows to show outflow channel 20220120
  dx.fraction <- c(0.1, 0.08, 0.08)
  for(i in 1:3) {
    x <- out4$metric_MG[i:(i+1)]
    y <- out4$metric_16S[i:(i+1)]
    lmxy <- lm(y ~ x)
    # Find coordinates for a fraction of the line length
    dx <- diff(x) * dx.fraction[i]
    xs <- c(x[1] + dx, x[2] - dx)
    ys <- predict(lmxy, data.frame(x = xs))
    arrows(xs[1], ys[1], xs[2], ys[2], length = 0.1)
  }
  text(-0.16, -0.212, "Bison Pool")
  text(-0.1905, -0.2172, "Hot spring source", cex = 0.9)
  text(-0.1875, -0.208, "Outflow channel", cex = 0.9)
  out5 <- MG16S("Mono_Lake")
  text(-0.182, -0.148, "Mono Lake water")
  # Add legend
  legend("topleft", rep("                                ", 3), pch = c(21, 21, 22), col = c(1, 1, 8), pt.bg = c(4, "white", "white"))
  ltxt <- as.expression(c(quote("Highest "*O[2]), "Metagenome", "Metatranscriptome"))
  legend("topleft", ltxt, pch = c(22, NA, NA), col = c(8, NA, NA), pt.bg = c(4, NA, NA), inset = c(0.025, 0), bty = "n")
  # Add title and figure label
  title("Various Environments", font.main = 1, cex.main = 1.1)
  label.figure("A", cex = 1.5, font = 2, xfrac = 0.04, yfrac = 0.96)

  ### Panels B-E: Comparisons with metagenomes analyzed in this study

  # Start plot B
  xlab <- quote(italic(Z)[C]~"from shotgun metagenome")
  plot(xylim, xylim, type = "n", xlab = xlab, ylab = ylab)
  lines(xylim, xylim, lty = 2, col = "gray40")
  out6 <- MG16S("Marcellus_Shale")
  dy <- c(-0.004, 0, 0, -0.005, -0.005)
  dx <- c(0.007, -0.008, 0.009, 0, 0)
  text(out6$metric_MG + dx, out6$metric_16S + dy, out6$Sample)
  # Add legend
  legend("topleft", c("Injected", "Flowback", "Produced"), pch = 23, pt.bg = c(4, 8, 2))
  # Add title and figure label
  title("Marcellus Shale Fluids", font.main = 1, cex.main = 1.1)
  label.figure("B", cex = 1.5, font = 2, xfrac = 0.04, yfrac = 0.96)

  # Use semi-transparent colors 20220122
  c1 <- addalpha(1, "80")
  c2 <- addalpha(2, "b0")
  c4 <- addalpha(4, "b0")
  c5 <- addalpha(5, "b0")
  c6 <- addalpha(6, "b0")
  c8 <- addalpha(8, "b0")

  # Start plot C
  xlab <- quote(italic(Z)[C]~"from shotgun metagenome")
  plot(xylim, xylim, type = "n", xlab = xlab, ylab = ylab)
  lines(xylim, xylim, lty = 2, col = "gray40")
  out7 <- MG16S("Manus_Basin", cex = 1.4)
  # Plot sample names and O2 concentrations (from Figure S5 of Meier et al., 2017)
  dx <- c(0.010, 0.009, -0.012, -0.008, 0.007)
  dy <- c(0, 0, 0.012, 0, -0.0018)
  samptxt <- substr(out7$Sample, 7, 13)
  text(out7$metric_MG + dx, out7$metric_16S - 0.004 + dy, samptxt, cex = 0.9)
  O2txt <- as.expression(c(quote(0.07~"mM"~O[2]), quote(0.14~"mM"~O[2]), quote(0.17~"mM"~O[2]), quote("ND"~O[2]), quote(0.2~"mM"~O[2])))
  text(out7$metric_MG + dx, out7$metric_16S - 0.009 + dy, O2txt, cex = 0.9)

  # Add Black Sea 20220115
  out8 <- MG16S("Black_Sea", cex = 0.9)
  # Add legends
  legend("topleft", c("Depth < 100 m", "Depth >= 100 m"), pch = c(24, 25), pt.bg = c(c4, c2), col = c1, pt.cex = 0.9, title = "Black Sea")
  legend <- as.expression(c(quote(italic(T)~"< 10 \u00B0C"), quote("10 \u00B0C <"~italic(T)~"< 50 \u00B0C"), quote(italic(T)~"> 50 \u00B0C")))
  legend("bottomright", legend = legend, pch = 22, pt.bg = c(c4, "white", c2), col = c1, pt.cex = 1.3, title = "Manus Basin")
  # Add title and figure label
  title("Manus Basin Vents and Black Sea", font.main = 1, cex.main = 1.1)
  label.figure("C", cex = 1.5, font = 2, xfrac = 0.04, yfrac = 0.96)

  ### Panels D-E: Plot ZC of 16S vs metagenomes for datasets used in Tax4Fun/PICRUSt papers 20211215

  # Start plot D
  xlab <- quote(italic(Z)[C]~"from shotgun metagenome")
  plot(xlimHMP, xylim, type = "n", xlab = xlab, ylab = ylab)
  lines(xylim, xylim, lty = 2, col = "gray40")
  out9 <- MG16S("HMP")
  legend("topleft", c("Skin", "Nasal cavity", "Oral cavity", "GI tract", "UG tract"), pch = c(24, 25, 21, 24, 25), pt.bg = c(c5, c4, c8, c2, c6), col = c1)
  title("Human Microbiome Project", font.main = 1, cex.main = 1.1)
  label.figure("D", cex = 1.5, font = 2, xfrac = 0.04, yfrac = 0.96)

  # Start plot E
  xlab <- quote(italic(Z)[C]~"from shotgun metagenome")
  plot(xylim, xylim, type = "n", xlab = xlab, ylab = ylab)
  lines(xylim, xylim, lty = 2, col = "gray40")
  out10 <- MG16S("Guts")
  out11 <- MG16S("Soils")
  legend("topleft", c("Soils", "Guts"), pch = 21, pt.bg = c(c4, c2), col = c1)
  title("Soils and Mammalian Guts", font.main = 1, cex.main = 1.1)
  label.figure("E", cex = 1.5, font = 2, xfrac = 0.04, yfrac = 0.96)

  if(pdf) dev.off()

  # Return values for Supplemental Table 20220125
  out <- list(out1, out2, out3, out4, out5, out6, out7, out8, out9, out10, out11)
  out <- data.frame(
    name = unlist(lapply(out, "[[", "which")),
    sample = unlist(lapply(out, "[[", "Sample")),
    ID_MG = unlist(lapply(out, "[[", "Metagenome")),
    ID_16S = unlist(lapply(out, "[[", "Amplicon")),
    ZC_MG = round(unlist(lapply(out, "[[", "metric_MG")), 6),
    ZC_16S = round(unlist(lapply(out, "[[", "metric_16S")), 6)
  )
  invisible(out)

}


# RefSeq and 16S rRNA data processing outline 20220104
geo16S_S1 <- function(pdf = FALSE) {

  if(pdf) pdf("geo16S_S1.pdf", width = 17.1, height = 5)

  par(mar = c(0.1, 0.1, 0.1, 0.1))
  plot(c(1, 10.2), c(1.9, 10.1), type = "n", axes = FALSE, xlab = "", ylab = "", xaxs = "i", yaxs = "i")
  text(1, 10, "Reference Sequence Processing", adj = c(0, 1), font = 2)
  text(1, 9, "1. Use taxids classified at species level", adj = 0)
  text(1, 8, "2. Sum AA of all sequences for each species", adj = 0)
  text(1, 7, "3. Divide by number of sequences to get mean AA for each species (AA[species])", adj = 0)
  text(1, 6, "4. Calculate mean AA[genus] from all AA[species] in each genus", adj = 0)
  text(1, 5, "5. Calculate mean AA[family] from all AA[species] in each family", adj = 0)
  text(1, 4, "6. Calculate mean AA[order] from all AA[species] in each order", adj = 0)
  text(1, 3, "7. Calculate mean AA[class] from all AA[species] in each class", adj = 0)
  text(1, 2, "8. Calculate mean AA[phylum] from all AA[species] in each phylum", adj = 0)
#  text(1, 1, "9. Calculate and record ZC and nH2O for all taxonomic groups at each level", adj = 0)

  text(3.35, 10, "(Data Source)", adj = c(0, 1), font = 2)
  text(3.35, 9, "(NCBI taxonomy files)", adj = 0, font = 2)
  text(3.35, 8, "(RefSeq protein sequences)", adj = 0, font = 2)

  text(5, 10, "Amino Acid Compositions\nfor Taxonomic Groups\n(Reference Proteomes)", adj = c(0, 1), font = 2)
  # Read file with precomputed metrics for taxa at different ranks
  metrics <- read.csv(system.file("extdata/chem16S/taxon_metrics.csv", package = "JMDplots"))
  ranks <- c("superkingdom", "phylum", "class", "order", "family", "genus")
  plural <- c("superkingdoms (*)", "phyla", "classes", "orders", "families", "genera")
  for(irank in 1:6) {
    nrank <- sum(metrics$rank == ranks[irank])
    text(5, 9-irank, paste(nrank, plural[irank]), adj = 0)
  }
  lines(c(3.87, 3.87), c(1.8, 6.2))
  lines(c(4.9, 4.9), c(7.2, 2.8))
  arrows(3.9, 4, 4.88, 5.00)
  lines(c(5.68, 5.68), c(7.2, 2.8))
  arrows(5.71, 5.00, 6.43, 5.00)

  text(6.5, 10, "16S rRNA Classification and Analysis", adj = c(0, 1), font = 2)
  text(6.5, 9, "1. Run RDP Classifier on public 16S rRNA datasets (see Tables 1 and S1)", adj = 0)
  text(6.5, 8, "2. Assemble counts of lowest-level classifications (genus to phylum)", adj = 0)
  text(6.5, 7, "3. Mapping step 1: Manually convert some RDP names to NCBI names (see Methods)", adj = 0)
  text(6.5, 6, "4. Mapping step 2: Automatic text match between RDP and NCBI names", adj = 0)
  text(6.5, 5, "5. Multiply classification counts by reference proteomes", adj = 0)
  text(6.5, 4, "6. Divide by number of classifications to get estimated community proteome (AA[community])", adj = 0)
  text(6.5, 3, "7. Use AA[community] to calculate ZC and nH2O (see Dick et al., 2020)", adj = 0)
  text(6.5, 2, "8. Visualize data", adj = 0)

  if(pdf) dev.off()

}

# Scatterplots of ZC and nH2O for genera vs higher taxonomic levels 20211130
geo16S_S2 <- function(pdf = FALSE) {
  if(pdf) pdf("geo16S_S2.pdf", width = 12, height = 6)

  # Read file with precomputed metrics for taxa at different ranks
  metrics <- read.csv(system.file("extdata/chem16S/taxon_metrics.csv", package = "JMDplots"))
  # Get tables of all taxonomic names and amino acid compositions
  names <- read.csv(system.file("extdata/refseq/taxid_names.csv.xz", package = "JMDplots"))
  # Take out viruses
  ivirus <- names$superkingdom == "Viruses"
  ivirus[is.na(ivirus)] <- TRUE
  names <- names[!ivirus, ]
  # Take out NA genera and dereplicate lineages 20220107
  names <- names[!is.na(names$genus), ]
  gfocp <- apply(names[, c("genus", "family", "order", "class", "phylum")], 1, paste, collapse = " ")
  names <- names[!duplicated(gfocp), ]

  # Initialize list of ZC and nH2O values for taxa in each rank
  NAvec <- rep(NA, nrow(names))
  nH2O <- ZC <- list(genus = NAvec, family = NAvec, order = NAvec, class = NAvec, phylum = NAvec, superkingdom = NAvec)

  for(rank in c("genus", "family", "order", "class", "phylum", "superkingdom")) {
    icol <- match(rank, colnames(names))
    isrank <- !is.na(names[, icol])
    # Use a copy of the metrics table with blanked-out names of taxa not at this rank
    thismet <- metrics
    thismet$group[metrics$rank != rank] <- ""
    # Lookup these taxa in the metrics table
    imetrics <- match(names[, icol][isrank], thismet$group)
    # Make sure all taxa are in this rank
    stopifnot(all(na.omit(unique(thismet$rank[imetrics])) %in% rank))
    # Store ZC and nH2O values
    ZC[[rank]][isrank] <- thismet$ZC[imetrics]
    names(ZC[[rank]])[isrank] <- names[, icol][isrank]
    nH2O[[rank]][isrank] <- thismet$nH2O[imetrics]
    names(nH2O[[rank]])[isrank] <- names[, icol][isrank]
  }

#  # Print number of genera
#  ngenus <- length(names(ZC[["genus"]]))
#  print(paste("There are", ngenus, "genera"))

  # Make plots
  par(mfrow = c(2, 5))
  par(mgp = c(2.9, 1, 0))
  # Identify superkingdoms
  iArc <- names(ZC[["superkingdom"]]) == "Archaea"
  iBac <- names(ZC[["superkingdom"]]) == "Bacteria"

  cex <- 2
  cex.an <- 1
  # Outer loop: ZC or nH2O
  for(metric in c("ZC", "nH2O")) {
    
    # Inner loop: higher-level ranks
    ranks <- c("family", "order", "class", "phylum", "superkingdom")
    for(rank in ranks) {
      if(metric == "ZC") { x <- ZC[["genus"]]; y <- ZC[[rank]] }
      if(metric == "nH2O") { x <- nH2O[["genus"]]; y <- nH2O[[rank]] }
      xylim <- range(na.omit(x))
      # Start plot and add 1:1 line
      if(metric == "ZC") plot(xylim, xylim, type = "n", xlab = quote("genus"~italic(Z)[C]), ylab = bquote(.(rank)~italic(Z)[C]))
      if(metric == "nH2O") plot(xylim, xylim, type = "n", xlab = quote("genus"~italic(n)[H[2]*O]), ylab = bquote(.(rank)~italic(n)[H[2]*O]))
      lines(xylim, xylim, lty = 2, col = "gray40")
      # Add points: Bacteria then Archaea
      points(x[iBac], y[iBac], pch = ".", col = "steelblue3", cex = cex)
      points(x[iArc], y[iArc], pch = ".", col = 2, cex = cex)
      # Add points and linear regression
      thislm <- lm(y ~ x)
      lines(xylim, predict(thislm, data.frame(x = xylim)), col = "#00000080")
      # Show R2 and slope
      R2 <- summary(thislm)$r.squared
      R2txt <- bquote(italic(R)^2 == .(formatC(R2, digits = 3, format = "f")))
      legend("topleft", legend = R2txt, bty = "n")
      Slope <- coef(thislm)["x"]
      Slopetxt <- paste("slope =", formatC(Slope, digits = 3, format = "f"))
      legend("bottomright", legend = Slopetxt, bty = "n")
      if(metric == "ZC") {
        # Add title with number of taxa at this rank
        plural <- switch(rank, "genus" = "genera", "family" = "families", "order" = "orders",
          "class" = "classes", "phylum" = "phyla", "superkingdom" = "superkingdoms")
        ntaxa <- length(unique(names(ZC[[rank]])))
        title(paste(ntaxa, plural), font.main = 1)
      }
    }
  }

  if(pdf) dev.off()

}

# nH2O-ZC plots for major phyla and their genera 20220114
geo16S_S3 <- function(pdf = FALSE) {

  if(pdf) pdf("geo16S_S3.pdf", width = 13, height = 11)

  mat <- matrix(c(1,1,2,2, 0,3,3,0), nrow = 2, byrow = TRUE)
  layout(mat)
  par(mar = c(4, 4, 2, 1))
  par(cex= 1.2)
  xlim <- c(-0.3, 0)
  ylim <- c(-0.85, -0.65)

  metrics <- read.csv(system.file("extdata/chem16S/taxon_metrics.csv", package = "JMDplots"))
  names <- read.csv(system.file("extdata/refseq/taxid_names.csv.xz", package = "JMDplots"))
  # Only keep taxa with non-NA genus and phylum
  names <- names[!(is.na(names$genus) | is.na(names$phylum)), ]

  # Plot phylum connected to all genera
  plotit <- function(phylum, col = 1) {
    genera <- unique(names$genus[names$phylum == phylum])
    # Lookup phylum and genera in metrics table
    iphylum <- match(phylum, metrics$group)
    igenera <- match(genera, metrics$group)
    # Use thicker lines and less transparency for phyla with fewer genera
    if(length(genera) < 50) lwd <- 2 else lwd <- 1
    if(length(genera) < 50) alpha <- "a0" else alpha <- "50"
    newcol <- addalpha(col, alpha)
    for(i in igenera) lines(c(metrics$ZC[iphylum], metrics$ZC[i]), c(metrics$nH2O[iphylum], metrics$nH2O[i]), col = newcol, lwd = lwd)
    lwd
  }

  # Get "majorcellular" phyla used in Figure 1
  phyla <- metrics[metrics$rank == "phylum" & metrics$parent != "Viruses", ]
  phyla <- phyla[phyla$ntaxa > 60, ]
  phyla <- phyla[order(phyla$ntaxa, decreasing = TRUE), ]
  # Swap Chloroflexi and Crenarchaeota so latter doesn't have same color as Euryarchaeota 20210527
  phyla[14:15, ] <- phyla[15:14, ]

  # Plot first 8 phyla
  plot(xlim, ylim, xlab = canprot::cplab$ZC, ylab = canprot::cplab$nH2O, type = "n", xaxs = "i", yaxs = "i")
  lwd <- lapply(1:8, function(i) {plotit(phyla$group[i], col = i)} )
  legend("topright", phyla$group[1:8], col = 1:8, lwd = lwd, cex = 0.8, bg = "white")
  label.figure("A", font = 2, cex = 1.5, xfrac = 0.03)

  # Plot second 8 phyla
  plot(xlim, ylim, xlab = canprot::cplab$ZC, ylab = canprot::cplab$nH2O, type = "n", xaxs = "i", yaxs = "i")
  lwd <- lapply(9:16, function(i) {plotit(phyla$group[i], col = i)} )
  legend("topright", phyla$group[9:16], col = 1:8, lwd = lwd, cex = 0.8, bg = "white")
  label.figure("B", font = 2, cex = 1.5, xfrac = 0.03)

  # Plot 8 more phyla 20220126
  phyla <- metrics[metrics$rank == "phylum" & metrics$parent != "Viruses", ]
  phyla <- phyla[phyla$ntaxa >= 18 & phyla$ntaxa <= 60, ]
  phyla <- phyla[order(phyla$ntaxa, decreasing = TRUE), ]
  plot(xlim, ylim, xlab = canprot::cplab$ZC, ylab = canprot::cplab$nH2O, type = "n", xaxs = "i", yaxs = "i")
  lwd <- lapply(1:8, function(i) {plotit(phyla$group[i], col = i)} )
  legend <- phyla$group[1:8]
  legend <- gsub("Candidatus Thermoplasmatota", "Candidatus", legend)
  legend <- c(legend, "Thermoplasmatota")
  legend("topright", legend, col = c(1:8, NA), lwd = lwd, cex = 0.8, bg = "white")
  label.figure("C", font = 2, cex = 1.5, xfrac = 0.03)

  if(pdf) dev.off()

}

# Venn diagrams for phylum and genus names in the RefSeq (NCBI), RDP, and SILVA taxonomies 20220117
geo16S_S4 <- function(pdf = FALSE) {

  # Read lists of RDP and SILVA names
  datadir <- system.file("extdata/geo16S/taxonomy", package = "JMDplots")
  RDPphyla <- readLines(file.path(datadir, "RDPphyla.txt"))
  RDPgenera <- readLines(file.path(datadir, "RDPgenera.txt"))
  SILVAphyla <- readLines(file.path(datadir, "SILVAphyla.txt"))
  SILVAgenera <- readLines(file.path(datadir, "SILVAgenera.txt"))
  # Read NCBI names
  NCBI <- read.csv(system.file("extdata/refseq/taxid_names.csv.xz", package = "JMDplots"))
  NCBIphyla <- unique(na.omit(NCBI$phylum))
  NCBIgenera <- unique(na.omit(NCBI$genus))

  pl <- list()

  # Make Venn diagrams
  fillsRDP <- list(fill = c("transparent", "slategray3"))
  fillsSILVA <- list(fill = c("transparent", "lightpink2"))

  A <- length(setdiff(NCBIphyla, RDPphyla))
  B <- length(setdiff(RDPphyla, NCBIphyla))
  AB <- length(intersect(NCBIphyla, RDPphyla))
  pl[[1]] <- plot( euler(c(A = A, B = B, "A&B" = AB)), legend = list(labels = c("RefSeq (NCBI)", "RDP"), side = "bottom"), quantities = TRUE, fills = fillsRDP )

  A <- length(setdiff(NCBIphyla, SILVAphyla))
  B <- length(setdiff(SILVAphyla, NCBIphyla))
  AB <- length(intersect(NCBIphyla, SILVAphyla))
  pl[[2]] <- plot( euler(c(A = A, B = B, "A&B" = AB)), legend = list(labels = c("RefSeq (NCBI)", "SILVA"), side = "bottom"), quantities = TRUE, fills = fillsSILVA )

  A <- length(setdiff(NCBIgenera, RDPgenera))
  B <- length(setdiff(RDPgenera, NCBIgenera))
  AB <- length(intersect(NCBIgenera, RDPgenera))
  pl[[3]] <- plot( euler(c(A = A, B = B, "A&B" = AB)), legend = list(labels = c("RefSeq (NCBI)", "RDP"), side = "bottom"), quantities = TRUE, fills = fillsRDP )

  A <- length(setdiff(NCBIgenera, SILVAgenera))
  B <- length(setdiff(SILVAgenera, NCBIgenera))
  AB <- length(intersect(NCBIgenera, SILVAgenera))
  pl[[4]] <- plot( euler(c(A = A, B = B, "A&B" = AB)), legend = list(labels = c("RefSeq (NCBI)", "SILVA"), side = "bottom"), quantities = TRUE, fills = fillsSILVA )

  # Make the plot
  if(pdf) pdf("geo16S_S4.pdf", width = 8, height = 6)

  lg <- gridExtra::tableGrob(c("A. phyla", "B. genera"), theme = gridExtra::ttheme_minimal(base_size = 18))
  rg <- gridExtra::arrangeGrob(grobs = pl, ncol = 2)
  grid.draw(cbind(lg, rg, size = "last"))
  if(pdf) dev.off()

  # Save RDP and SILVA names not in NCBI for Supplementary Table
  out <- list(
    RDP_phylum_name_not_in_RefSeq = sort(setdiff(RDPphyla, NCBIphyla)),
    RDP_genus_name_not_in_RefSeq = sort(setdiff(RDPgenera, NCBIgenera)),
    SILVA_phylum_name_not_in_RefSeq = sort(setdiff(SILVAphyla, NCBIphyla)),
    SILVA_genus_name_not_in_RefSeq = sort(setdiff(SILVAgenera, NCBIgenera))
  )
  # Pad each list element to the same length with NA
  len <- max(sapply(out, length))
  out <- lapply(out, "[", 1:len)
  out <- do.call(cbind, out) 
  invisible(out)

}

# Correlations between ZC estimated from MG and 16S 20220112
geo16S_S5 <- function(pdf = FALSE, H2O = FALSE) {

  if(pdf) pdf("geo16S_S5.pdf", width = 10, height = 6)
  par(mfrow = c(2, 3))
  par(mar = c(4, 4, 4, 1), mgp = c(2.8, 1, 0))
  if(H2O) xylim <- c(-0.8, -0.70) else xylim <- c(-0.22, -0.12)

  # Plot regression line and R2
  lmfun <- function(dat, xylim) {
    x <- unlist(lapply(dat, "[[", "metric_MG"))
    y <- unlist(lapply(dat, "[[", "metric_16S"))
    # Add points and linear regression
    thislm <- lm(y ~ x)
    lines(xylim, predict(thislm, data.frame(x = xylim)), col = "#00000080")
    # Add R2 value
    R2 <- summary(thislm)$r.squared
    R2txt <- bquote(italic(R)^2 == .(formatC(R2, digits = 3, format = "f")))
    legend("topleft", legend = R2txt, bty = "n")
  }

  for(row in 1:2) {

    if(row == 2) lineage = "genus" else lineage <- NULL
    # Not used: show the estimates using only phylum-level classifications
    if(row == 3) lowest.level <- "phylum" else lowest.level <- NULL

    ## Panel 1: Various Environments
    if(H2O) {
      xlab <- quote(italic(n)[H[2]*O]~"from shotgun metagenome or metatranscriptome")
      ylab <- quote(italic(n)[H[2]*O]~"estimated from 16S rRNA")
    } else {
      xlab <- quote(italic(Z)[C]~"from shotgun metagenome or metatranscriptome")
      ylab <- quote(italic(Z)[C]~"estimated from 16S rRNA")
    }
    plot(xylim, xylim, type = "n", xlab = xlab, ylab = ylab)
    lines(xylim, xylim, lty = 2, col = "gray40")
    cex <- c(1, 1, 1, 1.4, 1, 1, 1.4, 0.9)
    dat <- mapply(MG16S, which = c("Guerrero_Negro", "ETNP_MG", "ETNP_MT", "Bison_Pool", "Mono_Lake", "Marcellus_Shale", "Manus_Basin", "Black_Sea"),
                   cex = cex, plot.lines = FALSE, MoreArgs = list(lowest.level = lowest.level, lineage = lineage, H2O = H2O), SIMPLIFY = FALSE)
    lmfun(dat, xylim)
    # Add title and figure label
    title("Various Environments + Marcellus + Manus + Black Sea     ", font.main = 1, cex.main = 1.1, line = 1)

    ## Panel 2: Human Microbiome Project
    if(H2O) xlab <- quote(italic(n)[H[2]*O]~"from shotgun metagenome") else xlab <- quote(italic(Z)[C]~"from shotgun metagenome")
    if(H2O) xylimHMP <- c(-0.95, -0.70) else xylimHMP <- xylim
    plot(xylimHMP, xylimHMP, type = "n", xlab = xlab, ylab = ylab)
    lines(xylimHMP, xylimHMP, lty = 2, col = "gray40")
    dat <- lapply("HMP", MG16S, plot.lines = FALSE, lowest.level = lowest.level, lineage = lineage, rm.outliers = TRUE, H2O = H2O)
    lmfun(dat, xylimHMP)
    title("Human Microbiome Project", font.main = 1, cex.main = 1.1, line = 1)
    par(xpd = NA)
    if(row == 1) title("A. Lowest-level classifications (genus to phylum)", line = 2.5)
    if(row == 2) title("B. Genus-level classifications only (higher-level classifications discarded)", line = 2.5)
    if(row == 3) title("C. Phylum-level classifications (lineages truncated at phylum)", line = 2.5)
    par(xpd = FALSE)

    ## Panel 3: Soils and Mammamlian Guts
    plot(xylim, xylim, type = "n", xlab = xlab, ylab = ylab)
    lines(xylim, xylim, lty = 2, col = "gray40")
    dat <- lapply(c("Guts", "Soils"), MG16S, plot.lines = FALSE, lowest.level = lowest.level, lineage = lineage, H2O = H2O)
    lmfun(dat, xylim)
    title("Soils and Mammalian Guts", font.main = 1, cex.main = 1.1, line = 1)

  }

  if(pdf) dev.off()

}

# Correlation of ZC with GC content of metagenomic and 16S amplicon reads 20220126
geo16S_S6 <- function(pdf = FALSE) {

  if(pdf) pdf("geo16S_S6.pdf", width = 8, height = 7)
  par(mfrow = c(2, 2))
  par(mar = c(4, 4, 2, 1))
  par(mgp = c(2.8, 1, 0))

  # Change these to extract specific parts of the taxonomy
  lowest.level <- NULL
  lineage <- NULL

  # For MG datasets analyzed in this study (others are from gradox paper)
  ARASTdir <- system.file("extdata/geo16S/ARAST", package = "JMDplots")
  # Read data for paired metagenomes and amplicon sequences expanded from Tax4Fun paper (Aßhauer et al., 2015)
  AWDM15file <- system.file("extdata/geo16S/AWDM15.csv", package = "JMDplots")
  AWDM15 <- read.csv(AWDM15file)

  # Use semi-transparent colors 20220122
  c1 <- addalpha(1, "80")
  c2 <- addalpha(2, "b0")
  c4 <- addalpha(4, "b0")
  c5 <- addalpha(5, "b0")
  c6 <- addalpha(6, "b0")
  c8 <- addalpha(8, "b0")

  for(name in c("Marcellus", "HMP")) {

    if(name == "HMP") {
      # HMP 16S
      mdat <- getmdat("HMP12")
      met <- getmetrics("HMP12", mdat = mdat, lowest.level = lowest.level, lineage = lineage)
      # HMP metagenomes
      aa <- read.csv(file.path(ARASTdir, "HMP_AA.csv"))
      # Put data in same order
      dat <- AWDM15[AWDM15$Name == "HMP", ]
      imet <- match(dat$Amplicon, met$Run)
      met <- met[imet, ]
      iaa <- match(dat$Metagenome, aa$protein)
      aa <- aa[iaa, ]
      # Make sure the 16S and metagenomes are paired correctly
      stopifnot(all(na.omit(met$Run == dat$Amplicon)))
      stopifnot(all(aa$protein == dat$Metagenome))
      # Don't plot MG with low numbers of protein fragments 20220122
      ilow <- aa$chains < 50000
      if(any(ilow)) {
        mdat <- mdat[!ilow, ]
        met <- met[!ilow, ]
        aa <- aa[!ilow, ]
        dat <- dat[!ilow, ]
      }
      # Get ZC
      ZC_16S <- met$ZC
      ZC_MG <- ZCAA(aa)
      # Get GC
      GC_16S <- dat$GC_16S
      GC_MG <- dat$GC_MG
      # Colors: blue (Skin), green (Nasal cavity), gray (Oral cavity), red (GI tract), magenta (UG tract)
      col <- sapply(mdat$"Body site", switch, "Skin" = c5, "Nasal cavity" = c4, "Oral cavity" = c8, "GI tract" = c2, "UG tract" = c6)
      # Symbols: up triangle (skin, GI tract), circle (Oral cavity), down triangle (Nasal cavity, UG tract)
      pch <- sapply(mdat$"Body site", switch, "Skin" = 24, "Nasal cavity" = 25, "Oral cavity" = 21, "GI tract" = 24, "UG tract" = 25)
      legend.x <- "bottomright"
    }

    if(name == "Marcellus") {
      ## Marcellus metagenomes (Daly et al., 2016) 20211218
      aa <- read.csv(file.path(ARASTdir, "Marcellus_Shale_AA.csv"))
      ZC_MG <- ZCAA(aa)
      # Marcellus 16S (Cluff et al., 2014)
      mdat <- getmdat("CHM+14")
      dat_16S <- getmetrics("CHM+14", mdat = mdat, lowest.level = lowest.level, lineage = lineage)
      # Time points: input, T7, T13, T82, T328
      Sample = paste("Day", c(0, 7, 13, 82, 328))
      # List run IDs here
      Metagenome = c("SRR3111417", "SRR3111625", "SRR3111724", "SRR3111729", "SRR3111737")
      Amplicon = c("SRR1184016", "SRR1184060", "SRR1184062", "SRR1184081", "SRR1184083")
      # GC content from https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR******* 20220126
      GC_MG <- c(58.4, 37.7, 57, 39.6, 53.1)
      GC_16S <- c(53.9, 55, 54.8, 51.8, 53.2)
      # Make sure metagenomes are in correct order
      stopifnot(all(aa$protein == Metagenome))
      # Get 16S runs corresponding to metagenomes
      idat <- match(Amplicon, dat_16S$Run)
      dat_16S <- dat_16S[idat, ]
      ZC_16S <- dat_16S$ZC
      # Assign colors: open circle for injected fluid, gray for flowback, red for produced
      mdat <- mdat[idat, ]
      col <- rep(4, nrow(mdat))
      col[mdat$type == "flowback fluid"] <- 8
      col[mdat$type == "produced fluid"] <- 2
      pch <- 23
      legend.x <- "topleft"
    }

  #  ## Plot A: GC 16S vs GC MG
  #  plot(range(GC_MG), range(GC_16S), xlab = "GC content of shotgun metagenome reads (%)", ylab = "GC content of 16S amplicon reads (%)", type = "n")
  #  points(GC_MG, GC_16S, pch = pch, bg = col, col = c1)

    lmfun <- function(x, y) {
      # Add points and linear regression
      thislm <- lm(y ~ x)
      lines(range(x), predict(thislm, data.frame(x = range(x))), col = "#00000080")
      # Show R2 and slope
      R2 <- summary(thislm)$r.squared
      R2txt <- bquote(italic(R)^2 == .(formatC(R2, digits = 3, format = "f")))
      R2txt
    }

    ## Plot A: ZC MG vs GC MG
    plot(range(GC_MG), range(ZC_MG, na.rm = TRUE), xlab = "%GC content - metagenome", ylab = quote("Protein"~italic(Z)[C]~"- metagenome"), type = "n")
    points(GC_MG, ZC_MG, pch = pch, bg = col, col = c1)
    R2txt <- lmfun(GC_MG, ZC_MG)
    legend(legend.x, legend = R2txt, bty = "n")

    if(name == "Marcellus") legend("bottomright", c("Injected", "Flowback", "Produced"), pch = 23, pt.bg = c(4, 8, 2))

    ## Plot B: ZC 16S vs GC 16S
    plot(range(GC_16S), range(ZC_16S, na.rm = TRUE), xlab = "%GC content - 16S amplicon", ylab = quote("Protein"~italic(Z)[C]~"- 16S amplicon"), type = "n")
    points(GC_16S, ZC_16S, pch = pch, bg = col, col = c1)
    R2txt <- lmfun(GC_16S, ZC_16S)
    legend(legend.x, legend = R2txt, bty = "n")

    if(name == "HMP") legend("topleft", c("Skin", "Nasal cavity", "Oral cavity", "GI tract", "UG tract"), pch = c(24, 25, 21, 24, 25), pt.bg = c(c5, c4, c8, c2, c6), col = c1)

    # Add title
    par(xpd = NA)
    if(name == "Marcellus") title("A. Marcellus Shale                                                                                       ")
    if(name == "HMP") title("B. Human Microbiome Project                                                                                       ")
    par(xpd = FALSE)

  }

  if(pdf) dev.off()

}

### UNEXPORTED FUNCTION ###
addalpha <- function(col, alpha) {
  # Add transparency to given color
  x <- col2rgb(col)
  newcol <- rgb(x[1], x[2], x[3], maxColorValue = 255)
  newcol <- paste0(newcol, alpha)
  newcol
}
