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

  # Return Source Data 20210831
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

  p2 <- plotmet("SVH+19", title = FALSE, points = FALSE, ylim = c(-0.757, -0.737))
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

  p4 <- plotmet("MPB+17", title = FALSE, points = FALSE, ylim = c(-0.777, -0.725))
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

  p7 <- plotmet("HCW+13", title = FALSE, points = FALSE, ylim = c(-0.764, -0.7545))
#  title("Guerrero Negro mat layers\nHarris et al., 2013", font.main = 1)
  title("Guerrero Negro mat layers", font.main = 1)
  addhull(p7$ZC, p7$nH2O, "tan1", TRUE, lty = 2)
  pointfun(p7)
  text(c(-0.1518, -0.1577, -0.1578), c(-0.7547, -0.7602, -0.7637), c("0-1 mm", "1-2 mm", "2-3 mm"))
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

  # Return Source Data 20210831
  out <- rbind(p1, p2, p3, p4, p5, p6, p7, p8)
  out$ZC <- round(out$ZC, 6)
  out$nH2O <- round(out$nH2O, 6)
  invisible(out)

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
    # Setup plot metrics
    par(mgp = c(1.8, 0.5, 0), mar = c(3, 3, 3, 1))
  }

  # Identify datasets to plot
  study <- c("SVH+19", "MZG+20", "MZG+20", "GBL+15", "GBL+15",
             "HXZ+20", "HXZ+20", "HXZ+20", "HXZ+20",
             "BCA+21", "BCA+21", "BCA+21", "BCA+21")
  column <- c("study", "lake", "lake", "size", "size",
              "Station", "Station", "Station", "Station",
              "Month", "Month", "Month", "Month")
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

  # Make object to hold Source Data
  out <- list()
  # Loop over studies
  for(i in 1:length(study)) {
    # ZC range for plots
    if(study[i] == "BCA+21") ZClim <- c(-0.180, -0.145) else ZClim <- c(-0.172, -0.140)
    # Get the metadata and compositional metrics for this study
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

    if(plot.it) {
      # Reverse y-axis (depth)
      ylim <- rev(range(depth))
      # Visualize deeper O2 concentrations in Ursu Lake
      if(study[i] == "BCA+21") ylim <- c(11, 0)
      if(study[i] == "SVH+19") {
        # Plot 1000 and 2000 m samples closer to the others 20210608
        ylim <- c(700, 50)
        depth[match(c(1000, 2000), depth)] <- c(600, 700)
      }
      # Determine whether the title has changed
      newplot <- TRUE
      if(i > 1) if(titlesub[i]==titlesub[i-1]) newplot <- FALSE
      if(newplot) {
        if(study[i] == "SVH+19") {
          plot(ZC, depth, xlim = ZClim, ylim = ylim, xlab = axis.label("ZC"), ylab = "Depth (m)", type = "b", yaxt = "n")
          axis(2, at = seq(100, 700, 100), labels = c(100, 200, 300, 400, 500, 1000, 2000), gap.axis = 0)
          # Plot y-axis break 20210715
          par(xpd = NA)
          rect(-0.174, 557, -0.173, 542, col = "white", border = NA)
          text(-0.1734, 542, "/", srt = 90)
          text(-0.1734, 556, "/", srt = 90)
          par(xpd = FALSE)
        } else plot(ZC, depth, xlim = ZClim, ylim = ylim, xlab = axis.label("ZC"), ylab = "Depth (m)", type = "b")
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
          if(study[i] == "BCA+21") xlim <- c(0, 25) else xlim <- c(0, 220)
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
          text(44, 76, "0.2-\n1.6 \u00B5m")
          text(128, 138, "1.6-\n30 \u00B5m")
        }
        # Restore xlim for plotting ZC
        par(new = TRUE)
        plot(0, 0, type = "n", axes = FALSE, xlab = "", ylab = "", xlim = ZClim, ylim = ylim)
      }

      # Don't repeat Blue Hole in Source Data (there are loops for O2 and NO3-/NO2-; just use one)
      if(substr(title[i], 1, 1) == " ") next
      # Assemble the Source Data
      sd <- alldat[, c("study", "name", "Run", "sample", "depth")]
      sd <- cbind(sd, "O2 (umol kg-1)" = NA, "O2 (umol L-1)" = NA, "O2 (mg L-1)" = NA, "NO3- (umol L-1)" = NA, "NO2- (umol L-1)" = NA, ZC = NA, nH2O = NA)
      # Names of the columns with chemical concentrations
      cnames <- c("O2 (umol kg-1)", "O2 (umol L-1)", "O2 (mg L-1)", "NO3- (umol L-1)", "NO2- (umol L-1)")
      for(cname in cnames) if(cname %in% colnames(alldat)) sd[, cname] <- alldat[, cname]
      # Put ZC and nH2O values in correct place
      metrics <- metrics[imet, ]
      metrics <- metrics[!is.na(metrics$Run), ]
      isd <- match(metrics$Run, sd$Run)
      sd$ZC[isd] <- metrics$ZC
      sd$nH2O[isd] <- metrics$nH2O

      # Replace NA name with study name
      sd.name <- na.omit(sd$name)[1]
      sd$name[is.na(sd$name)] <- sd.name
      # Use NA instead of "" for missing Run and sample name
      sd$Run[sd$Run == ""] <- NA
      sd$sample[sd$sample == ""] <- NA
      # Place data frame into output list
      out[[i]] <- sd

    } # end if(plot.it)

  } # end loop

  if(plot.it) {
    if(pdf) dev.off()
  }

  # Return Source Data 20210831
  out <- do.call(rbind, out)
  out$ZC <- round(out$ZC, 6)
  out$nH2O <- round(out$nH2O, 6)
  invisible(out)

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
  p <- groupmet("MPB+17", "ZC", "domain", pch1 = 21, pch2 = 23, ylim = c(-0.23, -0.13), xlim = c(0, 100),
    xadj = c(Bacteria = 1), yadj = c(Bacteria = -5),
    mdat = mdat, RDP = RDP, map = map
  )
  title("Manus Basin")
  text(40, -0.161, "   < 50 \u00B0C", font = 2)
  text(10, -0.161, "T", font = 4)
  text(40, -0.204, "   > 50 \u00B0C", font = 2)
  text(10, -0.204, "T", font = 4)

  gg <- groupmet("MPB+17", "ZC", "phylum", pch1 = 21, pch2 = 23, ylim = c(-0.23, -0.13),
    xadj = c(Proteobacteria = 1, Bacteroidetes = -0.18, Campilobacterota = 0.42),
    yadj = c(Bacteroidetes = 1, Campilobacterota = -3),
    mdat = mdat, RDP = RDP, map = map
  )
  # Calculate total percentage of community represented by these taxa
  p <- sum(gg$Pboth[!is.na(gg$X1) & !is.na(gg$X2)])
  title(paste0("Phylum (", round(p), "% of total)"))

  gg <- groupmet("MPB+17", "ZC", "class", pch1 = 21, pch2 = 23, ylim = c(-0.23, -0.13),
    xadj = c(Flavobacteriia = -0.17, Gammaproteobacteria = 0.2, Campylobacteria = 0.45, Deltaproteobacteria = 0.1),
    yadj = c(Flavobacteriia = 1.2, Gammaproteobacteria = 1.7, Campylobacteria = -3, Deltaproteobacteria = -0.8),
    mdat = mdat, RDP = RDP, map = map
  )
  p <- sum(gg$Pboth[!is.na(gg$X1) & !is.na(gg$X2)])
  title(paste0("Class (", round(p), "% of total)"))

  gg <- groupmet("MPB+17", "ZC", "genus", pch1 = 21, pch2 = 23, ylim = c(-0.23, -0.13), minpercent = 1,
    xadj = c(Alteromonas = 0.1, Sulfurimonas = 1.05, Alcanivorax = -1, Halomonas = -0.65, Thiogranum = 0.17, Sulfurovum = -0.1, Pseudomonas = 0.1, Pseudoalteromonas = -0.25, Acinetobacter = 0.72),
    yadj = c(Sulfurimonas = 1.5, Pseudomonas = 1.8, Sulfurovum = -0.5, Thiogranum = 5.5, Marinimicrobia_genera_incertae_sedis = -0.8, Alteromonas = 1.8, Acinetobacter = 1.4),
    mdat = mdat, RDP = RDP, map = map
  )
  lines(c(0, 0), c(-0.1535, -0.1425))
  p <- sum(gg$Pboth[!is.na(gg$X1) & !is.na(gg$X2)])
  title(paste0("Genus (", round(p), "% of total)"))

  # Make plots for Baltic Sea
  mdat <- getmdat("HLA+16")
  RDP <- getRDP("HLA+16", mdat = mdat)
  map <- getmap("HLA+16", RDP = RDP)
  gg <- groupmet("HLA+16", "nH2O", "domain", pch1 = 21, pch2 = 24, ylim = c(-0.78, -0.71), xlim = c(0, 100),
    xadj = c(Bacteria = 1), yadj = c(Bacteria = 1.5),
    mdat = mdat, RDP = RDP, map = map
  )
  title("Baltic Sea")
  text(40, -0.743, "Salinity < 6", font = 2)
  text(40, -0.7534, "Salinity > 20", font = 2)

  gg <- groupmet("HLA+16", "nH2O", "phylum", pch1 = 21, pch2 = 24, ylim = c(-0.78, -0.71),
    xadj = c(Proteobacteria = -0.2, Planctomycetes = 0.1, "Cyanobacteria/Chloroplast" = 0.35, Bacteroidetes = -0.1),
    yadj = c(Planctomycetes = 1.8, "Cyanobacteria/Chloroplast" = 2),
    mdat = mdat, RDP = RDP, map = map
  )
  p <- sum(gg$Pboth[!is.na(gg$X1) & !is.na(gg$X2)])
  title(paste0("Phylum (", round(p), "% of total)"))

  gg <- groupmet("HLA+16", "nH2O", "class", pch1 = 21, pch2 = 24, ylim = c(-0.78, -0.71),
    xadj = c(Acidimicrobiia = 0.5, Gammaproteobacteria = 0.4, Flavobacteriia = -0.05, Verrucomicrobiae = -0.3, Betaproteobacteria = 0.55, Cyanobacteria = -0.07),
    yadj = c(Acidimicrobiia = -1, Alphaproteobacteria = -0.6, Gammaproteobacteria = -1, Verrucomicrobiae = 1.5, Flavobacteriia = 1.2, Betaproteobacteria = -1, Planctomycetacia = 1),
    mdat = mdat, RDP = RDP, map = map
  )
  p <- sum(gg$Pboth[!is.na(gg$X1) & !is.na(gg$X2)])
  title(paste0("Class (", round(p), "% of total)"))

  gg <- groupmet("HLA+16", "nH2O", "genus", pch1 = 21, pch2 = 24, ylim = c(-0.78, -0.71), minpercent = 1,
    xadj = c(Spartobacteria_genera_incertae_sedis = 1.02, `Candidatus Pelagibacter` = 0, GpIIa = 1.1),
    yadj = c(Spartobacteria_genera_incertae_sedis = -0.7, GpI = 0.2),
    mdat = mdat, RDP = RDP, map = map
  )
  p <- sum(gg$Pboth[!is.na(gg$X1) & !is.na(gg$X2)])
  title(paste0("Genus (", round(p), "% of total)"))

#  # Between-study comparisons
#  groupmet("MPB+17", "ZC", ylim = c(-0.02, 0.02), study2 = "HLA+16")
#  mtext("Manus Basin - Baltic Sea", line = 1)
#
#  groupmet("BGPF13", "ZC", ylim = c(-0.12, 0.02), study2 = "XDZ+17")
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
  plotmet("UKD+18.water_2014", xlim = xlim, ylim = ylim, title = FALSE)
  legend("topleft", c("MSA-", "MSA+"), pch = c(1, 21), pt.bg = c(1, 2), bg = "white", title = "NW PA water")
  label.figure("A", cex = 1.5, xfrac = 0.03, font = 2)

  ## Plot B: Comparison of different studies on Pennsylvania Streams 20210327

  studies <- c("UKD+18.water_2014", "UKD+18.sediment_2014", "MMA+20_spring", "MMA+20_fall")
  # Start plot
  plot(c(-0.148, -0.139), c(-0.745, -0.735), type = "n", xlab = cplab$ZC, ylab = cplab$nH2O)
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
    # Save values for Source Data 20210901
    outB[[i]] <- pm
  }
  # Add labels
  text(-0.1423, -0.7388, "NW PA\nwater")
  text(-0.1473, -0.742, "NW PA\nsediment")
  text(-0.1427, -0.7432, "PASF water (spring)")
  text(-0.143, -0.7447, "PASF water (fall)")
  # Add legend
  legend("topleft", c("Lowest disturbance", "Highest disturbance"), pch = c(21, 21), pt.bg = c("#ffffffa0", "#df536ba0"), pt.cex = c(1.4, 1.7), lwd = 2, lty = NA)
  label.figure("B", cex = 1.5, xfrac = 0.03, font = 2)

  ## Plots C-D: Comparison of different studies on produced water 20210330

  # Panel C: Cluff et al., 2014
  plotmet("CHM+14", title = FALSE)
  legend("topright", c("Injected fluids (day 0)", "Produced water (day 49 and after)"),
         pch = c(21, 21), pt.bg = c("white", 2), bg = "white", title = "Marcellus Shale")
  box()
  label.figure("C", cex = 1.5, xfrac = 0.03, font = 2)

  # Panel D: Multiple studies
  studies <- c("CHM+14", "HRR+18", "ZLF+19")
  # Start plot
  plot(c(-0.22, -0.14), c(-0.75, -0.71), type = "n", xlab = cplab$ZC, ylab = cplab$nH2O)
  pch <- 21:25
  xadj <- c(0.3, -1, 0.5)
  yadj <- c(-1.1, 0, -1)
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
    text(means$ZC1, means$nH2O1, n1, adj = c(xadj[i], yadj[i]))
    text(means$ZC2, means$nH2O2, n2, adj = c(xadj[i], yadj[i]))
    outD[[i]] <- pm
  }
  # Add labels
  text(-0.164, -0.726, "Marcellus Shale")
  text(-0.204, -0.729, "Denver-Julesburg Basin")
  text(-0.172, -0.742, "Duvernay Formation")
  # Add legend
  legend("topright", c("Source water\nor injected fluids", "Produced water"), pch = c(21, 21), pt.bg = c("#ffffffa0", "#df536ba0"), pt.cex = c(1.4, 1.7), lwd = 2, lty = NA)
  label.figure("D", cex = 1.5, xfrac = 0.03, font = 2)

  if(pdf) dev.off()

  # Return Source Data 20210901
  outB <- do.call(rbind, outB)
  outD <- do.call(rbind, outD)
  out <- rbind(outB, outD)
  out <- out[, !colnames(out) %in% c("pch", "col")]
  out$nH2O <- round(out$nH2O, 6)
  out$ZC <- round(out$ZC, 6)
  invisible(out)

}

# Abundance and ZC of classes in oxidizing and reducing conditions 20210609
geo16S_S1 <- function(pdf = FALSE) {
  if(pdf) pdf("geo16S_S1.pdf", width = 9, height = 12)
  par(mfrow = c(4, 3))
  par(mar = c(4, 4, 2, 1))
  par(mgp = c(2, 1, 0))

  study <- c(
    "GBL+15", "JHM+16", "MPB+17", "BCA+21",
    "SVH+19", "MZG+20", "HXZ+20",
    "UKD+18.water", "MMA+20_spring",
    "CHM+14", "HRR+18", "ZLF+19"
  )
  description <- c(
    "ETNP water", "Lake Fryxell mat", "Manus Basin vents", "Ursu Lake water",
    "Black Sea water", "Swiss lakes", "Sansha Yongle Blue Hole",
    "NW Pennsylvania water", "PASF Streams (spring)",
    "Marcellus Shale", "Denver-Julesburg Basin", "Duvernay Formation"
  )
  pch1 <- c(24, 24, 21, 24, 24, 24, 24, 1, 1, 1, 1, 1, 1)
  pch2 <- c(25, 25, 23, 25, 25, 25, 25, 21, 21, 21, 21, 21, 21)
  xadj <- list(
    c(Gammaproteobacteria = 1.1, Flavobacteriia = -0.9, Planctomycetacia = 0.5),
    c(Alphaproteobacteria = -0.25, Saprospiria = 0.5),
    c(Campylobacteria = 0.5, Deltaproteobacteria = 0.1),
    c(Saprospiria = -0.6, Gammaproteobacteria = -0.5, Balneolia = -0.65, Deltaproteobacteria = 1.3, Cytophagia = -0.35, Clostridia = 1.5, Verrucomicrobiae = -0.1, Planctomycetacia = -0.22),
    c(Planctomycetacia = -0.2, Cyanobacteria = -0.9, Alphaproteobacteria = -0.75, Verrucomicrobiae = -0.3, Flavobacteriia = -0.45, Deltaproteobacteria = 0.9, Gammaproteobacteria = -0.35),
    c(Phycisphaerae = -0.3, Planctomycetacia = -0.27, Alphaproteobacteria = -0.15, Nitrospira = -0.15, Cyanobacteria = -0.65, Gammaproteobacteria = 1),
    c(Acidimicrobiia = -0.05, Alphaproteobacteria = -0.05, Cyanobacteria = -0.8, Gammaproteobacteria = 1.05),
    c(Alphaproteobacteria = 1, Bacilli = -0.2, Actinobacteria = -0.05),
    c(Acidobacteria_Gp6 = -0.1, Acidobacteria_Gp1 = -0.1, Spartobacteria = -0.1, Subdivision3 = -0.1, Deltaproteobacteria = 0.3, Alphaproteobacteria = 0.87, Betaproteobacteria = 0.2),
    c(Gammaproteobacteria = -0.1, Clostridia = 1.2, Campylobacteria = -0.55),
    c(Gammaproteobacteria = -0.2),
    c(Actinobacteria = -0.05, Betaproteobacteria = -0.3, Alphaproteobacteria = 0.1, Flavobacteriia = -0.5, Clostridia = 1.7)
  )
  yadj <- list(
    c(Cyanobacteria = 1.5, Flavobacteriia = 0, Planctomycetacia = -0.8),
    c(Saprospiria = -0.5, Planctomycetacia = 1.5, Verrucomicrobiae = -0.4, Caldilineae = 2, Deltaproteobacteria = -0.2, Cytophagia = 1.4),
    c(Gammaproteobacteria = 1.2, Deltaproteobacteria = -0.6, Flavobacteriia = -0.4, Campylobacteria = -3.5),
    c(Saprospiria = 0, Balneolia = 0, Deltaproteobacteria = -1.7, Cytophagia = -0.6, Flavobacteriia = 1.5, Clostridia = 0, Verrucomicrobiae = 2),
    c(Alphaproteobacteria = -2.6, Verrucomicrobiae = 12, Flavobacteriia = -2, Deltaproteobacteria = 1.8),
    c(Phycisphaerae = -1, Planctomycetacia = -12, Nitrospira = 1.9, Deltaproteobacteria = -0.5, Gammaproteobacteria = -1.2),
    c(Acidimicrobiia = 1.3, Cyanobacteria = 1.5, Nitrospinia = 1.7, Flavobacteriia = -0.5),
    c(Betaproteobacteria = 1.3, Clostridia = -0.2),
    c(Acidobacteria_Gp3 = -1.3, Subdivision3 = 1, Acidobacteria_Gp1 = -0.5, Spartobacteria = 0.3, Deltaproteobacteria = 1.3, Betaproteobacteria = -0.8, Alphaproteobacteria = 1.9),
    c(Clostridia = 1.3),
    c(Gammaproteobacteria = 1.3),
    c(Actinobacteria = 1.2, Betaproteobacteria = 1, Alphaproteobacteria = 1.6)
  )

  n1 <- n2 <- DZC <- numeric()
  for(i in 1:length(study)) {
    gg <- groupmet(study[i], param = "ZC", rank = "class", pch1 = pch1[i], pch2 = pch2[i], xadj = xadj[[i]], yadj = yadj[[i]], scale100 = TRUE, minpercent = 2)
    title(description[i])
    # Assemble percent contribution by each taxonomic group
    P <- round(gg$DX / diff(gg$Xall) * 100)
    P <- as.data.frame(t(P))
    # Include study name and merge with other studies
    P <- cbind(study = study[i], P)
    if(i==1) percent <- P else percent <- merge(percent, P, all = TRUE)
    # Keep track of the overall ZC change
    DZC <- c(DZC, diff(gg$Xall))
    # Keep track of number of samples 20210902
    n1 <- c(n1, gg$n1)
    n2 <- c(n2, gg$n2)
  }
  if(pdf) dev.off()
  # Put studies in correct order (they are alphabetized by merge()) 20210902
  percent <- percent[match(study, percent$study), ]

  # Save results for geo16S6() 20210610
  out <- cbind(n1 = n1, n2 = n2, DZC = round(DZC, 6), percent)
  #write.csv(out, "geo16S_S1.csv", row.names = FALSE, quote = FALSE)
}

# Contributions of classes to overall ZC difference between oxidizing and reducing conditions 20210610
geo16S6 <- function(pdf = FALSE) {
  # Setup plot
  if(pdf) pdf("geo16S6.pdf", width = 14, height = 8)
  layout(matrix(c(1, 2), nrow = 1), widths = c(1, 8))
  
  # Read file created by geo16S_S1()
  file <- system.file("extdata/geo16S/geo16S_S1.csv", package = "JMDplots")
  dat <- read.csv(file, as.is = TRUE)
  DZC <- dat$DZC
  n1 <- dat$n1
  n2 <- dat$n2
  # Get the columns with the study key and percent change for classes
  percent <- dat[, -(1:3)]

  study <- c(
    "GBL+15", "JHM+16", "MPB+17", "BCA+21",
    "SVH+19", "MZG+20", "HXZ+20",
    "UKD+18.water", "MMA+20_spring",
    "CHM+14", "HRR+18", "ZLF+19"
  )
  description <- c(
    "ETNP water", "Lake Fryxell mat", "Manus Basin vents", "Ursu Lake water",
    "Black Sea water", "Swiss lakes", "Sansha Yongle Blue Hole",
    "NW Pennsylvania water", "PASF Streams (spring)",
    "Marcellus Shale", "Denver-Julesburg Basin", "Duvernay Formation"
  )
  # Include numbers of samples in condition text 20210902
  cond2 <- c( "> 100 m", "anoxic", "> 50 \u00B0C", "> 4 m", "euxinic", "deepest", "anoxic", "MSA+", "highest", "PW day 49+", "PW day 130+", "FW day 18")
  cond1 <- c("< 100 m", "oxic", "< 50 \u00B0C", "< 3 m", "oxic", "shallowest", "oxic", "MSA-", "lowest", "IF day 0", "SW day 0", "SW day 0")
  condition <- paste0(cond2, " (", n2, ") - ", cond1, " (", n1, ")")

  # Plot ZC 
  par(mar = c(8, 1, 1, 0.2))
  par(las = 1)
  plot(range(DZC), c(12.5, 0.5), ylim = c(12.5, 0.5), yaxs = "i", yaxt = "n", ylab = "", xlab = cplab$DZC, type = "n")
  for(i in 1:12) lines(c(DZC[i], DZC[i]), c(i-0.5, i+0.5), lwd = 2)

  # Make heatmap 20210609
  # Move "study" column to rownames
  rownames(percent) <- percent$study
  percent <- percent[, -1]
  # Make sure rows are in right order
  stopifnot(all(rownames(percent) == study))
  # Reorder columns to group classes by phyla
  classes <- c(
  "Alphaproteobacteria", "Betaproteobacteria", "Gammaproteobacteria", "Deltaproteobacteria",  # Proteobacteria
  "Campylobacteria",  # Campylobacterota
  "Acidobacteria_Gp1", "Acidobacteria_Gp3", "Acidobacteria_Gp6",  # Acidobacteria
  "Actinobacteria", "Acidimicrobiia",  # Actinobacteria
  "Clostridia", "Bacilli",  # Firmicutes
  "Balneolia", "Chitinophagia", "Cytophagia", "Flavobacteriia", "Saprospiria",  # Bacteroidetes
  "Verrucomicrobiae", "Spartobacteria", "Subdivision3",  # Verrucomicrobia
  "Planctomycetacia", "Phycisphaerae", # Planctomycetes 
  "Anaerolineae", "Caldilineae",  # Chloroflexi
  "Cyanobacteria", # Cyanobateria
  "Chlamydiia",  # Chlamydiae
  "Aquificae",   # Aquificae
  "Chlorobia",   # Chlorobi
  "Nitrospira",  # Nitrospirae
  "Nitrospinia", # Nitrospinae
  "Synergistia", # Synergistetes
  "Thermotogae"  # Thermotogae
  )
  icol <- match(classes, colnames(percent))
  # Put any unmatched columns at the end
  allcol <- 1:ncol(percent)
  icol <- c(icol, allcol[!allcol %in% icol])
  percent <- percent[, icol]
  # Save the names for making rotated labels
  labels <- colnames(percent)
  colnames(percent) <- rep("", ncol(percent))
  # Move names of single classes in phyla upwards
  isingle <- labels %in% c("Campylobacteria", "Planctomycetacia", "Cyanobacteria", "Chlamydiia",
    "Aquificae", "Chlorobia", "Nitrospira", "Nitrospinia", "Synergistia", "Thermotogae",
    "Actinobacteria", "Acidimicrobiia", "Balneolia", "Chitinophagia", "Cytophagia", "Flavobacteriia", "Saprospiria",
    "Planctomycetacia", "Phycisphaerae")
  labels[!isingle] <- paste0(labels[!isingle], "   ")

  # Make it a matrix
  percent <- as.matrix(percent)
  # Remove rownames - we will label the axis later
  rownames(percent) <- rep("", nrow(percent))
  # Truncate values to [-100, 100]
  Pall <- percent
  percent[percent < -100] <- -100
  percent[percent > 100] <- 100

  # Setup margins, colors, breaks
  par(cex = 1)
  par(mar = c(8, 12, 1, 4))
  par(tcl = -0.3)
  col <- function(n) hcl.colors(n, "RdYlBu")
  breaks <- c(-100, -50, -20, -10, 0, 10, 20, 50, 100)

  # Plot heatmap
  requireNamespace("plot.matrix")
  plot(percent, col = col, breaks = breaks, main = "", xlab = "", ylab = "")
  # Add triangles to indicate values > 100
  iup <- Pall > 100
  iup[is.na(iup)] <- FALSE
  xyup <- which(iup, arr.ind = TRUE)
  x <- xyup[, 2]
  y <- 13 - xyup[, 1]
  points(x, y, pch = 2, cex = 1, lwd = 2, col = "#ffffffc0")
  # We would plot down-pointing triange for values < -100, but there aren't any
  stopifnot(all(na.omit(as.vector(Pall >= -100))))

  # Make rotated labels (modified from https://www.r-bloggers.com/rotated-axis-labels-in-r-plots/)
  text(x = seq_along(labels) + 0.3, y = par()$usr[3] - 2 * strheight("A"), labels = labels, srt = 40, adj = 1, xpd = TRUE)

  # Add phylum labels
  par(xpd = NA)
  text(2.5, 0.1, "Proteobacteria", font = 2, cex = 0.7)
  text(7, 0.05, "Acidobacteria", font = 2, cex = 0.7)
  text(9.5, 0.23, "Actinobacteria", font = 2, cex = 0.7)
  text(11.5, 0.05, "Firmicutes", font = 2, cex = 0.7)
  text(15, 0.2, "Bacteroidetes", font = 2, cex = 0.7)
  lines(c(12.8, 13.5), c(0.2, 0.2))
  lines(c(16.5, 17.2), c(0.2, 0.2))
  text(19, 0.05, "Verrucomicrobia", font = 2, cex = 0.7)
  text(21.5, 0.23, "Planctomycetes", font = 2, cex = 0.7)
  text(23.5, 0.05, "Chloroflexi", font = 2, cex = 0.7)

  # Add outer legend label
  text(33.8, 12.4, expression(Delta * italic(Z)[C] * "%"))
  par(xpd = FALSE)

  # Add dataset description and conditions
  axis(2, at = 12:1 + 0.2, labels = description, tick = FALSE)
  axis(2, at = 12:1 - 0.2, labels = condition, tick = FALSE, cex.axis = 0.9)

  if(pdf) dev.off()
}

# Comparison of protein ZC from metagenomic or metatranscriptomic data with estimates from 16S and reference sequences 20211017
geo16S7 <- function(pdf = FALSE) {

  if(pdf) pdf("geo16S7.pdf", width = 7, height = 5)
  par(mar = c(4, 4, 1, 1), mgp = c(2.8, 1, 0))

  # Start plot
  xylim <- c(-0.22, -0.13)
  xlab <- quote(italic(Z)[C]~"from metagenome or metatranscriptome")
  ylab <- quote(italic(Z)[C]~"estimated from 16S and reference sequences")
  plot(xylim, xylim, type = "n", xlab = xlab, ylab = ylab)
  lines(xylim, xylim, lty = 2, col = "gray40")

  # Get ZC calculated from Guerrero Negro metagenome (Kunin et al., 2008)
  dat_MG <- mplot("Guerrero_Negro", "IMG_MGP", plot.it = FALSE)
  # Reverse the order because upper mat layers are plotted on the right
  ZC_MG <- rev(dat_MG$AA)
  # Get ZC calculated from Guerrero Negro 16S data (Harris et al., 2013)
  dat_16S <- getmetrics("HCW+13")
  ZC_16S <- dat_16S$ZC
  # Check that the sample names are the same
  stopifnot(all.equal(rev(rownames(dat_MG$meancomp)), gsub(".*_", "", dat_16S$sample)))
  # Add lines to plot
  lines(ZC_MG, ZC_16S, type = "b")
  # Fill symbol for most oxidized (surface) sample
  points(ZC_MG[1], ZC_16S[1], pch = 19)
  # Label points for upper 3 layers
  text(c(-0.135, -0.1382, -0.1422), c(-0.1472, -0.1598, -0.1618), c("1", "2", "3"), cex = 0.8)
  text(-0.135, -0.17, "Guerrero\nNegro")

  # Get ZC calculated from ETNP OMZ metagenome (Glass et al., 2015)
  dat_MG <- mplot("ETNP_OMZ", "SRA_MGP", plot.it = FALSE)
  # Reverse the order because smaller water depths are plotted on the right
  ZC_MG <- rev(dat_MG$AA)
  # Get ZC calculated from ETNP OMZ 16S data (Ganesh et al., 2015)
  dat_16S <- getmetrics("GBL+15")
  # Use smallest size fraction
  dat_16S <- subset(dat_16S, grepl("1.6micron", dat_16S$sample))
  ZC_16S <- dat_16S$ZC
  # Check that the sample names are the same
  all.equal(rev(rownames(dat_MG$meancomp)), gsub("_.*", "", dat_16S$sample))
  # Add lines to plot
  lines(ZC_MG, ZC_16S, type = "b")
  # Fill symbol for most oxidized (surface) sample
  points(ZC_MG[1], ZC_16S[1], pch = 19)
  text(-0.153, -0.177, "ETNP")

  # Get ZC calculated from Bison Pool metagenome (Havig et al., 2011)
  dat_MG <- mplot("Bison_Pool", "IMG_MGP", plot.it = FALSE)
  ZC_MG <- dat_MG$AA
  # Get ZC calculated from Bison Pool 16S data (Swingley et al., 2012)
  dat_16S <- getmetrics("SMS+12")
  ZC_16S <- dat_16S$ZC
  # Check that the sample names are the same
  stopifnot(all.equal(rownames(dat_MG$meancomp), dat_16S$sample))
  # Add lines to plot
  lines(ZC_MG, ZC_16S, type = "b")
  # Fill symbol for most oxidized (low-T) sample
  points(ZC_MG[5], ZC_16S[5], pch = 19)
  text(-0.18, -0.197, "Bison Pool")

  # Get ZC calculated from Mono Lake metatranscriptome (Edwardson and Hollibaugh, 2017)
  dat_MT <- mplot("Mono_Lake", "SRA_MTP", plot.it = FALSE)
  # NOTE: reverse the order because smaller water depths are plotted on the right
  ZC_MT <- rev(dat_MT$AA)
  # Get ZC calculated from Mono Lake 16S data (Edwardson and Hollibaugh, 2018)
  dat_16S <- getmetrics("EH18")
  ZC_16S <- dat_16S$ZC
  # Check that the sample names are the same
  stopifnot(all.equal(rev(rownames(dat_MT$meancomp)), gsub(".*_", "", dat_16S$sample)))
  # Add lines to plot
  lines(ZC_MT, ZC_16S, type = "b", pch = 0)
  # Fill symbol for most oxidized (surface) sample
  points(ZC_MT[1], ZC_16S[1], pch = 15)
  text(-0.172, -0.145, "Mono Lake")

  # Add legend
  legend("topleft", character(3), pch = c(1, 0, 19), bty = "n")
  legend("topleft", c("Metagenome", "Metatranscriptome", "Most oxidized", "(surface or low-T sample)"), pch = c(NA, NA, 15, NA), bty = "n", inset = c(0.02, 0))

  if(pdf) dev.off()


}

