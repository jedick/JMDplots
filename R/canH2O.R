# JMDplots/canH2O.R
# make plots for the paper:
# Water as a reactant in the differential expression of proteins in cancer
# 20191126 jmd first version
# 20200203 started moving functions to JMDplots
# 20200404 include liver cancer

# study overview 20191203
canH2O1 <- function(pdf = FALSE) {
  if(pdf) pdf("canH2O1.pdf", width = 7, height = 7)
  opar <- par(mar = c(0, 0, 0, 0))
  mat <- matrix(c(1, 1, 1, 1, 2, 1), nrow = 2, byrow = TRUE)
  layout(mat, widths = c(0.2, 1, 1), heights = c(1.5, 1))
  par(cex = 1)
  openplotmat(frame.plot = TRUE)
  # three rows with eleven positions
  npos <- c(11, 11, 11)
  pos <- coordinates(npos)
  # the positions and x-sizes of the boxes
  p1 <- 2; r1 <- 0.08
  p2 <- 5; r2 <- 0.1
  p3 <- 9; r3 <- 0.15
  p4 <- 19; r4 <- 0.115
  p5 <- 16; r5 <- 0.16
  p6 <- 13; r6 <- 0.07
  p7 <- 28; r7 <- 0.09
  pO <- p1 + 2*npos[1] # the position for the overview plot

  straightarrow(from = pos[p1, ], to = pos[p2, ], arr.pos = 0.56, endhead = TRUE)
  straightarrow(from = pos[p2, ], to = pos[p3, ], arr.pos = 0.53, endhead = TRUE)
  curvedarrow(from = pos[p3, ], to = pos[p4, ] + c(0, 0.08), curve = 0.5, arr.pos = 0.67, endhead = TRUE)
  curvedarrow(from = pos[p3, ] + c(0, -0.05), to = pos[p5, ] + c(0.08, 0.08), curve = 0.3, arr.pos = 0.725, endhead = TRUE)
  curvedarrow(from = pos[p4, ] + c(-0.06, -0.01), to = pos[pO, ] + c(0.27, 0.33), curve = -0.3, arr.pos = 0.40, endhead = TRUE)
  curvedarrow(from = pos[p6, ] + c(-0.07, 0), to = pos[pO, ] + c(0.01, 0.27), curve = 0.2, arr.pos = 0.75, endhead = TRUE)

  # hard-coded colors (palette.colors is not available in R < 4.0.0) 20200415
  #cols <- palette.colors(8, "Set 2")
  cols <- c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3")
  cex <- 1.3
  ry <- 0.08
  textrect(pos[p1, ], r1, ry, lab = c("Papers", "SI Tables"), cex = cex, box.col = cols[1])
  textplain(pos[p1, ] + c(0, 0.12), lab = "Literature Search", font = 2)
  textrect(pos[p2, ], r2, ry, lab = c("UniProt IDs", "Up   Down", "...      ...   "), cex = cex, box.col = cols[2])
  textplain(pos[p2, ] + c(0, 0.12), lab = c("Differential", "Expression"), font = 2, height = 0.04)
  Zctext <- quote(italic(Z)[C]~"oxidation state")
  nH2Otext <- quote(italic(n)[H[2]*O]~"hydration state")
  comptext <- as.expression(c(Zctext, nH2Otext))
  textrect(pos[p3, ], r3, ry, lab = comptext, cex = cex, box.col = cols[2])
  textplain(pos[p3, ] + c(0, 0.12), lab = c("Compositional", "Analysis"), font = 2, height = 0.04)

  # show numbers of datasets 
  cond1 <- c("hypoxia", "secreted", "osmotic_euk", "glucose", "3D")
  cond2 <- c("breast", "colorectal", "liver", "lung", "pancreatic", "prostate")
  vigout <- system.file("vignettes", package = "canprot")
  conddat <- function(cond) read.csv(paste0(vigout, "/", cond, "_", getOption("basis"), ".csv"), as.is = TRUE)
  culture <- lapply(cond1, conddat); names(culture) <- cond1
  cancer <- lapply(cond2, conddat); names(cancer) <- cond2
  total <- sum(sapply(culture, nrow), sapply(cancer, nrow))
  # calculate number of studies
  cultsets <- unlist(sapply(culture, "[[", "dataset"))
  cansets <- unlist(sapply(cancer, "[[", "dataset"))
  nstudies <- length(unique(sapply(strsplit(c(cultsets, cansets), "_"), "[", 1)))

  # add "Cancer vs normal" box
  textrect(pos[p4, ] + c(0.02, 0.065), r4, ry + 0.015, lab = "", cex = cex, box.col = cols[4])
  calab <- paste(sapply(cancer, nrow), cond2)
  textplain(pos[p4, ] + c(-0.087, 0.065), lab = calab, height = 0.105, cex = cex, adj = c(0, 0.5))
  textplain(pos[p4, ] + c(0.02, 0.2), lab = c("Cancer", "vs normal"), font = 2, height = 0.04)
  # add color legend 20200506
  #col2 <- palette.colors(8, "Classic Tableau")[c(7, 6, 2, 8, 5, 4)]
  col2 <- c("#E377C2", "#8C564B", "#FF7F0E", "#7F7F7F", "#9467BD", "#D62728")
  # function to get y-coordinates of text lines, extracted from diagram::textplain()
  ycoords <- function (mid, height = 0.1, lab = "", adj = c(0.5, 0.5), ...) {
    y <- numeric()
    if (length(lab) == 1) y <- mid[2]
    else {
      y1 <- mid[2] + height
      ddy <- 2 * height/(length(lab) + 1)
      for (i in 1:length(lab)) y <- c(y, y1 - ddy * i)
    }
    y
  }
  ycancer <- ycoords(pos[p4, ] + c(-0.08, 0.065), lab = calab, height = 0.105)
  xcancer <- rep(pos[p4, 1] + 0.117, length(ycancer))
  points(xcancer, ycancer, bg = col2, pch = 21, cex = 1.7)

  # add "Cell culture" box
  textrect(pos[p5, ] + c(-0.011, 0.08), r5, ry, lab = "", cex = cex, box.col = cols[3])
  # sum osmotic (salt) and glucose datasets for hyperosmotic stress 20200411
  culab <- sapply(culture, nrow)
  culab["osmotic_euk"] <- culab["osmotic_euk"] + culab["glucose"]
  culab <- culab[-4]
  # change "osmotic_euk" to "hyperosmotic" 20200414
  names(culab)[3] <- "hyperosmotic"
  # write "secreted in hypoxia" and "3D vs 2D culture" 20200118
  culab <- paste(culab, names(culab))
  shlab <- c(paste(culab[2], "in"), "     hypoxia")
  textplain(pos[p5, ] + c(-0.163, 0.095), lab = shlab, height = 0.038, cex = cex, adj = c(0, 0.5))
  culab <- c(culab[1], "", "", culab[3:4])
  culab[5] <- paste(culab[5], "vs 2D culture")
  textplain(pos[p5, ] + c(-0.163, 0.08), lab = culab, height = 0.09, cex = cex, adj = c(0, 0.5))
  textplain(pos[p5, ] + c(-0.02, 0.2), lab = c("Cell", "culture"), font = 2, height = 0.04)
  # add color legend 20200506
  #col1 <- palette.colors(9, "Okabe-Ito")[c(6, 4, 9, 3, 2)]
  col1 <- c(blue = "#0072B2", bluishgreen = "#009E73", gray = "#999999", skyblue = "#56B4E9", orange = "#E69F00")
  yculture <- ycoords(pos[p5, ] + c(-0.16, 0.08), lab = culab, height = 0.09)
  xculture <- rep(pos[p5, 1] + 0.13, length(yculture))
  # point salt-inducded hyperosmotic stress on same line as glucose
  points(xculture[-3], yculture[-3], bg = col1[-3], pch = 21, cex = 1.7)
  points(xculture[3] - 0.025, yculture[4], bg = col1[3], pch = 21, cex = 1.7)

  # add dataset summary, pan-cancer, and findings text
  textplain(pos[21, ] + c(-0.02, 0.1), lab = c(paste("Total:", total), "datasets", paste("from", nstudies), "studies"),
            adj = c(0, 0.5), font = 3, height = 0.07, cex = 1.2)

  textrect(pos[p6, ] + c(0, 0.08), r6, ry, lab = c("GEPIA2", "(TCGA /", "GTEx)", ""), cex = cex, box.col = cols[7])
  textplain(pos[p6, ] + c(-0.025, 0.025), r6, ry, lab = "HPA", cex = cex)
  textplain(pos[p6, ] + c(0, 0.2), lab = c("Pan-cancer", "datasets"), font = 2, height = 0.04)

  textplain(pos[25, ] + c(0.105, 0.295), lab = "Main finding", font = 2)
  cantext1 <- "Most cancer types have higher"
  cantext2 <- "hydration state of proteins"
  cantext3 <- "compared to normal tissue."
  cantext <- c(cantext1, cantext2, cantext3)
  textplain(pos[25, ] + c(0.105, 0.24), lab = cantext, height = 0.05, font = 4)

  textplain(pos[27, ] + c(0.16, 0.07), lab = "Other findings", font = 2)
  hyptext1 <- quote(italic("Hypoxia experiments show no consistent"))
  hyptext2 <- quote(italic("difference in oxidation state of proteins."))
  hyptext <- as.expression(c(hyptext1, hyptext2))
  textplain(pos[26, ] + c(0.265, 0.005), lab = hyptext, height = 0.04)

  hydtext1 <- quote(italic("Most hyperosmotic stress and 3D vs 2D culture"))
  hydtext2 <- quote(italic("experiments yield lower hydration state of proteins."))
  hydtext <- as.expression(c(hydtext1, hydtext2))
  textplain(pos[26, ] + c(0.21, -0.08), lab = hydtext, height = 0.04)

  # make overview plot
  plot.new()
  par(mar = c(2, 2, 0, 0), xaxs = "i", yaxs = "i")
  plot.window(c(-1, 1), c(-1, 1))
  axis(1, tck = 0, labels = FALSE)
  axis(2, tck = 0, labels = FALSE)
  mtext("oxidation state", 1, 0.5)
  mtext("hydration state", 2, 0.5)
  # draw arrows to mean values, muliplied by a constant to scale to (-1, 1) axis range
  for(i in 5:1) arrows(0, 0, 40*mean(culture[[i]]$Zc.diff), 40*mean(culture[[i]]$nH2O.diff), col = col1[i], lwd = 3, length = 0.15)
  for(i in 1:6) arrows(0, 0, 40*mean(cancer[[i]]$Zc.diff), 40*mean(cancer[[i]]$nH2O.diff), col = col2[i], lwd = 3, length = 0.15)

  # restore these defaults to be able to re-run this script with expected results
  par(xaxs = "r", yaxs = "r")
  if(pdf) {
    dev.off()
    addexif("canH2O1", "Study overview", "Dick (2021)")
  }
}

# median differences of nH2O and Zc for cell culture and cancer tissue 20191126
canH2O2 <- function(pdf = FALSE) {
  if(pdf) pdf("canH2O2.pdf", width = 8, height = 6)
  # layout with spaces between groups of plots
  # spaces
  p00 <- rep(0, 2); p0 <- rep(0, 15)
  # (A) cell culture
  p1 <- rep(1, 15); p2 <- rep(2, 15); p3 <- rep(3, 15); p4 <- rep(4, 15); p5 <- rep(5, 15)
  # (B) cancer
  p6 <- rep(6, 15); p7 <- rep(7, 15); p8 <- rep(8, 15)
  p9 <- rep(9, 15); p10 <- rep(10, 15); p11 <- rep(11, 15)
  # (C) credible regions
  p12 <- rep(12, 15); p13 <- rep(13, 15)
  # assemble columns (each one is 48 units high)
  c1 <- c(0, p1, p00, p6, p9)
  c2 <- c(0, p2, p00, p7, p10)
  c3 <- c(0, p3, p00, p8, p11)
  c4 <- c(0, p4, p00, p0, p0)
  c5 <- c(0, p4, p00, p12, p13)
  c6 <- c(0, p5, p00, p12, p13)
  mat <- matrix(c(c1, c2, c3, c4, c5, c6), nrow = 48)
  layout(mat, widths = c(1, 1, 1, 0.5, 0.5, 1))
  par(mar = c(3.2, 3.3, 1.5, 1), mgp = c(2.2, 1, 0))

  # read data for all conditions
  cond1 <- c("hypoxia", "secreted", "osmotic_euk", "glucose", "3D")
  cond2 <- c("breast", "colorectal", "liver", "lung", "pancreatic", "prostate")
  conds <- c(cond1, cond2)
  vigout <- system.file("vignettes", package = "canprot")
  alldat <- lapply(conds, function(cond) {
    file <- paste0(vigout, "/", cond, "_", getOption("basis"), ".csv")
    read.csv(file)
  })
  names(alldat) <- conds

  #col1 <- palette.colors(9, "Okabe-Ito")[c(6, 4, 9, 3, 2)]
  col1 <- c(blue = "#0072B2", bluishgreen = "#009E73", gray = "#999999", skyblue = "#56B4E9", orange = "#E69F00")
  #col2 <- palette.colors(8, "Classic Tableau")[c(7, 6, 2, 8, 5, 4)]
  col2 <- c("#E377C2", "#8C564B", "#FF7F0E", "#7F7F7F", "#9467BD", "#D62728")
  col <- c(col1, col2)

  # make scatter plots for nH2O and Zc in cell culture and cancer types
  lapply(1:length(alldat), function(i) {
    thisone <- names(alldat)[i]
    if(thisone %in% cond1) {
      # use open/filled symbols for cancer/non-cancer cells 20200103
      pch <- rep(21, nrow(alldat[[i]]))
      pch[grepl("cancer", alldat[[i]]$tags)] <- 19
      # use squares for yeast datasets 20200407
      pch[grepl("yeast", alldat[[i]]$tags)] <- 0
    } else {
      # use filled/open symbols for human/[mouse or rat] cancer 20200415
      pch <- rep(19, nrow(alldat[[i]]))
      pch[grepl("mouse", alldat[[i]]$tags)] <- 21
      pch[grepl("rat", alldat[[i]]$tags)] <- 21
    }
    dpargs <- list(comptab = alldat[i], pch = pch, pt.text = NA, col = col[i], cex = 1, labtext = NA,
                   xlim = c(-0.06, 0.06), ylim = c(-0.07, 0.07), axes = FALSE, frame.plot = TRUE)
    do.call(diffplot, dpargs)
    labels <- at <- seq(-0.06, 0.06, 0.02)
    labels[c(2, 3, 5, 6)] <- NA
    axis(1, at, labels)
    axis(2, at, labels)
    main <- names(alldat)[i]
    substr(main, 1, 1) <- toupper(substr(main, 1, 1))
    if(thisone == "secreted") main <- "hypoxia"
    if(thisone == "osmotic_euk") main <- "(salt)"
    if(thisone == "glucose") main <- "(glucose)"
    if(thisone == "3D") main <- "3D / 2D"
    title(main, font.main = 1)
    if(thisone == "secreted") title("Secreted in", font.main = 1, line = 1.4, xpd = NA)
    if(thisone %in% c("osmotic_euk", "glucose")) title("Hyperosmotic", font.main = 1, line = 1.4, xpd = NA)
    if(thisone == "hypoxia") label.figure("A", cex = 2, font = 2, yfrac = 1)
    if(thisone == "breast") label.figure("B", cex = 2, font = 2, yfrac = 1)
  })
  
  # compare density contours for cell culture and cancer types 20191126
  # CSV files generated by running the cancer and cell culture vignettes
  vigout <- system.file("vignettes", package = "canprot")
  conddat <- function(cond) read.csv(paste0(vigout, "/", cond, "_", getOption("basis"), ".csv"), as.is = TRUE)
  culture <- lapply(cond1, conddat); names(culture) <- cond1
  cancer <- lapply(cond2, conddat); names(cancer) <- cond2
  names <- names(culture)
  names[3] <- "salt"
  contplot(culture, "Cell culture", col1,
           dx = c(0.038, 0.030, -0.03, 0.0042, -0.026), dy = c(0.025, -0.054, 0.025, -0.036, -0.02), names = names)
  lines(c(0.0147, 0.024), c(0.0035, 0.016))
  lines(c(-0.0255, -0.0158), c(0.0143, 0.013))
  label.figure("C", cex = 2, font = 2, xfrac = 0.02, yfrac = 0.9)
  contplot(cancer, "Cancer tissue", col2,
           dx = c(-0.019, NA, 0.022, -0.022, 0.030, -0.025), dy = c(-0.028, NA, 0.027, 0.038, -0.035, 0.023))
  text(0.009, -0.021, "CRC", col = col2[2])

  if(pdf) {
    dev.off()
    addexif("canH2O2", "Median differences of protein length, nH2O and Zc for cell culture and cancer tissue", "Dick (2021)")
  }
}

# nH2O-Zc plots for TCGA and HPA datasets 20191126
canH2O3 <- function(pdf = FALSE) {
  vigout2 <- system.file("vignettes", package = "JMDplots")
  HPA <- read.csv(file.path(vigout2, paste0("HPA_", getOption("basis"), ".csv")), as.is = TRUE)
  TCGA <- read.csv(file.path(vigout2, paste0("TCGA_", getOption("basis"), ".csv")), as.is = TRUE)
  TCGA_labels <- TCGA$description
  HPA_labels <- HPA$description
  HPA_labels <- sapply(strsplit(HPA_labels, " "), "[", 1)
  HPA_labels[grepl("head", HPA_labels)] <- "head and neck"

  # get colors for six cancers in paper 20191208
  cond2 <- c("breast", "colorectal", "liver", "lung", "pancreatic", "prostate")
  #col2 <- palette.colors(8, "Classic Tableau")[c(7, 6, 2, 8, 5, 4)]
  col2 <- c("#E377C2", "#8C564B", "#FF7F0E", "#7F7F7F", "#9467BD", "#D62728")
  jHPA <- match(cond2, sapply(strsplit(HPA$description, " "), "[", 1))
  colHPA <- rep("darkslategray", nrow(HPA))
  colHPA[jHPA] <- col2
  shapeHPA <- rep(15, nrow(HPA))
  shapeHPA[jHPA] <- 19
  # now do TCGA
  TCGAnames <- names(HTmap)[match(cond2, sapply(strsplit(HTmap, " "), "[", 1))]
  jTCGA <- match(TCGAnames, TCGA$description)
  colTCGA <- rep("darkslategray", nrow(TCGA))
  colTCGA[jTCGA] <- col2
  shapeTCGA <- rep(15, nrow(TCGA))
  shapeTCGA[jTCGA] <- 19

  # calculate 2D density and 50% probability level 20200318
  densfun <- function(x, y) {
    n <- 200
    dens <- kde2d(x, y, n = n)
    # find 50% probability level
    # https://stackoverflow.com/questions/16225530/contours-of-percentiles-on-level-plot
    # (snippet from emdbook::HPDregionplot from @benbolker)
    dx <- diff(dens$x[1:2])
    dy <- diff(dens$y[1:2])
    sz <- sort(dens$z)
    c1 <- cumsum(sz) * dx * dy
    probs <- 0.5
    levels <- sapply(probs, function(x) {
      approx(c1, sz, xout = 1 - x)$y
    })
    # turn density into data frame for ggplot
    # get x- and y- values
    xr <- range(x)
    xs <- seq(xr[1], xr[2], length.out = n)
    yr <- range(y)
    ys <- seq(yr[1], yr[2], length.out = n)
    dat <- expand.grid(xs, ys)
    colnames(dat) <- c("x", "y")
    dat <- cbind(dat, z = as.vector(dens$z))
    list(levels = levels, dat = dat)
  }

  TCGAstuff <- densfun(TCGA$Zc.diff, TCGA$nH2O.diff)
  TCGAdens <- TCGAstuff$dat
  TCGAlevels <- TCGAstuff$levels
  HPAstuff <- densfun(HPA$Zc.diff, HPA$nH2O.diff)
  HPAdens <- HPAstuff$dat
  HPAlevels <- HPAstuff$levels

  # median differences of nH2O-Zc for HPA and TCGA datasets
  # common elements for both plots
  pl1.common <- list(
    theme_bw(),
    xlab(canprot::cplab$DZc),
    ylab(canprot::cplab$DnH2O),
    geom_hline(yintercept = 0, linetype = 3, colour = "gray30"),
    geom_vline(xintercept = 0, linetype = 3, colour = "gray30"),
    theme(plot.tag = element_text(size = 20), plot.title = element_text(hjust = 0.5))
  )
  # create plots
  nudge_x <- ifelse(TCGA_labels %in% c("SKCM"), 0.001, 0)
  # workaround for "no visible binding for global variable ‘Zc.diff’" etc. in R CMD check 20200317
  Zc.diff <- nH2O.diff <- x <- y <- z <- NULL
  pl1 <- list(
    ggplot(TCGA, aes(Zc.diff, nH2O.diff, label = TCGA_labels)) + 
      pl1.common + geom_point(color = colTCGA, size = 1.5, shape = shapeTCGA, stroke = 1.5) + 
      geom_contour(data = TCGAdens, aes(x, y, z = z), breaks = TCGAlevels, inherit.aes = FALSE, color = "slateblue4", lty = 2) +
      ggrepel::geom_text_repel(size = 2.5, nudge_x = nudge_x, seed = 42, box.padding = 0.12, point.padding = 0.1) +
      ggtitle("TCGA / GTEx") + labs(tag = expression(bold(A))),
    ggplot(HPA, aes(Zc.diff, nH2O.diff, label = HPA_labels)) + 
      pl1.common + geom_point(color = colHPA, size = 1.5, shape = shapeHPA, stroke = 1.5) + 
      geom_contour(data = HPAdens, aes(x, y, z = z), breaks = HPAlevels, inherit.aes = FALSE, color = "darkslategray4", lty = 2) +
      ggrepel::geom_text_repel(size = 3, seed = 42) +
      ggtitle("HPA") + labs(tag = expression(bold(B)))
  )

  # put together the figure
  mat <- matrix(c(1, 2), byrow = TRUE, nrow = 1)
  ml <- gridExtra::marrangeGrob(pl1, layout_matrix = mat, top = NULL)
  if(pdf) {
    ggsave("canH2O3.pdf", ml, width = 10, height = 4)
    addexif("canH2O3", "nH2O-Zc and plots for TCGA and HPA datasets", "Dick (2021)")
  } else ml
}

# differentially expressed genes in aneuploid and osmotically shocked yeast cells 20200505
canH2O4 <- function(pdf = FALSE) {
  if(pdf) pdf("canH2O4.pdf", width = 5, height = 4)
  # aneuploid yeast cells (Tsai et al., 2019)
  pdat <- pdat_aneuploidy("TNC+19")
  layout(matrix(c(1, 2, 4, 1, 3, 5), nrow = 3), heights = c(0.2, 1, 1))
  par(mar = c(0, 0, 0, 0), mgp = c(2.2, 0.8, 0))
  plot.new()
  uptxt <- paste("Proteins coded by", sum(pdat$up2), "up-regulated genes")
  dntxt <- paste("Proteins coded by", sum(!pdat$up2), "down-regulated genes")
  legend("center", c(dntxt, uptxt, ""), lty = c(1, 2, NA), col = c(1, 2, NA), bty = "n")
  par(mar = c(3.5, 3.5, 0.5, 0.5))
  qdist(pdat, "Zc")
  label.figure("A", xfrac = 0.04, yfrac = 1, font = 2, cex = 1.5)
  qdist(pdat, "nH2O")
  label.figure("B", xfrac = 0.04, yfrac = 1, font = 2, cex = 1.5)
  # hyper- and hypo-osmotic experiments (Gasch et al., 2000)
  vigout2 <- system.file("vignettes", package = "JMDplots")
  yeast_stress <- read.csv(file.path(vigout2, paste0("yeast_stress_", getOption("basis"), ".csv")), as.is = TRUE)
  ihyper <- grepl("sorbitol", yeast_stress$dataset)
  ihypo <- grepl("Hypo", yeast_stress$dataset)
  # plot Zc
  plot(c(1, 7), extendrange(yeast_stress$Zc.diff), xlab = "Time (minutes)", ylab = cplab$DZc, xaxt = "n", type = "n")
  abline(h = 0, lty = 3, col = "gray30")
  axis(1, at = 1:7, labels = c(5, 15, 30, 45, 60, 90, 120))
  lines(1:7, yeast_stress$Zc.diff[ihyper], pch = 1, type = "b")
  lines(1:5, yeast_stress$Zc.diff[ihypo], pch = 0, type = "b")
  text(c(2.5, 2.7), c(-0.006, 0.0135), c("hyperosmotic", "hypoosmotic"))
  label.figure("C", xfrac = 0.04, yfrac = 1, font = 2, cex = 1.5)
  # plot nH2O
  plot(c(1, 7), extendrange(yeast_stress$nH2O.diff), xlab = "Time (minutes)", ylab = cplab$DnH2O, xaxt = "n", type = "n")
  abline(h = 0, lty = 3, col = "gray30")
  axis(1, at = 1:7, labels = c(5, 15, 30, 45, 60, 90, 120))
  lines(1:7, yeast_stress$nH2O.diff[ihyper], pch = 1, type = "b")
  lines(1:5, yeast_stress$nH2O.diff[ihypo], pch = 0, type = "b")
  text(c(2.2, 2.8), c(-0.036, 0.024), c("hyperosmotic", "hypoosmotic"))
  label.figure("D", xfrac = 0.04, yfrac = 1, font = 2, cex = 1.5)
  if(pdf) {
    dev.off()
    addexif("canH2O4", "Differentially expressed genes in aneuploid and osmotically shocked yeast cells", "Dick (2021)")
  }
}

###############
### TABLE 2 ###
###############

# mean differences and p-values for all datasets in each condition 20200125
canH2OT2 <- function() {
  cond1 <- c("hypoxia", "secreted", "osmotic_euk", "glucose", "3D")
  cond2 <- c("breast", "colorectal", "liver", "lung", "pancreatic", "prostate")
  cond3 <- c("TCGA", "HPA")
  vigout <- system.file("vignettes", package = "canprot")
  conddat <- function(cond) read.csv(paste0(vigout, "/", cond, "_", getOption("basis"), ".csv"), as.is = TRUE)
  culture <- lapply(cond1, conddat); names(culture) <- cond1
  cancer <- lapply(cond2, conddat); names(cancer) <- cond2
  vigout2 <- system.file("vignettes", package = "JMDplots")
  conddat <- function(cond) read.csv(paste0(vigout2, "/", cond, "_", getOption("basis"), ".csv"))
  pancan <- lapply(cond3, conddat); names(pancan) <- cond3
  # make data frames for mean differences and p-values
  # include comparison of up- and down-regulated proteins for "Secreted" in hypoxia compared to "Hypoxia" (whole-cell)
  pdat <- mdat <- data.frame(condition = c(names(culture), names(cancer), names(pancan), "up", "down"))
  for(property in c("Zc", "nH2O")) {
    pvals <- mvals <- numeric()
    idn <- grep(paste0(property, ".down"), colnames(culture[[1]]))
    iup <- grep(paste0(property, ".up"), colnames(culture[[1]]))
    for(type in c("culture", "cancer", "pancan", "up", "down")) {
      if(type %in% c("up", "down")) {
        if(type=="up") idiff <- iup
        if(type=="down") idiff <- idn
        hypoxia <- na.omit(culture$hypoxia[, idiff])
        secreted <- na.omit(culture$secreted[, idiff])
        mvals <- c(mvals, mean(secreted) - mean(hypoxia))
        pvals <- c(pvals, t.test(hypoxia, secreted)$p.value)
      } else {
        thisdat <- get(type)
        for(i in 1:length(thisdat)) {
          dn <- thisdat[[i]][, idn]
          up <- thisdat[[i]][, iup]
          isNA <- is.na(dn) | is.na(up)
          dn <- dn[!isNA]
          up <- up[!isNA]
          mvals <- c(mvals, mean(up) - mean(dn))
          pvals <- c(pvals, t.test(dn, up, paired = TRUE)$p.value)
        }
      }
    }
    mcol <- data.frame(X = mvals)
    pcol <- data.frame(X = round(log10(pvals), 1))
    colnames(mcol)[1] <- property
    colnames(pcol)[1] <- property
    mdat <- cbind(mdat, mcol)
    pdat <- cbind(pdat, pcol)
  }
  # format values and put p-values in parentheses (bold for values below 0.05 significance threshold)
  fmts <- c(NA, "%.3f", "%.3f")
  for(icol in 2:3) {
    for(irow in 1:nrow(mdat)) {
      pround <- round(pdat[irow, icol], 1)
      ptxt <- sprintf("%.1f", pdat[irow, icol])
      if(pround < -1.3) ptxt <- paste0("(**", ptxt, "**)")
      else ptxt <- paste0("(", ptxt, ")")
      mdat[irow, icol] <- paste(sprintf(fmts[icol], as.numeric(mdat[irow, icol])), ptxt)
    }
  }
  # adjustments for pretty kable output
  rownames(mdat) <- mdat[, 1]
  rownames(mdat)[rownames(mdat)=="osmotic_euk"] <- "salt"
  mdat <- mdat[, -1]
  colnames(mdat) <- c("&Delta;*Z*~C~", "&Delta;*n*~H2O~")
  mdat
}

#############################
### SI TABLES AND FIGURES ###
#############################

# stoichiometric matrix for amino acids with QEC basis species 20200104
canH2OT1 <- function() {
  basis("QEC")
  species(aminoacids(""))
  out <- species()[, c(9, 1:5)]
  # adjustments for pretty kable output
  rownames(out) <- out[, 1]
  out <- out[, -1]
  colnames(out) <- gsub("([[:digit:]])", "~\\1~", colnames(out))
  out
}

# nO2-Zc and nH2O-Zc correlations using QEC basis species 20201016
canH2OS1 <- function(pdf = FALSE) {
  # set up figure
  if(pdf) pdf("canH2OS1.pdf", width = 7, height = 3.5)
  par(mfrow = c(1, 2))
  par(mar = c(3.5, 3.5, 1, 1))
  par(mgp = c(2.5, 0.7, 0))
  par(las = 1)
  par(cex.lab = 1.2)

  # define axis labels
  nH2Olab <- expression(italic(n)[H[2] * O])
  nO2lab <- expression(italic(n)[O[2]])
  Zclab <- expression(italic(Z)[C])

  # function to plot values for amino acids
  aaplot <- function(x, y, xlab, ylab, legend.x, lmlim = c(-1, 1)) {
    plot(x, y, type = "p", pch = aminoacids(1), xlab = xlab, ylab = NA)
    mtext(ylab, side = 2, line = 2.4, las = 0, cex = 1.2)
    lmfun(x, y, legend.x, lmlim)
  }

  # set up amino acid compositions (one amino acid per composition)
  AAcomp <- as.data.frame(diag(20))
  names(AAcomp) <- aminoacids(3)

  # plot 1: nO2-Zc of amino acids
  aaplot(Zc(AAcomp), nO2(AAcomp, "QEC"), Zclab, nO2lab, "bottomright")
  label.figure("A", cex = 1.7, font = 2)

  # plot 2: nH2O-Zc of amino acids
  aaplot(Zc(AAcomp), nH2O(AAcomp, "QEC"), Zclab, nH2Olab, "bottomright")
  label.figure("B", cex = 1.7, font = 2)

  if(pdf) {
    dev.off()
    addexif("canH2OS1", "nO2-Zc and nH2O-Zc correlations using QEC basis species", "Dick (2021)")
  }
}

# HPA-TCGA scatterplots for Zc and nH2O 20200127
canH2OS2 <- function(pdf = FALSE) {
  vigout2 <- system.file("vignettes", package = "JMDplots")
  HPA <- read.csv(file.path(vigout2, paste0("HPA_", getOption("basis"), ".csv")), as.is = TRUE)
  TCGA <- read.csv(file.path(vigout2, paste0("TCGA_", getOption("basis"), ".csv")), as.is = TRUE)
  TCGA_labels <- TCGA$description
  HPA_labels <- HPA$description
  HPA_labels <- sapply(strsplit(HPA_labels, " "), "[", 1)
  HPA_labels[grepl("head", HPA_labels)] <- "head and neck"

  # get colors for 6 cancers in paper 20191208
  cond2 <- c("breast", "colorectal", "liver", "lung", "pancreatic", "prostate")
  #col2 <- palette.colors(8, "Classic Tableau")[c(7, 6, 2, 8, 5, 4)]
  col2 <- c("#E377C2", "#8C564B", "#FF7F0E", "#7F7F7F", "#9467BD", "#D62728")
  jHPA <- match(cond2, sapply(strsplit(HPA$description, " "), "[", 1))
  colHPA <- rep("slateblue4", nrow(HPA))
  colHPA[jHPA] <- col2
  sizeHPA <- rep(1.5, nrow(HPA))
  sizeHPA[jHPA] <- 2
  shapeHPA <- rep(15, nrow(HPA))
  shapeHPA[jHPA] <- 19

  # HPA-TCGA mappings
  iHPA <- match(HTmap, HPA$description)
  iTCGA <- match(names(HTmap), TCGA$description)

  # HPA-TCGA scatterplots for Zc and nH2O
  Zc <- data.frame(TCGA = TCGA$Zc.diff[iTCGA], HPA = HPA$Zc.diff[iHPA])
  nH2O <- data.frame(TCGA = TCGA$nH2O.diff[iTCGA], HPA = HPA$nH2O.diff[iHPA])

  labels <- TCGA_labels[iTCGA]
  col <- "slateblue4"
  size <- 1.5
  shape <- 15

  r.squared.Zc <- format(summary(lm(HPA ~ TCGA, Zc))$r.squared, digits = 2)
  Zc.title <- paste0("italic(R)^2 == '", r.squared.Zc, "'")
  pl1 <- ggplot(Zc, aes(x = TCGA, y = HPA, label = labels)) +
    theme_classic() + geom_smooth(method = "lm") +
    annotate("text", -Inf, Inf, label = Zc.title, parse = TRUE, hjust = -0.2, vjust = 1.5) +
    xlab(quote(Delta*italic(Z)[C]*" (TCGA / GTEx)")) +
    ylab(quote(Delta*italic(Z)[C]*" (HPA)")) +
    geom_hline(yintercept = 0, linetype = 3, colour = "gray30") +
    geom_vline(xintercept = 0, linetype = 3, colour = "gray30") +
    geom_point(shape = shapeHPA, size = sizeHPA, col = colHPA, stroke = 1.5) +
    ggrepel::geom_text_repel(size = 3, seed = 42) +
    labs(tag = expression(bold(A))) +
    theme(plot.tag = element_text(size = 20), plot.title = element_text(hjust = 0.5),
          panel.border = element_rect(colour = "black", fill=NA))
  pl1 <- list(pl1)

  r.squared.nH2O <- format(summary(lm(HPA ~ TCGA, nH2O))$r.squared, digits = 2)
  nH2O.title <- paste0("italic(R)^2 == '", r.squared.nH2O, "'")
  pl2 <- ggplot(nH2O, aes(x = TCGA, y = HPA, label = labels)) +
    theme_classic() + geom_smooth(method = "lm") +
    annotate("text", Inf, -Inf, label = nH2O.title, parse = TRUE, hjust = 1.1, vjust = -0.6) +
    xlab(quote(Delta*italic(n)[H[2]*O]*" (TCGA / GTEx)")) +
    ylab(quote(Delta*italic(n)[H[2]*O]*" (HPA)")) +
    geom_hline(yintercept = 0, linetype = 3, colour = "gray30") +
    geom_vline(xintercept = 0, linetype = 3, colour = "gray30") +
    geom_point(shape = shapeHPA, size = sizeHPA, col = colHPA, stroke = 1.5) +
    ggrepel::geom_text_repel(size = 3, seed = 42) +
    labs(tag = expression(bold(B))) +
    theme(plot.tag = element_text(size = 20), plot.title = element_text(hjust = 0.5),
          panel.border = element_rect(colour = "black", fill=NA))
  pl2 <- list(pl2)

  # put together the figure
  mat <- matrix(1:2, nrow = 1)
  ml <- gridExtra::marrangeGrob(c(pl1, pl2), layout_matrix = mat, top = NULL)
  if(pdf) {
    ggsave("canH2OS2.pdf", ml, width = 8, height = 4)
    addexif("canH2OS2", "HPA-TCGA scatterplots for Zc and nH2O", "Dick (2021)")
  }
  else ml
}

#########################
### UNEXPORTED OBJECT ###
#########################

# mapping between HPA and TCGA names
HTmap <- c(
  BRCA = "breast cancer / breast",
  CESC = "cervical cancer / cervix, uterine",
  COAD = "colorectal cancer / colon",
  UCEC = "endometrial cancer / endometrium 1",
  GBM = "glioma / cerebral cortex",
  HNSC = "head and neck cancer / salivary gland",
  LIHC = "liver cancer / liver",
  LUAD = "lung cancer / lung",
  DLBC = "lymphoma / lymph node",
  SKCM = "melanoma / skin 1",
  OV = "ovarian cancer / ovary",
  PAAD = "pancreatic cancer / pancreas",
  PRAD = "prostate cancer / prostate",
  KIRC = "renal cancer / kidney",
  STAD = "stomach cancer / stomach 1",
  TGCT = "testis cancer / testis",
  THCA = "thyroid cancer / thyroid gland",
  BLCA = "urothelial cancer / urinary bladder"
)

############################
### UNEXPORTED FUNCTIONS ###
############################

# function to plot 50 percentile contours / extracted from canH2O2 20191201
contplot <- function(dat, main, col, xvar = "Zc", yvar = "nH2O", dx = NULL, dy = NULL,
                     labtext = NULL, names = NULL) {
  xlab <- canprot::cplab[[paste0("D", strsplit(xvar, "_")[[1]][1])]]
  ylab <- canprot::cplab[[paste0("D", strsplit(yvar, "_")[[1]][1])]]
  if(!is.null(labtext)) {
    xlab <- substitute(xlab ~ "("*labtext*")", list(xlab = xlab[[1]], labtext = labtext))
    ylab <- substitute(ylab ~ "("*labtext*")", list(ylab = ylab[[1]], labtext = labtext))
  }
  plot(c(-0.06, 0.06), c(-0.06, 0.06), xlab = xlab, ylab = ylab, type = "n")
  abline(h = 0, v = 0, lty = 3, col = "grey30")
  if(is.null(dx)) dx <- rep(0, length(dat))
  if(is.null(dy)) dy <- rep(0, length(dat))
  if(is.null(names)) names <- names(dat)
  lapply(1:length(dat), function(i) {
    diffplot(dat[[i]], c(xvar, yvar), pch = NA, pt.text = NA, col.contour = col[i], add = TRUE)
    x <- mean(dat[[i]][, paste0(xvar, ".diff")], na.rm = TRUE)
    y <- mean(dat[[i]][, paste0(yvar, ".diff")], na.rm = TRUE)
    text(x + dx[i], y + dy[i], names[i], col = col[i])
  })
  title(main, font.main = 1)
}

# function to plot linear model
# extracted from canH2OS1 20200224
lmfun <- function(x, y, legend.x = NULL, lmlim = NULL, ...) {
  mylm <- lm(y ~ x)
  if(is.null(lmlim)) lmlim <- range(x)
  lines(lmlim, predict(mylm, data.frame(x = lmlim)), ...)
  # add R-squared text
  if(!is.null(legend.x)) {
    R2 <- format(round(summary(mylm)$r.squared, 2), nsmall = 2)
    R2txt <- substitute(bolditalic(R)^bold("2")*bold(" = "*R2), list(R2 = R2))
    legend(legend.x, legend = R2txt, bty = "n")
  }
  invisible(round(residuals(mylm), 3))
}
