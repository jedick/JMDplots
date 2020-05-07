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
  mat <- matrix(c(1, 1, 1, 1, 1, 1, 3, 1, 2), nrow = 3, byrow = TRUE)
  layout(mat, widths = c(1, 0.5, 1), heights = c(1, 1, 1.5))
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
  pP <- p3 + 2*npos[1] # the position for prostate cancer plot --> changed to phylostratigraphy box
  pO <- p1 + 2*npos[1] # the position for the overview plot

  straightarrow(from = pos[p1, ], to = pos[p2, ], arr.pos = 0.56, endhead = TRUE)
  straightarrow(from = pos[p2, ], to = pos[p3, ], arr.pos = 0.53, endhead = TRUE)
  curvedarrow(from = pos[p3, ], to = pos[p4, ] + c(0, 0.08), curve = 0.5, arr.pos = 0.67, endhead = TRUE)
  curvedarrow(from = pos[p3, ] + c(0, -0.05), to = pos[p5, ] + c(0.08, 0.08), curve = 0.3, arr.pos = 0.725, endhead = TRUE)
  curvedarrow(from = pos[p4, ] + c(-0.105, -0.01), to = pos[pP, ] + c(0, 0.25), curve = 0.2, arr.pos = 0.52, endhead = TRUE)
  curvedarrow(from = pos[p4, ] + c(-0.106, -0.01), to = pos[pO, ] + c(0.05, 0.27), curve = -0.2, arr.pos = 0.40, endhead = TRUE)
  curvedarrow(from = pos[p6, ] + c(-0.07, 0), to = pos[pO, ] + c(-0.02, 0.23), curve = 0.3, arr.pos = 0.55, endhead = TRUE)

  # hard-coded colors (palette.colors is not available in R < 4.0.0) 20200415
  #cols <- palette.colors(8, "Set 2")
  cols <- c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3")
  cex <- 1.3
  ry <- 0.08
  textrect(pos[p1, ], r1, ry, lab = c("Papers", "SI Tables"), cex = cex, box.col = cols[1])
  textplain(pos[p1, ] + c(0, 0.12), lab = "Literature Search", font = 2)
  textrect(pos[p2, ], r2, ry, lab = c("UniProt IDs", "Up   Down", "...      ...   "), cex = cex, box.col = cols[2])
  textplain(pos[p2, ] + c(0, 0.12), lab = c("Differential", "Expression"), font = 2, height = 0.04)
  ZCtext <- quote(italic(Z)[C]~"oxidation state")
  nH2Otext <- quote(italic(n)[H[2]*O]~"hydration state")
  nAAtext <- quote(italic(n)[AA]~"protein length")
  comptext <- as.expression(c(ZCtext, nH2Otext, nAAtext))
  textrect(pos[p3, ], r3, ry, lab = comptext, cex = cex, box.col = cols[2])
  textplain(pos[p3, ] + c(0, 0.12), lab = c("Compositional", "Analysis"), font = 2, height = 0.04)

  # show numbers of datasets 
  cond1 <- c("hypoxia", "secreted", "osmotic_euk", "glucose", "3D")
  cond2 <- c("colorectal", "pancreatic", "breast", "lung", "prostate", "liver")
  vigout <- system.file("vignettes", package = "canprot")
  conddat <- function(cond) read.csv(paste0(vigout, "/", cond, ".csv"), as.is = TRUE)
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
  #col2 <- palette.colors(8, "Classic Tableau")[c(6, 5, 7, 8, 4, 2)]
  col2 <- c("#8C564B", "#9467BD", "#E377C2", "#7F7F7F", "#D62728", "#FF7F0E")
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
  #col1 <- palette.colors(8, "Okabe-Ito")[c(2, 4, 3, 1, 6)]
  col1 <- c(orange = "#E69F00", bluishgreen = "#009E73", skyblue = "#56B4E9", black = "#000000", blue = "#0072B2")
  yculture <- ycoords(pos[p5, ] + c(-0.16, 0.08), lab = culab, height = 0.09)
  xculture <- rep(pos[p5, 1] + 0.13, length(yculture))
  # point salt-inducded hyperosmotic stress on same line as glucose
  points(xculture[-3], yculture[-3], bg = col1[-3], pch = 21, cex = 1.7)
  points(xculture[3] - 0.025, yculture[4], bg = col1[3], pch = 21, cex = 1.7)


  textrect(pos[p6, ] + c(0, 0.08), r6, ry, lab = c("GEPIA2", "(TCGA/", "GTEx)", ""), cex = cex, box.col = cols[1])
  textplain(pos[p6, ] + c(-0.025, 0.025), r6, ry, lab = "HPA", cex = cex)
  textplain(pos[p6, ] + c(0, 0.2), lab = c("Pan-cancer", "comparison"), font = 2, height = 0.04)

  textplain(pos[21, ] + c(-0.02, 0.1), lab = c(paste("Total:", total), "datasets", paste("from", nstudies), "studies"),
            adj = c(0, 0.5), font = 3, height = 0.07, cex = 1.2)

  textrect(pos[pP, ] + c(0.02, 0.095), r7 + 0.015, ry, lab = c("Evolutionary", "trends of", "composition"), cex = cex, box.col = cols[1])
  textplain(pos[pP, ] + c(0.02, 0.215), lab = c("Phylostratigraphic", "analysis"), font = 2, height = 0.04)

  textplain(pos[24, ] + c(0.12, 0.29), lab = "Main finding", font = 2)
  cantext1 <- "Most cancer types have higher"
  cantext2 <- "hydration state of proteins"
  cantext3 <- "compared to normal tissue."
  cantext <- c(cantext1, cantext2, cantext3)
  textplain(pos[24, ] + c(0.12, 0.235), lab = cantext, height = 0.05, font = 4)

  textplain(pos[25, ] + c(0.21, 0.085), lab = "Other findings", font = 2)
  hyptext1 <- quote(italic("Hypoxia experiments show no consistent"))
  hyptext2 <- quote(italic("difference in oxidation state of proteins."))
  hyptext <- as.expression(c(hyptext1, hyptext2))
  textplain(pos[25, ] + c(0.21, 0.03), lab = hyptext, height = 0.04)

  hydtext1 <- quote(italic("Most high-glucose (hyperosmotic) and 3D culture"))
  hydtext2 <- quote(italic("experiments yield lower hydration state of proteins."))
  hydtext <- as.expression(c(hydtext1, hydtext2))
  textplain(pos[25, ] + c(0.25, -0.045), lab = hydtext, height = 0.04)

  ## make diffplot for prostate cancer
  # no longer - just make empty plot to fill space 20200404
  plot.new()

  # make overview plot
  plot.new()
  par(mar = c(3, 2, 1, 1), xaxs = "i", yaxs = "i")
  #plot.window(c(-0.05, 0.05), c(-0.05, 0.05))
  plot.window(c(-1, 1), c(-1, 1))
  axis(1, tck = 0, labels = FALSE)
  axis(2, tck = 0, labels = FALSE)
  mtext("oxidation state", 1, 0.5)
  mtext("hydration state", 2, 0.5)
  # draw arrows to mean values
  # offset x = -0.1 to make room for text
  for(i in 5:1) arrows(-0.1, 0, 35*mean(culture[[i]]$ZC.diff) - 0.1, 35*mean(culture[[i]]$nH2O_rQEC.diff), col = col1[i], lwd = 3, length = 0.15)
  for(i in 1:6) arrows(-0.1, 0, 35*mean(cancer[[i]]$ZC.diff) - 0.1, 35*mean(cancer[[i]]$nH2O_rQEC.diff), col = col2[i], lwd = 3, length = 0.15)
  

  # restore these defaults to be able to re-run this script with expected results
  par(xaxs = "r", yaxs = "r")
  if(pdf) {
    dev.off()
    addexif("canH2O1", "Study overview", "Dick (2020) (preprint)")
  }
}

# median differences of nH2O and ZC for cell culture and cancer tissue 20191126
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
  p12 <- rep(12, 16); p13 <- rep(13, 16)
  # assemble columns (each one is 48 units high)
  c1 <- c(0, p1, p00, p6, p9)
  c2 <- c(0, p2, p00, p7, p10)
  c3 <- c(0, p3, p00, p8, p11)
  c4 <- c(0, p4, p00, p0, p0)
  c5 <- c(0, p4, p12, p13)
  c6 <- c(0, p5, p12, p13)
  mat <- matrix(c(c1, c2, c3, c4, c5, c6), nrow = 48)
  layout(mat, widths = c(1, 1, 1, 0.5, 0.5, 1))
  par(mar = c(4, 4, 1.5, 1), mgp = c(2, 1, 0))

  # read data for all conditions
  cond1 <- c("hypoxia", "secreted", "osmotic_euk", "glucose", "3D")
  cond2 <- c("colorectal", "pancreatic", "breast", "lung", "prostate", "liver")
  conds <- c(cond1, cond2)
  vigout <- system.file("vignettes", package = "canprot")
  alldat <- lapply(conds, function(cond) {
    file <- paste0(vigout, "/", cond, ".csv")
    read.csv(file)
  })
  names(alldat) <- conds

  #col1 <- palette.colors(8, "Okabe-Ito")[c(2, 4, 3, 1, 6)]
  col1 <- c(orange = "#E69F00", bluishgreen = "#009E73", skyblue = "#56B4E9", black = "#000000", blue = "#0072B2")
  #col2 <- palette.colors(8, "Classic Tableau")[c(6, 5, 7, 8, 4, 2)]
  col2 <- c("#8C564B", "#9467BD", "#E377C2", "#7F7F7F", "#D62728", "#FF7F0E")
  col <- c(col1, col2)

  # make scatter plots for nH2O and ZC in cell culture and cancer types
  par(mgp = c(2.5, 1, 0))
  lapply(1:length(alldat), function(i) {
    thisone <- names(alldat)[i]
    if(thisone %in% cond1) {
      # use open/filled symbols for cancer/non-cancer cells 20200103
      pch <- rep(21, nrow(alldat[[i]]))
      pch[grepl("cancer", alldat[[i]]$tags)] <- 19
      # use squares for microbial (yeast and bacteria) datasets 20200407
      pch[grepl("microbial", alldat[[i]]$tags)] <- 0
      pch[grepl("yeast", alldat[[i]]$tags)] <- 0
    } else {
      # use filled/open symbols for human/[mouse or rat] cancer 20200415
      pch <- rep(19, nrow(alldat[[i]]))
      pch[grepl("mouse", alldat[[i]]$tags)] <- 21
      pch[grepl("rat", alldat[[i]]$tags)] <- 21
    }
    diffplot(alldat[i], pch = pch, pt.text = NA, col = col[i], cex = 1, labtext = NA)
    main <- names(alldat)[i]
    if(thisone == "secreted") main <- "hypoxia"
    if(thisone == "osmotic_euk") main <- "(salt)"
    if(thisone == "glucose") main <- "(glucose)"
    title(main, font.main = 1)
    if(thisone == "secreted") title("secreted in", font.main = 1, line = 1.4, xpd = NA)
    if(thisone %in% c("osmotic_euk", "glucose")) title("hyperosmotic", font.main = 1, line = 1.4, xpd = NA)
    if(thisone == "hypoxia") label.figure("A", cex = 2, font = 2, yfrac = 1)
    if(thisone == "colorectal") label.figure("B", cex = 2, font = 2, yfrac = 1)
  })
  
  # compare density contours for cell culture and cancer types 20191126
  # CSV files generated by running the cancer and cell culture vignettes
  vigout <- system.file("vignettes", package = "canprot")
  conddat <- function(cond) read.csv(paste0(vigout, "/", cond, ".csv"), as.is = TRUE)
  culture <- lapply(cond1, conddat); names(culture) <- cond1
  cancer <- lapply(cond2, conddat); names(cancer) <- cond2
  names <- names(culture)
  names[1] <- "hyp-\noxia"
  names[3] <- "salt"
  contplot(culture, "Cell culture", col1, ylim = c(-0.06, 0.04), dx = c(0.027, 0.01, 0.022, 0, -0.02), dy = c(-0.005, 0.03, -0.037, -0.025, -0.03), names = names)
  label.figure("C", cex = 2, font = 2)
  contplot(cancer, "Cancer tissue", col2, ylim = c(-0.04, 0.06), dx = c(NA, 0.012, -0.005, -0.022, -0.006, -0.008), dy = c(NA, -0.035, -0.035, 0.035, -0.02, 0.035))
  text(0.025, -0.008, "CRC", col = col2[1])

  if(pdf) {
    dev.off()
    addexif("canH2O2", "Median differences of protein length, nH2O and ZC for cell culture and cancer tissue", "Dick (2020) (preprint)")
  }
}

# nH2O-ZC and phylostrata plots for TCGA and HPA datasets 20191126
canH2O3 <- function(pdf = FALSE) {
  vigout2 <- system.file("vignettes", package = "JMDplots")
  HPA <- read.csv(file.path(vigout2, "HPA.csv"), as.is = TRUE)
  TCGA <- read.csv(file.path(vigout2, "TCGA.csv"), as.is = TRUE)
  TCGA_labels <- TCGA$description
  HPA_labels <- HPA$description
  HPA_labels <- sapply(strsplit(HPA_labels, " "), "[", 1)
  HPA_labels[grepl("head", HPA_labels)] <- "head and neck"

  # get colors for six cancers in paper 20191208
  cond2 <- c("colorectal", "pancreatic", "breast", "lung", "prostate", "liver")
  #col2 <- palette.colors(8, "Classic Tableau")[c(6, 5, 7, 8, 4, 2)]
  col2 <- c("#8C564B", "#9467BD", "#E377C2", "#7F7F7F", "#D62728", "#FF7F0E")
  jHPA <- match(cond2, sapply(strsplit(HPA$description, " "), "[", 1))
  colHPA <- rep("darkslategray", nrow(HPA))
  colHPA[jHPA] <- col2
  shapeHPA <- rep(15, nrow(HPA))
  shapeHPA[jHPA] <- 1
  # now do TCGA
  TCGAnames <- names(HTmap)[match(cond2, sapply(strsplit(HTmap, " "), "[", 1))]
  jTCGA <- match(TCGAnames, TCGA$description)
  colTCGA <- rep("darkslategray", nrow(TCGA))
  colTCGA[jTCGA] <- col2
  shapeTCGA <- rep(15, nrow(TCGA))
  shapeTCGA[jTCGA] <- 1

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

  TCGAstuff <- densfun(TCGA$ZC.diff, TCGA$nH2O_rQEC.diff)
  TCGAdens <- TCGAstuff$dat
  TCGAlevels <- TCGAstuff$levels
  HPAstuff <- densfun(HPA$ZC.diff, HPA$nH2O_rQEC.diff)
  HPAdens <- HPAstuff$dat
  HPAlevels <- HPAstuff$levels

  # median differences of nH2O-ZC for HPA and TCGA datasets
  # common elements for both plots
  pl1.common <- list(
    theme_bw(),
    xlab(canprot::cplab$DZC),
    ylab(canprot::cplab$DnH2O),
    geom_hline(yintercept = 0, linetype = 3, colour = "gray30"),
    geom_vline(xintercept = 0, linetype = 3, colour = "gray30"),
    theme(plot.tag = element_text(size = 20), plot.title = element_text(hjust = 0.5))
  )
  # create plots
  nudge_x <- ifelse(TCGA_labels %in% c("SKCM"), 0.001, 0)
  # workaround for "no visible binding for global variable ‘ZC.diff’" etc. in R CMD check 20200317
  ZC.diff <- nH2O_rQEC.diff <- z <- NULL
  pl1 <- list(
    ggplot(TCGA, aes(ZC.diff, nH2O_rQEC.diff, label = TCGA_labels)) + 
      pl1.common + geom_point(color = colTCGA, size = 1.5, shape = shapeTCGA, stroke = 1.5) + 
      geom_contour(data = TCGAdens, aes(x, y, z = z), breaks = TCGAlevels, inherit.aes = FALSE, color = "slateblue4", lty = 2) +
      ggrepel::geom_text_repel(size = 2.5, nudge_x = nudge_x, seed = 42, box.padding = 0.12, point.padding = 0.1) +
      ggtitle("TCGA/GTEx") + labs(tag = expression(bold(A))),
    ggplot(HPA, aes(ZC.diff, nH2O_rQEC.diff, label = HPA_labels)) + 
      pl1.common + geom_point(color = colHPA, size = 1.5, shape = shapeHPA, stroke = 1.5) + 
      geom_contour(data = HPAdens, aes(x, y, z = z), breaks = HPAlevels, inherit.aes = FALSE, color = "darkslategray4", lty = 2) +
      ggrepel::geom_text_repel(size = 3, seed = 42) +
      ggtitle("HPA") + labs(tag = expression(bold(B)))
  )

  # table of HPA-TCGA mappings
  iHPA <- match(HTmap, HPA$description)
  iTCGA <- match(names(HTmap), TCGA$description)
  df <- data.frame(x = rep(c(0.8, 2.2), each = 9), y = rep(8:0, 2))
  x <- y <- NULL
  pl2 <- ggplot(df, aes(x, y, label = TCGA_labels[iTCGA])) + xlim(-0.2, 2.8) + ylim(-1.5, 9.5) +
           theme_void() + geom_text(hjust = 1, nudge_x = -0.03) +
           annotate("text", label = HPA_labels[iHPA], x = df$x, y = df$y, hjust = 0) +
           annotate("text", label = "TCGA - HPA pairs", x = 1.5, y = 9, vjust = 0, size = 5) +
           theme(plot.margin = unit(c(0, 0, 0, 5.5), "pt"), plot.tag = element_text(size = 20))
  pl2 <- ggplot_gtable(ggplot_build(pl2))
  pl2$layout$clip[pl2$layout$name == "panel"] <- "off"
  pl2 <- list(pl2)

  # HPA-TCGA scatterplot for PS
  colname <- "PS_TPPG17.diff"
  dat <- data.frame(TCGA = TCGA[iTCGA, colname], HPA = HPA[iHPA, colname])
  # use different symbols for 5 cancers in this paper
  TCGAnames <- names(HTmap)[match(cond2, sapply(strsplit(HTmap, " "), "[", 1))]
  kTCGA <- match(TCGAnames, TCGA$description[iTCGA])
  col <- rep("darkslategray", nrow(dat))
  col[kTCGA] <- col2
  shape <- rep(15, nrow(dat))
  shape[kTCGA] <- 21
  # use bold labels for cancers studied in Trigos et al., 2017
  labels <- TCGA_labels[iTCGA]
  fontface <- ifelse(labels %in% c("LUAD", "LUSC", "BRCA", "PRAD", "LIHC", "COAD", "STAD"), "bold.italic", "plain")

  pl3 <- ggplot(dat, aes(x = TCGA, y = HPA, label = labels)) +
    theme_classic() +
    xlab(quote(Delta*"PS (TCGA/GTEx)")) +
    ylab(quote(Delta*"PS (HPA)")) +
    geom_hline(yintercept = 0, linetype = 3, colour = "gray30") +
    geom_vline(xintercept = 0, linetype = 3, colour = "gray30") +
    geom_point(shape = shape, size = 1.5, col = col, stroke = 1.5) +
    ggrepel::geom_text_repel(size = 3, fontface = fontface, seed = 42) +
    labs(tag = expression(bold(C))) +
    theme(plot.tag = element_text(size = 20), plot.title = element_text(hjust = 0.5),
          panel.border = element_rect(colour = "black", fill=NA))
  pl3 <- list(pl3)

  # put together the figure
  mat <- matrix(c(1, 1, 2, 2, 3, 3, 4, 4), byrow = TRUE, nrow = 2)
  ml <- gridExtra::marrangeGrob(c(pl1, pl2, pl3), layout_matrix = mat, top = NULL, heights = c(2.5, 2))
  if(pdf) {
    ggsave("canH2O3.pdf", ml, width = 10, height = 6.2)
    addexif("canH2O3", "nH2O-ZC and phylostrata plots for TCGA and HPA datasets", "Dick (2020) (preprint)")
  } else ml
}

# compositional analysis of phylostrata; phylostrata differences in cancer 20191201; 20200507
canH2O4 <- function(pdf = FALSE, SI = FALSE) {
  if(pdf) {
    if(SI) pdf("canH2OS4.pdf", width = 9, height = 5)
    else pdf("canH2O4.pdf", width = 9, height = 5)
  }
  par(mar = c(4, 4, 1, 1), mgp = c(2.5, 1, 0))
  mat <- matrix(c(1, 1, 2, 5, 3, 6, 4, 7), nrow = 2)
  layout(mat, widths = c(0.5, 1, 1, 1))

  # make legend for phylostrata
  plot.new()
  par(xpd = NA)
  if(SI) {
    # (Liebeskind et al., 2016)
    ys <- seq(0.85, 0.15, length.out = 8)
    text1 <- c("Cellular_organisms", "Euk_Archaea", "Euk+Bac", "Eukaryota", "Opisthokonta", "Eumetazoa", "Vertebrata", "Mammalia")
    text(-0.6, ys, 1:8, adj = c(1, 0.5))
    text(-0.52, ys, text1, adj = c(0, 0.5))
  } else {
    # (Trigos et al., 2017)
    ys <- seq(0.9, 0.1, length.out = 16)
    text1 <- c("Cellular organisms", "Eukaryota", "Opisthokonta", "Metazoa", "Eumetazoa", "Bilateria", "Chordata", "Euteleostomi")
    text2 <- c("Ammiota", "Mammalia", "Theria", "Eutheria", "Euarchontoglires", "Catarrhini", "Homininae", "")
    text(-0.6, ys, 1:16, adj = c(1, 0.5))
    text(-0.52, ys, c(text1, text2), adj = c(0, 0.5))
    text(-0.52, ys[16], "Homo sapiens", font = 3, adj = c(0, 0.5))
  }
  par(xpd = FALSE)

  # choose PS_source based on plot number
  if(SI) PS_source <- "LMM16" else PS_source <- "TPPG17"

  # plots 1-3: nAA, ZC, nH2O for phylostrata
  memo <- plotphylo("nAA", PS_source = PS_source)
  label.figure("A", font = 2, cex = 1.7, yfrac = 0.96, xfrac = 0.05)
  plotphylo("ZC", PS_source = PS_source, memo = memo)
  plotphylo("nH2O", PS_source = PS_source, memo = memo)

  # get cancer data
  cond2 <- c("colorectal", "pancreatic", "breast", "lung", "prostate", "liver")
  vigout <- system.file("vignettes", package = "canprot")
  conddat <- function(cond) read.csv(paste0(vigout, "/", cond, ".csv"), as.is = TRUE)
  cancer <- lapply(cond2, conddat)
  names(cancer) <- cond2
  #col2 <- palette.colors(8, "Classic Tableau")[c(6, 5, 7, 8, 4, 2)]
  col2 <- c("#8C564B", "#9467BD", "#E377C2", "#7F7F7F", "#D62728", "#FF7F0E")

  # plot 4: 50 percentile contours for nAA-PS of proteomics datasets 20191201
  if(SI) {
    names <- names(cancer)
    names[1] <- "colo-\nrectal"
    contplot(cancer, "", col2,
             xvar = "PS_LMM16", yvar = "nAA", xlim = c(-2, 2), ylim = c(-250, 200),
             dx = c(-0.2, 1.25, -1, -0.65, -0.6, 0.72), dy = c(120, 50, 120, -140, 170, -110), names = names)
  } else {
    contplot(cancer, "", col2,
             xvar = "PS_TPPG17", yvar = "nAA", xlim = c(-4, 5), ylim = c(-250, 200),
             dx = c(1.1, 2.5, -1.8, 0.5, -1.2, -1), dy = c(-210, 60, 140, -140, -185, 130))
  }
  label.figure("B", cex = 1.7, font = 2, yfrac = 0.96, xfrac = 0.05)
  # plot 5: ZC-PS
  if(SI) {
    contplot(cancer, "", col2, xvar = "PS_LMM16", yvar = "ZC", xlim = c(-2, 2), ylim = c(-0.05, 0.05))
  } else {
    contplot(cancer, "", col2, xvar = "PS_TPPG17", yvar = "ZC", xlim = c(-4, 4), ylim = c(-0.05, 0.05))
  }
  # plot 6: nH2O-PS
  if(SI) {
    contplot(cancer, "", col2, xvar = "PS_LMM16", yvar = "nH2O_rQEC", xlim = c(-2, 2), ylim = c(-0.02, 0.05))
  } else {
    contplot(cancer, "", col2, xvar = "PS_TPPG17", yvar = "nH2O_rQEC", xlim = c(-4, 4), ylim = c(-0.02, 0.05))
  }

  if(pdf) {
    dev.off()
    if(SI) addexif("canH2OS4", "Compositional analysis of phylostrata representing Liebeskind et al. gene ages", "Dick (2020) (preprint)")
    else addexif("canH2O4", "Compositional analysis of phylostrata; phylostrata differences in cancer", "Dick (2020) (preprint)")
  }
}


# quantile distributions for proteins coded by differentially expressed genes in aneuploid yeast cells 20200505
canH2O5 <- function(pdf = FALSE) {
  if(pdf) pdf("canH2O5.pdf", width = 6, height = 3)
  pd <- pdat_aneuploidy("TNC+19")
  layout(matrix(c(1, 2, 1, 3), nrow = 2), heights = c(0.2, 1))
  par(mar = c(0, 0, 0, 0), mgp = c(2.5, 1, 0))
  plot.new()
  uptxt <- paste("Proteins coded by", sum(pd$up2), "up-regulated genes")
  dntxt <- paste("Proteins coded by", sum(!pd$up2), "down-regulated genes")
  legend("center", c(dntxt, uptxt), lty = c(1, 2), col = c(1, 2), bty = "n")
  par(mar = c(3.5, 3.5, 0.5, 0.5))
  qdist(pd, "ZC")
  label.figure("A", xfrac = 0.04, yfrac = 1.05, font = 2, cex = 1.5)
  qdist(pd, "nH2O")
  label.figure("B", xfrac = 0.95, yfrac = 1.05, font = 2, cex = 1.5)
  if(pdf) {
    dev.off()
    addexif("canH2O5", "Quantile distributions for proteins coded by differentially expressed genes in aneuploid yeast cells", "Dick (2020) (preprint)")
  }
}

###############
### TABLE 2 ###
###############

# mean differences and p-values across all datasets 20200125
canH2OT2 <- function() {
  cond1 <- c("hypoxia", "secreted", "osmotic_euk", "glucose", "3D")
  cond2 <- c("breast", "colorectal", "liver", "lung", "pancreatic", "prostate")
  cond3 <- c("HPA", "TCGA")
  vigout <- system.file("vignettes", package = "canprot")
  conddat <- function(cond) read.csv(paste0(vigout, "/", cond, ".csv"), as.is = TRUE)
  culture <- lapply(cond1, conddat); names(culture) <- cond1
  cancer <- lapply(cond2, conddat); names(cancer) <- cond2
  vigout2 <- system.file("vignettes", package = "JMDplots")
  conddat <- function(cond) read.csv(paste0(vigout2, "/", cond, ".csv"))
  pancan <- lapply(cond3, conddat); names(pancan) <- cond3
  # make data frames for mean differences and p-values
  # include comparison of up- and down-regulated proteins for "Secreted" in hypoxia compared to "Hypoxia" (whole-cell)
  pdat <- mdat <- data.frame(condition = c(names(culture), names(cancer), names(pancan), "up", "down"))
  for(property in c("ZC", "nH2O_rQEC", "PS_TPPG17", "PS_LMM16", "nAA")) {
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
  # round and combine values
  mdat[, 2:3] <- round(mdat[, 2:3], 3)
  for(icol in 2:3) mdat[, icol] <- paste0(sprintf("%.3f", mdat[, icol]), " (", sprintf("%.1f", pdat[, icol]), ")")
  mdat[, 4:5] <- round(mdat[, 4:5], 2)
  for(icol in 4:5) mdat[, icol] <- paste0(sprintf("%.2f", mdat[, icol]), " (", sprintf("%.1f", pdat[, icol]), ")")
  mdat[, 6] <- round(mdat[, 6], 1)
  for(icol in 6) mdat[, icol] <- paste0(sprintf("%.1f", mdat[, icol]), " (", sprintf("%.1f", pdat[, icol]), ")")
  # adjustments for pretty kable output
  rownames(mdat) <- mdat[, 1]
  mdat <- mdat[, -1]
  colnames(mdat) <- c("&Delta;*Z*~C~", "&Delta;*n*~H2O~", "&Delta;PS (TPPG17)", "&Delta;PS (LMM16)", "&Delta;*n*~AA~")
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
  # subtract 1 H2O to make residues
  out$H2O <- out$H2O - 1
  # adjustments for pretty kable output
  rownames(out) <- out[, 1]
  out <- out[, -1]
  colnames(out) <- gsub("([[:digit:]])", "~\\1~", colnames(out))
  out
}

# rQEC derivation 20190713
canH2OS1 <- function(pdf = FALSE) {
  # set up figure
  if(pdf) pdf("canH2OS1.pdf", width = 6, height = 3)
  par(mfrow = c(1, 2))
  par(mar = c(3.2, 3.2, 1, 1))
  par(mgp = c(2, 0.7, 0))
  par(las = 1)
  par(cex.lab = 1.2)

  # define axis labels
  nH2Olab.QEC <- expression(italic(n)[H[2] * O]~"(QEC)")
  nH2Olab.rQEC <- expression(italic(n)[H[2] * O]~"(rQEC)")
  ZClab <- expression(italic(Z)[C])

  # function to plot values for amino acids
  aaplot <- function(x, y, xlab, ylab, legend.x, lmlim = c(-1, 1)) {
    plot(x, y, type = "p", pch = aminoacids(1), xlab = xlab, ylab = NA)
    mtext(ylab, side = 2, line = 1.8, las = 0)
    lmfun(x, y, legend.x, lmlim)
  }

  # set up amino acid compositions to get compositional values for residues
  AAcomp <- as.data.frame(diag(20))
  names(AAcomp) <- aminoacids(3)

  # plot 1: nH2O-ZC of amino acid residues (QEC)
  # subtract 1 to make residues
  aaplot(ZCAA(AAcomp), H2OAA(AAcomp, "QEC") - 1, ZClab, nH2Olab.QEC, "bottomright")
  label.figure("A", cex = 1.7, font = 2)

  # plot 2: nH2O-ZC of amino acid residues (rQEC)
  par(mgp = c(2, 0.7, 0))
  aaplot(ZCAA(AAcomp), H2OAA(AAcomp, "rQEC") - 1, ZClab, nH2Olab.rQEC, "bottomright")
  label.figure("B", cex = 1.7, font = 2)

  if(pdf) {
    dev.off()
    addexif("canH2OS1", "rQEC derivation", "Dick (2020) (preprint)")
  }
}

# HPA-TCGA scatterplots for ZC and nH2O 20200127
canH2OS2 <- function(pdf = FALSE) {
  vigout2 <- system.file("vignettes", package = "JMDplots")
  HPA <- read.csv(file.path(vigout2, "HPA.csv"), as.is = TRUE)
  TCGA <- read.csv(file.path(vigout2, "TCGA.csv"), as.is = TRUE)
  TCGA_labels <- TCGA$description
  HPA_labels <- HPA$description
  HPA_labels <- sapply(strsplit(HPA_labels, " "), "[", 1)
  HPA_labels[grepl("head", HPA_labels)] <- "head and neck"

  # get colors for 5 cancers in paper 20191208
  cond2 <- c("colorectal", "pancreatic", "breast", "lung", "prostate", "liver")
  #col2 <- palette.colors(8, "Classic Tableau")[c(6, 5, 7, 8, 4, 2)]
  col2 <- c("#8C564B", "#9467BD", "#E377C2", "#7F7F7F", "#D62728", "#FF7F0E")
  jHPA <- match(cond2, sapply(strsplit(HPA$description, " "), "[", 1))
  colHPA <- rep("slateblue4", nrow(HPA))
  colHPA[jHPA] <- col2
  sizeHPA <- rep(1.5, nrow(HPA))
  sizeHPA[jHPA] <- 2
  shapeHPA <- rep(15, nrow(HPA))
  shapeHPA[jHPA] <- 1
  # now do TCGA
  TCGAnames <- names(HTmap)[match(cond2, sapply(strsplit(HTmap, " "), "[", 1))]
  jTCGA <- match(TCGAnames, TCGA$description)
  colTCGA <- rep("slateblue4", nrow(TCGA))
  colTCGA[jTCGA] <- col2
  sizeTCGA <- rep(1.5, nrow(TCGA))
  sizeTCGA[jTCGA] <- 2
  shapeTCGA <- rep(15, nrow(TCGA))
  shapeTCGA[jTCGA] <- 1

  # HPA-TCGA mappings
  iHPA <- match(HTmap, HPA$description)
  iTCGA <- match(names(HTmap), TCGA$description)

  # HPA-TCGA scatterplots for ZC and nH2O
  ZC <- data.frame(TCGA = TCGA$ZC.diff[iTCGA], HPA = HPA$ZC.diff[iHPA])
  nH2O <- data.frame(TCGA = TCGA$nH2O_rQEC.diff[iTCGA], HPA = HPA$nH2O_rQEC.diff[iHPA])

  labels <- TCGA_labels[iTCGA]
  col <- "slateblue4"
  size <- 1.5
  shape <- 15

  r.squared.ZC <- format(summary(lm(HPA ~ TCGA, ZC))$r.squared, digits = 2)
  ZC.title <- paste0("italic(R)^2 == '", r.squared.ZC, "'")
  pl1 <- ggplot(ZC, aes(x = TCGA, y = HPA, label = labels)) +
    theme_classic() + geom_smooth(method = "lm") +
    annotate("text", -Inf, Inf, label = ZC.title, parse = TRUE, hjust = -0.2, vjust = 1.5) +
    xlab(quote(Delta*italic(Z)[C]*" (TCGA/GTEx)")) +
    ylab(quote(Delta*italic(Z)[C]*" (HPA)")) +
    geom_hline(yintercept = 0, linetype = 3, colour = "gray30") +
    geom_vline(xintercept = 0, linetype = 3, colour = "gray30") +
    geom_point(shape = shape, size = size, col = col, stroke = 1.5) +
    ggrepel::geom_text_repel(size = 3, seed = 42) +
    labs(tag = expression(bold(A))) +
    theme(plot.tag = element_text(size = 20), plot.title = element_text(hjust = 0.5),
          panel.border = element_rect(colour = "black", fill=NA))
  pl1 <- list(pl1)

  r.squared.nH2O <- format(summary(lm(HPA ~ TCGA, nH2O))$r.squared, digits = 2)
  nH2O.title <- paste0("italic(R)^2 == '", r.squared.nH2O, "'")
  pl2 <- ggplot(nH2O, aes(x = TCGA, y = HPA, label = labels)) +
    theme_classic() + geom_smooth(method = "lm") +
    annotate("text", -Inf, Inf, label = nH2O.title, parse = TRUE, hjust = -0.2, vjust = 1.5) +
    xlab(quote(Delta*italic(n)[H[2]*O]*" (TCGA/GTEx)")) +
    ylab(quote(Delta*italic(n)[H[2]*O]*" (HPA)")) +
    geom_hline(yintercept = 0, linetype = 3, colour = "gray30") +
    geom_vline(xintercept = 0, linetype = 3, colour = "gray30") +
    geom_point(shape = shape, size = size, col = col, stroke = 1.5) +
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
    addexif("canH2OS2", "HPA-TCGA scatterplots for ZC and nH2O", "Dick (2020) (preprint)")
  }
  else ml
}

# scatterplots of ZC and nH2O vs phylostrata for TCGA and HPA datasets 20200505
canH2OS3 <- function(pdf = FALSE) {
  # read files and get labels
  vigout2 <- system.file("vignettes", package = "JMDplots")
  HPA <- read.csv(file.path(vigout2, "HPA.csv"), as.is = TRUE)
  TCGA <- read.csv(file.path(vigout2, "TCGA.csv"), as.is = TRUE)
  TCGA_labels <- TCGA$description
  HPA_labels <- HPA$description
  HPA_labels <- sapply(strsplit(HPA_labels, " "), "[", 1)
  HPA_labels[grepl("head", HPA_labels)] <- "head and neck"

  # ZC, nH2O and PS values
  HPAdat <- data.frame(ZC = HPA$ZC.diff, nH2O = HPA$nH2O_rQEC.diff, PS = HPA$PS_TPPG17.diff)
  TCGAdat <- data.frame(ZC = TCGA$ZC.diff, nH2O = TCGA$nH2O_rQEC.diff, PS = TCGA$PS_TPPG17.diff)

  # get colors for 5 cancers in paper 20191208
  cond2 <- c("colorectal", "pancreatic", "breast", "lung", "prostate", "liver")
  #col2 <- palette.colors(8, "Classic Tableau")[c(6, 5, 7, 8, 4, 2)]
  col2 <- c("#8C564B", "#9467BD", "#E377C2", "#7F7F7F", "#D62728", "#FF7F0E")

  # workaround for "no visible binding for global variable" in R CMD check 20200505
  PS <- ZC <- nH2O <- NULL

  # TCGA plots
  TCGAnames <- c("COAD", "PAAD", "BRCA", "LUAD", "PRAD", "LIHC")
  jTCGA <- match(TCGAnames, TCGA$description)
  colTCGA <- rep("slateblue4", nrow(TCGA))
  colTCGA[jTCGA] <- col2
  sizeTCGA <- rep(1.5, nrow(TCGA))
  sizeTCGA[jTCGA] <- 2
  shapeTCGA <- rep(15, nrow(TCGA))
  shapeTCGA[jTCGA] <- 1

  r.squared.ZC <- format(summary(lm(ZC ~ PS, TCGAdat))$r.squared, digits = 2)
  ZC.title <- paste0("italic(R)^2 == '", r.squared.ZC, "'")
  pl1 <- ggplot(TCGAdat, aes(x = PS, y = ZC, label = TCGA_labels)) +
    theme_classic() + geom_smooth(method = "lm") +
    annotate("text", -Inf, -Inf, label = ZC.title, parse = TRUE, hjust = -0.2, vjust = -1.5) +
    xlab(quote(Delta*PS*" (TCGA)")) +
    ylab(quote(Delta*italic(Z)[C]*" (TCGA)")) +
    geom_hline(yintercept = 0, linetype = 3, colour = "gray30") +
    geom_vline(xintercept = 0, linetype = 3, colour = "gray30") +
    geom_point(shape = shapeTCGA, size = sizeTCGA, col = colTCGA, stroke = 1.5) +
    ggrepel::geom_text_repel(size = 3, seed = 42) +
    labs(tag = expression(bold(A))) +
    theme(plot.tag = element_text(size = 20), plot.title = element_text(hjust = 0.5),
          panel.border = element_rect(colour = "black", fill=NA))
  pl1 <- list(pl1)

  r.squared.nH2O <- format(summary(lm(nH2O ~ PS, TCGAdat))$r.squared, digits = 2)
  nH2O.title <- paste0("italic(R)^2 == '", r.squared.nH2O, "'")
  pl2 <- ggplot(TCGAdat, aes(x = PS, y = nH2O, label = TCGA_labels)) +
    theme_classic() + geom_smooth(method = "lm") +
    annotate("text", -Inf, -Inf, label = nH2O.title, parse = TRUE, hjust = -0.2, vjust = -1.5) +
    xlab(quote(Delta*PS*" (TCGA)")) +
    ylab(quote(Delta*italic(n)[H[2]*O]*" (TCGA)")) +
    geom_hline(yintercept = 0, linetype = 3, colour = "gray30") +
    geom_vline(xintercept = 0, linetype = 3, colour = "gray30") +
    geom_point(shape = shapeTCGA, size = sizeTCGA, col = colTCGA, stroke = 1.5) +
    ggrepel::geom_text_repel(size = 3, seed = 42) +
    labs(tag = expression(bold(B))) +
    theme(plot.tag = element_text(size = 20), plot.title = element_text(hjust = 0.5),
          panel.border = element_rect(colour = "black", fill=NA))
  pl2 <- list(pl2)

  # HPA plots
  jHPA <- match(cond2, sapply(strsplit(HPA$description, " "), "[", 1))
  colHPA <- rep("slateblue4", nrow(HPA))
  colHPA[jHPA] <- col2
  sizeHPA <- rep(1.5, nrow(HPA))
  sizeHPA[jHPA] <- 2
  shapeHPA <- rep(15, nrow(HPA))
  shapeHPA[jHPA] <- 1

  r.squared.ZC <- format(summary(lm(ZC ~ PS, HPAdat))$r.squared, digits = 2)
  ZC.title <- paste0("italic(R)^2 == '", r.squared.ZC, "'")
  pl3 <- ggplot(HPAdat, aes(x = PS, y = ZC, label = HPA_labels)) +
    theme_classic() + geom_smooth(method = "lm") +
    annotate("text", -Inf, Inf, label = ZC.title, parse = TRUE, hjust = -0.2, vjust = 1.5) +
    xlab(quote(Delta*PS*" (HPA)")) +
    ylab(quote(Delta*italic(Z)[C]*" (HPA)")) +
    geom_hline(yintercept = 0, linetype = 3, colour = "gray30") +
    geom_vline(xintercept = 0, linetype = 3, colour = "gray30") +
    geom_point(shape = shapeHPA, size = sizeHPA, col = colHPA, stroke = 1.5) +
    ggrepel::geom_text_repel(size = 3, seed = 42) +
    labs(tag = expression(bold(C))) +
    theme(plot.tag = element_text(size = 20), plot.title = element_text(hjust = 0.5),
          panel.border = element_rect(colour = "black", fill=NA))
  pl3 <- list(pl3)

  r.squared.nH2O <- format(summary(lm(nH2O ~ PS, HPAdat))$r.squared, digits = 2)
  nH2O.title <- paste0("italic(R)^2 == '", r.squared.nH2O, "'")
  pl4 <- ggplot(HPAdat, aes(x = PS, y = nH2O, label = HPA_labels)) +
    theme_classic() + geom_smooth(method = "lm") +
    annotate("text", -Inf, Inf, label = nH2O.title, parse = TRUE, hjust = -0.2, vjust = 1.5) +
    xlab(quote(Delta*PS*" (HPA)")) +
    ylab(quote(Delta*italic(n)[H[2]*O]*" (HPA)")) +
    geom_hline(yintercept = 0, linetype = 3, colour = "gray30") +
    geom_vline(xintercept = 0, linetype = 3, colour = "gray30") +
    geom_point(shape = shapeHPA, size = sizeHPA, col = colHPA, stroke = 1.5) +
    ggrepel::geom_text_repel(size = 3, seed = 42) +
    labs(tag = expression(bold(D))) +
    theme(plot.tag = element_text(size = 20), plot.title = element_text(hjust = 0.5),
          panel.border = element_rect(colour = "black", fill=NA))
  pl4 <- list(pl4)

  # put together the figure
  mat <- matrix(1:4, nrow = 2, byrow = TRUE)
  ml <- gridExtra::marrangeGrob(c(pl1, pl2, pl3, pl4), layout_matrix = mat, top = NULL)
  if(pdf) {
    ggsave("canH2OS3.pdf", ml, width = 8, height = 8)
    addexif("canH2OS3", "Scatterplots of ZC and nH2O vs phylostrata for TCGA and HPA datasets", "Dick (2020) (preprint)")
  }
  else ml
}

# Compositional analysis of phylostrata representing Liebeskind et al. gene ages 20200507
canH2OS4 <- function(pdf = FALSE) {
  canH2O4(pdf = pdf, SI = TRUE)
}

# Scatterplots of hypoxia scores and ZC or nH2O for TCGA or HPA datasets 20200224
canH2OS5 <- function(pdf = FALSE) {
  if(pdf) pdf("canH2OS5.pdf", width = 6, height = 6)
  # get TCGA and HPA data
  vigout2 <- system.file("vignettes", package = "JMDplots")
  HPA <- read.csv(file.path(vigout2, "HPA.csv"), as.is = TRUE)
  TCGA <- read.csv(file.path(vigout2, "TCGA.csv"), as.is = TRUE)
  # map HPA datasets to TCGA abbreviations
  iHPA <- match(HTmap, HPA$description)
  iTCGA <- match(names(HTmap), TCGA$description)
  HPA$description[iHPA] <- TCGA$description[iTCGA]
  # read hypoxia scores
  file <- system.file("extdata/canH2O/BHL+19_Fig1a.csv", package = "JMDplots")
  scores <- read.csv(file, as.is = TRUE)

  # function to plot values and regression line
  myplot <- function(x, y, xlab, ylab, main = NULL, legend.x = NULL, lmlim = NULL) {
    plot(x, y, type = "p", xlab = xlab, ylab = ylab, main = main)
    lmfun(x, y, legend.x, lmlim)
  }
  # scatterplots of TCGA and HPA metrics vs hypoxia score
  par(mfrow = c(2, 2))
  par(mar = c(4, 4, 2, 1), mgp = c(2.3, 1, 0))
  iTCGA <- match(scores$cancer, TCGA$description)
  myplot(scores$score, TCGA$ZC.diff[iTCGA], xlab = "hypoxia score", ylab = canprot::cplab$DZC, legend.x = "bottomright")
  label.figure("A", cex = 1.7, font = 2, yfrac = 0.96, xfrac = 0.05)
  myplot(scores$score, TCGA$nH2O_rQEC.diff[iTCGA], xlab = "hypoxia score", ylab = canprot::cplab$DnH2O, legend.x = "bottomleft")
  label.figure("B", cex = 1.7, font = 2, yfrac = 0.96, xfrac = 0.05)
  iHPA <- match(scores$cancer, HPA$description)
  myplot(scores$score, HPA$ZC.diff[iHPA], xlab = "hypoxia score", ylab = canprot::cplab$DZC, legend.x = "topleft")
  label.figure("C", cex = 1.7, font = 2, yfrac = 0.96, xfrac = 0.05)
  myplot(scores$score, HPA$nH2O_rQEC.diff[iHPA], xlab = "hypoxia score", ylab = canprot::cplab$DnH2O, legend.x = "topright")
  label.figure("D", cex = 1.7, font = 2, yfrac = 0.96, xfrac = 0.05)
  if(pdf) {
    dev.off()
    addexif("canH2OS5", "Scatterplots of hypoxia scores and ZC or nH2O for TCGA or HPA datasets", "Dick (2020) (preprint)")
  }
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
contplot <- function(dat, main, col, xvar = "ZC", yvar = "nH2O_rQEC", xlim = c(-0.04, 0.04), ylim = c(-0.05, 0.05), dx = NULL, dy = NULL,
                     labtext = NULL, names = NULL) {
  xlab <- canprot::cplab[[paste0("D", strsplit(xvar, "_")[[1]][1])]]
  ylab <- canprot::cplab[[paste0("D", strsplit(yvar, "_")[[1]][1])]]
  if(!is.null(labtext)) {
    xlab <- substitute(xlab ~ "("*labtext*")", list(xlab = xlab[[1]], labtext = labtext))
    ylab <- substitute(ylab ~ "("*labtext*")", list(ylab = ylab[[1]], labtext = labtext))
  }
  plot(xlim, ylim, xlab = xlab, ylab = ylab, type = "n")
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

# mean ZC and nH2O of phylostrata 20191122
# extracted from canH2O4 20191205
plotphylo <- function(vars = c("ZC", "nH2O"), basis = "rQEC", PS_source = "TPPG17", memo = NULL) {
  if(is.null(memo)) {
    dat <- read.csv(system.file("extdata/phylostrata/TPPG17.csv.xz", package = "canprot"), as.is = TRUE)
    if(PS_source == "LMM16") {
      dat <- read.csv(system.file("extdata/phylostrata/LMM16.csv.xz", package = "canprot"), as.is = TRUE)
      colnames(dat)[c(1,3)] <- c("Entry", "Phylostrata")
      # remove entries that have ENSP instead of UniProt IDs
      dat <- dat[!grepl("^ENSP", dat$Entry), ]
    }
    dat <- check_IDs(dat, "Entry")
    dat <- cleanup(dat, "Entry")
    pcomp <- protcomp(dat$Entry, basis = basis)
  } else {
    dat <- memo$dat
    pcomp <- memo$pcomp
  }
  nH2O <- pcomp$residue.basis[, "H2O"]
  ZC <- pcomp$ZC
  nAA <- pcomp$protein.length
  # get mean ZC and nH2O for each phylostratum
  Phylostratum <- sort(unique(dat$Phylostrata))
  cum.ZC <- cum.nH2O <- cum.nAA <- mean.ZC <- mean.nH2O <- mean.nAA <- numeric()
  for(p in Phylostratum) {
    # point mean
    mean.ZC <- c(mean.ZC, mean(ZC[dat$Phylostrata == p]))
    mean.nH2O <- c(mean.nH2O, mean(nH2O[dat$Phylostrata == p]))
    mean.nAA <- c(mean.nAA, mean(nAA[dat$Phylostrata == p]))
    # cumulative mean
    cum.ZC <- c(cum.ZC, mean(ZC[dat$Phylostrata <= p]))
    cum.nH2O <- c(cum.nH2O, mean(nH2O[dat$Phylostrata <= p]))
    cum.nAA <- c(cum.nAA, mean(nAA[dat$Phylostrata <= p]))
  }
  if("ZC" %in% vars) {
    plot(Phylostratum, mean.ZC, type = "b", ylab = expression(italic(Z)[C]))
    lines(Phylostratum, cum.ZC, col = 2)
  }
  if("nH2O" %in% vars) {
    ylab <- expression(italic(n)[H[2] * O])
    plot(Phylostratum, mean.nH2O, type = "b", ylab = ylab)
    lines(Phylostratum, cum.nH2O, col = 2)
  }
  if("nAA" %in% vars) {
    plot(Phylostratum, mean.nAA, type = "b", ylab = expression(italic(n)[AA]))
    lines(Phylostratum, cum.nAA, col = 2)
  }
  # return the dat and pcomp for memoization 20191211
  invisible(list(dat = dat, pcomp = pcomp))
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
