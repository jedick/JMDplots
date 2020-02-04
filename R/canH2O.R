# JMDplots/canH2O.R
# make plots for the paper:
# Water as a reactant in the differential expression of proteins in cancer
# 20191126 jmd first version
# 20200203 started moving functions to JMDplots

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
  p5 <- 16; r5 <- 0.125
  p6 <- 13; r6 <- 0.07
  p7 <- 28; r7 <- 0.09
  pP <- p3 + 2*npos[1] # the position for prostate cancer plot
  pO <- p1 + 2*npos[1] # the position for the overview plot

  straightarrow(from = pos[p1, ], to = pos[p2, ], arr.pos = 0.56, endhead = TRUE)
  straightarrow(from = pos[p2, ], to = pos[p3, ], arr.pos = 0.53, endhead = TRUE)
  curvedarrow(from = pos[p3, ], to = pos[p4, ] + c(0, 0.08), curve = 0.5, arr.pos = 0.67, endhead = TRUE)
  curvedarrow(from = pos[p3, ] + c(0, -0.05), to = pos[p5, ] + c(0.08, 0.08), curve = 0.3, arr.pos = 0.725, endhead = TRUE)
  curvedarrow(from = pos[p4, ] + c(-0.129, -0.01), to = pos[pP, ] + c(0, 0.25), curve = 0.2, arr.pos = 0.55, endhead = TRUE)
  curvedarrow(from = pos[p4, ] + c(-0.13, -0.01), to = pos[pO, ] + c(0, 0.27), curve = -0.2, arr.pos = 0.35, endhead = TRUE)
  curvedarrow(from = pos[p6, ] + c(-0.07, 0), to = pos[pO, ] + c(0, 0.23), curve = 0.3, arr.pos = 0.55, endhead = TRUE)
  curvedarrow(from = pos[pO, ] + c(0.16, 0.19), to = pos[p7, ] + c(-0.03, 0.23), curve = 0.3, arr.pos = 0.4, endhead = TRUE)

  cols <- palette.colors(8, "Set 2")
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
  cond1 <- c("hypoxia", "hyperosmotic", "secreted", "3D")
  cond2 <- c("colorectal", "pancreatic", "breast", "lung", "prostate")
  vigout <- system.file("extdata/vignette_output", package = "canprot")
  conddat <- function(cond) read.csv(paste0(vigout, "/", cond, ".csv"), as.is = TRUE)
  culture <- lapply(cond1, conddat); names(culture) <- cond1
  cancer <- lapply(cond2, conddat); names(cancer) <- cond2
  total <- sum(sapply(culture, nrow), sapply(cancer, nrow))
  # calculate number of studies
  cultsets <- unlist(sapply(culture, "[[", "dataset"))
  cansets <- unlist(sapply(cancer, "[[", "dataset"))
  nstudies <- length(unique(sapply(strsplit(c(cultsets, cansets), "_"), "[", 1)))

  textrect(pos[p4, ] + c(0, 0.08), r4, ry, lab = "", cex = cex, box.col = cols[4])
  calab <- paste(sapply(cancer, nrow), cond2)
  textplain(pos[p4, ] + c(0, 0.08), lab = calab, height = 0.09, cex = cex)
  textplain(pos[p4, ] + c(0, 0.2), lab = c("Cancer", "vs normal"), font = 2, height = 0.04)

  textrect(pos[p5, ] + c(0, 0.08), r5, ry, lab = "", cex = cex, box.col = cols[3])
  # write "secreted in hypoxia" and "3D growth" 20200118
  culab <- paste(sapply(culture, nrow), cond1)
  shlab <- c(paste(culab[3], "in"), "hypoxia")
  textplain(pos[p5, ] + c(0, 0.065), lab = shlab, height = 0.038, cex = cex)
  culab <- c(culab[1:2], "", "", culab[4])
  culab[5] <- paste(culab[5], "growth")
  textplain(pos[p5, ] + c(0, 0.08), lab = culab, height = 0.09, cex = cex)
  textplain(pos[p5, ] + c(0, 0.2), lab = c("Cell", "culture"), font = 2, height = 0.04)

  textrect(pos[p6, ] + c(0, 0.08), r6, ry, lab = c("TCGA", "GTEx", "GEPIA", "HPA"), cex = cex, box.col = cols[1])
  textplain(pos[p6, ] + c(0, 0.2), lab = c("Pan-cancer", "comparison"), font = 2, height = 0.04)

  textrect(pos[p7, ] + c(0.01, 0.035), r7, ry, lab = c("Evolutionary", "trends of", "composition"), cex = cex, box.col = cols[1])
  textplain(pos[p7, ] + c(0.01, 0.155), lab = c("Phylostratigraphic", "analysis"), font = 2, height = 0.04)

  textplain(pos[21, ] + c(-0.04, 0.08), lab = c(paste("Total", total), "datasets", paste("from", nstudies), "studies"), adj = c(0, 0.5), font = 3, height = 0.06)
  textplain(pos[pP, ] + c(0.07, 0.25), lab = c("Example: 23 datasets", "for prostate cancer"), font = 3, height = 0.035)

  cantext1 <- quote(italic("Up-regulated proteins in most"))
  cantext2 <- quote(italic("cancers are relatively hydrated."))
  cantext <- as.expression(c(cantext1, cantext2))
  textplain(pos[24, ] + c(0.12, 0.23), lab = cantext, height = 0.04)
  textplain(pos[24, ] + c(0.12, 0.28), lab = "Main finding", font = 2)

  celtext1 <- quote(italic("Hydration state of"))
  celtext2 <- quote(italic("proteins decreases"))
  celtext3 <- quote(italic("in hyperosmotic"))
  celtext4 <- quote(italic("and 3D culture."))
  celtext <- as.expression(c(celtext1, celtext2, celtext3, celtext4))
  textplain(pos[25, ] + c(0.08, -0.03), lab = celtext, height = 0.06)

  # make diffplot for prostate cancer
  par(mar = c(4, 4, 1.5, 1), mgp = c(2, 1, 0))
  dat <- read.csv(paste0(vigout, "/prostate.csv"), as.is = TRUE)
  pch <- rep(19, nrow(dat))
  # use triangle for transcriptomes 20191202
  pch[grepl("transcriptome", dat$tags)] <- 17
  col <- palette.colors(8, "Classic Tableau")[4]
  diffplot(dat, pch = pch, pt.text = NA, col = col, cex = 1)
  #title("Prostate cancer", font.main = 1)

  # make overview plot
  plot.new()
  par(mar = c(3, 2, 1, 1), xaxs = "i", yaxs = "i")
  #plot.window(c(-0.05, 0.05), c(-0.05, 0.05))
  plot.window(c(-1, 1), c(-1, 1))
  axis(1, tck = 0, labels = FALSE)
  axis(2, tck = 0, labels = FALSE)
  mtext("oxidation state", 1, 0.5)
  mtext("hydration state", 2, 0.5)
  col1 <- palette.colors(8, "Okabe-Ito")[c(2:4, 6)]
  col2 <- palette.colors(8, "Classic Tableau")[c(6, 5, 7, 8, 4)]
  # draw arrows to mean values
  # offset x = -0.1 to make room for text
  for(i in 1:4) arrows(-0.1, 0, 35*mean(culture[[i]]$ZC.diff) - 0.1, 35*mean(culture[[i]]$nH2O_rQEC.diff), col = col1[i], lwd = 3, length = 0.15)
  for(i in 1:5) arrows(-0.1, 0, 35*mean(cancer[[i]]$ZC.diff) - 0.1, 35*mean(cancer[[i]]$nH2O_rQEC.diff), col = col2[i], lwd = 3, length = 0.15)

  # restore these defaults to be able to re-run this script with expected results
  par(xaxs = "r", yaxs = "r")
  if(pdf) {
    dev.off()
    addexif("canH2O1", "Study overview", "Dick (2020) (preprint)")
  }
}

# median differences of protein length, nH2O and ZC for cell culture and cancer 20191126
canH2O2 <- function(pdf = FALSE) {
  if(pdf) pdf("canH2O2.pdf", width = 9, height = 6)
  # layout with spaces between groups of plots
  p0 <- rep(0, 3)
  p1 <- rep(1, 15); p2 <- rep(2, 15); p3 <- rep(3, 15)
  p4 <- rep(4, 15); p5 <- rep(5, 15); p6 <- rep(6, 15)
  p7 <- rep(7, 15); p8 <- rep(8, 15); p9 <- rep(9, 15); p10 <- rep(10, 15)
  p11 <- rep(11, 16); p12 <- rep(12, 16); p13 <- rep(13, 16)
  mat <- matrix(c(p1, p0, p3, p7, p1, p0, p4, p8, p2, p0, p5, p9, p2, p0, p6, p10, rep(0, 48), p11, p12, p13), nrow = 48)
  layout(mat, widths = c(1, 1, 1, 1, 0.2, 1.3))
  par(mar = c(4, 4, 1.5, 1), mgp = c(2, 1, 0))

  # read data for all conditions
  cond1 <- c("hypoxia", "hyperosmotic", "secreted", "3D")
  cond2 <- c("colorectal", "pancreatic", "breast", "lung", "prostate")
  conds <- c(cond1, cond2)
  vigout <- system.file("extdata/vignette_output", package = "canprot")
  alldat <- lapply(conds, function(cond) {
    file <- paste0(vigout, "/", cond, ".csv")
    read.csv(file)
  })
  names(alldat) <- conds

  col1 <- palette.colors(8, "Okabe-Ito")[c(2:4, 6)]
  col2 <- palette.colors(8, "Classic Tableau")[c(6, 5, 7, 8, 4)]
  col <- c(col1, col2)

  # make beanplots for protein length in cell culture and cancer datasets
  nAA <- lapply(alldat, "[[", "nAA.diff")
  beanplot_mod(nAA[1:4], main = "Cell culture", ylab = "protein length\n(median difference)", col = as.list(col1), xaxt = "n", font.main = 1)
  ## https://www.r-bloggers.com/rotated-axis-labels-in-r-plots/
  axis(1, at = 1:4, labels = FALSE)
  # write "secreted in hypoxia" 20200118
  labels <- names(nAA)[1:4]
  labels[3] <- "secreted in\nhypoxia"
  text(x = 1:4, y = par()$usr[3]-0.15*(par()$usr[4]-par()$usr[3]), labels = labels, srt = 30, adj = 1, xpd = TRUE)
  label.figure("A", cex = 2, font = 2, xfrac = 0.03)
  beanplot_mod(nAA[5:9], main = "Cancer tissue", ylab = "protein length\n(median difference)", col = as.list(col2), xaxt = "n", font.main = 1)
  axis(1, at = 1:5, labels = FALSE)
  text(x = 1:5, y = par()$usr[3]-0.15*(par()$usr[4]-par()$usr[3]), labels = names(nAA)[5:9], srt = 30, adj = 1, xpd = TRUE)

  # make scatter plots for nH2O and ZC in cell culture and cancer types (except prostate which is shown in Fig. 1)
  par(mgp = c(2.5, 1, 0))
  lapply(1:8, function(i) {
    if(i %in% 1:4) {
      # use open/filled symbols for cancer/non-cancer cells 20200103
      pch <- rep(21, nrow(alldat[[i]]))
      pch[grepl("cancer", alldat[[i]]$tags)] <- 19
      # use triangles for transcriptomes 20191202
      pch[grepl("transcriptome", alldat[[i]]$tags)] <- 2
      pch[grepl("transcriptome", alldat[[i]]$tags) & grepl("cancer", alldat[[i]]$tags)] <- 17
    } else {
      pch <- rep(19, nrow(alldat[[i]]))
      pch[grepl("transcriptome", alldat[[i]]$tags)] <- 17
    }
    # use gray color for transcriptomes 20200121
    thiscol <- rep(col[i], nrow(alldat[[i]]))
    thiscol[grepl("transcriptome", alldat[[i]]$tags)] <- "darkslategray"
    diffplot(alldat[i], pch = pch, pt.text = NA, col = thiscol, cex = 1, labtext = NA)
    main <- names(alldat)[i]
    if(i==3) main <- "hypoxia"
    title(main, font.main = 1)
    if(i==3) title("secreted in", font.main = 1, line = 1.2, xpd = NA)
    if(i==1) label.figure("B", cex = 2, font = 2, yfrac = 1)
  })
  
  # compare density contours for cell culture, cancer types, and TCGA and HPA data 20191126
  cond3 <- c("HPA", "TCGA")
  # CSV files generated by running the cancer and cell culture vignettes
  vigout <- system.file("extdata/vignette_output", package = "canprot")
  conddat <- function(cond) read.csv(paste0(vigout, "/", cond, ".csv"), as.is = TRUE)
  culture <- lapply(cond1, conddat); names(culture) <- cond1
  cancer <- lapply(cond2, conddat); names(cancer) <- cond2
  # CSV files generated by running the TCGA and HPA vignettes
  vigout2 <- system.file("extdata/vignette_output", package = "JMDplots")
  conddat <- function(cond) read.csv(paste0(vigout2, "/", cond, ".csv"))
  pancan <- lapply(cond3, conddat); names(pancan) <- cond3
  names <- names(culture)
  names[1] <- "hyp-\noxia"
  contplot(culture, "Cell culture", col1, ylim = c(-0.06, 0.04), dx = c(0.023, 0.022, 0.015, -0.012), dy = c(-0.005, -0.025, 0.028, -0.01), names = names)
  label.figure("C", cex = 2, font = 2)
  contplot(cancer, "Cancer tissue", col2, ylim = c(-0.04, 0.06), dx = c(0.015, 0.012, -0.005, -0.025, -0.006), dy = c(0.025, -0.035, -0.035, 0.035, -0.02))
  contplot(pancan, "Pan-cancer", c("darkslategray4", "slateblue4"), dx = c(0.015, -0.02), dy = c(0.012, 0.012))

  if(pdf) {
    dev.off()
    addexif("canH2O2", "Median differences of protein length, nH2O and ZC for cell culture and cancer", "Dick (2020) (preprint)")
  }
}

# nH2O-ZC and phylostrata plots for TCGA and HPA datasets 20191126
canH2O3 <- function(pdf = FALSE) {
  vigout2 <- system.file("extdata/vignette_output", package = "JMDplots")
  HPA <- read.csv(file.path(vigout2, "HPA.csv"), as.is = TRUE)
  TCGA <- read.csv(file.path(vigout2, "TCGA.csv"), as.is = TRUE)
  TCGA_labels <- TCGA$description
  HPA_labels <- HPA$description
  HPA_labels <- sapply(strsplit(HPA_labels, " "), "[", 1)
  HPA_labels[grepl("head", HPA_labels)] <- "head and neck"

  # get colors for five cancers in paper 20191208
  cond2 <- c("colorectal", "pancreatic", "breast", "lung", "prostate")
  col2 <- palette.colors(8, "Classic Tableau")[c(6, 5, 7, 8, 4)]
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

  # median differences of nH2O-ZC for HPA and TCGA datasets
  # common elements for both plots
  pl1.common <- list(
    theme_bw(),
    xlab(canprot::cplab$DZC),
    ylab(canprot::cplab$DnH2O),
#    coord_cartesian(xlim = range(HPA$ZC.diff, TCGA$ZC.diff),
#                    ylim = range(HPA$nH2O_rQEC.diff, TCGA$nH2O_rQEC.diff)),
    geom_hline(yintercept = 0, linetype = 3, colour = "gray30"),
    geom_vline(xintercept = 0, linetype = 3, colour = "gray30"),
    theme(plot.tag = element_text(size = 20), plot.title = element_text(hjust = 0.5))
  )
  # create plots
  nudge_x <- ifelse(TCGA_labels %in% c("SKCM"), 0.001, 0)
  pl1 <- list(
    ggplot(TCGA, aes(TCGA$ZC.diff, TCGA$nH2O_rQEC.diff, label = TCGA_labels)) + ggrepel::geom_text_repel(size = 2.5, nudge_x = nudge_x) +
      pl1.common + geom_point(color = colTCGA, size = sizeTCGA, shape = shapeTCGA, stroke = 1.5) + ggtitle("TCGA/GTEx") + labs(tag = expression(bold(A))),
    ggplot(HPA, aes(HPA$ZC.diff, HPA$nH2O_rQEC.diff, label = HPA_labels)) + ggrepel::geom_text_repel(size = 3) +
      pl1.common + geom_point(color = colHPA, size = sizeHPA, shape = shapeHPA, stroke = 1.5) + ggtitle("HPA") + labs(tag = expression(bold(B)))
  )

  # table of HPA-TCGA mappings
  iHPA <- match(HTmap, HPA$description)
  iTCGA <- match(names(HTmap), TCGA$description)
  df <- data.frame(x = rep(c(0.8, 2.2), each = 9), y = rep(8:0, 2))
  pl2 <- ggplot(df, aes(df$x, df$y, label = TCGA_labels[iTCGA])) + xlim(-0.2, 2.8) + ylim(-1.5, 9.5) +
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
  col <- rep("slateblue4", nrow(dat))
  col[kTCGA] <- col2
  size <- rep(1.5, nrow(dat))
  size[kTCGA] <- 2
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
    geom_point(shape = shape, size = size, col = col, stroke = 1.5) +
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

# compositional analysis of phylostrata and phylostrata-nAA plots for cancer tissue 20191201
canH2O4 <- function(pdf = FALSE) {
  if(pdf) pdf("canH2O4.pdf", width = 9, height = 5)
  par(mar = c(4, 3, 1, 1), mgp = c(2, 1, 0))
  mat <- matrix(1:10, byrow = TRUE, nrow = 2)
  layout(mat, widths = c(1, 1, 1, 1.3, 1.5))

  # plots 1-3: ZC, nH2O, nAA for Trigos et al. phylostrata
  memo <- plotphylo("ZC", PS_source = "TPPG17")
  label.figure("A", font = 2, cex = 1.7, yfrac = 0.96, xfrac = 0.05)
  plotphylo("nH2O", PS_source = "TPPG17", memo = memo)
  label.figure("B", font = 2, cex = 1.7, yfrac = 0.96, xfrac = 0.05)
  plotphylo("nAA", PS_source = "TPPG17", memo = memo)
  label.figure("C", font = 2, cex = 1.7, yfrac = 0.96, xfrac = 0.05)

  # legend for phylostrata (Trigos et al., 2017)
  plot.new()
  par(xpd = NA)
  ys <- seq(0.8, 0, length.out = 8)
  text1 <- c("Cellular organisms", "Eukaryota", "Opisthokonta", "Metazoa", "Eumetazoa", "Bilateria", "Chordata", "Euteleostomi")
  text2 <- c("Ammiota", "Mammalia", "Theria", "Eutheria", "Euarchontoglires", "Catarrhini", "Homininae", "")
  text(-0.25, ys, 1:8, adj = c(1, 0.5))
  text(-0.2, ys, text1, adj = c(0, 0.5))
  text(0.55, ys, 9:16, adj = c(1, 0.5))
  text(0.6, ys, text2, adj = c(0, 0.5))
  text(0.6, ys[8], "Homo sapiens", font = 3, adj = c(0, 0.5))
  text(0.35, 1, "Phylostrata\n(Trigos et al., 2017)", cex = 1.2)
  par(xpd = FALSE)

  # plot 4: 50 percentile contours for nAA-PS of proteomics datasets 20191201
  cond2 <- c("colorectal", "pancreatic", "breast", "lung", "prostate")
  vigout <- system.file("extdata/vignette_output", package = "canprot")
  conddat <- function(cond) read.csv(paste0(vigout, "/", cond, ".csv"), as.is = TRUE)
  cancer <- lapply(cond2, conddat)
  names(cancer) <- cond2
  col2 <- palette.colors(8, "Classic Tableau")[c(6, 5, 7, 8, 4)]
  opar <- par(mar = c(4, 4, 1, 1), mgp = c(2.3, 1, 0), las = 1)
  contplot(cancer, "", col2,
           xvar = "nAA", yvar = "PS_TPPG17", xlim = c(-250, 200), ylim = c(-5, 5),
           dx = c(-130, 60, 140, -130, -190), dy = c(1.7, 2.5, -1.9, 0.7, -1.3))
  label.figure("D", cex = 1.7, font = 2, yfrac = 0.96, xfrac = 0.05)
  par(opar)

  # plots 5-7: ZC, nH2O, nAA for Liebeskin et al. phylostrata
  memo <- plotphylo("ZC", PS_source = "LMM16")
  label.figure("E", font = 2, cex = 1.7, yfrac = 0.96, xfrac = 0.05)
  plotphylo("nH2O", PS_source = "LMM16", memo = memo)
  label.figure("F", font = 2, cex = 1.7, yfrac = 0.96, xfrac = 0.05)
  plotphylo("nAA", PS_source = "LMM16", memo = memo)
  label.figure("G", font = 2, cex = 1.7, yfrac = 0.96, xfrac = 0.05)

  # legend for phylostrata (Liebeskind et al., 2016)
  plot.new()
  par(xpd = NA)
  ys <- seq(0.8, 0, length.out = 8)
  text1 <- c("Cellular_organisms", "Euk_Archaea", "Euk+Bac", "Eukaryota", "Opisthokonta", "Eumetazoa", "Vertebrata", "Mammalia")
  text(0.15, ys, 1:8, adj = c(1, 0.5))
  text(0.2, ys, text1, adj = c(0, 0.5))
  text(0.35, 1, "Phylostrata\n(Liebeskind et al., 2016)", cex = 1.2)
  par(xpd = FALSE)

  # plot 8: 50 percentile contours for nAA-PS of proteomics datasets 20191201
  names <- names(cancer)
  names[1] <- "colo-\nrectal"
  opar <- par(mar = c(4, 4, 1, 1), mgp = c(2.3, 1, 0), las = 1)
  contplot(cancer, "", col2,
           xvar = "nAA", yvar = "PS_LMM16", xlim = c(-250, 200), ylim = c(-2, 2),
           dx = c(130, 50, 10, -90, 120), dy = c(-0.7, 1.25, 0.3, 0.85, -0.7), names = names)
  par(opar)
  label.figure("H", cex = 1.7, font = 2, yfrac = 0.96, xfrac = 0.05)

  if(pdf) {
    dev.off()
    addexif("canH2O4", "Compositional analysis of phylostrata and phylostrata-nAA plots for cancer tissue", "Dick (2020) (preprint)")
  }
}

# plots with nH2O and nO2 of amino acid biosynthesis reactions 20191205
canH2O5 <- function(pdf = FALSE) {
  if(pdf) pdf("canH2O5.pdf", width = 6, height = 6)
  par(mar = c(4, 4.1, 1.5, 1), mgp = c(2.8, 1, 0), las = 1)
  par(mfrow = c(2, 2))

  # plot nO2(biosynth) vs ZC for TCGA data
  ZClab <- quote(Delta*italic(Z)[C]~"(carbon oxidation state)")
  O2lab <- quote(Delta*italic(n)[O[2]]~"(AA biosynthesis)")
  dat <- read.csv(system.file("extdata/vignette_output/TCGA.csv", package = "JMDplots"), as.is = TRUE)
  plot(dat$ZC.diff, dat$nO2_biosynth.diff, xlab = ZClab, ylab = O2lab)
  abline(h = 0, v = 0, lty = 3, col = "grey30")
  title("TCGA/GTEx", font.main = 1)
  label.figure("A", font = 2, cex = 1.7, yfrac = 0.97, xfrac = 0.03)

  # plot nH2O(biosynth) vs nH2O(rQEC) for TCGA data
  rQEClab <- quote(Delta*italic(n)[H[2] * O]~"(stoichiometric hydration state)    ")
  biolab <- quote(Delta*italic(n)[H[2] * O]~"(AA biosynthesis)")
  plot(dat$nH2O_rQEC.diff, dat$nH2O_biosynth.diff, xlab = rQEClab, ylab = biolab)
  abline(h = 0, v = 0, lty = 3, col = "grey30")
  title("TCGA/GTEx", font.main = 1)
  label.figure("B", font = 2, cex = 1.7, yfrac = 0.97, xfrac = 0.03)

  vigout <- system.file("extdata/vignette_output", package = "canprot")
  conddat <- function(cond) read.csv(paste0(vigout, "/", cond, ".csv"), as.is = TRUE)

  cond2 <- c("colorectal", "pancreatic", "breast", "lung", "prostate")
  cancer <- lapply(cond2, conddat); names(cancer) <- cond2
  col2 <- palette.colors(8, "Classic Tableau")[c(6, 5, 7, 8, 4)]
  contplot(cancer, "Cancer", col2, "nO2_biosynth", "nH2O_biosynth", xlim = c(-0.15, 0.1), ylim = c(-0.1, 0.1),
           dx = c(-0.05, 0.05, 0.052, -0.06, 0.045), dy = c(-0.03, -0.065, 0.022, -0.012, 0.05), labtext = "AA biosynthesis")
  label.figure("C", font = 2, cex = 1.7, yfrac = 0.97, xfrac = 0.03)

  cond1 <- c("hypoxia", "hyperosmotic", "secreted", "3D")
  culture <- lapply(cond1, conddat); names(culture) <- cond1
  col1 <- palette.colors(8, "Okabe-Ito")[c(2:4, 6)]
  names <- names(culture)
  names[1] <- "hyp-\noxia"
  names[2] <- "hyper-\nosmotic"
  contplot(culture, "Cell culture", col1, "nO2_biosynth", "nH2O_biosynth", xlim = c(-0.15, 0.1), ylim = c(-0.1, 0.1),
           dx = c(-0.05, 0.05, -0.1, 0.004), dy = c(0.072, -0.06, -0.025, 0.01), labtext = "AA biosynthesis", names = names)
  label.figure("D", font = 2, cex = 1.7, yfrac = 0.97, xfrac = 0.03)

  if(pdf) {
    dev.off()
    addexif("canH2O5", "Plots with nH2O and nO2 of amino acid biosynthesis reactions", "Dick (2020) (preprint)")
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


# modified beanplot.R from package 'beanplot' version 1.2
# modifications by Jeffrey M. Dick (jmd)
# 20190401 modified beanplot function to plot zero line, not overall average
# 20191203 also modified to accept a font.main argument

#`beanplot` <-
`beanplot_mod` <-
function (..., bw = "SJ-dpi", kernel = "gaussian", cut = 3, cutmin = -Inf, 
    cutmax = Inf, grownage = 10, what = c(TRUE, TRUE, TRUE, TRUE), 
    add = FALSE, col, axes = TRUE, log = "auto", handlelog = NA, 
    ll = 0.16, wd = NA, maxwidth = 0.8, maxstripline = 0.96, 
    method = "stack", names, overallline = "mean", beanlines = overallline, 
    horizontal = FALSE, side = "no", jitter = NULL, beanlinewd = 2, 
    frame.plot = axes, border = NULL, innerborder = NA, at = NULL, 
    boxwex = 1, ylim = NULL, xlim = NULL, show.names = NA) 
{
    #internal functions (later on, mlog and mexp will be defined)
    mdensityxy <- function(x) {
        if (length(x) > 0) {
            from <- max(cutmin, (min(mlog(x)) - cut * bw))
            to <- min(cutmax, max(mlog(x)) + cut * bw)
            stats::density(mlog(x), bw = bw, kernel = kernel, from = from, 
                to = to)[c("x", "y")]
        }
        else list(x = numeric(), y = numeric())
    }

    #get and store function arguments
    args <- match.call()
    mcall <- as.list(args)

    #settings with multiple options
    method <- pmatch(method, c("overplot", "stack", "jitter"))
    if (is.na(method) || method == 0) 
        stop("invalid plotting method")
    beanlines <- pmatch(beanlines, c("mean", "median", "quantiles"))
    if (is.na(beanlines) || beanlines == 0) 
        stop("invalid beanlines")
    overallline <- pmatch(overallline, c("mean", "median"))
    if (is.na(overallline) || overallline == 0) 
        stop("invalid overallline")
    side <- pmatch(side, c("no", "first", "second", "both"))
    if (is.na(side) || side == 0) 
        stop("invalid side")

    #get the groups dataset, we will generate one bean(-side) per group
    #and set the name and position settings
    groups <- getgroupsfromarguments(args)
    groups <- lapply(groups, stats::na.omit)
    n <- length(groups)
    displayn <- if (side == 4) 
        ceiling(n/2)
    else n
    if (n == 0) 
        stop("no data found to beanplot")
    if (missing(names)) {
        if (is.null(base::names(groups))) 
            attr(groups, "names") = 1:displayn
        names <- base::names(groups)
    }
    else {
        attr(groups, "names") <- names
        if (is.na(show.names)) 
            show.names <- TRUE
    }
    if (is.null(at)) {
        at <- 1:displayn
    }
    if ((side == 4) && (length(names) > length(at))) {
        for (i in 1:length(at)) {
            names[i] <- makecombinedname(names[i * 2 - 1], names[i * 
                2])
        }
        length(names) <- length(at)
    }

    #color settings
    combinedpolygons <- ((side == 4) && (length(border) < 2) && 
        (n%%2 == 0))
    if (missing(col)) 
        col <- graphics::par("fg")
    if (!is.list(col)) 
        col <- list(col)
    else combinedpolygons <- FALSE
    col <- lapply(col, fixcolorvector)
    col <- rep(col, length.out = n)
	if (!is.null(border))
        border <- rep(border, length.out = n)

    #set the logarithm handling settings
    if (!add && log == "auto") {
        if (seemslog(groups)) {
            log <- "y"
            message("log=\"y\" selected")
        }
        else log <- ""
    }
    if (is.na(handlelog)) 
        if (add && ((horizontal & graphics::par()$xlog) || (!horizontal & 
            graphics::par()$ylog))) 
            handlelog <- TRUE
        else if (!add && (log != "")) 
            handlelog <- TRUE
        else handlelog <- FALSE
    if (handlelog) {
        mlog <- base::log
        mexp <- base::exp
    }
    else {
        mlog <- function(x) {
            x
        }
        mexp <- mlog
    }

    #generate the necessary data for the density shapes from the group data
    if (!is.numeric(bw)) {
        bw <- mean(sapply(groups, function(x) {
            ifelse(length(x) > 1, stats::density(mlog(x), kernel = kernel, 
                bw = bw)$bw, NA)
        }), na.rm = TRUE)
        if (is.nan(bw)) 
            bw <- 0.5
    }
    dens <- sapply(groups, mdensityxy)
    for (i in 1:n) dens[["y", i]] <- dens[["y", i]] * min(1, 
        length(groups[[i]])/grownage)
    if (is.na(wd)) 
        wd <- maxwidth/max(unlist(dens["y", ]))
    wd2 <- wd * boxwex/2

    #plot windows and axes
    axespars <- lapply(mcall[base::names(mcall) %in% c("xaxt", 
        "yaxt", "las", "cex.axis", "col.axis", "format", "tick", 
        "xaxp", "yaxp")], eval, parent.frame())
    if (!add) {
        if (!is.numeric(xlim)) {
            if (side == 2) 
                xlim <- c(0, displayn)
            else if (side == 3) 
                xlim <- c(1, displayn + 1)
            else xlim <- c(0.5, displayn + 0.5)
        }
        if (!is.numeric(ylim)) 
            ylim <- range(groups, mexp(unlist(dens["x", ])))
        graphics::plot.new()
        windowpars <- lapply(mcall[base::names(mcall) %in% c("yaxs", 
            "xaxs")], eval)
        if (horizontal) {
            names(windowpars)[names(windowpars) %in% c("xaxs", 
                "yaxs")] <- rev(names(windowpars)[names(windowpars) %in% 
                c("xaxs", "yaxs")])
            if (log == "y") 
                log <- "x"
            do.call("plot.window", c(list(xlim = ylim, ylim = xlim, 
                log = log), windowpars))
        }
        else {
            do.call("plot.window", c(list(xlim = xlim, ylim = ylim, 
                log = log), windowpars))
        }
        if (frame.plot) 
            graphics::box()
        if (axes) 
            do.call("axis", c(list(side = 2 - horizontal), axespars))
    }
    if (axes) {
        if (is.na(show.names)) 
            show.names <- (n > 1)
        if (show.names) 
            do.call("axis", c(list(1 + horizontal, at = at, labels = names), 
                axespars))
    }

    #plot the window contents
    if (what[1]) {
        if (overallline == 2) 
            val <- mexp(stats::median(mlog(unlist(groups))))
        else val <- mexp(mean(mlog(unlist(groups))))
# just plot a horizontal line at 0 / jmd 20190401
#        if (horizontal) 
#            abline(v = val, lty = 3)
#        else abline(h = val, lty = 3)
        graphics::abline(h = 0, lty = 3, lwd = 1.5)
    } else {
		val = "value not calculated, overall line was omitted"
	}
    if (what[2]) {
        beanplotpolyshapes(side, dens, at, wd2, combinedpolygons, 
            displayn, n, col, border, horizontal, mlog, mexp)
    }
    if (what[3]) {
        stats = beanplotbeanlines(groups, side, beanlines, beanlinewd, 
            at, boxwex, n, col, horizontal, mlog, mexp)
    } else {
		stats = "not calculated, beanlines were omitted"
	}
    if (what[4]) {
        beanplotscatters(groups, side, method, jitter, dens, 
            at, wd2, boxwex, n, ll, maxstripline, col, horizontal, 
            mlog, mexp)
    }
    if (any(!is.na(innerborder))) {
        beanplotinnerborders(innerborder, at, dens, side, displayn, 
            n, horizontal, mexp)
    }

    #finally, prints labels
    titlepars <- lapply(mcall[base::names(mcall) %in% c("main", 
        "sub", "xlab", "ylab", "cex.main", "col.main", "cex.lab", 
        "col.lab", "cex.sub", "col.sub", "font.main")], eval, parent.frame())
    do.call("title", titlepars)

    #return generated data that can be used for subsequent calls
    invisible(list(bw = bw, wd = wd, names = names, stats = stats, overall = val, log = log))
}

# 20191121 https://stackoverflow.com/questions/24331690/modify-package-function
environment(beanplot_mod) <- asNamespace("beanplot")

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
    if(basis=="rQEC") ylab <- expression(italic(n)[H[2] * O])
    if(basis=="biosynth") ylab <- expression(italic(n)[H[2] * O]~"(AA biosynthetic reactions)")
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

