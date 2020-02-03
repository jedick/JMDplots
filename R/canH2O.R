# JMDplots/canH2O.R
# make plots for the paper:
# Water as a reactant in the differential expression of proteins in cancer
# 20191126 jmd first version

# study overview plot 20191203
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
  if(pdf) dev.off()
}

