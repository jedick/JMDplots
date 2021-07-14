# JMDplots/demo/evdevH2O.R
# Graphical abstract for research profile on website 20210203
# Moved to JMDplots 20210715
require("JMDplots")
# png("evdevH2O.png", width = 800, height = 600)

mat <- matrix(c(1,1,0, 1,1,0, 1,1,0, 1,2,2, 1,2,2, 0,2,2, 0,2,2), ncol = 3, byrow = TRUE)
layout(mat)
par(mar = c(3, 3.1, 3, 1), mgp = c(2.2, 0.5, 0), bg = "transparent")
par(cex = 1.5)
logH2Olab <- quote(phantom(.)~~~~~~~~~~~bold(log)*bolditalic(a)[bold(H[2]*O)])
logO2lab <- quote(bold(log)*bolditalic(f)[bold(O[2])]~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~phantom(.))

# Get mean amino acid compositions
gpa <- JMDplots:::getphyloaa("TPPG17")
PS <- gpa$aa$protein
# Load PS model proteins 
ipPS <- add.protein(gpa$aa)
# Set up system
basis("QEC")
O2 <- c(-72, -67)
H2O <- c(-3, 2)

# Plot B: 16 PS model proteins and n = 200 sample of human proteins
set.seed(3)
iind <- sample(1:nrow(gpa$pcomp$aa), 200)
ipind <- add.protein(gpa$pcomp$aa[iind, ])
a <- affinity(O2 = c(O2, 128), H2O = c(H2O, 128), iprotein = c(ipPS, ipind))
e <- equilibrate(a, as.residue = TRUE, loga.balance = 0)
d <- diagram(e, plot.it = FALSE)
par.orig <- JMDplots:::my.filled.contour(e$vals$O2, e$vals$H2O, d$predominant.values, xlab = logO2lab, ylab = NA,
  nlevels = 50, col = hcl.colors(75, "YlGnBu")[20:75], plot.key = FALSE,
  # use plot.axes to label the contour plot (see ?filled.contour)
  plot.axes = {
    mtext(logH2Olab, side = 2, line = 2, las = 0, cex = par("cex"))
    names <- sapply(strsplit(d$species$name, "\\|"), "[", 2)
    #diagram(e, add = TRUE, names = names, format.names = FALSE)
    # Use a higher resolution for making the lines
    a <- affinity(O2 = c(O2, 400), H2O = c(H2O, 400), iprotein = c(ipPS, ipind))
    diagram(a, as.residue = TRUE, add = TRUE, names = FALSE, format.names = FALSE)
    opar <- par(tcl = 0.3)
    thermo.axis()
    axis(1)
    axis(2)
    par(opar)
    # Show location of maximum activity for each PS model protein
    par(xpd = TRUE)
    for(i in 1:16) {
      imax <- arrayInd(which.max(e$loga.equil[[i]]), dim(e$loga.equil[[i]]))
      optO2 <- e$vals$O2[imax[1]]
      optH2O <- e$vals$H2O[imax[2]]
      points(optO2, optH2O, pch = 21, bg = 7, cex = 1.5)
    }
    par(xpd = FALSE)
    title("Maximum activities of phylostrata\non a background of human proteins", font.main = 1)
  },
  key.axes = {
    opar <- par(tcl = 0)
    axis(4, at = par("usr")[3:4], labels = round(par("usr")[3:4], 2))
    title(quote(bold(log)*bolditalic(a)[bold(protein)]), cex.main = 1, line = 1)
    par(opar)
  },
  add2 = TRUE
)

## Plot 2: Effective Eh for Liebeskind gene ages
PS_source <- "LMM16"
# Read results
datadir <- system.file("extdata/evdevH2O", package = "JMDplots")
H2Ofile <- file.path(datadir, paste0(PS_source, "_H2O_Hsa.csv"))
O2file <- file.path(datadir, paste0(PS_source, "_O2_Hsa.csv"))
H2O <- read.csv(H2Ofile, as.is = TRUE, check.names = FALSE)
O2 <- read.csv(O2file, as.is = TRUE, check.names = FALSE)
# Get phylostrata
iPS <- 2:ncol(H2O)
PS <- as.numeric(colnames(H2O)[iPS])
# Get mean values of logaH2O and logfO2
meanH2O <- colMeans(H2O[, iPS])
meanO2 <- colMeans(O2[, iPS])

# Calculate Eh
# H2O = 0.5 O2 + 2 H+ + 2 e- 
# logK = 0.5 logfO2 - 2 pH - 2 pe - logaH2O
# pe = 0.25 logfO2 - pH - 0.5 logaH2O - 0.5 logK
logK <- subcrt(c("H2O", "oxygen", "H+", "e-"), c(-1, 0.5, 2, 2), T = 25)$out$logK
pH <- 7
pe <- 0.25 * meanO2 - pH - 0.5 * meanH2O - 0.5 * logK
Eh <- convert(pe, "Eh")
# Another way to calculate using CHNOSZ::convert
Eh2 <- convert(meanO2, "E0", pH = pH, logaH2O = meanH2O)
stopifnot(all.equal(Eh, Eh2))
mV <- Eh * 1000

# Make Eh plot
xlim <- range(PS)
ylim <- c(-325, -150)
plot.new()
par(mar = c(4.9, 4.4, 0.8, 0.8), mgp = c(2.5, 1, 0))
plot.window(xlim, ylim, xaxs = "i")
# Draw white rectangle over figure region
X <- grconvertX(c(0, 1), "nfc", "user")
Y <- grconvertY(c(0, 1), "nfc", "user")
rect(X[1], Y[1]+1, X[2]-0.03, Y[2], col = "white", xpd = NA, border = "black")
box()
mtext("Eh (mV)", side = 2, line = 3, las = 0, cex = par("cex"), font = 2)
# Make rotated labels (modified from https://www.r-bloggers.com/rotated-axis-labels-in-r-plots/)
labels <- c("Cellular\norganisms", "Euk_Archaea", "Euk+Bac", "Eukaryota", "Opisthokonta", "Eumetazoa", "Vertebrata", "Mammalia")
text(x = 1:8, y = par()$usr[3] - 1.5 * strheight("A"), labels = labels, srt = 45, adj = 1, xpd = TRUE)
axis(1, at = 1:8, labels = NA)
axis(2, at = seq(-300, -150, 50))
# PS 1 (Cellular_organisms), 4 (Eukaryota), 6 (Eumetazoa), 8 (Mammalia)
abline(v = c(1.02, 4, 6, 7.98), col = 5, lwd = 3)
xtext <- 5.5
# Eh = -150 mV (plasma GSH/GSSG) Jones and Sies, 2015
# Eh = -199 mV (erythrocyte GSH/GSSG) van 't Erve et al., 2013
# Eh = -241 mV (cytosolic NADH/NAD+) Jones and Sies, 2015
# Eh = -318 mV (mitochondrial NADH/NAD+) Jones and Sies, 2015
abline(h = c(-150, -199, -241, -318), lty = 2, lwd = 1.5, col = "slategray4")
lines(PS, mV, lwd = 2, col = 2) 
text(2.85, -150, "Plasma GSH/GSSG", adj = c(0.5, 1.3))
text(2.85, -199, "Erythrocyte GSH/GSSG", adj = c(0.5, -0.3))
text(5.5, -241, "Cytosolic NADH/NAD+", adj = c(0.5, 1.3))
text(5.5, -318, "Mitochondrial NADH/NAD+", adj = c(0.5, -0.3))

# Run these to close and compress the PNG file
if(FALSE) {
  dev.off()
  system("pngquant evdevH2O.png")
  system("mv evdevH2O-fs8.png evdevH2O.png")
}
