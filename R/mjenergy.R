# JMDplots/mjenergy.R
# Make plots for mjenergy paper
# 20201207 jmd first version
# 20210205 moved to JMDplots

# Affinities for methanogenesis and amino acid synthesis 20201207
# Modified from Rainbow example from CHNOSZ/anintro.Rmd
mjenergy1 <- function(pdf = FALSE) {

  # Setup figure
  if(pdf) pdf("mjenergy1.pdf", width = 8, height = 5)
  layout(matrix(1:2, nrow = 1), widths = c(1, 1.5))
  par(mar = c(3, 3.5, 1, 1))

  # Read Shock and Canovas (2010) modelling results
  file <- system.file("extdata/cpetc/SC10_Rainbow.csv", package = "CHNOSZ")
  rb <- read.csv(file, check.names = FALSE)
  # Round highest value so we get a 350 °C label
  rb$T[1] <- round(rb$T[1])
  # Set basis species
  basis(c("CO2", "H2", "NH4+", "H2O", "H2S", "H+"))

  # Calculation for methanogenesis
  species("CH4", -3)
  a <- affinity(T = rb$T, CO2 = rb$CO2, H2 = rb$H2, `NH4+` = rb$`NH4+`, H2S = rb$H2S, pH = rb$pH)
  # Convert dimensionless affinities to kJ/mol
  T <- convert(a$vals[[1]], "K")
  a$values <- lapply(a$values, convert, "G", T)
  a$values <- lapply(a$values, convert, "J")
  a$values <- lapply(a$values, `*`, -0.001)
  # Make plot
  E.units("J")  # only affects axis label, not conversion above
  diagram(a, balance = 1, ylim = c(-200, 400), ylab = axis.label("A", prefix = "k"), lwd = 2, names = NA)
  abline(h = 0, lty = 3, lwd = 2, col = "gray40")
  label.figure("a", cex = 1.5)

  # Calculation for amino acids
  # Stay below 250 °C
  rb <- rb[rb$T < 250, ]
  species(c("serine", "asparagine", "glycine", "aspartic acid", "valine", "leucine", "tryptophan", "phenylalanine"), -3)
  a <- affinity(T = rb$T, CO2 = rb$CO2, H2 = rb$H2, `NH4+` = rb$`NH4+`, H2S = rb$H2S, pH = rb$pH)
  # Convert dimensionless affinities to kJ/mol
  T <- convert(a$vals[[1]], "K")
  a$values <- lapply(a$values, convert, "G", T)
  a$values <- lapply(a$values, convert, "J")
  a$values <- lapply(a$values, `*`, -0.001)
  # Make plot
  E.units("J")  # only affects axis label, not conversion above
  col <- c(1, 1, 2, 2, 3, 3, 4, 4)
  lty <- c(1, 2, 1, 2, 1, 2, 1, 2)
  diagram(a, balance = 1, ylim = c(-200, 400), ylab = axis.label("A", prefix = "k"), lwd = 2, col = col, lty = lty, names = NA)
  abline(h = 0, lty = 3, lwd = 2, col = "gray40")
  text(c(104, 99, 151, 33, 33, 33, 33, 33), c(-52, -79, -45, 17, 201, 271, 320, 344), c("S", "N", "G", "D", "V", "L", "W", "F"), adj = 0.5)
  label.figure("b", cex = 1.5)

  # FIXME: Reset units so they don't interfere with protein ionization calculations 20201210
  reset()
  if(pdf) dev.off()
}

# ZC of amino acids vs frequency in Mj proteome 20201207
mjenergy2 <- function(pdf = FALSE) {

  # Setup figure
  if(pdf) pdf("mjenergy2.pdf", width = 8, height = 4)
  par(mfrow = c(1, 2))
  par(mar = c(3.5, 3.5, 1, 1), mgp = c(2.45, 1, 0), las = 1)

  # Read amino acid compositions of Mj proteins
  #faafile <- c("UP000000805_243232.fasta.gz")
  #aa <- read.fasta(faafile)
  aa <- read.csv(system.file("extdata/organisms/UP000000805_243232.csv.gz", package = "JMDplots"), as.is = TRUE)
  # Calculate frequency (percentage) of each amino acid
  AAsum <- colSums(aa[, 6:25])
  AAperc <- 100 * AAsum / sum(AAsum)
  # Calculate ZC of amino acids
  ZC <- ZC(info(aminoacids("")))
  # Swap N and D so D is on top (plotted later)
  AAlab <- aminoacids()
  AAlab <- AAlab[c(1,2,12,4,5,6,7,8,9,10,11,3,13,14,15,16,17,18,19,20)]
  AAperc <- AAperc[c(1,2,12,4,5,6,7,8,9,10,11,3,13,14,15,16,17,18,19,20)]
  ZC <- ZC[c(1,2,12,4,5,6,7,8,9,10,11,3,13,14,15,16,17,18,19,20)]
  # Make plot
  plot(ZC, AAperc, xlab = axis.label("ZC"), ylab = quote("Percentage in "*italic("M. jannaschii")*" proteome"), cex = 2.5, pch = 21, bg = "white")
  # Add amino acid labels, nudge C and N
  ZC[2] <- ZC[2] + 0.175
  ZC[3] <- ZC[3] - 0.175
  text(ZC, AAperc, AAlab)
  # Add lines for C and N
  lines(c(ZC[2] - 0.05, ZC[2] - 0.1), c(AAperc[2], AAperc[2]))
  lines(c(ZC[3] + 0.05, ZC[3] + 0.1), c(AAperc[3], AAperc[3]))
  label.figure("a", cex = 1.2, xfrac = 0.02, yfrac = 0.97)

  # Calculate ZC of Mj proteins
  ZC <- ZC(protein.formula(aa))
  # Make barplot
  b <- seq(-0.55, 0, by = 0.01)
  hist(ZC, b, xlab = axis.label("ZC"), ylab = "Number of proteins", main = NULL)
  # Calculate the 1st and 3rd quartiles
  quartiles <- quantile(ZC, c(1,3)/4)
  abline(v = quartiles, lty = 2, lwd = 2, col = "gray40")
  text(-0.1, 140, "Upper\nquartile", col = 4)
  text(-0.4, 140, "Lower\nquartile", col = 2)
  label.figure("b", cex = 1.2, xfrac = 0.02, yfrac = 0.97)

  if(pdf) dev.off()
}

# Affinities of overall synthesis of proteins in Mj proteome 20201208
mjenergy3 <- function(pdf = FALSE) {

  # Setup figure
  if(pdf) pdf("mjenergy3.pdf", width = 8, height = 4)
  par(mfrow = c(1, 2))
  par(mar = c(3.5, 3.5, 1, 1), mgp = c(2.5, 1, 0), las = 1)

  # Read amino acid compositions of Mj proteins
  #faafile <- c("UP000000805_243232.fasta.gz")
  #aa <- read.fasta(faafile)
  aa <- read.csv(system.file("extdata/organisms/UP000000805_243232.csv.gz", package = "JMDplots"), as.is = TRUE)
  # Calculate ZC of Mj proteins
  ZC <- ZC(protein.formula(aa))
  # Calculate the 1st and 3rd quartiles
  quartiles <- quantile(ZC, c(1,3)/4)

  # Read Shock and Canovas (2010) modelling results
  file <- system.file("extdata/cpetc/SC10_Rainbow.csv", package = "CHNOSZ")
  rb <- read.csv(file, check.names = FALSE)
  # Stay below 250 °C
  rb <- rb[rb$T < 250, ]
  # Set basis species
  basis(c("CO2", "H2", "NH4+", "H2O", "H2S", "H+"))

  # Calculate affinities of protein synthesis
  ip <- add.protein(aa)
  a <- affinity(T = rb$T, CO2 = rb$CO2, H2 = rb$H2, `NH4+` = rb$`NH4+`, H2S = rb$H2S, pH = rb$pH, iprotein = ip, loga.protein = -3)
  # Convert dimensionless affinities to MJ/mol
  T <- convert(a$vals[[1]], "K")
  a$values <- lapply(a$values, convert, "G", T)
  a$values <- lapply(a$values, convert, "J")
  a$values <- lapply(a$values, `*`, -1e-6)

  # Make plots
  E.units("J")  # only affects axis label, not conversion above
  # Use red for lower quartile, black for interquartile range, and blue for highest quartile
  col <- rep(1, length(ZC))
  col[ZC < quartiles[1]] <- 2
  col[ZC > quartiles[2]] <- 4
  ylab <- axis.label("A", prefix = "M")
  ylab[[3]][[3]][[2]] <- "(mol protein)"
  diagram(a, balance = 1, ylim = c(-50, 150), ylab = ylab, lty = 1, lwd = 0.2, names = NA, col = col)
  abline(h = 0, lty = 3, lwd = 2, col = "gray40")
  label.figure("a", cex = 1.5)
  # Second plot: divide by protein length
  pl <- protein.length(ip)
  # Use kJ/mol
  a$values <- lapply(a$values, `*`, 1e3)
  ylab <- axis.label("A", prefix = "k")
  ylab[[3]][[3]][[2]] <- "(mol residue)"
  diagram(a, balance = pl, ylim = c(-50, 150), ylab = ylab, lty = 1, lwd = 0.2, names = NA, col = col)
  abline(h = 0, lty = 3, lwd = 2, col = "gray40")
  legend("topright", c("Lower quartile", "Middle 50%", "Upper quartile"), title = as.expression(quote(italic(Z)[C]~values)), lty = 1, col = c(2, 1, 4), bty = "n")
  label.figure("b", cex = 1.5)

  # FIXME: Reset units so they don't interfere with protein ionization calculations 20201210
  reset()
  if(pdf) dev.off()
}

# Read Shock and Canovas (2010) modelling results
# and interpolate on 1 °C increments 20201210
get_rainbow <- function(plot.it = FALSE) {
  file <- system.file("extdata/cpetc/SC10_Rainbow.csv", package = "CHNOSZ")
  dat <- read.csv(file, check.names = FALSE)
  T <- 5:350
  rb <- data.frame(
    T = T,
    CO2 = splinefun(dat$T, dat$CO2, method = "monoH.FC")(T),
    H2 = splinefun(dat$T, dat$H2, method = "monoH.FC")(T),
    pH = splinefun(dat$T, dat$pH, method = "monoH.FC")(T),
    `NH4+` = splinefun(dat$T, dat$`NH4+`, method = "monoH.FC")(T),
    H2S = splinefun(dat$T, dat$H2S, method = "monoH.FC")(T),
    CH4 = splinefun(dat$T, dat$CH4, method = "monoH.FC")(T),
    check.names = FALSE
  )
  if(plot.it) {
    par(mfrow = c(2, 3))
    plot(dat$T, dat$CO2); lines(rb$T, rb$CO2)
    plot(dat$T, dat$H2); lines(rb$T, rb$H2, col = 2)
    plot(dat$T, dat$pH); lines(rb$T, rb$pH)
    plot(dat$T, dat$`NH4+`); lines(rb$T, rb$`NH4+`)
    plot(dat$T, dat$H2S); lines(rb$T, rb$H2S)
    plot(dat$T, dat$CH4); lines(rb$T, rb$CH4)
  }
  rb
}

# Calculate affinity for amino acid synthesis and polymerization 20201213
# Add 'T' argument to set temperature and 'protein' argument to choose protein 20210125
# NOTE: CSG_METJA does not include signal peptide (1-28)
calc_affinity <- function(T = 85, protein = "CSG_METJA") {

  # Read Shock and Canovas (2010) modelling results
  rb <- get_rainbow()
  # Use values for selected temperature
  rb <- rb[rb$T==T, ]
  # Set basis species
  basis(c("CO2", "H2", "NH4+", "H2O", "H2S", "H+"), c(rb$CO2, rb$H2, rb$`NH4+`, 0, rb$H2S, -rb$pH))

  # Get amino acid composition of protein
  aa <- pinfo(pinfo(protein))
  AA <- as.numeric(aa[, 6:25])
  # Calculate affinity of synthesis of all amino acids in the protein
  E.units("J")
  sAA <- suppressMessages(subcrt(aminoacids(""), AA, logact = rep(-3, 20), T = T)$out)
  AAMJ <- formatC(sAA$A / 1e6, format = "f", digits = 1)
  print(paste("Affinity for synthesis of", sum(AA), "amino acids:", AAMJ, "MJ/mol"))

  # Calculate affinity of polymerization
  spoly <- suppressMessages(subcrt(c(aminoacids(""), protein), c(-AA, 1), logact = rep(-3, 21), T = T)$out)
  polyMJ <- formatC(spoly$A / 1e6, format = "f", digits = 1)
  print(paste0("Affinity for polymerization to form ", protein, ": ", polyMJ, " MJ/mol"))

  # Calculate affinity of protein synthesis (non-ionized)
  sprot <- suppressMessages(subcrt(protein, 1, logact = -3, T = T)$out)
  protMJ <- formatC(sprot$A / 1e6, format = "f", digits = 1)
  print(paste("Affinity for synthesis of", protein, "(non-ionized):", protMJ, "MJ/mol"))

  # Calculate affinity of protein synthesis (ionized)
  E.units("cal")
  species(protein)
  a <- suppressMessages(affinity(T = T))
  protJ.ion <- convert(-convert(a$values[[1]], "G", T = convert(T, "K")), "J")
  protMJ.ion <- formatC(protJ.ion / 1e6, format = "f", digits = 1)
  print(paste("Affinity for synthesis of", protein, "(ionized):", protMJ.ion, "MJ/mol"))

}
