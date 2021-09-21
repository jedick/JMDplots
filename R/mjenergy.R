# JMDplots/mjenergy.R
# Make plots for mjenergy paper
# 20201207 jmd first version
# 20210205 moved to JMDplots

# Affinities for methanogenesis and amino acid synthesis 20201207
# Modified from Rainbow example from CHNOSZ/anintro.Rmd
mjenergy1 <- function(pdf = FALSE) {

  ## Setup figure
  if(pdf) pdf("mjenergy1.pdf", width = 9, height = 9)
  mat <- matrix(c(1,1,1,1, 2,2,2,2,2,2, 0,0,3,3,3,3,3,3,0,0), byrow = TRUE, nrow = 2)
  layout(mat, heights = c(1.5, 1))
  par(mar = c(3, 3.5, 1, 1))
  par(cex = 1.2)

  ## Plot (a): Methanogenesis

  # Calculate affinities for Rainbow and Endeavour vents 20210815
  for(vent in c("Rainbow", "Endeavour")) {

    # Read Shock and Canovas (2010) modelling results
    file <- system.file(paste0("extdata/mjenergy/SC10_", vent, ".csv"), package = "JMDplots")
    SC10 <- read.csv(file, check.names = FALSE)
    # Set basis species
    basis(c("CO2", "H2", "NH4+", "H2O", "H2S", "H+"))

    # Calculate affinity for methanogenesis
    species("CH4", -3)
    a <- affinity(T = SC10$T, CO2 = SC10$CO2, H2 = SC10$H2, `NH4+` = SC10$`NH4+`, H2S = SC10$H2S, pH = SC10$pH)
    # Convert dimensionless affinities to kJ/mol
    TK <- convert(SC10$T, "K")
    a$values <- lapply(a$values, convert, "G", TK)
    a$values <- lapply(a$values, convert, "J")
    a$values <- lapply(a$values, `*`, -0.001)
    if(vent == "Rainbow") {
      # Make plot
      E.units("J")  # only affects axis label, not conversion above
      diagram(a, balance = 1, xlim = c(0, 350), ylim = c(-200, 400), ylab = axis.label("A", prefix = "k"), lwd = 2, names = NA)
      abline(h = 0, lty = 3, lwd = 2, col = "gray40")
    }
    if(vent == "Endeavour") {
      # Add line
      diagram(a, balance = 1, lwd = 2, col = 2, names = NA, add = TRUE)
    }

    # Add dashed line for variable activity of CH4  20210814
    a_vals <- numeric()
    for(i in 1:nrow(SC10)) {
      basis("CO2", SC10$CO2[i])
      basis("H2", SC10$H2[i])
      species("CH4", SC10$CH4[i])
      a <- suppressMessages(affinity(T = SC10$T[i]))
      a_vals <- c(a_vals, a$values[[1]])
    }
    # Convert dimensionless affinities to kJ/mol
    a_vals <- convert(a_vals, "G", T = TK)
    a_vals <- convert(a_vals, "J")
    a_vals <- - a_vals / 1000
    col <- 1
    if(vent == "Endeavour") col <- 2
    lines(SC10$T, a_vals, lty = 2, col = col)

  }
  # Add legend
  ltxt <- as.expression(c(quote(italic(a)[CH[4]] == 10^-3), quote(italic(a)[CH[4]]~from~mixing~model), quote("(Shock and Canovas, 2010)")))
  legend("topright", legend = ltxt, bty = "n", lty = c(1, 2, NA), lwd = c(1.7, 1.2), cex = 0.8)
  # Add labels
  text(175, 250, "Methanogenesis", font = 2)
  text(250, 110, "Rainbow", font = 2)
  text(250, -30, "Endeavour", font = 2)
  label.figure("a", cex = 1.5)

  # Plots (b) and (c): Amino acids
  for(vent in c("Rainbow", "Endeavour")) {
    # Read Shock and Canovas (2010) modelling results
    file <- system.file(paste0("extdata/mjenergy/SC10_", vent, ".csv"), package = "JMDplots")
    SC10 <- read.csv(file, check.names = FALSE)
    # Stay below 250 °C
    SC10 <- SC10[SC10$T < 250, ]
    species(aminoacids(""), -3)
    a <- affinity(T = SC10$T, CO2 = SC10$CO2, H2 = SC10$H2, `NH4+` = SC10$`NH4+`, H2S = SC10$H2S, pH = SC10$pH)
    # Convert dimensionless affinities to kJ/mol
    T <- convert(a$vals[[1]], "K")
    a$values <- lapply(a$values, convert, "G", T)
    a$values <- lapply(a$values, convert, "J")
    a$values <- lapply(a$values, `*`, -0.001)
    # Make plot
    E.units("J")  # only affects axis label, not conversion above
    # Use colors for ZC 20210813
    zc <- ZC(info(aminoacids(""), "aq"))
    col <- rep(1, length(zc))
    col[zc < -0.3] <- 2
    col[zc > 0.3] <- 4
    # Make plot
    if(vent == "Rainbow") ylim <- c(-200, 400)
    if(vent == "Endeavour") ylim <- c(-400, 100)
    diagram(a, balance = 1, xlim = c(0, 200), ylim = ylim, ylab = axis.label("A", prefix = "k"), lty = 1, lwd = 2, col = col, names = NA)
    abline(h = 0, lty = 3, lwd = 2, col = "gray40")
    if(vent == "Rainbow") {
      # Add labels for all amino acids 20210813
      lT <- c(A = 7, C = 5, D = 7, E = 7, F = 7, G = 7, H = 7, I = 7, K = 7, L = 7,
              M = 7, N = 4, P = 7, Q = 5, R = 7, S = 4, T = 7, V = 7, W = 7, Y = 7)
      lx <- SC10$T[lT]
      dy <- c(A = 8.5, C = 8.5, D = 8.5, E = 8.5, F = 8.5, G = 8.5, H = -8.5, I = -8.5, K = 8.5, L = 8.5,
              M = 8.5, N = -9.5, P = 8.5, Q = 8.5, R = 8.5, S = 9.5, T = 8.5, V = 8.5, W = 8.5, Y = 8.5)
      ly <- mapply("[", a$values, lT) + dy
      text(lx, ly, aminoacids(), cex = 0.7)
      # Add legend
      ltxt <- as.expression(c(quote(italic(Z)[C] < -0.3), quote("-0.3 <"~italic(Z)[C] < 0.3), quote(italic(Z)[C] > 0.3)))
      legend("topright", legend = ltxt, lty = 1, col = c(2, 1, 4), lwd = 2, bty = "n", cex = 0.8)
      text(150, 250, "Amino acids\nRainbow", font = 2)
      label.figure("b", cex = 1.5)
    }
    if(vent == "Endeavour") {
      text(30, -100, "Amino acids\nEndeavour", font = 2)
      label.figure("c", cex = 1.5)
    }
  }

  # Reset units so they don't interfere with protein ionization calculations 20201210
  # (bug fixed in CHNOSZ 1.4.0, but we leave this here for compatibility with earlier versions)
  reset()
  if(pdf) dev.off()
}

# ZC of amino acids vs frequency in the Mj proteome 20201207
mjenergy2 <- function(pdf = FALSE) {

  # Setup figure
  if(pdf) pdf("mjenergy2.pdf", width = 8, height = 4)
  par(mfrow = c(1, 2))
  par(mar = c(3.5, 3.5, 1, 1), mgp = c(2.45, 1, 0), las = 1)

  # Read amino acid compositions of Mj proteins
  aa <- read.csv(system.file("extdata/organisms/UP000000805_243232.csv.xz", package = "JMDplots"), as.is = TRUE)
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

# Affinities of overall synthesis of proteins in the Mj proteome 20201208
mjenergy3 <- function(pdf = FALSE, write.csv = FALSE) {

  # Setup figure
  if(pdf) pdf("mjenergy3.pdf", width = 8, height = 7)
  par(mfrow = c(2, 2))
  par(mar = c(3.5, 3.5, 1, 1), mgp = c(2.5, 1, 0), las = 1)
  par(cex = 1.1)

  # Read amino acid compositions of Mj proteins
  aa <- read.csv(system.file("extdata/organisms/UP000000805_243232.csv.xz", package = "JMDplots"), as.is = TRUE)
  # Calculate ZC of Mj proteins
  ZC <- ZC(protein.formula(aa))
  # Calculate the 1st and 3rd quartiles
  quartiles <- quantile(ZC, c(1,3)/4)

  # Save affinities to write CSV files  20210816
  out <- list(Rainbow = NULL, Endeavour = NULL)
  for(vent in c("Rainbow", "Endeavour")) {

    # Read Shock and Canovas (2010) modelling results
    file <- system.file(paste0("extdata/mjenergy/SC10_", vent, ".csv"), package = "JMDplots")
    SC10 <- read.csv(file, check.names = FALSE)
    # Stay below 250 °C
    SC10 <- SC10[SC10$T < 250, ]
    # Set basis species
    basis(c("CO2", "H2", "NH4+", "H2O", "H2S", "H+"))

    # Calculate affinities of protein synthesis
    ip <- add.protein(aa)
    a <- affinity(T = SC10$T, CO2 = SC10$CO2, H2 = SC10$H2, `NH4+` = SC10$`NH4+`, H2S = SC10$H2S, pH = SC10$pH, iprotein = ip, loga.protein = -3)
    # Convert dimensionless affinities to MJ/mol
    T <- convert(a$vals[[1]], "K")
    a$values <- lapply(a$values, convert, "G", T)
    a$values <- lapply(a$values, convert, "J")
    a$values <- lapply(a$values, `*`, -1e-6)
    # Save values with temperature as column names
    outvals <- do.call(rbind, a$values)
    colnames(outvals) <- SC10$T
    out[[vent]] <- outvals

    # First plot: whole proteins
    if(vent == "Rainbow") ylim <- c(-50, 150)
    if(vent == "Endeavour") ylim <- c(-200, 0)
    # Use red for lower quartile, black for interquartile range, and blue for highest quartile
    col <- rep(1, length(ZC))
    col[ZC < quartiles[1]] <- 2
    col[ZC > quartiles[2]] <- 4
    # Get axis label
    E.units("J")  # only affects axis label, not energy conversion above
    ylab <- axis.label("A", prefix = "M")
    ylab[[3]][[3]][[2]] <- "(mol protein)"
    diagram(a, balance = 1, ylim = ylim, ylab = ylab, lty = 1, lwd = 0.2, names = NA, col = col)
    abline(h = 0, lty = 3, lwd = 2, col = "gray40")
    if(vent == "Rainbow") {
      legend <- c("Lower quartile", "Middle 50%", "Upper quartile")
      legend("topright", legend, title = as.expression(quote(italic(Z)[C]~values)), lty = 1, col = c(2, 1, 4), bty = "n", cex = 0.8, inset = c(0.02, 0.05))
    }
    if(vent == "Rainbow") label.figure("a", cex = 1.5) else label.figure("c", cex = 1.5)

    # Second plot: divide by protein length
    pl <- protein.length(ip)
    # Use kJ/mol
    a$values <- lapply(a$values, `*`, 1e3)
    ylab <- axis.label("A", prefix = "k")
    ylab[[3]][[3]][[2]] <- "(mol residue)"
    diagram(a, balance = pl, xlim = c(0, 200), ylim = ylim, ylab = ylab, lty = 1, lwd = 0.2, names = NA, col = col)
    abline(h = 0, lty = 3, lwd = 2, col = "gray40")
    if(vent == "Rainbow") {
      text(125, 125, "Rainbow", font = 2)
      label.figure("b", cex = 1.5)
    }
    if(vent == "Endeavour") {
      text(125, -25, "Endeavour", font = 2)
      label.figure("d", cex = 1.5)
    }

  }

  # FIXME: Reset units so they don't interfere with protein ionization calculations 20201210
  reset()
  if(pdf) dev.off()

  # Return and save affinities 20210816
  out$Rainbow <- round(out$Rainbow, 4)
  out$Endeavour <- round(out$Endeavour, 4)
  rownames(out$Rainbow) <- aa$protein
  rownames(out$Endeavour) <- aa$protein
  if(write.csv) {
    write.csv(out$Rainbow, "Dataset_S2.csv", quote = FALSE)
    write.csv(out$Endeavour, "Dataset_S3.csv", quote = FALSE)
  }
  invisible(out)

}

# Read Shock and Canovas (2010) modelling results
# and interpolate on 1 °C increments 20201210
get_SC10 <- function(vent = "Rainbow", plot.it = FALSE) {
  file <- system.file(paste0("extdata/mjenergy/SC10_", vent, ".csv"), package = "JMDplots")
  dat <- read.csv(file, check.names = FALSE)
  T <- 5:350
  SC10 <- data.frame(
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
    plot(dat$T, dat$CO2); lines(SC10$T, SC10$CO2)
    plot(dat$T, dat$H2); lines(SC10$T, SC10$H2, col = 2)
    plot(dat$T, dat$pH); lines(SC10$T, SC10$pH)
    plot(dat$T, dat$`NH4+`); lines(SC10$T, SC10$`NH4+`)
    plot(dat$T, dat$H2S); lines(SC10$T, SC10$H2S)
    plot(dat$T, dat$CH4); lines(SC10$T, SC10$CH4)
  }
  SC10
}

# Calculate affinity for amino acid synthesis and polymerization 20201213
# Add 'T' argument to set temperature and 'protein' argument to choose protein 20210125
# NOTE: CSG_METJA does not include signal peptide (1-28)
calc_affinity <- function(T = 85, protein = "CSG_METJA") {

  # Read Shock and Canovas (2010) modelling results
  SC10 <- get_SC10()
  # Use values for selected temperature
  SC10 <- SC10[SC10$T==T, ]
  # Set basis species
  basis(c("CO2", "H2", "NH4+", "H2O", "H2S", "H+"), c(SC10$CO2, SC10$H2, SC10$`NH4+`, 0, SC10$H2S, -SC10$pH))

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

# Functions to make Appendix (SI) files 20210812

# Table S3: Overall synthesis reactions of amino acids
mjenergy_Table_S3 <- function() {
  # Set basis species
  basis(c("CO2", "NH4+", "H2S", "H2O", "H2", "H+"))
  # Get names of amino acids
  AA <- sort(aminoacids(""))
  # Start table
  out <- data.frame("Amino acid" = AA, Reaction = NA)
  # Loop over amino acids
  for(i in seq_along(AA)) {
    # Use subcrt() to balance formation reaction for each amino acid
    reaction <- suppressMessages(subcrt(AA[i], 1, T = 25))$reaction
    # Where to store text for reactants and products
    rtxt <- ptxt <- character()
    # Loop over species in reaction
    for(j in 1:nrow(reaction)) {
      # Get absolute value of coefficient
      cabs <- abs(reaction$coeff[j])
      # Paste coefficient (if it's != 1) and chemical formula
      ccf <- reaction$formula[j]
      if(cabs != 1) ccf <- paste(cabs, ccf)
      ## Append name of amino acid
      #if(j == 1) ccf <- paste0(ccf, " (", reaction$name[j], ")")
      # Assign text to reactants or products
      if(reaction$coeff[j] > 0) ptxt <- c(ptxt, ccf)
      if(reaction$coeff[j] < 0) rtxt <- c(rtxt, ccf)
    }
    # Collapse species text for reactants and products
    rtxt <- paste(rtxt, collapse = " + ")
    ptxt <- paste(ptxt, collapse = " + ")
    # Write reaction
    reaction_txt <- paste(rtxt, ptxt, sep = " = ")
    # Put reaction into table
    out$Reaction[i] <- reaction_txt
  }
  out
}

# Data Set S1: Chemical formulas of proteins 20210813
mjenergy_Dataset_S1 <- function(write.csv = FALSE) {
  # Read amino acid compositions of Mj proteins
  aa <- read.csv(system.file("extdata/organisms/UP000000805_243232.csv.xz", package = "JMDplots"), as.is = TRUE)
  # Calculate protein formulas and length
  pf <- protein.formula(aa)
  pl <- protein.length(aa)
  # Calculate protein ZC (round to 5 decimal places for table)
  zc <- round(ZC(pf), 5)
  # Make data frame
  out <- data.frame(Protein = rownames(pf), pf, length = pl, ZC = zc)
  # Command to save CSV file
  if(write.csv) write.csv(out, "Dataset_S1.csv", row.names = FALSE, quote = FALSE) else out
}

