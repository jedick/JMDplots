# This is adapted from the supporting R code for
# "Growth and stability of stratiform carrollite (CuCo2S4) in the Tenke-Fungurume ore district, Central African Copperbelt"
# by Bjorn von der Heyden et al.
# R script written by Jeffrey Dick with input from Bjorn von der Heyden

# History:
# 20220524 to 20221214
#   Revisions 01 to 12
# 20230216 v13_RH95
#   Use data for cattierite and linnaeite from Robie & Hemingway (1995)
#   Put all functions into one file
# 20230217 v14_carrollite_V
#   Add V for carrollite
#   Use same number of S atoms for both reactions in Figure S5
# 20230815 v15_Cu0.92Co2.07S4
#   Use formula and properties for carrollite listed by Gibson et al. (2023)
# 20230822 v16_PMW87
#   Revert to cattierite and linnaeite from Pankratz et al. (1987)
# 20230831 v17_uDSC7
#   Use transition enthalpy and fit post-transition Cp from uDSC7
# 20240205 Moved to JMDplots with modifications:
#   Thermodynamic data for sulfides (including carrollite) are now in OBIGT instead of added via functions
#   Add "pdf" option to functions (defaults to FALSE for compiling carrollite.Rmd)

# Temperature (K) and Cp (J) from uDSC7 measurements for carrollite (Cu0.92Co2.07S4)
uDSC7 <- structure(list(Temp..K = c(300.03, 305.03, 310.02, 315.02, 320.01,
  325.01, 330, 335.09, 340.08, 345.08, 350.07, 355.07, 360.06,
  364.08, 365.06, 366.04, 367.02, 368, 369.07, 370.05, 371.03,
  372.01, 373.08, 374.07, 375.05, 376.03, 377, 378.08, 379.061,
  380.04, 380.5, 381.02, 382, 383.07, 384.06, 385.04, 386.02, 387.09,
  388.07, 389.05, 390.03, 391.01, 391.99), Cp..J.mol.K = c(158.76,
  160.29, 160.9, 162.52, 160.99, 163.37, 162.48, 164.68, 163.82,
  165.58, 165.33, 166.76, 166.88, 166.19, 169.19, 169.41, 168.03,
  166.32, 167.86, 168.32, 170.2, 171.99, 174.63, 179.51, 184.92,
  191.48, 191.84, 179.85, 166.57, 165.26, 172.87, 164.55, 165.33,
  164.71, 165.16, 165.26, 165.16, 164.51, 164.5, 163.31, 162.98,
  163.89, 163.25), Stage = c("Pre-transition", "Pre-transition",
  "Pre-transition", "Pre-transition", "Pre-transition", "Pre-transition",
  "Pre-transition", "Pre-transition", "Pre-transition", "Pre-transition",
  "Pre-transition", "Pre-transition", "Pre-transition", "Pre-transition",
  "Transition", "Transition", "Transition", "Transition", "Transition",
  "Transition", "Transition", "Transition", "Transition", "Transition",
  "Transition", "Transition", "Transition", "Transition", "Transition",
  "Transition", "Transition", "Post-transition", "Post-transition",
  "Post-transition", "Post-transition", "Post-transition", "Post-transition",
  "Post-transition", "Post-transition", "Post-transition", "Post-transition",
  "Post-transition", "Post-transition")), class = "data.frame", row.names = c(NA,
  -43L)
)

# Update OBIGT with fits to formation constants from
# Migdisov et al. (2011) (https://doi.org/10.1016/j.gca.2011.05.003)
# 20220721 jmd
add_Co_aqueous <- function() {

  # Temperature and formation constants from Table 3 of Migdisov et al. (2011)
  # (From spectroscopic data)
  T_table3 <- c(200, 250, 300)
  species.Cl <- list(
    Cl1 = c("Co+2", "Cl-", "CoCl+"),
    Cl2 = c("Co+2", "Cl-", "CoCl2"),
    Cl3 <- c("Co+2", "Cl-", "CoCl3-"),
    Cl4 <- c("Co+2", "Cl-", "CoCl4-2")
  )
  coeff.Cl <- list(
    Cl1 = c(-1, -1, 1),
    Cl2 <- c(-1, -2, 1),
    Cl3 <- c(-1, -3, 1),
    Cl4 <- c(-1, -4, 1)
  )
  logB_table3 <- list(
    Cl1 = c(1.37, 1.82, 2.42),
    Cl2 = c(2.71, 3.63, 4.94),
    Cl3 = c(3.72, 5.06, 6.44),
    Cl4 = c(5.12, 6.49, 8.33)
  )

  # Temperature and formation constants from Table 4 of Migdisov et al. (2011)
  # (From solubility data)
  T_table4 <- c(120, 150, 200, 250, 300)
  species.HS <- list(
    HS = c("Co+2", "HS-", "CoHS+"),
    H2S = c("Co+2", "H2S", "H+", "CoHS+")
  )
  coeff.HS <- list(
    HS = c(-1, -1, 1),
    H2S <- c(-1, -1, 1, 1)
  )
  logB_table4 <- list(
    HS = c(6.24, 6.02, 5.84, 5.97, 6.52),
    H2S = c(-0.23, -0.48, -0.85, -1.05, -1.04),
    # Use solubility-determined formation constants of Cl species for plots only (not in logB.to.OBIGT)
    Cl1 = c(NA, NA, NA, NA, NA),
    Cl2 = c(NA, 1.26, 2.42, 3.85, 5.36),
    Cl3 = c(NA, 2.89, 3.50, NA, NA),
    Cl4 = c(NA, NA, 4.86, 6.62, 8.39)
  )

  # Fit spectroscopic data for Cl complexes
  tolerance <- c(0.05, 0.05, 0.1, 0.05)
  for(i in 1:4) {
    # Don't try to fit NA values
    ina <- is.na(logB_table3[[i]])
    logB.to.OBIGT(logB_table3[[i]][!ina], species.Cl[[i]], coeff.Cl[[i]], T = T_table3[!ina], P = "Psat", npar = 2, tolerance = tolerance[i])
  }
  # Fit solubility data for HS complexes
  ## Use HS- reaction only
  for(i in 1:1) {
    ina <- is.na(logB_table4[[i]])
    logB.to.OBIGT(logB_table4[[i]][!ina], species.HS[[i]], coeff.HS[[i]], T = T_table4[!ina], P = "Psat", npar = 2, tolerance = 0.1)
  }

  # Return values for making Figure S4 20240206
  list(T_table3 = T_table3, species.Cl = species.Cl, coeff.Cl = coeff.Cl, logB_table3 = logB_table3,
       T_table4 = T_table4, species.HS = species.HS, coeff.HS = coeff.HS, logB_table4 = logB_table4)

}

# Add carrollite GHS and fitted Cp coefficients
# 20220724 First version
# 20230831 Calculate parameters for second polymorph
# 20240205 Change species name to carrollite_test to not clash with carrollite in OBIGT
calc_carrollite <- function() {

  MK_coeffs <- lapply(c("Pre-transition", "Post-transition"), function(Stage) {
    # Get uDSC7 data for pre-transition or post-transition
    dat <- uDSC7[uDSC7$Stage == Stage, ]
    TK <- dat$Temp..K
    Cp <- dat$Cp..J.mol.K
    # Calculate values of TK^-2
    TKn2 <- TK^-2
    data <- data.frame(Cp, TK, TKn2)
    # Regress Maier-Kelley parameters for pre-transition
    if(Stage == "Pre-transition") Cplm <- lm(Cp ~ TK + TKn2, data = data)
    # Regress linear equation for post-transition
    if(Stage == "Post-transition") Cplm <- lm(Cp ~ TK, data = data)
    Cplm$coefficients
  })

  ##
  ## The remaining steps are adapted from the CHNOSZ FAQ 20230831
  ##

  # The formula of the new mineral
  formula <- "Cu0.92Co2.07S4"
  # Use temperature in Kelvin for the calculations below
  old.T.units <- T.units()
  T.units("K")
  # Thermodynamic properties of polymorph 1 at 25 °C (298.15 K)
  # GHS and 25 °C Cp from Gibson et al. (2023) doi:10.1016/j.jct.2023.107096
  G1 <- -331140
  H1 <- -344460
  S1 <- 176.33
  Cp1 <- 158.48
  # Heat capacity coefficients for polymorph 1
  a1 <- MK_coeffs[[1]][1]
  b1 <- MK_coeffs[[1]][2]
  c1 <- MK_coeffs[[1]][3]
  # V from molar mass and calculated density (https://www.mindat.org/min-911.html)
  # 308.694 g/mol / 4.83 g/cm3 = 63.91 cm3/mol
  V1 <- V2 <- 63.91
  # Transition temperature (at maximum Cp of uDSC7 peak)
  Ttr <- 377
  # Transition enthalpy (J/mol)
  DHtr <- 162
  # Heat capacity coefficients for polymorph 2
  a2 <- MK_coeffs[[2]][1]
  b2 <- MK_coeffs[[2]][2]
  # Maximum temperature of polymorph 2
  T2 <- 650

  # Use the temperature (Ttr) and enthalpy of transition (DHtr) to calculate the entropy of transition (DStr)
  DGtr <- 0  # DON'T CHANGE THIS
  TDStr <- DHtr - DGtr  # TΔS° = ΔH° - ΔG°
  DStr <- TDStr / Ttr

  # Start new database entries that include basic information, volume, and heat capacity coefficients for each polymorph
  mod.OBIGT("carrollite_test", formula = formula, state = "cr",
    E_units = "J", G = 0, H = 0, S = 0, V = V1, Cp = Cp1,
    a = a1, b = b1, c = c1, d = 0, e = 0, f = 0, lambda = 0, T = Ttr)
  mod.OBIGT("carrollite_test", formula = formula, state = "cr2",
    E_units = "J", G = 0, H = 0, S = 0, V = V2,
    a = a2, b = b2, c = 0, d = 0, e = 0, f = 0, lambda = 0, T = T2)

  # Calculate changes of entropy from 298.15 K to Ttr, then the difference between polymorphs at 298.15 K
  DS1 <- subcrt("carrollite_test", "cr", P = 1, T = Ttr, use.polymorphs = FALSE)$out[[1]]$S
  DS2 <- subcrt("carrollite_test", "cr2", P = 1, T = Ttr)$out[[1]]$S
  DS298 <- DS1 + DStr - DS2

  # Put the values of S° at 298.15 into OBIGT, then calculate the changes of all thermodynamic properties of each polymorph between 298.15 K and Ttr
  mod.OBIGT("carrollite_test", state = "cr", S = S1)
  mod.OBIGT("carrollite_test", state = "cr2", S = S1 + DS298)
  D1 <- subcrt("carrollite_test", "cr", P = 1, T = Ttr, use.polymorphs = FALSE)$out[[1]]
  D2 <- subcrt("carrollite_test", "cr2", P = 1, T = Ttr)$out[[1]]

  # It’s a good idea to check that the entropy of transition is calculated correctly
  stopifnot(all.equal(D2$S - D1$S, DStr))

  # Calculate differences of ΔG° and ΔH° between the polymorphs at 298.15 K
  DG298 <- D1$G + DGtr - D2$G
  DH298 <- D1$H + DHtr - D2$H
  mod.OBIGT("carrollite_test", state = "cr", G = G1, H = H1)
  mod.OBIGT("carrollite_test", state = "cr2", G = G1 + DG298, H = H1 + DH298)

  # Check that the values of G, H, and S at 25 °C for a given polymorph are consistent with each other
  # (to within 33 J/mol)
  cr <- info(info("carrollite_test", "cr"))
  cr2 <- info(info("carrollite_test", "cr2"))
  stopifnot(abs(check.GHS(cr)) < 33)
  stopifnot(abs(check.GHS(cr2)) < 33)

  # Reset T units
  T.units(old.T.units)

  # Return the parameter values
  list(cr = cr, cr2 = cr2)

  # chk <- check_carrollite()
  # all.equal(chk$cr[, c(10:22)], info(info("carrollite", "cr"))[, c(10:22)], check.attributes = FALSE, tolerance = 1e-4)
  # all.equal(chk$cr2[, c(10:22)], info(info("carrollite", "cr2"))[, c(10:22)], check.attributes = FALSE, tolerance = 1e-4)

}

# logfO2-pH diagram for Cu and Co minerals with solubility contours
# 20220524 Extrapolate high-T Cp values for carrollite
# 20220715 Start solubility calculations (derived from CHNOSZ/demo/minsol.R)
# 20220721 Update aqueous Co complexes (fits to Migdisov et al.)
# 20220724 Use stack_mosaic() in CHNOSZ devel version
# 20220825 Suppress linnaeite
# 20230209 Don't suppress linnaeite
carrollite_5 <- function(res = 500, pdf = FALSE) {

  if(pdf) pdf("Figure_5.pdf", width = 9, height = 4)

  # System variables held constant
  P <- "Psat"
  pH <- c(2, 11, res)
  # This gets us close to total S = 0.003 m
  Stot <- 3e-3
  ## Mass fraction NaCl in saturated solution at 100 degC, from CRC handbook
  #wNaCl <- 0.2805  
  # 15 % NaCl (Vasyukova and Williams-Jones, 2022)
  wNaCl <- 0.15
  # Molality of NaCl
  mNaCl <- 1000 * wNaCl / (mass("NaCl") * (1 - wNaCl))

  add_solubility <- function(metal = "Co", T, pH, O2) {

    # Set up basis species
    basis(c(metal, "H2S", "Cl-", "oxygen", "H2O", "H+"))
    basis("H2S", log10(Stot))
    # Estimate ionic strength and molality of Cl-
    NaCl <- NaCl(m_NaCl = mNaCl, T = T, P = P)
    basis("Cl-", log10(NaCl$m_Clminus))

    # Add minerals and aqueous species
    icr <- retrieve(metal, c("H", "O", "S", "Cl"), state = "cr")
    iaq <- retrieve(metal, c("H", "O", "S", "Cl"), state = "aq")
    # Swap through these basis species to make mosaic diagram
    bases <- c("H2S", "HS-", "HSO4-", "SO4-2")

    # Calculate minimum solubility among all the minerals 20201008
    # (i.e. saturation condition for the solution)
    # Use solubility() 20210303
    species(icr)
    sout <- solubility(iaq, bases = bases, pH = pH, O2 = O2, T = T, P = P, IS = NaCl$IS, in.terms.of = metal)
    # Convert to ppm
    sout <- convert(sout, "ppm")
    # Specify contour levels
    levels <- c(1, 10, 100)
    # Color for solubility contours
    if(metal == "Cu") scol <- 2
    if(metal == "Co") scol <- "blue2"
    diagram(sout, levels = levels, contour.method = "flattest", add = TRUE, col = scol, lwd = 1, cex = 0.8)

  }

  make_stack <- function(T, pH, O2) {

    # System: Cu-Co-O-S
    # NOTE: this must include the first species listed in each of bases, species1, and species2 below
    basis(c("Cu", "Co", "Cl-", "H2S", "H2O", "oxygen", "H+"))
    basis("H2S", log10(Stot))
    # Estimate ionic strength and molality of Cl-
    NaCl <- NaCl(m_NaCl = mNaCl, T = T)
    basis("Cl-", log10(NaCl$m_Clminus))

    # Speciate aqueous sulfur
    bases <- c("H2S", "HS-", "HSO4-", "SO4-2")
    # Cu-bearing minerals
    names1 <- species1 <- c("copper", "cuprite", "tenorite", "chalcocite", "covellite")
    # Co-bearing and CoCu-bearing minerals
    names2 <- species2 <- c("cobalt", "cobalt monoxide", "guite", "cattierite", "linnaeite", "Co-pentlandite")
    names12 <- species12 <- c("carrollite")
    ## Uncomment these to use mineral formulas instead of names
    #names1 <- info(info(species1))$formula
    #names2 <- info(info(species2))$formula
    #names12 <- info(info(species12))$formula
    names2[names2 == "cobalt monoxide"] <- "CoO"
    # Define plotting parameters
    lty <- list(2, 0, 0)
    lwd <- list(1.5, 0, 0)
    col <- list(8, 4, 4)
    col.names <- list("#888888", c(1, 1, 1, 1, "white", 1), 6)
    fill <- list(NA, c("white", "#aaacac", "#e0e2e2", "#e8eca7", "#5e5f60", "#8f9091"), adjustcolor(6, alpha.f = 0.312))
    names <- list(names1, names2, names12)

    if(T < 200) {
      dx <- list(c(-4.8, -3.5, 0, -1.6, 0.6), c(0, 2.5, 3.7, 2.2, 0, 0.6), 2.2)
      dy <- list(c(3.5, 0, 0, 1.5, -0.8), c(0, -2, 0, -1.7, 0, 0), 0)
      srt <- list(0, c(0, 0, 0, 32, 38, 0), 38)
    } else {
      dx <- list(c(-3.8, -1.5, 0, -2.5, 0), c(0, -0.7, 0, -0.8, 2.53, 1.3), 2.15)
      dy <- list(c(3, -1, 0, 8, 0), c(0, 0, 0, 0.8, 1.7, 1), 0.14)
      srt <- list(0, c(0, 0, 0, 44, 44, 0), 44)
    }

    # Create mosaic stack (Cu in species1, Co in species2)
    sm <- stack_mosaic(bases, species1, species2, species12, names = names, col = col, col.names = col.names, fill = fill,
      dx = dx, dy = dy, srt = srt, lwd = lwd, lty = lty, pH = pH, O2 = O2, T = T, P = P, IS = NaCl$IS)
    # Because solid fill of Co fields covers the Cu lines, replot them 20220815
    diagram(sm[[1]], add = TRUE, lty = lty[[1]], lwd = lwd[[1]], col = col[[1]], col.names = col.names[[1]],
      fill = fill[[1]], dx = dx[[1]], dy = dy[[1]], srt = srt[[1]])
    # Add water stability line 20220825
    water.lines(sm[[1]], lty = 3)
    # Also replot tick marks and borders of plot
    thermo.axis()
    box()
    # Add title
    main <- bquote("Cu-Co-O-S, "*.(T)*"\u00b0C")
    title(main, font.main = 1)

  }

  layout(matrix(1:3, nrow = 1), widths = c(2, 2, 1))

  # Define temperature (degrees C) and logfO2 ranges
  T <- 150
  O2 <- c(-57, -32, res)
  make_stack(T = T, pH = pH, O2 = O2)
  add_solubility("Cu", T = T, pH = pH, O2 = O2)
  add_solubility("Co", T = T, pH = pH, O2 = O2)
  label.figure("a", font = 2, cex = 2)

  # Increase temperature
  T <- 250
  O2 <- c(-43, -23, res)
  make_stack(T = T, pH = pH, O2 = O2)
  add_solubility("Cu", T = T, pH = pH, O2 = O2)
  add_solubility("Co", T = T, pH = pH, O2 = O2)
  label.figure("b", font = 2, cex = 2)

  # Add legend
  par(mar = c(4, 0, 4, 0))
  plot.new()
  l <- c(
    lP(P),
    lNaCl(mNaCl),
    lS(Stot)
  )
  legend("topleft", legend = lex(l), bty = "n", cex = 1.2, xpd = NA)
  legend <- as.expression(c("Cu solubility", "Co solubility", quote(H[2]*O~"stability limit")))
  legend("bottomleft", legend = legend, col = c(2, "blue2", 1), lty = c(1, 1, 3), lwd = 1.2, bty = "n", cex = 1.2, xpd = NA)
  #legend("left", "Formation of linnaeite\nis NOT suppressed", text.font = 3, bty = "n")

  if(pdf) dev.off()

}

# Comparison of Cu-Co and Fe-Cu diagrams
# 20220823 Adapted from ../v9_summary/Figure_1.R
# 20220825 Renamed from Fe-Cu-Co.R
carrollite_8 <- function(res, pdf = FALSE) {

  if(pdf) pdf("Figure_8.pdf", width = 9, height = 4)

  # System variables held constant
  P <- "Psat"
  pH <- c(2, 11, res)
  # This gets us close to total S = 0.003 m
  Stot <- 3e-3
  ## Mass fraction NaCl in saturated solution at 100 degC, from CRC handbook
  #wNaCl <- 0.2805  
  # 15 % NaCl (Vasyukova and Williams-Jones, 2022)
  wNaCl <- 0.15
  # Molality of NaCl
  mNaCl <- 1000 * wNaCl / (mass("NaCl") * (1 - wNaCl))

  stack_CuCo <- function(T, pH, O2) {

    # System: Cu-Co-O-S
    # NOTE: this must include the first species listed in each of bases, species1, and species2 below
    basis(c("Cu", "Co", "Cl-", "H2S", "H2O", "oxygen", "H+"))
    basis("H2S", log10(Stot))
    # Estimate ionic strength and molality of Cl-
    NaCl <- NaCl(m_NaCl = mNaCl, T = T)
    basis("Cl-", log10(NaCl$m_Clminus))

    # Speciate aqueous sulfur
    bases <- c("H2S", "HS-", "HSO4-", "SO4-2")
    # Cu-bearing minerals
    names1 <- species1 <- c("copper", "cuprite", "tenorite", "chalcocite", "covellite")
    # Co-bearing and CoCu-bearing minerals
    names2 <- species2 <- c("cobalt", "cobalt monoxide", "guite", "cattierite", "Co-pentlandite")
    names12 <- species12 <- c("carrollite")
    ## Uncomment these to use mineral formulas instead of names
    #names1 <- info(info(species1))$formula
    #names2 <- info(info(species2))$formula
    #names12 <- info(info(species12))$formula
    names2[names2 == "cobalt monoxide"] <- "CoO"
    names <- list(names1, names2, names12)
    # Define plotting parameters
    lty <- list(2, 0, 0)
    lwd <- list(1.5, 0, 0)
    col <- list(8, 4, 4)
    col.names <- list("#888888", c(1, 1, 1, 1, "white", 1), 6)
    fill <- list(NA, c("white", "#aaacac", "#e0e2e2", "#e8eca7", "#8f9091"), adjustcolor(6, alpha.f = 0.312))

    dx <- list(c(0, -3.5, 0, -2.3, 0), c(0, 2.5, 3.7, 2, 0.3), 2.2)
    dy <- list(c(0, 0, 0, 0, 0), c(0, -2, 0, -1.5, -1), 0)
    srt <- list(0, c(0, 0, 0, 44, 44, 44), 44)

    # Create mosaic stack (Cu in species1, Co in species2)
    sm <- stack_mosaic(bases, species1, species2, species12, names = names, col = col, col.names = col.names, fill = fill,
      dx = dx, dy = dy, srt = srt, lwd = lwd, lty = lty, pH = pH, O2 = O2, T = T, P = P, IS = NaCl$IS)
    # Because solid fill of Co fields covers the Cu lines, replot them 20220815
    diagram(sm[[1]], add = TRUE, lty = lty[[1]], lwd = lwd[[1]], col = col[[1]], col.names = col.names[[1]],
      fill = fill[[1]], dx = dx[[1]], dy = dy[[1]], srt = srt[[1]])
    # Add water stability line 20220825
    water.lines(sm[[1]], lty = 3)
    # Also replot tick marks
    thermo.axis()
    # Add title
    main <- bquote("Cu-Co-O-S, "*.(T)*"\u00b0C")
    title(main, font.main = 1)

  }

  stack_FeCu <- function(T, pH, O2) {

    # System: Fe-Cu-O-S
    # NOTE: this must include the first species listed in each of bases, species1, and species2 below
    basis(c("Fe", "Cu", "Cl-", "H2S", "H2O", "oxygen", "H+"))
    basis("H2S", log10(Stot))
    # Estimate ionic strength and molality of Cl-
    NaCl <- NaCl(m_NaCl = mNaCl, T = T)
    basis("Cl-", log10(NaCl$m_Clminus))

    # Speciate aqueous sulfur
    bases <- c("H2S", "HS-", "HSO4-", "SO4-2")
    # Fe-bearing minerals
    names1 <- species1 <- c("iron", "ferrous-oxide", "pyrite", "pyrrhotite", "magnetite", "hematite")
    # Cu-bearing and FeCu-bearing minerals
    names2 <- species2 <- c("copper", "cuprite", "tenorite", "chalcocite", "covellite")
    names12 <- species12 <- c("bornite", "chalcopyrite")
    ## Uncomment these to use mineral formulas instead of names
    #names1 <- info(info(species1))$formula
    #names2 <- info(info(species2))$formula
    #names12 <- info(info(species12))$formula
    names <- list(names1, names2, names12)
    # Define plotting parameters
    lty <- list(2, 1, 1)
    lwd <- list(1, 1, 1)
    col <- list(2, 8, 8)
    col.names <- list(2, "#888888", "#888888")
    fill <- list(NA, NA, NA)

    dx <- list(c(0, 0, 2.5, 0.1, 0, 1), c(0, 0, 0, 0, 0), c(0, 0))
    dy <- list(c(0, 0, 0, 0.4, 0, -4), c(0, 0, 0, 0, 0), c(0, 0))
    srt <- list(c(0, 0, 0, 0, 0, 0), c(0, 0, 0, 0, 0), c(0, 0))

    # Create mosaic stack (Fe in species1, Cu in species2)
    sm <- stack_mosaic(bases, species1, species2, species12, names = names, col = col, col.names = col.names, fill = fill,
      dx = dx, dy = dy, srt = srt, lwd = lwd, lty = lty, pH = pH, O2 = O2, T = T, P = P, IS = NaCl$IS)
    # Because solid fill of Cu fields covers the Fe lines, replot them 20220815
    diagram(sm[[1]], add = TRUE, lty = lty[[1]], lwd = lwd[[1]], col = col[[1]], col.names = col.names[[1]],
      fill = fill[[1]], dx = dx[[1]], dy = dy[[1]], srt = srt[[1]])
    # Add water stability line 20220825
    water.lines(sm[[1]], lty = 3)
    # Also replot tick marks
    thermo.axis()
    # Add title
    main <- bquote("Fe-Cu-O-S, "*.(T)*"\u00b0C")
    title(main, font.main = 1)

  }

  layout(matrix(1:3, nrow = 1), widths = c(2, 2, 1))

  # Define temperature (degrees C) and logfO2 ranges
  T <- 150
  O2 <- c(-55, -35, res)

  # Cu-Co diagram
  stack_CuCo(T = T, pH = pH, O2 = O2)
  label.figure("a", font = 2, cex = 2)
  # Highlight bornite-chalcocite boundary
  lines(c(3.616049, 9.330938), c(-41.71679, -47.47867), col = 4, lwd = 2)
  # Highlight chalcopyrite-bornite boundary
  lines(c(3.495, 8.9586), c(-41.9348, -47.4475), col = 3, lwd = 2)

  # Cu-Fe diagram
  stack_FeCu(T = T, pH = pH, O2 = O2)
  label.figure("b", font = 2, cex = 2)
  # Highlight bornite-chalcocite boundary
  lines(c(3.616049, 9.330938), c(-41.71679, -47.47867), col = 4, lwd = 2)
  # Highlight chalcopyrite-bornite boundary
  lines(c(3.495, 8.9586), c(-41.9348, -47.4475), col = 3, lwd = 2)

  # Add legend
  par(mar = c(4, 0, 4, 0))
  plot.new()
  l <- c(
    lP(P),
    lNaCl(mNaCl),
    lS(Stot)
  )
  legend("topright", legend = lex(l), bty = "n", cex = 1.4, xpd = NA)

  legend <- as.expression(c("Bn-Cct", "Ccp-Bn", quote(H[2]*O~"stability limit")))
  legend("bottomleft", legend = legend, col = c(4, 3, 1), lty = c(1, 1, 3), lwd = 1.5, bty = "n", cex = 1.4, xpd = NA)

  if(pdf) dev.off()

}

# Experimental and fitted Cp for carrollite
# 20220816 jmd
carrollite_S3 <- function(pdf = FALSE) {

  if(pdf) pdf("Figure_S3.pdf", width = 5, height = 5)
  par(mar = c(4.1, 4.1, 1, 1))
  par(mgp = c(2.9, 1, 0))
  par(las = 1)

  ## Cp plot adapted from ../v10_Cp_and_Fe/specific_Cp.R
  ## 20220825
  # Define temperature (°C) values for calculating and plotting Cp
  TC_vals <- seq(25, 250, 1)
  # Define y-axis limits
  ylim <- c(140, 200)
  # Start plot
  plot(range(TC_vals), ylim, type = "n", xlab = axis.label("T"), ylab = axis.label("Cp"))

  # Assemble uDSC7 data
  TK <- uDSC7$Temp..K
  TC <- convert(TK, "C")
  Cp <- uDSC7$Cp..J.mol.K
  # Blue points for pre-transition
  blue <- adjustcolor(4, 0.8)
  ipre <- uDSC7$Stage == "Pre-transition"
  points(TC[ipre], Cp[ipre], pch = 19, cex = 0.5, col = blue)
  # Black points for transition
  black <- adjustcolor(1, 0.8)
  itran <- uDSC7$Stage == "Transition"
  points(TC[itran], Cp[itran], pch = 19, cex = 0.5, col = black)
  # Red points for post-transition
  red <- adjustcolor(2, 0.8)
  ipost <- uDSC7$Stage == "Post-transition"
  points(TC[ipost], Cp[ipost], pch = 19, cex = 0.5, col = red)

  # Calculate heat capacity
  # Plot line for pre-transition
  TC_pre <- TC_vals[TC_vals < (377 - 273.15)]
  spre <- subcrt("carrollite", T = TC_pre, P = 1)$out[[1]]
  lines(TC_pre, spre$Cp)
  # Plot line for post-transition
  TC_post <- TC_vals[TC_vals > (377 - 273.15)]
  spost <- subcrt("carrollite", T = TC_post, P = 1)$out[[1]]
  lines(TC_post, spost$Cp, lty = 2)

  legend <- as.expression(c(
    quote(mu*"DSC7 pre-transition"),
    quote(mu*"DSC7 transition"),
    quote(mu*"DSC7 post-transition"),
    "Pre-transition M-K fit",
    "Post-transition linear fit"
  ))
  legend("topright", legend = legend, pch = c(19, 19, 19, NA, NA), pt.cex = c(0.5, 0.5, 0.5, NA, NA),
    lty = c(NA, NA, NA, 1, 2), col = c(blue, black, red, 1, 1), cex = 0.8)

  if(pdf) dev.off()

}

# Experimental and fitted formation constants for Cl- complexes
# 20220721 jmd
# 'complexes' is result from add_Co_aqueous()
carrollite_S4 <- function(complexes, pdf = FALSE) {

  if(pdf) pdf("Figure_S4.pdf", width = 5, height = 5)
  par(mar = c(4, 4, 2, 1))
  par(mgp = c(2.5, 1, 0))

  # Plot 1: Cl
  plot(c(100, 300), c(-1, 9), xlab = axis.label("T"), ylab = quote(log~beta), type = "n")
  T <- seq(100, 300, 5)
  dy <- c(-0.3, 0.3, 0.3, 0.3)
  srt <- c(0, 12, 13, 14)
  for(i in 1:4) {
    logB.calc <- subcrt(complexes$species.Cl[[i]], complexes$coeff.Cl[[i]], T = T, P = "Psat")$out$logK
    lines(T, logB.calc, col = i)
    points(complexes$T_table3, complexes$logB_table3[[i]], col = i)
    if(i > 1) points(complexes$T_table4, complexes$logB_table4[[i+2]], pch = 0, col = i)
    sptxt <- expr.species(tail(complexes$species.Cl[[i]], 1))
    text(125, logB.calc[6] + dy[i], sptxt, srt = srt[i])
  }
  legend("topleft", c("Fit to spectroscopic data", "Spectroscopic data", "Solubility data"), lty = c(1, 0, 0), pch = c(NA, 1, 0), bty = "n")
  title("Co-Cl complexes", font.main = 1)

  if(pdf) dev.off()

}

# logK of reactions showing temperatures of linnaeite in and carrollite out
# 20230215 jmd
carrollite_S5 <- function(pdf = FALSE) {

  # Temperature and pressure for logK calculations
  T <- seq(100, 600)
  P <- 1
  # Start with empty plot
  if(pdf) pdf("Figure_S5.pdf", width = 7, height = 5)
  par(mar = c(4, 4, 1, 1))
  plot(range(T), c(-3, 1), type = "n", xlab = axis.label("T"), ylab = axis.label("logK"))
  abline(h = 0, lty = 2, col = 8)

  # Temporarily adjust Tmax for carrollite in order to calculate logK at higher temperatures
  Tmax <- info(info("carrollite", "cr2"))$T
  mod.OBIGT("carrollite", state = "cr2", T = 1000)

  # Balance reaction to form linnaeite
  basis(c("carrollite", "Co-pentlandite", "chalcocite"))
  # Set coefficient to get 10 S on each side of the reaction 20230217
  rxn_linnaeite <- subcrt("linnaeite", 2.2794, T = T, P = P)
  # Find the temperature where linnaeite becomes stable (where logK = 0)
  ilinn <- which.min(abs(rxn_linnaeite$out$logK))
  # Plot the logK values and temperature of linnaeite in
  lines(T, rxn_linnaeite$out$logK, col = 4)
  linntxt <- bquote("linnaeite in at"~.(T[ilinn])~degree*C)
  text(T[ilinn], rxn_linnaeite$out$logK[ilinn], linntxt, adj = c(0, 1.5))
  # Label the reaction
  rlinn <- rxn_linnaeite$reaction
  rlinn$coeff <- round(rlinn$coeff, 2)
  rlinn <- rlinn[rlinn$coeff != 0, ]
  text(185, -0.2, describe.reaction(rlinn), col = 4, srt = 50, cex = 0.7)

  # Balance reaction to react carrollite
  basis(c("cattierite", "linnaeite", "chalcocite"))
  rxn_carrollite <- subcrt("carrollite", -10/4, T = T, P = P)
  ilinn <- which.min(abs(rxn_carrollite$out$logK))
  lines(T, rxn_carrollite$out$logK, col = 2)
  linntxt <- bquote("carrollite out at"~.(T[ilinn])~degree*C)
  text(T[ilinn], rxn_carrollite$out$logK[ilinn], linntxt, adj = c(0.05, 1.5))
  # Label the reaction
  rcarr <- rxn_carrollite$reaction
  rcarr$coeff <- round(rcarr$coeff, 2)
  text(320, -1.2, describe.reaction(rcarr), col = 2, srt = 45, cex = 0.7)

  # Add legend for pressure 20230217
  legend("bottomright", legend = describe.property("P", P), bty = "n")

  # Restore Tmax
  mod.OBIGT("carrollite", state = "cr2", T = Tmax)

  if(pdf) dev.off()

}

# Figure_S6.R (renamed from carrollite_solubility_25.R)
# Compare single solubility contour for:
#   - Mosaic stack with carrollite
#   - solubility() function for only Cu or Co (no carrollite)
# 20220715 jmd first version (derived from CHNOSZ/demo/minsol.R)
# 20220721 Update aqueous Co complexes
# 20220816 Overlay contour from solubility()
carrollite_S6 <- function(res = 500, pdf = FALSE) {

  if(pdf) pdf("Figure_S6.pdf", width = 8, height = 3)
  logm_metal <- -5

  # To setup system and define parameters 20220715
  setup <- function(metal = "Cu", metal2 = NULL) {

    # System variables
    T <- 150
    P <- "Psat"
    Stot <- 3e-3
    pH <- c(2, 12, res)
    O2 <- c(-57, -32, res)
    # 15 % NaCl (Vasyukova and Williams-Jones, 2022)
    wNaCl <- 0.15

    # Set up basis species
    basis(c(metal, "H2S", "Cl-", "oxygen", "H2O", "H+"))
    basis("H2S", log10(Stot))
    # Molality of NaCl
    mNaCl <- 1000 * wNaCl / (mass("NaCl") * (1 - wNaCl))
    # Estimate ionic strength and molality of Cl-
    NaCl <- NaCl(m_NaCl = mNaCl, T = T, P = P)
    basis("Cl-", log10(NaCl$m_Clminus))

    # Add minerals and aqueous species
    icr <- retrieve(metal, c("H", "O", "S", "Cl"), state = "cr")
    iaq <- retrieve(metal, c("H", "O", "S", "Cl"), state = "aq")
    # Swap through these basis species to make mosaic diagram
    bases <- c("H2S", "HS-", "HSO4-", "SO4-2")
    # Set color for Cu and Co 20230210
    if(metal == "Cu") col <- 2 else col <- 4

    if(is.null(metal2)) {

      # Return list of parameters for use by other functions 20220715
      list(metal = metal, icr = icr, iaq = iaq, bases = bases, T = T, P = P, Stot = Stot, pH = pH, O2 = O2, mNaCl = mNaCl, NaCl = NaCl, col = col)

    } else {

      # For second metal 20220724
      # NOTE: Both metal2 and metal are needed to retrieve carrollite (it has both Co and Cu)
      icr2 <- retrieve(metal2, c("H", "O", "S", "Cl", metal), state = "cr")
      iaq2 <- retrieve(metal2, c("H", "O", "S", "Cl"), state = "aq")
      # The basis needs to have the first mineral in icr (metal 1)
      basis(c(metal2, names(icr)[1], "H2S", "Cl-", "oxygen", "H2O", "H+"))
      # Don't forget to also set S and Cl- activity!
      basis("H2S", log10(Stot))
      basis("Cl-", log10(NaCl$m_Cl))
      # Return list of parameters for use by other functions 20220715
      list(metal = metal, icr = icr, iaq = iaq, bases = bases, T = T, P = P, Stot = Stot, pH = pH, O2 = O2,
           mNaCl = mNaCl, NaCl = NaCl, metal2 = metal2, icr2 = icr2, iaq2 = iaq2, col = col)

    }

  }

  # To make stacked (two-metal) mineral-aqueous species diagram 20220724
  aqueous_stack <- function(metal, icr, iaq, bases, T, P, Stot, pH, O2, mNaCl, NaCl, metal2, icr2, iaq2, col) {

    species(icr)
    species(iaq, logm_metal, add = TRUE)
    mout <- mosaic(bases, pH = pH, O2 = O2, T = T, P = P, IS = NaCl$IS)
    d1 <- diagram(mout$A.species, plot.it = FALSE)

    species(icr2)
    species(iaq2, logm_metal, add = TRUE)
    m12 <- mosaic(list(bases, c(names(icr), names(iaq))), pH = pH, O2 = O2, T = T, IS = NaCl$IS, stable = list(NULL, d1$predominant), loga_aq = c(NA, logm_metal))
    names <- names(c(icr2, iaq2))
    dy <- dx <- srt <- rep(0, length(names))
    srt[names %in% c("Cu0.92Co2.07S4", "Co9S8")] <- 37
    dy[names %in% c("Cu0.92Co2.07S4", "Co9S8")] <- -0.5
    dy[names == "CoCl4-2"] <- 3
    dx[names == "Cu2S"] <- -2.8
    dy[names == "Cu2S"] <- 3.2
    srt[names %in% c("CoO", "CoO2-2", "Cu(OH)2-")] <- 90
    diagram(m12$A.species, srt = srt, dx = dx, dy = dy)
    #title(paste0(metal, "-", metal2, " stack"), font.main = 1)
    title(paste0(metal2, " solubility"), font.main = 1)

  }

  # To overlay single solubility contour
  overlay_solubility <- function(metal, icr, iaq, bases, T, P, Stot, pH, O2, mNaCl, NaCl, col) {

    species(icr)
    sout <- solubility(iaq, bases = bases, pH = pH, O2 = O2, T = T, P = P, IS = NaCl$IS, in.terms.of = metal)
    diagram(sout, add = TRUE, col = col, levels = logm_metal)

  }

  # Setup plot
  layout(matrix(1:3, nrow = 1), widths = c(1, 1, 0.8))

  # Aqueous species diagram: Cu-Co stack
  CuCo <- setup("Cu", "Co")
  do.call(aqueous_stack, CuCo)
  Co <- setup("Co")
  # Overlay solubility contour for Co minerals only
  do.call(overlay_solubility, Co)
  label.figure("a", font = 2, cex = 2)

  # Aqueous species diagram: Co-Cu stack
  CoCu <- setup("Co", "Cu")
  do.call(aqueous_stack, CoCu)
  Cu <- setup("Cu")
  # Overlay solubility contour for Cu minerals only
  do.call(overlay_solubility, Cu)
  label.figure("b", font = 2, cex = 2)

  # Add legend
  par(mar = c(4, 1, 4, 0))
  plot.new()

  l <- c(
    lT(Co$T),
    lP(Co$P),
    lNaCl(Co$mNaCl),
    lS(Co$Stot)
  )
  legend("topleft", legend = lex(l), bty = "n")
  legend("left", c("Co minerals only", "Cu minerals only", "Co or Cu minerals incl. carrollite"),
    lty = c(1, 1, 1), col = c(4, 2, 1), bty = "n")

  Culeg <- bquote(log ~ italic(m)*"Cu(aq) species" == .(logm_metal))
  Coleg <- bquote(log ~ italic(m)*"Co(aq) species" == .(logm_metal))
  legend <- as.expression(c(Culeg, Coleg))
  legend("bottomleft", legend = legend, bty = "n")

  # Done!
  if(pdf) dev.off()

}
