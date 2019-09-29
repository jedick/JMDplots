# project/filename: chnszten/plot.R
# Code for the paper "CHNOSZ: Thermodynamic calculations and diagrams for geochemistry"
# By: Jeffrey M. Dick
# Available at: https://doi.org/10.5281/zenodo.2648521

# History:
# 20170918 first version for early draft of manuscript
# 20190422 manuscript submission (Zenodo Version 1)
# 20190607 manuscript revision (Zenodo Version 2)
# 20190929 added to JMDplots package

# This code depends on CHNOSZ version 1.3.2 (https://cran.r-project.org/package=CHNOSZ)
# and was run under R version 3.6.0 on Linux x64 (https://www.r-project.org/).
# Additional system tools may be needed to support the cairo_pdf device, used in chnszten7() and chnsztenS7().
# chnszten1() requires CRAN packages timevis and shiny.

##########################
### Figure 1: timeline ###
##########################

# make timeline of CHNOSZ development 20170920 / updated 20190417
chnszten1 <- function() {
  # uses timevis package
  # timeline data
  timedata <- data.frame(matrix(c(
    "2008-03-06", "Website launched", # date from server log
    "2008-10-03", "First paper about CHNOSZ", # date from paper (doi: 10.1186/1467-4866-9-10)
    "2009-04-23", "Released on CRAN", # date from CRAN
    "2009-11-30", "Multivariable transects", # date from ONEWS
    "2010-09-30", "Vignette: An introduction", # date from anintro.Rmd
    "2011-08-23", "Vignette: Hot-spring proteins", # date from hotspring.Rnw
    "2012-06-16", "Development on R-Forge", # date from SVN log (Initial import)
    "2012-09-30", "Vignette: Equilibrium", # date from equilibrium.Rnw
    "2013-09-02", "Gibbs energy of transformation", # published date of Dick and Shock, 2013 paper in PLoS ONE
    "2015-03-01", "Variable-pressure standard state", # date from SVN log (r76: add thermo$opt$varP to use variable-pressure standard state for gases)
    "2014-12-20", "Mosaic diagrams", # date from SVN log
    "2017-02-27", "Vignette: Thermodynamic data", # date from SVN log (r177: add obigt.Rmd)
    "2017-10-02", "Berman, DEW, and B-dot equations", # date from SVN log (r237: add modifications to Berman data (1990 to 1992))
    "2018-10-31", "Solubility calculations", # date from SVN log
    "2019-02-26", "Akinfiev-Diamond model", # date from version 1.3.0 submitted to CRAN
    "2008-03-08", "0.6", # date from ONEWS
    "2008-09-14", "0.7", # date from ONEWS
    "2009-04-23", "0.8", # all remaining dates from CRAN
    "2009-12-01", "0.9", "2010-07-25", "0.9-1", "2011-08-23", "0.9-7",
    "2013-03-29", "1.0.0", "2014-01-12", "1.0.3", "2015-05-20", "1.0.5", "2016-05-29", "1.0.8",
    # note: 1.1.0 was really released on 2017-05-04; move it forward a little to make room for 1.1.3
    "2017-04-20", "1.1.0", "2017-11-13", "1.1.3", "2019-04-20", "1.3.2"
    ), byrow=TRUE, ncol=2
  ))
  colnames(timedata) <- c("start", "content")
  # group 1 is versions, group 2 is updates
  group <- as.numeric(!grepl("[0-9]", timedata$content)) + 1
  # use different styles for the groups and Vignettes
  className <- rep("updates", nrow(timedata))
  className[group==1] <- "versions"
  className[grepl("Vignette", timedata$content)] <- "vignettes"
  # assemble the timeline data
  timedata <- cbind(timedata, group=group, className=className)
  groups <- data.frame(id = 1:2, content = c("versions", "updates"))
  # simple page
  #timevis(timedata, groups=groups, showZoom=FALSE, options = list(showCurrentTime = FALSE))
  # use shiny app to include custom CSS
  # rather than read from a named file, use an anonymous file (adapted from ?connections)
  Tfile <- file()
  cat('.vignettes.vis-item { border-color: green; background-color: lightgreen; }
       .versions.vis-item { border-color: #F991A3; background-color: pink; }
       .versions.vis-item.vis-line { border-width: 0px; }
       .versions.vis-item.vis-dot { border-width: 0px; border-radius: 0px; }\n',
       file = Tfile)
  ui <- fluidPage(
    includeCSS(Tfile),
    timevisOutput("timeline")
  )
  close(Tfile)
  server <- function(input, output, session) {
    output$timeline <- renderTimevis({
      timevis(timedata, groups=groups, showZoom=FALSE, options = list(showCurrentTime = FALSE))
    })
  }
  shinyApp(ui = ui, server = server)
}

################################
### Figure 4: mosaic diagram ###
################################

chnszten4 <- function() {
  pdf("chnszten4.pdf", width=8, height=3)
  par(mfrow=c(1, 2), cex=1.3)
  mar <- c(2.5, 3, 1, 1)
  mgp <- c(1.5, 0.3, 0)
  basis(c("Cu", "H2S", "Cl-", "H2O", "H+", "e-"))
  basis("H2S", -6)
  basis("Cl-", -0.7)
  species(c("copper", "cuprite", "tenorite", "chalcocite", "covellite"))
  species(c("CuCl", "CuCl2-", "CuCl3-2", "CuCl+", "CuCl2", "CuCl3-", "CuCl4-2"))
  T <- 200
  res <- 500
  bases <- c("H2S", "HS-", "HSO4-", "SO4-2")
  m <- mosaic(bases, blend = TRUE, pH = c(0, 12, res), Eh=c(-1, 1, res), T=T)
  # first plot: S basis species
  diagram(m$A.bases, col = "blue", col.names = "blue", lty = 2, fill = NA, mar=mar, mgp=mgp, cex.names=0.9)
  legend("topright", legend=c(expression(italic(T)==200~degree*C), expression(italic(P)==15.5~bar)), bty="n", cex=0.8)
  label.figure("A", cex=1.3)
  # second plot: Cu species
  names <- species()$name
  names[c(1, 2, 4)] <- ""
  diagram(m$A.species, lwd = 2, fill = NA, mar=mar, mgp=mgp, cex.names=0.9, names=names)
  diagram(m$A.bases, add = TRUE, col = "blue", names = "", lty = 2)
  legend("topright", legend=c(expression(italic(T)==200~degree*C), expression(italic(P)==15.5~bar)), bty="n", cex=0.8)
  legend(-0.5, -0.4, legend=c(expression(sum(S)==10^-6~M), expression(sum(Cl)==0.2~M)), bty="n", cex=0.8)
  text(4.4, -0.32, "chalcocite", srt=-22, cex=0.9)
  text(9.5, -0.62, "copper", srt=-22, cex=0.9)
  text(9, -0.33, "cuprite", srt=-21, cex=0.9)
  label.figure("B", cex=1.3)
  invisible(dev.off())
}

####################################
### Figure 5: solubility diagram ###
####################################

chnszten5 <- function() {
  pdf("chnszten5.pdf", width=4, height=3)
  mar <- c(2.5, 3, 1, 1)
  mgp <- c(1.5, 0.3, 0)
  add.obigt("SLOP98")
  basis(c("corundum", "H2O", "H+", "O2"))
  species(c("Al+3", "AlO2-", "AlOH+2", "AlO+", "HAlO2"))
  a <- affinity(pH = c(0, 10), IS = 0)
  s <- solubility(a, in.terms.of = "Al+3")
  diagram(s, type = "loga.balance", ylim = c(-10, 0), lwd = 3, col = "green3", mar=mar, mgp=mgp)
  diagram(s, add = TRUE, adj = c(0, 1, 2.5, 0, -1.5), dy = c(0, 0, 4, 0, 0), names = c(-3, -4))
  legend("topright", c("25 \u00b0C", "1 bar"), bty = "n")
  # add labels with custom placement 20190606
  text(1, -0.7, expr.species("AlOH+2"))
  text(1, -4, expr.species("AlO+"))
  ## show neutral pH (pure water)
  logK <- subcrt(c("H2O", "H+", "OH-"), c(-1, 1, 1), T = 25)$out$logK
  #abline(v = -logK/2, col = "gray50")
  # show pH of solution (corundum in water)
  # calculate charge balance
  mol <- sapply(s$loga.equil, function(x) 10 ^ x)
  # include H+ and OH-
  pH <- a$vals$pH
  pOH <- -logK - pH
  molH <- 10 ^ -pH
  molOH <- 10 ^ -pOH
  mol <- cbind(mol, molH, molOH)
  # get charges of all species
  Z <- sapply(makeup(species()$ispecies), "[", "Z")
  Z[is.na(Z)] <- 0
  # include H+ and OH-
  Z <- c(Z, 1, -1)
  # calculate total charge in the solution
  Ztot <- colSums(t(mol) * Z)
  # where is it closest to zero (i.e. electroneutrality)
  imin <- which.min(abs(Ztot))
  # a line at the pH of electroneutrality
  pH0 <- pH[imin]
  abline(v = pH0, lty = 2, col = "gray50")
  # add labels
  #text(7.15, -3, "water", srt = 90, cex = 0.7, col = "gray50")
  text(6.55, -3, "equilibrium pH", srt = 90, cex = 0.7, col = "gray50")
  invisible(dev.off())
}

###########################
### Figure 6: DEW model ###
###########################

chnszten6 <- function() {
  pdf("chnszten6.pdf", width=5, height=5)
  par(cex = 1.2)
  # conditions:
  # T = 600, 700, 800, 900, 1000 degC
  # P = 5.0GPa (50000 bar)
  # fO2 = QFM - 2
  # pH set by jadeite + kyanite + coesite (approximated here as constant)
  # output from EQ3NR calculations (SSH14 Supporting Information):
  # dissolved carbon: 0.03, 0.2, 1, 4, 20 molal
  # true ionic strength: 0.39, 0.57, 0.88, 1.45, 2.49
  # pH: 3.80, 3.99, 4.14, 4.25, 4.33
  ## use DEW model
  water("DEW")
  # add species data for DEW
  inorganics <- c("methane", "CO2", "HCO3-", "CO3-2")
  organics <- c("formic acid", "formate", "acetic acid", "acetate", "propanoic acid", "propanoate")
  add.obigt("DEW", c(inorganics, organics))
  ## set basis species
  basis(c("Fe", "SiO2", "CO3-2", "H2O", "oxygen", "H+"))
  ## calculate logfO2 in QFM buffer
  basis("O2", "QFM")
  T <- seq(600, 1000, 100)
  buf <- affinity(T = T, P = 50000, return.buffer = TRUE)
  ## add species
  species(c(inorganics, organics))
  ## values of IS, pH, and molC at every 100 degC
  IS <- c(0.39, 0.57, 0.88, 1.45, 2.49)
  pH <- c(3.80, 3.99, 4.14, 4.25, 4.33)
  molC <- c(0.03, 0.2, 1, 4, 20)
  ## use Debye-Huckel equation with B-dot set to zero
  nonideal("Bdot0")
  ## calculate affinities on the T-logfO2-pH-IS transect
  a <- affinity(T = T, O2 = buf$O2 - 2, IS = IS, pH = pH, P = 50000)
  ## calculate metastable equilibrium activities using the total
  ## carbon molality as an approximation of total activity
  e <- equilibrate(a, loga.balance = log10(molC))
  ## make the diagram; don't plot names of low-abundance species
  names <- c(inorganics, organics)
  names[c(4, 5, 7, 9)] <- ""
  ## also exclude label for HCO3- - it will be manually positioned
  names[c(3, 8)] <- ""
  col <- rep("black", length(names))
  col[c(1, 3, 6, 8, 10)] <- c("red", "darkgreen", "purple", "orange", "navyblue")
  diagram(e, alpha = "balance", names = names, col = col, ylim = c(0, 0.8), lty=1, lwd=2,
          mar = c(3, 3, 1, 1), ylab="carbon fraction", spline.method="natural")
  ## add HCO3- label
  text(720, 0.03, expr.species("HCO3-"))

  ## add legend and title
  ltxt1 <- quote(italic(P) == 50000~bar)
  ltxt2 <- substitute(logfO2=="QFM-2", list(logfO2 = axis.label("O2")))
  legend("left", legend = as.expression(c(ltxt1, ltxt2)), bty = "n")

  ## clear settings for next calculation
  reset()
  invisible(dev.off())
}

#####################################
### Figure 7: Al-bearing minerals ###
#####################################

# reactions of Al-bearing minerals
chnszten7 <- function() {
  # use cairo_pdf for better handling of symbols (e.g. reaction double arrow)
  cairo_pdf("chnszten7.pdf", width = 7.2, height = 7.2)
  # the code for the figure is very close to demo("aluminum"),
  # but the line color for SUPCRT92 in the legend of (A) needs to be changed to blue,
  # so we include the entire code here 20190607
  #demo("aluminum", ask = FALSE)
  ## set up plotting area
  opar <- par(mfrow = c(2, 2))

  ###########
  ### plot 1: boehmite - kaolinite equilibrium
  ###########
  # After Zhu and Lu, 2009 (doi:10.1016/j.gca.2009.03.015)
  # experimental data from Table 1 of Hemley et al., 1980 (doi:10.2113/gsecongeo.75.2.210)
  xT <- c(200, 200, 200, 200, 250, 250, 265, 300, 300, 300, 300)
  xlogaSiO2 <- -c(2.54, 2.59, 2.65, 2.77, 2.21, 2.32, 2.12, 1.90, 1.95, 1.94, 1.90)
  ## set up basis species so that axis.label shows activity of SiO2
  basis(c("Al2O3","SiO2", "H2O", "O2"))
  T <- 125:350
  thermo.plot.new(xlim = range(T), ylim = c(-3.5, -1.5), xlab = axis.label("T"), ylab = axis.label("SiO2"))
  points(xT, xlogaSiO2)
  basis(delete = TRUE)
  ## first calculation: as in SUPCRT92
  add.obigt("SUPCRT92") # gets kaolinite and boehmite from HDNB78
  r1 <- subcrt(c("boehmite", "H2O", "SiO2", "kaolinite"), c(-1, -0.5, -1, 0.5), T = T, P = 1000, exceed.Ttr = TRUE) 
  # we need exceed.Ttr = TRUE because the T limit for boehmite is 500 K (Helgeson et al., 1978)
  ## second calculation: CHNOSZ default
  # kaolinite from Berman, 1988
  # boehmite from Hemingway et al., 1991
  # SiO2 from Apps and Spycher, 2004
  reset()
  r2 <- subcrt(c("boehmite", "H2O", "SiO2", "kaolinite"), c(-1, -0.5, -1, 0.5), T = T, P = 1000, exceed.Ttr = TRUE) 
  ## third calculation: get SiO2(aq) from SHS89
  add.obigt("AS04")
  r3 <- subcrt(c("boehmite", "H2O", "SiO2", "kaolinite"), c(-1, -0.5, -1, 0.5), T = T, P = 1000, exceed.Ttr = TRUE) 
  ## log activity of SiO2 is -ve logK
  lines(T, -r1$out$logK, col = "blue1", lty = 2)
  lines(T, -r2$out$logK, lwd = 1.5)
  lines(T, -r3$out$logK, col = "red", lty = 2)
  ## add points calculated using the SUPCRTBL package
  points(seq(125, 350, 25), -c(3.489, 3.217, 2.967, 2.734, 2.517, 2.314, 2.124, 1.946, 1.781, 1.628), pch = 4, col = "red")
  ## add legend and title
  title(main = describe.reaction(r1$reaction), cex.main = 1.1)
  legend("bottomright", lty = c(0, 2, 0, 1, 2), pch = c(1, NA, 4, NA, NA), lwd = c(1, 1, 1, 1.5, 1),
         col = c("black", "blue", "red", "black", "red"), bty = "n", cex = 0.9,
         legend = c("Hemley et al., 1980", "SUPCRT92", "SUPCRTBL", "CHNOSZ", 'add.obigt("AS04")'))
  legend("topleft", c("Boehmite - Kaolinite", "After Zhu and Lu, 2009 Fig. A1"), bty = "n")
  reset()
  # Helgeson et al., 1978 (HDNB78): http://www.worldcat.org/oclc/13594862
  # Shock et al., 1989 (SHS89): doi:10.1016/0016-7037(89)90341-4
  # Berman, 1988 (Ber88): doi:10.1093/petrology/29.2.445
  # Holland and Powell, 2011 (HP11): 10.1111/j.1525-1314.2010.00923.x
  # Hemingway et al., 1991 (HRA91): http://pubs.er.usgs.gov/publication/70016664
  # Apps and Spycher, 2004 (AS04): Bechtel SAIC Company, LLC ANL-NBS-HS-000043 REV 00 (DOC.20041118.0004)

  ###########
  ### plot 2: dawsonite solubility
  ###########
  # After Zimmer et al., 2016 (doi:10.1016/j.cageo.2016.02.013)
  # experimental data from Benezeth et al., 2007 Table 5 (doi:10.1016/j.gca.2007.07.003)
  # (averages for each temperature in a single run)
  T <- c(100.1, 100.1, 150.1, 100.1, 150.1, 99.8, 99.8, 200.7, 99.8, 50.1, 75.1, 100.3, 150.1)
  logK <- -c(14.825, 14.735, 13.625, 14.79, 13.665, 14.725, 14.1775, 12.74, 14.4925, 16.8625, 15.61, 14.51, 13.455)
  thermo.plot.new(c(25, 250), c(-18, -10), axis.label("T"), axis.label("logK"))
  points(T, logK)
  # calculation 1: CHNOSZ default
  T <- 0:250
  species <- c("dawsonite", "H2O", "Al(OH)4-", "HCO3-", "Na+", "H+")
  coeffs <- c(-1, -2, 1, 1, 1, 1)
  Daw1 <- subcrt(species, coeffs, T = T)
  lines(T, Daw1$out$logK, lwd = 1.5)
  # calculation 2: dawsonite with Cp = 0
  mod.obigt("dawsonite", Cp = 0, a = 0, b = 0, c = 0)
  Daw2 <- subcrt(species, coeffs, T = T)
  lines(T, Daw2$out$logK, col = "red", lty = 2)
  ## add points calculated using the SUPCRTBL package
  #points(seq(25, 250, 25), c(-17.829, -16.523, -15.402, -14.425, -13.568, -12.815, -12.154, -11.581, -11.094, -10.699), pch=4, col="red")
  ## 20190417: recalculated using the SUPCRTBL package (timestamp: 20190309)
  ##   with a locally updated data file that includes heat capacity coefficients of dawsonite
  ##   from Robie and Hemingway, 1995, with typos corrected in Tutolo et al., 2014
  points(seq(25, 250, 25), c(-17.829, -16.546, -15.485, -14.599, -13.856, -13.236, -12.724, -12.312, -11.997, -11.782), pch=4, col="red")
  ## add legend and title
  title(main = describe.reaction(Daw1$reaction), cex.main = 0.95)
  legend("bottomright", lty = c(0, 0, 0, 1, 2), pch = c(1, 4, NA, NA, NA), col = c("black", "red", NA, "black", "red"), lwd = c(1, 1, 0, 1.5, 1),
         bty = "n", cex = 0.9, legend = c("Ben\u00e9z\u00e9th et al., 2007", "SUPCRTBL with Cp", "  coefficients for dawsonite", "CHNOSZ", "Cp(dawsonite) = 0"))
  legend("topleft", c("Dawsonite solubility", "After Zimmer et al., 2016 Fig. 2"), bty = "n")
  reset()

  ###########
  ### plot 3: kaolinite solubility
  ###########
  # After Tutolo et al., 2014, Fig. 2 (doi:10.1016/j.gca.2014.02.036)
  dat <- read.csv(system.file("extdata/cpetc/TKSS14_Fig2.csv", package = "CHNOSZ"))
  thermo.plot.new(c(3.5, 1.5), c(-2, 14), quote(1000 / italic(T)*"(K)"), quote(p*italic(K)))
  points(dat)
  # plot line: default database
  invTK <- seq(3.5, 1.6, -0.02)
  T <- 1000/invTK - 273.15
  sres <- subcrt(c("kaolinite", "OH-", "H2O", "Al(OH)4-", "SiO2"), c(-1, -2, -1, 2, 2), T = T)
  pK <- -sres$out$logK
  lines(invTK, pK, lwd = 1.5)
  # plot line: SiO2 from Apps and Spycher, 2004
  add.obigt("AS04")
  sres <- subcrt(c("kaolinite", "OH-", "H2O", "Al(OH)4-", "SiO2"), c(-1, -2, -1, 2, 2), T = T)
  pK <- -sres$out$logK
  lines(invTK, pK, col = "red", lty = 2)
  reset()
  # plot line: SUPCRT92
  add.obigt("SUPCRT92")
  sres <- subcrt(c("kaolinite", "OH-", "H2O", "Al(OH)4-", "SiO2"), c(-1, -2, -1, 2, 2), T = T)
  pK <- -sres$out$logK
  lines(invTK, pK, col = "blue", lty = 2)
  # add points calculated using the SUPCRTBL package
  T <- seq(25, 300, 25)
  invTK <- 1000/(T + 273.15)
  points(invTK, c(12.621, 11.441, 10.383, 9.402, 8.477, 7.597, 6.756, 5.948, 5.171, 4.422, 3.703, 3.023), pch = 4, col = "red")
  # add title and legend
  par(xpd = NA)
  title(main = describe.reaction(sres$reaction), cex.main = 1.1)
  par(xpd = FALSE)
  legend("topright", c("Kaolinite solubility", "After Tutolo et al., 2014 Fig. 2"), bty = "n")
  legend("bottomleft", lty = c(0, 0, 2, 0, 1, 2), pch = c(1, NA, NA, 4, NA, NA), lwd = c(1, 1, 1, 1, 1.5, 1), col = c("black", "black", "blue", "red", "black", "red"),
         legend = c("Various sources \u2013", "  see Tutolo et al., 2014", "SUPCRT92", "SUPCRTBL", "CHNOSZ", 'add.obigt("AS04")'), bty = "n", cex = 0.9)
  reset()

  ###########
  ### plot 4: albite - K-feldspar exchange
  ###########
  # After Tutolo et al., 2014, Fig. 5 (doi:10.1016/j.gca.2014.02.036)
  # experimental data from Merino, 1975, Table 4 (doi:10.1016/0016-7037(75)90085-X)
  # plot line calculated using default database
  basis(c("Al2O3", "SiO2", "K+", "Na+", "O2", "H2O", "H+"))
  species(c("albite", "K-feldspar"))
  T <- 100
  P <- 150
  a <- affinity("K+" = c(4, 7), "Na+" = c(6, 9), T = T, P = P)
  diagram(a, lwd = 1.5, xlab = ratlab("K+"), ylab = ratlab("Na+"), names = NULL)
  # plot experimental data
  dat <- read.csv(system.file("extdata/cpetc/Mer75_Table4.csv", package = "CHNOSZ"))
  points(dat$log.aK..aH.., dat$log.aNa..aH..)
  # plot line calculated using SUPCRT92 data
  add.obigt("SUPCRT92")
  a <- affinity("K+" = c(4, 7), "Na+" = c(6, 9), T = 100, P = 150)
  diagram(a, col = "blue", lty = 2, add = TRUE, names = NULL)
  # add SUPCRTBL calculation
  logK_BL <- 2.092
  logaK <- seq(4, 7, 0.5)
  logaNa <- logaK + logK_BL
  points(logaK, logaNa, pch = 4, col = "red")
  # add title and legend
  sres <- subcrt(c("albite", "K+", "K-feldspar", "Na+"), c(-1, -1, 1, 1))
  title(main = describe.reaction(sres$reaction), cex.main = 1.1)
  legend("topleft", c("Albite - K-feldspar", "After Tutolo et al., 2014 Fig. 5"), bty = "n", cex = 0.9)
  legend("bottomright", lty = c(0, 2, 0, 1), pch = c(1, NA, 4, NA), lwd = c(1, 1, 1, 1.5), col = c("black", "blue", "red", "black"),
         legend = c("Merino, 1975", "SUPCRT92", "SUPCRTBL", "CHNOSZ"), bty = "n", cex = 0.9)
  legend("left", describe.property(c("T", "P"), c(T, P)), bty = "n")
  reset()

  par(opar)
  # add plot labels
  par(xpd = NA)
  text(3.5, 9.2, "A", cex = 1.5)
  text(5.4, 9.2, "B", cex = 1.5)
  text(3.5, 7.17, "C", cex = 1.5)
  text(5.4, 7.17, "D", cex = 1.5)
  par(xpd = FALSE)
  invisible(dev.off())
}

############################
### Supplemental Figures ###
############################

# calculations and plot for comparing the logK and maximum affinity methods 20170927
chnsztenS1 <- function() {
  # calcualate logK of reactions and logfO2 for activities = 10^-3
  # CO2 - CH4 ... use "oxygen" for the gas!
  K_1 <- suppressMessages(subcrt(c("H2O", "CO2", "CH4", "oxygen"), c(-2, -1, 1, 2), T=25, P=1)$out$logK)
  O_1 <- K_1 / 2
  print(paste("CO2 - CH4 logK", round(K_1, 2)))
  print(paste("CO2 - CH4 logfO2", round(O_1, 2)))
  # CH4 - acetic acid (unnumbered reaction in text)
  K_0 <- suppressMessages(subcrt(c("oxygen", "CH4", "acetic acid", "H2O"), c(-2, -2, 1, 2), T=25, P=1)$out$logK)
  O_0 <- (-K_0 -3 +6) / 2
  print(paste("CH4 - CH3COOH logK", round(K_0, 2)))
  print(paste("CH4 - CH3COOH logK", round(O_0, 2)))
  # CO2 - acetic acid
  K_2 <- suppressMessages(subcrt(c("H2O", "CO2", "acetic acid", "oxygen"), c(-2, -2, 1, 2), T=25, P=1)$out$logK)
  O_2 <- (K_2 +3 -6) / 2
  print(paste("CO2 - CH3COOH logK", round(K_2, 2)))
  print(paste("CO2 - CH3COOH logK", round(O_2, 2)))
  # affinities of R1 and R2 at activities of CH4 and acetic acid equal to 10^-3 (logfO2 = O_0)
  Q_1 <- -3 +2*O_0 +3
  A_1 <- K_1 - Q_1
  Q_2 <- -3 +2*O_0 +6
  A_2 <- K_2 - Q_2
  print(paste("CH4 affinity", round(A_1, 2)))
  print(paste("CH3COOH affinity", round(A_2, 2)))
  print(paste("CH3COOH affinity divided by 2 is", round(A_2 / 2, 2)))
  ## now calculate affinities with CHNOSZ
  basis("CHNOS")
  basis("O2", O_0)
  species(c("CH4", "acetic acid", "CO2"), -3)
  A <- unlist(affinity()$values)
  paste("affinities of CH4, CH3COOH and CO2 are", paste(round(A, 2), collapse=", "))
  # calculate and plot affinites as a function of logfO2
  a <- affinity(O2=c(-78, -68, 21))
  # divide acetic acid affinity by 2
  a$values[[2]] <- a$values[[2]] / 2
  # set up plot
  pdf("chnsztenS1.pdf", width=6.5, height=5)
  par(mar=c(3.5, 3.5, 0.5, 0.5), mgp=c(2.5, 1, 0), cex=1.2)
  plot(range(a$vals[[1]]), range(unlist(a$values)), type="n", xlab=axis.label("O2"), ylab="")
  mtext(expression(italic(A) / italic(n)[CO[2]]), side=2, line=2, cex=1.2)
  # balance-divided affinity lines
  lines(a$vals[[1]], a$values[[1]])
  points(a$vals[[1]], a$values[[1]], pch=19)
  text(-74.4, 4.9, expression(CH[4]), srt=-35)
  lines(a$vals[[1]], a$values[[2]])
  points(a$vals[[1]], a$values[[2]], pch=19)
  text(-73.5, -1.75, expression(CH[3]*COOH), srt=-20)
  lines(a$vals[[1]], a$values[[3]])
  points(a$vals[[1]], a$values[[3]], pch=19)
  text(-70, 0.6, expression(CO[2]))
  # show logK-calculated lines
  dx <- 0.2
  abline(v=O_1, lty=2)
  text(O_1-dx, -7, expression(CH[4]), srt=90, cex=0.9)
  text(O_1+dx, -7, expression(CO[2]), srt=90, cex=0.9)
  abline(v=O_0, lty=2)
  text(O_0-dx, -4.5, expression(CH[4]), srt=90, cex=0.9)
  text(O_0+dx, -4.5, expression(CH[3]*COOH), srt=90, cex=0.9)
  abline(v=O_2, lty=2)
  text(O_2-dx, -5, expression(CH[3]*COOH), srt=90, cex=0.9)
  text(O_2+dx, -5, expression(CO[2]), srt=90, cex=0.9)
  legend(-71.8, 12, legend=c(expression(italic(T)==25~degree*C), expression(italic(P)==1~bar)), bty="n", cex=0.9)
  invisible(dev.off())
}

# Figure4B modified to reproduce Fig. 5A of Caporuscio et al. (2017)
chnsztenS2 <- function() {
  pdf("chnsztenS2.pdf", width=4, height=4)
  mar <- c(2.5, 3, 1, 1)
  mgp <- c(1.5, 0.3, 0)
  add.obigt("SLOP98")
  basis(c("Cu", "H2S", "Cl-", "H2O", "H+", "e-"))
  basis("H2S", -6)
  basis("Cl-", -0.7)
  species(c("copper", "cuprite", "tenorite", "chalcocite", "covellite"))
  species(c("CuCl", "CuCl2-", "CuCl3-2", "CuCl+", "CuCl2", "CuCl3-", "CuCl4-2"))
  T <- 200
  res <- 500
  bases <- c("H2S", "HS-", "HSO4-", "SO4-2")
  m <- mosaic(bases, blend = TRUE, pH = c(0, 12, res), Eh=c(-1, 1, res), T=T)
  diagram(m$A.species, lwd = 2, fill = NA, mar=mar, mgp=mgp, cex.names=0.9, balance = 1, limit.water=FALSE)
  diagram(m$A.bases, add = TRUE, col = "blue", names = "", limit.water=FALSE)
  legend("topright", legend=c(expression(italic(T)==200~degree*C), expression(italic(P)==15.5~bar)), bty="n", cex=0.8)
  legend(-0.5, -0.4, legend=c(expression(sum(S)==10^-6~M), expression(sum(Cl)==0.2~M)), bty="n", cex=0.8)
  reset()
  invisible(dev.off())
}

# calculate Gibbs energy of transformation for an assemblage of n-alkanes 20190604
chnsztenS3 <- function() {
  # set up plot
  pdf("chnsztenS3.pdf", width = 7, height = 4)
  par(mfrow = c(1, 2), mar = c(3, 4, 1, 1))
  ## transforming an equilibrium assemblage of n-alkanes
  basis(c("CH4", "H2"), c("gas", "gas"))
  species(c("methane", "ethane", "propane", "butane"), "aq")
  # set logact to -1 so total activity is 1 (total number of C is 10)
  species(1:4, -1)
  # calculate equilibrium assemblages over a range of logaH2
  a <- affinity(H2 = c(-10, -6, 101), exceed.Ttr = TRUE)
  e <- equilibrate(a)
  # plot the equilibrium values and reference state
  diagram(e, names = "")
  legend("bottomleft", species()$name, lty = 1:4, bty = "n")
  abline(v = -8, lty = 2)
  text(-8.15, -3, "reference state", srt = 90)
  label.figure("A", cex = 1.5)
  # take a reference equilibrium distribution at logfH2 = -8
  loga1 <- list2array(e$loga.equil)[51, ]
  Astar <- list2array(e$Astar)[51, ]
  # calculate the DGtr compared to the equilibrium assemblage at logfH2 = -8
  DGtr.out <- numeric()
  for(i in 1:length(a$vals[[1]])) {
    loga2 <- list2array(e$loga.equil)[i, ]
    DGtr.out <- c(DGtr.out, DGtr(t(loga1), loga2, t(Astar)))
  }
  # plot the DGtr values
  thermo.plot.new(xlim = range(a$vals[[1]]), xlab = axis.label("H2"),
    ylim = range(DGtr.out), ylab = expr.property("DGtr/(2.303RT)"), mgp = c(2, 0.3, 0))
  abline(v = -8, lty = 2)
  text(-8.15, 0.1, "reference state", srt = 90)
  lines(a$vals[[1]], DGtr.out)
  label.figure("B", cex = 1.5)
  invisible(dev.off())
}

# findit() calculations for sulfur species 20190604
chnsztenS4 <- function() {
  # set up plot
  pdf("chnsztenS4.pdf", width = 8, height = 8)
  par(mfrow = c(2, 2))
  # the ranges of the parameters we will optimize
  vars <- list(O2 = c(-80, -40), pH = c(4, 14), T = c(0, 200))
  # set constant values of pH and T for the 1-D optimization
  # (and T for the 2-D optimization)
  pH <- 7
  T <- 25
  # set up chemical system
  basis("CHNOS+")
  basis("pH", pH)
  species(c("H2S", "S2-2", "S3-2", "S2O3-2", "S2O4-2", "S3O6-2",
    "S5O6-2", "S2O6-2", "HSO3-", "SO2", "HSO4-"))
  # objective function: standard deviations of the logarithms of activity the species
  objective <- "SD"
  # optimize logfO2 at constant T and pH
  f1 <- findit(vars[1], objective, T = T, niter = 4)
  title("1-D optimization", font.main = 1)
  legend("bottomright", c(paste("T =", T, "degC"), paste("pH =", pH)), bg = "white")
  label.figure("A", cex = 1.5)
  # optimize logfO2 and pH at constant T
  f2 <- findit(vars[1:2], objective, T = T, res = 20, niter = 5)
  title("2-D optimization", font.main = 1)
  legend("bottomright", c(paste("T =", T, "degC")), bg = "white")
  label.figure("B", cex = 1.5)
  # optimize logfO2, pH and T (at constant P ...)
  f3 <- findit(vars, objective, res = 20, niter = 5)
  title("3-D optimization", font.main = 1)
  label.figure("C", cex = 1.5)
  # the results
  f1.out <- sapply(f1$value, tail, 1)
  f2.out <- sapply(f2$value, tail, 1)
  f3.out <- sapply(f3$value, tail, 1)
  f1.out <- data.frame(nd = 1, t(f1.out))
  f2.out <- data.frame(nd = 2, t(f2.out))
  f3.out <- data.frame(nd = 3, t(f3.out))
  out <- merge(merge(f1.out, f2.out, all = TRUE), f3.out, all = TRUE)
  # show the results in a legend
  ltxt <- c("n", out$nd,
            "logfO2", formatC(out$O2, digits = 1, format = "f"),
            "pH", formatC(out$pH, digits = 1, format = "f"),
            "T", formatC(out$T, digits = 1, format = "f"),
            "SD", formatC(out$SD, digits = 3, format = "f"))
  ltxt[ltxt == "NA"] <- "-"
  opar <- par(mar = c(0, 0, 0, 0))
  plot.new()
  legend("center", ltxt, ncol = 5)
  text(0.5, 0.7, "Optimized parameters")
  par(opar)
  invisible(dev.off())
}

# Debye-Hückel extended term parameter extrapolated from plots of Manning et al., 2013
chnsztenS5 <- function() {
  pdf("chnsztenS5.pdf", width = 6, height = 6)
  bgamma(showsplines = "T")
  invisible(dev.off())
}

# Figure 6 modified to exclude DEW data for acetate
chnsztenS6A <- function() {
  pdf("chnsztenS6A.pdf", width=5, height=5)
  # conditions:
  # T = 600, 700, 800, 900, 1000 degC
  # P = 5.0GPa (50000 bar)
  # fO2 = QFM - 2
  # pH set by jadeite + kyanite + coesite (approximated here as constant)
  # output from EQ3NR calculations (SSH14 Supporting Information):
  # dissolved carbon: 0.03, 0.2, 1, 4, 20 molal
  # true ionic strength: 0.39, 0.57, 0.88, 1.45, 2.49
  # pH: 3.80, 3.99, 4.14, 4.25, 4.33
  ## use DEW model
  water("DEW")
  # add species data for DEW
  inorganics <- c("methane", "CO2", "HCO3-", "CO3-2")
  organics <- c("formic acid", "formate", "acetic acid", "acetate", "propanoic acid", "propanoate")
  # skip updating acetate because the new data from the DEW spreadsheet give different logK
  add.obigt("DEW", c(inorganics, organics[-4]))
  ## set basis species
  basis(c("Fe", "SiO2", "CO3-2", "H2O", "oxygen", "H+"))
  ## calculate logfO2 in QFM buffer
  basis("O2", "QFM")
  T <- seq(600, 1000, 100)
  buf <- affinity(T = T, P = 50000, return.buffer = TRUE)
  ## add species
  species(c(inorganics, organics))
  ## values of IS, pH, and molC at every 100 degC
  IS <- c(0.39, 0.57, 0.88, 1.45, 2.49)
  pH <- c(3.80, 3.99, 4.14, 4.25, 4.33)
  molC <- c(0.03, 0.2, 1, 4, 20)
  ## use Debye-Huckel equation with B-dot set to zero
  nonideal("Bdot0")
  ## calculate affinities on the T-logfO2-pH-IS transect
  a <- affinity(T = T, O2 = buf$O2 - 2, IS = IS, pH = pH, P = 50000)
  ## calculate metastable equilibrium activities using the total
  ## carbon molality as an approximation of total activity
  e <- equilibrate(a, loga.balance = log10(molC))
  ## make the diagram; don't plot names of low-abundance species
  names <- c(inorganics, organics)
  names[c(4, 5, 7, 9)] <- ""
  ## also exclude labels for HCO3- and acetate - they will be manually positioned
  names[c(3, 8)] <- ""
  col <- rep("black", length(names))
  col[c(1, 3, 6, 8, 10)] <- c("red", "darkgreen", "purple", "orange", "navyblue")
  if(packageVersion("CHNOSZ") > "1.1.3") {
    diagram(e, alpha = "balance", names = names, col = col, ylim = c(0, 0.8), lty=1, lwd=2,
            mar = c(3, 3, 1, 1), ylab="carbon fraction", spline.method="natural")
  } else {
    diagram(e, alpha = "balance", names = names, col = col, ylim = c(0, 0.8), lty=1, lwd=2,
            mar = c(3, 3, 1, 1), ylab="carbon fraction")
  }
  ## add HCO3- and acetate labels
  text(795, 0.035, "acetate")
  lines(c(767, 752), c(0.035, 0.0125))
  text(720, 0.04, expr.species("HCO3-"))
  lines(c(697, 680), c(0.045, 0.022))

  ## add legend and title
  ltxt1 <- quote(italic(P) == 50000~bar)
  ltxt2 <- substitute(logfO2=="QFM-2", list(logfO2 = axis.label("O2")))
  legend("left", legend = as.expression(c(ltxt1, ltxt2)), bty = "n")
  label.figure("A", xfrac = 0.02, yfrac = 0.97, cex = 1.7)

  ## check that we're within 0.1 of the QFM-2 values used by SSH14
  stopifnot(maxdiff((buf$O2-2), c(-17.0, -14.5, -12.5, -10.8, -9.4)) < 0.1)

  # Here are the logKs of aqueous species dissociation reactions at 600 degC and 50000 bar,
  # taken from the Supporting Information of the paper (p. 103-109):
  inorganic.logK <- c(24.4765, -9.0784, -5.3468, 0)
  organic.logK <- c(1.7878, 2.5648, 15.3182, 16.9743, 30.4088, 28.9185)
  # calculate equilibrium constants of the reactions in CHNOSZ; use a negative sign to change from formation to dissociation
  logK.calc <- -unlist(affinity(T = 600, P = 50000, property = "logK")$values)
  logK.calc - c(inorganic.logK, organic.logK)
  ## check that we're within 0.021 of the logK values used by SSH14
  stopifnot(maxdiff(logK.calc, c(inorganic.logK, organic.logK)) < 0.021)

  ## check that we get similar activity coefficients
  # activity coefficients for monovalent species from EQ3NR output
  loggamma <- c(-0.15, -0.18, -0.22, -0.26, -0.31)
  # activity coefficients calculated in CHNOSZ
  sres <- subcrt("propanoate", T = seq(600, 1000, 100), P = 50000, IS = c(0.39, 0.57, 0.88, 1.45, 2.49))
  stopifnot(maxdiff(sres$out[[1]]$loggam, loggamma) < 0.023)
  
  ## clear settings for next calculation
  reset()
  invisible(dev.off())
}

# Figure S6A modified to use default bgamma equation (non-zero extended term parameter extrapolated from Manning et al., 2013)
chnsztenS6B <- function() {
  pdf("chnsztenS6B.pdf", width=5, height=5)
  # conditions:
  # T = 600, 700, 800, 900, 1000 degC
  # P = 5.0GPa (50000 bar)
  # fO2 = QFM - 2
  # pH set by jadeite + kyanite + coesite (approximated here as constant)
  # output from EQ3NR calculations (SSH14 Supporting Information):
  # dissolved carbon: 0.03, 0.2, 1, 4, 20 molal
  # true ionic strength: 0.39, 0.57, 0.88, 1.45, 2.49
  # pH: 3.80, 3.99, 4.14, 4.25, 4.33
  ## use DEW model
  water("DEW")
  # add species data for DEW
  inorganics <- c("methane", "CO2", "HCO3-", "CO3-2")
  organics <- c("formic acid", "formate", "acetic acid", "acetate", "propanoic acid", "propanoate")
  # skip updating acetate because the new data from the DEW spreadsheet give different logK
  add.obigt("DEW", c(inorganics, organics[-4]))
  ## set basis species
  basis(c("Fe", "SiO2", "CO3-2", "H2O", "oxygen", "H+"))
  ## calculate logfO2 in QFM buffer
  basis("O2", "QFM")
  T <- seq(600, 1000, 100)
  buf <- affinity(T = T, P = 50000, return.buffer = TRUE)
  ## add species
  species(c(inorganics, organics))
  ## values of IS, pH, and molC at every 100 degC
  IS <- c(0.39, 0.57, 0.88, 1.45, 2.49)
  pH <- c(3.80, 3.99, 4.14, 4.25, 4.33)
  molC <- c(0.03, 0.2, 1, 4, 20)
  ## use Debye-Huckel equation with B-dot set to zero
  nonideal("bgamma")
  ## calculate affinities on the T-logfO2-pH-IS transect
  a <- affinity(T = T, O2 = buf$O2 - 2, IS = IS, pH = pH, P = 50000)
  ## calculate metastable equilibrium activities using the total
  ## carbon molality as an approximation of total activity
  e <- equilibrate(a, loga.balance = log10(molC))
  ## make the diagram; don't plot names of low-abundance species
  names <- c(inorganics, organics)
  names[c(4, 5, 7, 9)] <- ""
  ## also exclude label for HCO3- - it will be manually positioned
  names[c(3, 8)] <- ""
  col <- rep("black", length(names))
  col[c(1, 3, 6, 8, 10)] <- c("red", "darkgreen", "purple", "orange", "navyblue")
  diagram(e, alpha = "balance", names = names, col = col, ylim = c(0, 0.8), lty=1, lwd=2,
          mar = c(3, 3, 1, 1), ylab="carbon fraction", spline.method="natural")
  ## add HCO3- and acetate labels
  text(780, 0.027, "acetate")
  lines(c(750, 737), c(0.028, 0.0125))
  text(705, 0.04, expr.species("HCO3-"))
  lines(c(682, 665), c(0.045, 0.022))

  ## add legend and title
  ltxt1 <- quote(italic(P) == 50000~bar)
  ltxt2 <- substitute(logfO2=="QFM-2", list(logfO2 = axis.label("O2")))
  legend("left", legend = as.expression(c(ltxt1, ltxt2)), bty = "n")
  label.figure("B", xfrac = 0.02, yfrac = 0.97, cex = 1.7)

  ## clear settings for next calculation
  reset()
  invisible(dev.off())
}

# logK of NaCl dissociation
chnsztenS7 <- function() {
  # use cairo_pdf for better handling of symbols (e.g. reaction double arrow)
  cairo_pdf("chnsztenS7.pdf", width = 7, height = 7)
  demo("NaCl", ask = FALSE)
  invisible(dev.off())
}

# calcite solubility: comparison with Manning et al., 2013
chnsztenS8 <- function() {
  pdf("chnsztenS8.pdf", width=7, height=7)
  par(mfrow = c(1, 2))
  ## set pH range and resolution, constant temperature and ionic strength
  pH <- c(0, 14)
  res <- 100
  T <- 25
  IS <- 0

  ## now do calcite (a dissociation reaction)
  calfun <- function(dissociation = NULL) {
    basis(c("calcite", "Ca+2", "H2O", "O2", "H+"))
    species(c("CO2", "HCO3-", "CO3-2"))
    a <- affinity(pH = c(pH, res), T = T, IS = IS)
    s <- solubility(a, dissociation = dissociation)
    diagram(s, ylim = c(-10, 4), type = "loga.balance", lwd = 4, col = "green2")
    diagram(s, add = TRUE, dy = c(1, 0.7, 0.2))
    lexpr <- as.expression(c("total", expr.species("CO2", state = "aq"),
      expr.species("HCO3-"), expr.species("CO3-2")))
    legend("topright", lty = c(1, 1:3), lwd = c(4, 2, 2, 2),
      col = c("green2", rep("black", 3)), legend = lexpr)
    title(main = substitute("Solubility of"~what~"at"~T~degree*"C",
      list(what = "calcite", T = T)))
  }

  # considering the dissociation reactions with a shared ion
  calfun()
  mtext("all reactions give total activity of Ca+2")
  label.figure("A", cex = 1.7)
  # considering the dissociation reactions individually
  calfun(2)
  mtext("reactions considered individually")
  label.figure("B", cex = 1.7)
  invisible(dev.off())
}

# compare gold solubility in HCh and CHNOSZ - hematite-magnetite buffer
chnsztenS9 <- function() {
  # use Helgeson et al., 1978 minerals here
  add.obigt("SUPCRT92")
  # set up plot
  pdf("chnsztenS9.pdf", width = 9, height = 7)
  par(mfrow = c(2, 2))
  # log(m_Au) from Fig. 2A Williams-Jones et al., 2009 (doi:10.2113/gselements.5.5.281)
  WBM09_Fig2A <- data.frame(
    T = c(158.8, 185.8, 212.8, 238.5, 264.9, 290.5, 315.5, 342.6, 369.6, 393.2, 416.9, 441.2, 465.5, 490.5, 514.9, 540.5, 567.6),
    logm = c(-5.12, -5.23, -5.36, -5.54, -5.73, -5.96, -6.21, -6.32, -6.23, -5.93, -5.57, -5.25, -4.95, -4.64, -4.39, -4.18, -3.98)
  )
  # log(ppm Au) (magnetite-hematite buffer) from Fig. 8a of Pokrovski et al., 2009
  PAB14_Fig8a_MH <- data.frame(
    T = c(150, 200, 250, 300, 350, 400, 450, 500, 550),
    logppm = c(0.6, 0.72, 0.65, 0.31, -0.21, -0.62, -0.61, 0.3, 0.99)
  )
  # keep the solubility calculated in CHNOSZ to add to the HCh plot
  sol.out <- list()
  for(i in 1:3) {
    if(i==1) {
      molS <- 0.03  # high sulfide, like Pokrovski et al., 2014
      logH <- "QMK" # buffered pH
      m_NaCl = 1.5
      m_KCl = 0.5
      legend.x <- 320
      names.x <- c(245, 300, 520, 495)
      names.y <- c(-5.5, -6.3, -6.3, -5.7)
    }
    if(i==2) {
      molS <- 0.01  # low sulfide, like Williams-Jones et al., 2009
      logH <- "QMK" # buffered pH
      m_NaCl = 1.5
      m_KCl = 0.5
      legend.x <- 250
      names.x <- c(225, 250, 520, 495)
      names.y <- c(-6.3, -6.9, -6.3, -5.7)
    }
    if(i==3) {
      molS <- 0.03
      logH <- -5    # constant pH, like Pokrovski et al., 2014
      m_NaCl = 1.7
      m_KCl = 0
      legend.x <- 330
      names.x <- c(275, 300, 520, 510)
      names.y <- c(-5.5, -6.3, -6.3, -5.7)
    }
    ## plot Au molality calculated as in CHNOSZ/demo/gold.R
    # define colors for Au(HS)2-, Au(HS), AuOH, AuCl2-
    col <- c("#ED4037", "#F58645", "#0F9DE2", "#22CC88")
    # set up system
    # use H2S here: it's the predominant species at the pH of the QMK buffer -- see sulfur() in demo("gold")
    basis(c("Al2O3", "quartz", "Fe", "Au", "K+", "Cl-", "H2S", "H2O", "oxygen", "H+"))
    # set molality of K+ in completely dissociated 0.5 molal KCl
    # NOTE: This value is used only for making the legend;
    # activities corrected for ionic strength are computed below
    basis("K+", log10(0.5))
    # create a pH buffer
    mod.buffer("QMK", c("quartz", "muscovite", "K-feldspar"), "cr", 0)
    # log(m_Au)-T diagram like Fig. 2A of Williams-Jones et al., 2009 and Fig. 8a of Pokrovski et al., 2014
    species(c("Au(HS)2-", "AuHS", "AuOH", "AuCl2-"))
    basis("H2S", log10(molS))
    basis("H+", logH)
    # apply HM buffer for fO2
    basis("O2", "HM")
    # estimate solution composition for given amounts of NaCl and KCl
    chl <- chloride(T = seq(100, 550, 10), P = 1000, m_NaCl = m_NaCl, m_KCl = m_KCl)
    # calculate affinity and solubility
    if(i==3) a <- affinity(T = seq(100, 550, 10), `Cl-` = log10(chl$m_Cl), P = 1000, IS = chl$IS)
    else a <- affinity(T = seq(100, 550, 10), `Cl-` = log10(chl$m_Cl), `K+` = log10(chl$m_K), P = 1000, IS = chl$IS)
    s <- solubility(a)
    sol.out[[i]] <- s
    # make diagram and show total log molality
    diagram(s, ylim = c(-9, -4), col = col, lwd = 2, lty = 1, names = "")
    diagram(s, add = TRUE, type = "loga.balance", lwd = 3)
    # label lines
    names <- as.expression(sapply(s$species$name, expr.species))
    text(names.x, names.y, names)
    # make legend and title
    dS <- substitute(sum(S) == molS~italic(m), list(molS = molS))
    dpH <- describe.basis(ibasis = 10)
    dNaCl <- substitute(m_NaCl~mol~NaCl, list(m_NaCl = m_NaCl))
    dKCl <- substitute(m_KCl~mol~KCl, list(m_KCl = m_KCl))
    legend(legend.x, -4, c(dS, dpH, dNaCl, dKCl), bty = "n")
    # save the legend for the HCh plot
    if(i==2) l1.HCh <- c(dS, dpH, dNaCl, dKCl)
    dP <- describe.property("P", 1000)
    dO2 <- describe.basis(ibasis = 9)
    legend("bottomleft", c(dP, dO2), bty = "n")
    if(i==2) l2.HCh <- c(dP, dO2)
    if(i==1) title(main="CHNOSZ: high S, buffered pH", font.main=1)
    if(i==2) title(main="CHNOSZ: low S, buffered pH", font.main=1)
    if(i==3) title(main="CHNOSZ: high S, constant pH", font.main=1)
    ## add lines from WBM09 or PAB+14
    if(i==2) {
      lines(WBM09_Fig2A, lty=3, lwd=3)
      legend("bottomright", c("Williams-Jones et al., 2009", "CHNOSZ"), lty=c(3, 1), lwd=2, bg = "white")
    }
    if(i %in% c(1, 3)) {
      logm <- logppm2logm(PAB14_Fig8a_MH$logppm)
      lines(PAB14_Fig8a_MH$T, logm, lty=3, lwd=3)
      legend("bottomright", c("Pokrovski et al., 2014", "CHNOSZ"), lty=c(3, 1), lwd=2, bg = "white")
    }
    label.figure(LETTERS[i], cex=1.7, xfrac=0.03)
  }

  # results from HCh
  HCh <- data.frame(T = c(100L, 125L, 150L, 175L, 200L, 225L, 250L, 275L, 300L,
                          325L, 350L, 375L, 400L, 425L, 450L, 475L, 500L, 525L, 550L),
    Cl. = c(1.77, 1.74, 1.7, 1.66, 1.61, 1.57, 1.52, 1.48, 1.45,
            1.37, 1.29, 1.19, 1.1, 0.99, 0.89, 0.77, 0.66, 0.53, 0.43),
    AuCl2. = c(7.83e-18, 1.52e-16, 2.16e-15, 2.4e-15, 2.18e-13, 1.66e-12, 1.08e-11, 6.08e-11, 2.96e-10,
               1.42e-09, 6.23e-09, 2.48e-08, 9.08e-08, 3.07e-07, 9.61e-07, 2.79e-06, 7.55e-06, 1.88e-05, 4.22e-05),
    AuOH = c(1.91e-12, 1.03e-11, 4.42e-11, 1.55e-10, 4.66e-10, 1.21e-09, 2.83e-09, 5.97e-09, 1.15e-08,
             2.08e-08, 3.52e-08, 5.62e-08, 8.55e-08, 1.24e-07, 1.75e-07, 2.37e-07, 3.14e-07, 4.04e-07, 5.1e-07),
    AuHS = c(8.5e-09, 2.52e-08, 5.73e-08, 1.03e-07, 1.53e-07, 1.94e-07, 2.13e-07, 2.07e-07, 1.78e-07,
             1.57e-07, 1.29e-07, 9.79e-08, 6.67e-08, 3.94e-08, 1.96e-08, 8.17e-09, 2.93e-09, 9.56e-10, 2.79e-10),
    Au.HS.2. = c(8.06e-07, 1.34e-06, 1.72e-06, 1.71e-06, 1.36e-06, 9.13e-07, 5.3e-07, 2.75e-07, 1.28e-07,
                 5.79e-08, 2.43e-08, 9.18e-09, 2.94e-09, 7.4e-10, 1.36e-10, 1.81e-11, 1.82e-12, 1.51e-13, 1.11e-14),
    I = c(1.792, 1.75, 1.71, 1.667, 1.622, 1.576, 1.53, 1.49, 1.466,
          1.386, 1.3, 1.208, 1.11, 1.008, 0.9, 0.787, 0.669, 0.545, 0.445),
    pH = c(5.113, 4.987, 4.898, 4.832, 4.784, 4.752, 4.736, 4.738, 4.762,
           4.764, 4.777, 4.802, 4.841, 4.894, 4.964, 5.053, 5.164, 5.298, 5.494)
  )

  # plot Au molalities calculated in HCh 20181108
  thermo.plot.new(range(HCh$T), c(-9, -4), xlab=axis.label("T"), ylab=quote(log~italic(m)))
  lines(HCh$T, log10(HCh$Au.HS.2.), col=col[1], lwd=2)
  lines(HCh$T, log10(HCh$AuHS), col=col[2], lwd=2)
  lines(HCh$T, log10(HCh$AuOH), col=col[3], lwd=2)
  lines(HCh$T, log10(HCh$AuCl2.), col=col[4], lwd=2)
  Autot <- rowSums(HCh[, 3:6])
  lines(HCh$T, log10(Autot), lwd=3, lty=2)
  #title(main="HCh version 3.7 with updated KCl(aq), NaCl(aq) and Au complexes", font.main=1)
  title(main="HCh: low S, buffered pH", font.main=1)
  legend(250, -4, l1.HCh, bty = "n")
  legend("bottomleft", l2.HCh, bty = "n")
  ## add lines from CHNOSZ and WBM09
  lines(sol.out[[2]]$vals$T, sol.out[[2]]$loga.balance, lwd = 3)
  lines(WBM09_Fig2A, lty=3, lwd=3)
  legend("bottomright", c("Williams-Jones et al., 2009", "CHNOSZ", "HCh"), lty=c(3, 1, 2), lwd=2, bg = "white")
  label.figure("D", cex=1.7, xfrac=0.03)
  # label lines
  names <- as.expression(sapply(sol.out[[2]]$species$name, expr.species))
  names.x <- c(200, 230, 515, 495)
  names.y <- c(-6.3, -6.9, -6.2, -5.7)
  text(names.x, names.y, names)

  # close plot and reset obigt
  obigt()
  invisible(dev.off())
}

# compare gold solubility in HCh and CHNOSZ - pyrite-pyrrhotite-magnetite buffer
chnsztenS10 <- function() {
  # use Helgeson et al., 1978 minerals here
  add.obigt("SUPCRT92")
  # set up plot
  pdf("chnsztenS10.pdf", width = 9, height = 7)
  par(mfrow = c(2, 2))
  # log(m_Au) from Fig. 2B Williams-Jones et al., 2009 (doi:10.2113/gselements.5.5.281)
  WBM09_Fig2B <- data.frame(
    T = c(157.5, 177.2, 196.9, 218.6, 238.9, 261.3, 283.6, 305.3, 328.3, 350.7, 373.7, 396.1, 418.4, 441.4, 463.8, 488.8),
    logm = c(-9.31, -8.87, -8.47, -8.09, -7.72, -7.35, -7, -6.67, -6.33, -6.01, -5.68, -5.34, -4.99, -4.65, -4.36, -4.09)
  )
  # log(ppm Au) (pyrite-pyrrhotite-magnetite buffer) from Fig. 8a of Pokrovski et al., 2009
  PAB14_Fig8a_PPM <- data.frame(
    T = c(250, 300, 350, 400, 450, 500, 550),
    logppm = c(-3.12, -2.23, -1.56, -1.06, -0.6, -0.1, 0.5)
  )
  # keep the solubility calculated in CHNOSZ to add to the HCh plot
  sol.out <- list()
  for(i in 1:2) {
    if(i==1) {
      logH <- "QMK" # buffered pH
      m_NaCl = 1.5
      m_KCl = 0.5
      names.x <- c(510, 520, 450, 385)
      names.y <- c(-8.4, -4.5, -8, -9)
    }
    if(i==2) {
      logH <- -5 # constant pH
      m_NaCl = 1.7
      m_KCl = 0
      names.x <- c(510, 520, 455, 400)
      names.y <- c(-8.4, -4.5, -7.9, -9)
    }
    ## plot Au molality calculated as in CHNOSZ/demo/gold.R
    # define colors for Au(HS)2-, Au(HS), AuOH, AuCl2-
    col <- c("#ED4037", "#F58645", "#0F9DE2", "#22CC88")
    # set up system
    # use H2S here: it's the predominant species at the pH of the QMK buffer -- see sulfur() in demo("gold")
    basis(c("Al2O3", "quartz", "Fe", "Au", "K+", "Cl-", "H2S", "H2O", "oxygen", "H+"))
    # set molality of K+ in completely dissociated 0.5 molal KCl
    # NOTE: This value is used only for making the legend;
    # activities corrected for ionic strength are computed below
    basis("K+", log10(0.5))
    # create a pH buffer
    mod.buffer("QMK", c("quartz", "muscovite", "K-feldspar"), "cr", 0)
    # log(m_Au)-T diagram like Fig. 2A of Williams-Jones et al., 2009 and Fig. 8a of Pokrovski et al., 2014
    species(c("Au(HS)2-", "AuHS", "AuOH", "AuCl2-"))
    basis("H+", logH)
    # apply PPM buffer for fO2 and aH2S
    basis("O2", "PPM")
    basis("H2S", "PPM")
    # estimate solution composition for given amounts of NaCl and KCl
    chl <- chloride(T = seq(100, 550, 10), P = 1000, m_NaCl = m_NaCl, m_KCl = m_KCl)
    # calculate affinity and solubility
    if(i==2) a <- affinity(T = seq(100, 550, 10), `Cl-` = log10(chl$m_Cl), P = 1000, IS = chl$IS)
    else a <- affinity(T = seq(100, 550, 10), `Cl-` = log10(chl$m_Cl), `K+` = log10(chl$m_K), P = 1000, IS = chl$IS)
    s <- solubility(a)
    sol.out[[i]] <- s
    # make diagram and show total log molality
    diagram(s, ylim = c(-12, -4), col = col, lwd = 2, lty = 1, names = "")
    diagram(s, add = TRUE, type = "loga.balance", lwd = 3)
    # label lines
    names <- as.expression(sapply(s$species$name, expr.species))
    text(names.x, names.y, names)
    if(i==1) lines(c(500, 500, NA, 520, 520), c(-8.1, -6.6, NA, -4.7, -6))
    if(i==2) lines(c(500, 500, NA, 520, 520), c(-8.1, -6.5, NA, -4.7, -6))
    # make legend and title
    dpH <- describe.basis(ibasis = 10)
    dNaCl <- substitute(m_NaCl~mol~NaCl, list(m_NaCl = m_NaCl))
    dKCl <- substitute(m_KCl~mol~KCl, list(m_KCl = m_KCl))
    legend(230, -4, c(dpH, dNaCl, dKCl), bty = "n")
    # save the legend for the HCh plot
    if(i==1) l1.HCh <- c(dpH, dNaCl, dKCl)
    dP <- describe.property("P", 1000)
    dO2 <- describe.basis(ibasis = 9)
    dS <- describe.basis(ibasis = 7)
    legend("topleft", c(dP, dO2, dS), bty = "n")
    if(i==1) l2.HCh <- c(dP, dO2, dS)
    if(i==1) title(main="CHNOSZ: buffered pH", font.main=1)
    if(i==2) title(main="CHNOSZ: constant pH", font.main=1)
    ## add lines from WBM09 or PAB+14
    if(i==1) {
      lines(WBM09_Fig2B, lty=3, lwd=3)
      legend("bottomright", c("Williams-Jones et al., 2009", "CHNOSZ"), lty=c(3, 1), lwd=2, bg = "white")
    }
    if(i==2) {
      logm <- logppm2logm(PAB14_Fig8a_PPM$logppm)
      lines(PAB14_Fig8a_PPM$T, logm, lty=3, lwd=3)
      legend("bottomright", c("Pokrovski et al., 2014", "CHNOSZ"), lty=c(3, 1), lwd=2, bg = "white")
    }
    label.figure(LETTERS[i], cex=1.7, xfrac=0.03)
  }

  # results from HCh
  HCh <- data.frame(T = c(100L, 125L, 150L, 175L, 200L, 225L, 250L, 275L, 300L,
                          325L, 350L, 375L, 400L, 425L, 450L, 475L, 500L, 525L, 550L),
    Cl. = c(1.75, 1.73, 1.7, 1.66, 1.61, 1.57, 1.52, 1.48, 1.45, 
            1.37, 1.29, 1.19, 1.1, 1, 0.89, 0.78, 0.66, 0.54, 0.45),
    AuCl2. = c(4.27e-19, 1.16e-17, 2.08e-16, 2.73e-15, 2.85e-14, 2.45e-13, 1.78e-12, 1.1e-11, 5.93e-11,
               3.12e-10, 1.49e-09, 6.43e-09, 2.53e-08, 9.13e-08, 3.03e-07, 9.29e-07, 2.63e-06, 6.87e-06, 1.63e-05),
    AuOH = c(1.38e-13, 8.92e-13, 4.46e-12, 1.81e-11, 6.17e-11, 1.81e-10, 4.71e-10, 1.1e-09, 2.35e-09,
             4.63e-09, 8.53e-09, 1.47e-08, 2.41e-08, 3.74e-08, 5.57e-08, 7.95e-08, 1.09e-07, 1.46e-07, 1.9e-07),
    AuHS = c(3.73e-14, 5.28e-13, 5e-12, 3.41e-11, 1.76e-10, 7.23e-10, 2.43e-09, 6.95e-09, 1.71e-08,
             3.74e-08, 7.35e-08, 1.31e-07, 2.16e-07, 3.3e-07, 4.75e-07, 6.44e-07, 8.32e-07, 1.02e-06, 1.22e-06),
    Au.HS.2. = c(2.81e-16, 7.61e-16, 1.35e-13, 1.62e-12, 1.37e-11, 8.61e-11, 4.21e-10, 1.69e-09, 5.92e-09,
                 1.5e-08, 3.29e-08, 6.38e-08, 1.1e-07, 1.75e-07, 2.54e-07, 3.4e-07, 4.22e-07, 4.83e-07, 5.47e-07),
    I = c(1.813, 1.756, 1.71, 1.665, 1.619, 1.571, 1.525, 1.484, 1.458,
          1.378, 1.292, 1.2, 1.108, 1, 0.894, 0.783, 0.66, 0.551, 0.475),
    pH = c(5.228, 5.034, 4.916, 4.841, 4.79, 4.757, 4.742, 4.744, 4.768, 
           4.77, 4.784, 4.809, 4.847, 4.9, 4.97, 5.059, 5.169, 5.299, 5.465)
  )

  # plot Au molalities calculated in HCh 20181108
  thermo.plot.new(range(HCh$T), c(-12, -4), xlab=axis.label("T"), ylab=quote(log~italic(m)))
  lines(HCh$T, log10(HCh$Au.HS.2.), col=col[1], lwd=2)
  lines(HCh$T, log10(HCh$AuHS), col=col[2], lwd=2)
  lines(HCh$T, log10(HCh$AuOH), col=col[3], lwd=2)
  lines(HCh$T, log10(HCh$AuCl2.), col=col[4], lwd=2)
  Autot <- rowSums(HCh[, 3:6])
  lines(HCh$T, log10(Autot), lwd=3, lty=2)
  #title(main="HCh version 3.7 with updated KCl(aq), NaCl(aq) and Au complexes", font.main=1)
  title(main="HCh: buffered pH", font.main=1)
  legend(230, -4, l1.HCh, bty = "n")
  legend("topleft", l2.HCh, bty = "n")
  ## add lines from CHNOSZ and WBM09
  lines(sol.out[[2]]$vals$T, sol.out[[2]]$loga.balance, lwd = 3)
  lines(WBM09_Fig2B, lty=3, lwd=3)
  legend("bottomright", c("Williams-Jones et al., 2009", "CHNOSZ", "HCh"), lty=c(3, 1, 2), lwd=2, bg = "white")
  label.figure("C", cex=1.7, xfrac=0.03)
  # label lines
  names <- as.expression(sapply(sol.out[[2]]$species$name, expr.species))
  names.x <- c(510, 520, 455, 385)
  names.y <- c(-8.3, -4.5, -7.7, -8.8)
  text(names.x, names.y, names)
  lines(c(500, 500, NA, 520, 520), c(-8, -6.4, NA, -4.7, -6))

  # close plot and reset obigt
  obigt()
  invisible(dev.off())
}

############################
### UNEXPORTED FUNCTIONS ###
############################

## Function used in chnsztenS9 and chnsztenS10
# estimate the Cl- molality and ionic strength for a hypothetical 
# NaCl solution with total chloride equal to specified NaCl + KCl solution,
# then estimate the molality of K+ in that solution 20181109
chloride <- function(T, P, m_NaCl, m_KCl) {
  NaCl <- NaCl(T = T, P = P, m_tot = m_NaCl + m_KCl)
  # calculate logK of K+ + Cl- = KCl, adjusted for ionic strength
  logKadj <- subcrt(c("K+", "Cl-", "KCl"), c(-1, -1, 1), T = T, P = P, IS = NaCl$IS)$out$logK
  # what is the molality of K+ from 0.5 mol/kg KCl, assuming total chloride from above
  m_K <- m_KCl / (10^logKadj * NaCl$m_Cl + 1)
  list(IS = NaCl$IS, m_Cl = NaCl$m_Cl, m_K = m_K)
}

## Function used in chnsztenS9 and chnsztenS10
# convert log(ppm) to log(m) 20190418
logppm2logm <- function(logppm) {
  ppm <- 10^logppm
  # how many grams of Au per kg of H2O?
  # 1 ppm = 1 mg Au / kg H2O
  gAu <- ppm / 1000
  # how many moles of Au per kg of H2O
  massAu <- mass("Au")
  mAu <- gAu / massAu
  # return log(m)
  log10(mAu)
}

