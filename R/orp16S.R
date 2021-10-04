# JMDplots/orp16S.R
# Plots for paper on ZC-ORP correlations 20210827

## To run these functions, chem16Sdir needs to be changed 20211003
#orp16Sdir = system.file("extdata/orp16S", package = "JMDplots")
#options(chem16Sdir = orp16Sdir)

# Group studies by environment types 20210828
envirotype <- list(
  "River & Seawater" = c("MLL+18", "SVH+19", "HXZ+20", "KLY+20", "WHL+21", "LXH+20", "JVW+20", "ZZL+21", "GZL21"),
  "Lake & Pond" = c("SAR+13", "ECS+18", "LLC+19", "BCA+21", "HLZ+18", "CNA+20", "BWD+19", "RSJ+21", "BOEM21", "FAV+21"),
  "Groundwater" = c("KLM+16", "SDH+19", "SRM+19", "SKP+21", "YHK+20", "JDP+20", "GWS+19", "SRM+21", "ZCZ+21"),
  # NOTE: Keep Hot Spring at #4 to get red color 20210904
  "Hot Spring" = c("SMS+12", "BMJ+19", "LMG+20", "GWSS21", "GWS+20", "PBU+20", "MWY+21", "OFY+19"),
  "Alkaline Spring" = c("SBP+20", "RMB+17", "CTS+17", "KSR+21", "NTB+21"),
  "Sediment" = c("JHL+12", "GFE+16", "ZML+17", "BYB+17", "BSPD17", "ABT+17", "HDZ+19", "SCM+18",
                 "CLS+19", "ZDA+20", "VMB+19", "HSF+19", "MCS+21", "LMBA21", "ZZLL21", "WFB+21"),
  "Soil" = c("SBW+17", "MLL+19", "BMOB18", "ZZZ+18", "PMM+20", "ZHZ+19", "CWC+20", "PSG+20", "XLD+20", "LJC+20", "DTJ+20", "LLL+21", "DLS21")
)
# Turn the list into a data frame for easier lookup 20210904
envirodat <- do.call(rbind, lapply(seq_along(envirotype), function(i) data.frame(study = envirotype[[i]], groupnum = i)))
envirodat <- cbind(envirodat, group = names(envirotype)[envirodat$groupnum])

# Set the palette for R colors (numeric index of 'col') 20210914
orp16Scol <- palette.colors(n = length(envirotype), palette = "Classic Tableau", alpha = 0.75)

# Figure 1: Geobiochemical predictive framework 20210830
orp16S1 <- function(pdf = FALSE) {

  if(pdf) pdf("Figure1.pdf", width = 7, height = 5)

  # Define colors; adapted from
  # https://cdn.elifesciences.org/author-guide/tables-colour.pdf
  Purple <- "#9E86C9"
  Purple80 <- "#9E86C980"
  PurpleText <- "#9C27B0"
  Blue80 <- "#90CAF980"
  BlueText <- "#366BFB"
  Red80 <- "#E5737380"
  RedText <- "#D50000"
  Orange80 <- "#FFB74D80"
  OrangeText <- "#f7831c"
  Gray <- "#E6E6E6"

  # Setup plot
  par(mar = c(1, 1, 2, 1))
  plot(c(5, 95), c(1, 80), type = "n", axes = FALSE, xlab = "", ylab = "")
  # Uncomment this as a guide for making the 'grid' one 20210927
  #box(lwd = 2)
  grid.roundrect(0.5, 0.48, 0.94, 0.88, gp = gpar(fill = Gray))

  # Plot bottom box and text for correlative framework
  # Uncomment this as a guide for making the 'grid' one 20211002
  #rect(5, 2, 95, 40, col = Purple80, border = Purple, lty = 2)
  grid.roundrect(0.5, 0.28, 0.88, 0.4, gp = gpar(col = Purple, fill = Purple80, lty = 2))
  text(50, 32, "Conventional\nCorrelative Framework", col = PurpleText, cex = 1.2)

  # Plot lines to go behind shapes
  lines(c(20, 50), c(20, 60))
  lines(c(20, 50), c(60, 60))
  lines(c(50, 80), c(60, 60))
  # Add arrows along lines 20210927
  arrows(20, 20, 20*0.35 + 50*0.65, 20*0.35 + 60*0.65, length = 0.1)
  arrows(20, 60, 37, 60, length = 0.1)
  arrows(50, 60, 66, 60, length = 0.1)

  # Plot shapes and text for biological methods
  text(20, 79, "Biological Methods", col = RedText, font = 2)
  for(bg in c("white", Red80)) points(20, 60, pch = 21, cex = 18, bg = bg)
  text(20, 66, "RefSeq", font = 2, col = RedText)
  text(20, 57, "Reference\nproteomes\nof taxa")
  for(bg in c("white", Red80)) points(20, 20, pch = 21, cex = 18, bg = bg)
  text(20, 25, "16S + RDP", font = 2, col = RedText)
  text(20, 18, "Taxonomic\nabundances")

  # Plot shapes and text for chemical methods
  text(80, 79, "Chemical Methods", col = BlueText, font = 2)
  for(bg in c("white", Blue80)) points(80, 60, pch = 22, cex = 18, bg = bg)
  text(80, 68, quote(bolditalic(Z)[bold(C)]), col = BlueText, cex = 1.2)
  text(80, 58.5, "Carbon\noxidation\nstate")
  # Show multiple physicochemical variables 20210927
  # Function to draw rectangle at x,y with width and height w,h
  myrect <- function(x, y, w, h, ...) rect(x - w/2, y - h/2, x + w/2, y + h/2, ...)
  # T, Eh, pH
  for(col in c("white", Blue80)) myrect(73, 28, 5, 6, col = col)
  text(73, 28, "T", col = BlueText)
  for(col in c("white", Blue80)) myrect(80, 30, 7, 7, col = col, lwd = 2)
  text(80, 30, "Eh", font = 2, cex = 1.2, col = BlueText)
  for(col in c("white", Blue80)) myrect(87.5, 28, 6, 6, col = col)
  text(87.5, 28, "pH", col = BlueText)
  # Salinity, Metals
  for(col in c("white", Blue80)) myrect(75, 20, 11, 6, col = col)
  text(75, 20, "Salinity", col = BlueText)
  for(col in c("white", Blue80)) myrect(86.5, 20, 10, 6, col = col)
  text(86.5, 20, "Metals", col = BlueText)
  # O2, CO2, SO4
  for(col in c("white", Blue80)) myrect(72, 12, 5, 6, col = col)
  text(72, 12, quote(O[2]), col = BlueText)
  for(col in c("white", Blue80)) myrect(79, 12, 7, 6, col = col)
  text(79, 12, quote(CO[2]), col = BlueText)
  for(col in c("white", Blue80)) myrect(87.5, 12, 8, 6, col = col)
  text(87.5, 12, quote(SO[4]^"2-"), col = BlueText)

  # Plot inner rectangle and text
  third <- 100/3
  # Uncomment this to make the original rectangle (as a guide for the 'grid' one)
  #for(bg in c("white", Orange80)) rect(1.2*third, 50, 1.8*third, 70, col = bg)
  # Use this to get rounded corners 20210927
  for(fill in c("white", Orange80)) grid.roundrect(0.5, 0.685, 0.19, 0.21, gp = gpar(fill = fill))
  text(50, 67, quote(bold(C[bolditalic(c)]*H[bolditalic(h)]*N[bolditalic(n)]*O[bolditalic(o)]*S[bolditalic(s)])), col = OrangeText)
  text(50, 58, "Estimated\nCommunity\nProteomes")

  # Plot arrows and text labels
  arrows(1*third - 1, 20, 2*third + 1, 20, code = 3, lty = 2, length = 0.1)
  text(50, 20, "Multivariate\ncorrelations")
  arrows(80, 35, 80, 51, code = 3, lwd = 1.5, length = 0.1)
  text(80, 44, "Univariate  prediction", font = 2)

  # Add title
  title("Geobiochemical Predictive Framework", font.main = 1)

  if(pdf) dev.off()

}

# Figure 2: Depth profiles in Winogradsky Columns 20210829
orp16S2 <- function(pdf = FALSE) {

  if(pdf) pdf("Figure2.pdf", width = 7, height = 5)

  layout(t(matrix(1:2)), widths = c(3, 4))
  par(mar = c(4, 4, 1, 1), mgp = c(2.5, 1, 0))

  # ORP measurements from Diez-Ercilla et al. (2019)
  depth <- c(25, 22, 18, 14, 12, 9, 5)
  ORP <- c(-229, -241, -229, -201, 302, 501, 641)
  pH <- c(5.5, 5.5, 5.5, 4.7, 3.0, 2.7, 2.6)
  SWIdepth <- 11.3
  # Shift depth values so SWI is at zero
  MODdepth <- depth - SWIdepth
  plot(ORP, MODdepth, ylim = c(14.7, -6.3), type = "b", lty = 2, xlab = "ORP or Eh7 (mV)", ylab = "Depth (cm)", xlim = c(-400, 800), las = 1)
  abline(h = 0, lty = 2, col = "darkgray", lwd = 1.5)
  # Calculate and plot Eh7 20210926
  Eh7 <- ORP + -59.16 * (7 - pH)
  lines(Eh7, MODdepth, type = "b", pch = 19)
  text(50, -2.5, "Eh7")
  text(700, -2.5, "ORP")
  label.figure("A", font = 2, cex = 1.5)

  # Get 16S metadata and chemical metrics for Rundell et al. (2014) experiments
  mdat <- getmdat("RBW+14")
  metrics <- getmetrics("RBW+14")
  # Get ZC values for each layer
  layers <- c("12 cm", "8 cm", "4 cm", "SWI", "Top")
  ZC <- lapply(layers, function(layer) metrics$ZC[mdat$layer == layer])
  # Make boxplots
  boxplot(ZC, horizontal = TRUE, show.names = FALSE, xlab = axis.label("ZC"), ylim = c(-0.18, -0.145), yaxs = "i")
  axis(2, 1:5, labels = layers, las = 1)
  # Add sample sizes
  par(xpd = NA)
  for(i in 1:5) {
    label <- bquote(italic(n) == .(length(ZC[[i]])))
    text(-0.182, i - 0.25, label, adj = 1)
  }
  par(xpd = FALSE)
  label.figure("B", font = 2, cex = 1.5)

  if(pdf) dev.off()

}

# Figure 3: Sample locations on world map
orp16S3 <- function(pdf = FALSE) {

  if(pdf) pdf("Figure3.pdf", width = 26, height = 14)

  # Coordinates for orp16S datasets
  file <- tempfile()
  # Write spaces here (but don't save them in the file) to make this easier to read
  writeLines(con = file, text = gsub(" ", "", c(
    # Column names
    "study, latitude, longitude",
    # River & Seawater (comments indicate source of coordinates from paper or SRA metadata)
    "MLL+18, 22.20, 113.09", # Fig. 1
    "SVH+19, 42.90, 30.68", # Materials and methods
    "HXZ+20, 16.52, 111.77", # Table 1
    "KLY+20, 37.60189, 126.8128", # Materials and methods
    "WHL+21, 30.76, 115.36", # SAMN13242327
    "LXH+20, 25.42, 99.34", # SAMN15090995
    "JVW+20, 45.432025, 12.260878", # SAMN15796698
    "ZZL+21, 22.77, 113.79", # SAMN16964962
#    "GZL21, 29.568, 106.668", # SAMN19460485
    # Lake & Pond
    "SAR+13, 1.96, -157.33", # Materials and methods
    "ECS+18, -33.65, -70.117", # Materials and methods
    "LLC+19, 24.82, 118.15", # SAMN04549101
    "BCA+21, 46.3615, 25.0509", # SAMN07409474
    "HLZ+18, 24.795, 118.138", # SAMN07638080
    "CNA+20, 39.441, -77.371", # Web search for geographic center of Maryland --> https://sos.maryland.gov/mdkids/Pages/Geography.aspx
    "BWD+19, 47.120571, -88.545425", # SAMN09980099
    "RSJ+21, 61.833, 24.283", # Materials and methods
    "BOEM21, 43.051389, -75.965", # Materials and methods
    "FAV+21, 0.757, 36.372", # SAMN19267646
    # Hot Spring
    "SMS+12, 44.6, -110.9", # JGI IMG/M sample name 1_050719N
    "BMJ+19, 14.089567, 40.348583", # SAMN11581539
    "LMG+20, -37.855, -71.158", # Table 1
    "GWSS21, 24.86, 98.33", # SAMN16802401
    "GWS+20, 30.12, 101.94", # SAMN13430433
    "PBU+20, 54.4395, 160.144194", # SAMN14538724
#    "MWY+21, 30.5848, 81.5806", # Table 1
    "OFY+19, 31.228, 130.613", # Online Resource 1
    # Alkaline Spring
    "SBP+20, 38.862, -122.414", # SAMN03850954
    "RMB+17, 22.9052, 58.6606", # SAMN05981641
    "CTS+17, 10.94323, -85.63485", # SAMN06226041
    "KSR+21, 44.264340, 8.46442", # SAMN17101425
    "NTB+21, 22.881, 58.701", # SAMN19998441
    # Soil - put this group before Groundwater and Sediment for clearer visualization in GBA 20210927
    "SBW+17, 28.25, 116.92", # Materials and methods  ### Laboratory
    "MLL+19, 26.1, 112.5", # Materials and methods
    "BMOB18, 40.60842, -74.19258", # SAMN07828017  ### Laboratory
    "ZZZ+18, 21.816, 112.464", # Materials and methods  ### Laboratory
    "PMM+20, 43.397, -80.311",  # Web search for Cambridge, ON, Canada  ### Laboratory
    "ZHZ+19, 23.130, 113.671", # Materials and methods
    "CWC+20, 28.226, 116.898", # Materials and methods  ### Laboratory
    "PSG+20, 36.61, -119.53", # Web search for Parlier, CA   ### Mesocosm
    "XLD+20, 27.35, 112.05", # Materials and methods   ### Laboratory
#    "LJC+20, 35.19, 118.69", # SAMN14149974
    "DTJ+20, 26.45, 111.52", # SAMN14332759   ### Laboratory
    "LLL+21, 27.78, 113.13", # Materials and methods   ### Laboratory
    "DLS21, 39.39, -75.44", # SAMN17245435  ### Mesocosm
    # Groundwater
    "KLM+16, 42.99, -82.30", # SAMN04423023
    "SDH+19, 23.03, 113.38", # SAMN07692244
    "SRM+19, 12.67417, 101.3889", # Materials and methods
    "SKP+21, 16.263306, 100.647778", # SAMN11191517
    "YHK+20, 51.209467, 10.791968", # SAMEA5714424
    "JDP+20, 44.8883, 110.1353", # SAMN12236980
    "GWS+19, 36.31, 94.81", # SAMN07765433
    "SRM+21, 14.83, 99.35", # SAMN14829351
    "ZCZ+21, 45.21, 9.57", # Table 1
    # Sediment
    "JHL+12, 73.566167, 8.1585", # SAMN00744894
    "GFE+16, -36.69, -73.07", # Table 1
    "ZML+17, 22.494, 114.029", # Table 1
    "BYB+17, 41.423, -112.087", # Table 1
    "BSPD17, 57.89297, 16.5855", # SAMN05163191   ### Laboratory
    "ABT+17, 43.42, -2.7", # Fig. 1
    "HDZ+19, 29.901, 113.52435", # SAMN05990289
    "SCM+18, 25.25, -97.23", # Table 1
    "CLS+19, 32.22, 118.83", # SAMN08683376
    "ZDA+20, -21.42996, -70.05874", # SAMEA4858706
    "VMB+19, 52.11, 79.17", # SAMN08987150
    "HSF+19, 47.803, 16.709", # methods
    "MCS+21, -32.15, -71.1", # Materials and methods
    "LMBA17, 43.42, -2.7", # Fig. 1
    "ZZLL21, 22.68, 113.97", # SAMN16964887
    "WFB+21, 56.440893, -2.863194", # methods   ### Mesocosm
    "NA, NA, NA"
  )))

  # This reads the data
  coords <- read.csv(file, as.is = TRUE)

  # Start map
  # https://cran.r-project.org/web/packages/oce/vignettes/map_projections.html
  par(mar = c(2, 0.5, 0, 0.5))
  # We don't need data(coastlineWorld) ... it's the default map 20211003
  mapPlot(col = "lightgray", projection = "+proj=wintri", border = "white", drawBox = FALSE)

  # Add Great Lakes
  # https://www.sciencebase.gov/catalog/item/530f8a0ee4b0e7e46bd300dd
  # https://gis.stackexchange.com/questions/296170/r-shapefile-transform-longitude-latitude-coordinates
  hydro_pdir <- system.file("extdata/orp16S/hydro_p", package = "JMDplots")
  for(lake in c("Erie", "Huron", "Michigan", "Ontario", "Superior", "StClair")) {
    file <- file.path(hydro_pdir, paste0("hydro_p_Lake", lake, ".shp"))
    tx <- st_read(file, quiet = TRUE)
    # https://mhweber.github.io/AWRA_2020_R_Spatial/coordinate-reference-systems.html
    tx_ll <- st_transform(tx, "+proj=longlat +ellps=GRS80 +datum=NAD83")
    gl.coords <- st_coordinates(tx_ll)
    mapPolygon(gl.coords[, "X"], gl.coords[, "Y"], col = "white", border = NA)
  }

  # Plot transects 20210929
  # Coordinates for East Asia Paddy Soil dataset are from
  # Sourcedata.xlsx from https://doi.org/10.6084/m9.figshare.12622829
  dat <- getmdat("LJC+20")
  mapPoints(dat$longitude, dat$latitude, col = orp16Scol[7], lwd = 2)
  # Coordinates for Southern Tibetan Plateau dataset are from Table 1 of MWY+21
  dat <- getmdat("MWY+21")
  mapPoints(dat$longitude, dat$latitude, col = orp16Scol[4], lwd = 2)
  # Coordinates for Three Gorges Reservoir are from Table S1 of GZL21
  dat <- getmdat("GZL21")
  dat <- dat[!is.na(dat$"ORP (mV)"), ]
  latlon <- paste(dat$Latitude, dat$Longitude)
  isuniq <- !duplicated(latlon)
  dat <- dat[isuniq, ]
  mapPoints(dat$Longitude, dat$Latitude, col = 1, lwd = 1)

  # Get colors for studies
  icol <- envirodat$groupnum[match(coords$study, envirodat$study)]
  # Identify studies that use samples from laboratory or mesocosm experiments
  lab <- c(
    "BSPD17", "WFB+21", # Sediment
    "SBW+17", "BMOB18", "ZZZ+18", "PMM+20", "ZHZ+19", "CWC+20", "PSG+20", "XLD+20", "DTJ+20", "LLL+21", "DLS21" # Soil
  )
  pch <- ifelse(coords$study %in% lab, 15, 19)
  # Use smaller points for high-density regions 20210915
  cex <- ifelse(coords$study %in% c(
    "MLL+19", "XLD+20", "LLL+21", "DTJ+20", # Hunan
    "ZZL+21", "MLL+18", "SDH+19", "ZML+17", "ZZLL21", "ZZZ+18", "ZHZ+19" # GD-HK-MO GBA
#    "JVW+20", "ZCZ+21", "KSR+21"  # Northern Italy
  ), 1.5, 2.5)
  # Plot sample locations
  mapPoints(coords$longitude, coords$latitude, pch = pch, col = orp16Scol[icol], cex = cex)

  # Add legend
  ienv = c(1, 2, 4, 5, 3, 6, 7)
  par(xpd = NA)
  legend("bottomleft", names(envirotype)[ienv], pch = 19, col = orp16Scol[ienv], bty = "n", cex = 2, inset = c(0, -0.03))
  ltext <- c("Field samples", "Laboratory or mesocosm", "Smaller symbols for", "densely sampled areas", "Open symbols for transects")
  legend("bottomright", ltext, pch = c(19, 15, 20, NA, 1), bty = "n", cex = 2, pt.cex = c(2, 2, 2, 2, 1), inset = c(0, -0.03))
  par(xpd = FALSE)

  if(pdf) dev.off()

}

# Figure 4: Selected plots for each environment type 20211003
orp16S4 <- function(pdf = FALSE) {

  if(pdf) pdf("Figure4.pdf", width = 12, height = 6)
  par(mfrow = c(2, 4))
  par(mar = c(4, 4, 3, 1))
  par(cex.lab = 1.1)
  par(mgp = c(2.5, 1, 0))

  # Black Sea (River & Seawater)
  plotZC("SVH+19", "Bacteria", groupby = "Type", groups = c("Oxic", "Suboxic", "Euxinic"), title.line = 1.5)
  # Fayetteville Green Lake (Lake & Pond)
  plotZC("BOEM21", "Bacteria", groupby = "Stratum", groups = c("Upper", "Chemocline", "Lower"), title.line = 1.5)
  # Samail Ophiolite (Alkaline Spring)
  plotZC("RMB+17", "Bacteria", groupby = "pH Group", groups = c("< 10", "> 10"), title.line = 1.5)
  # Eastern Tibetan Plateau (Hot Spring)
  plotZC("GWS+20", "Bacteria", groupby = "Field", groups = c("Batang", "Litang", "Kangding"), title.line = 1.5)
  # Hainich Critical Zone (Groundwater)
  plotZC("YHK+20", "Bacteria", groupby = "Location", groups = c("Upper Hillslope", "Middle Slope", "Lower Footslope"), title.line = 1.5)
  # Po Plain (Groundwater)
  plotZC("ZCZ+21", "Bacteria", groupby = "Location", groups = c("LO", "CR1", "MN", "VA", "BS", "CR2"), legend.x = "topright", title.line = 1.5)
  # Bay of Biscay (Sediment)
  plotZC("LMBA21", "Bacteria", groupby = "Year", groups = c(2013, 2015, 2016, 2017, 2018, 2019), legend.x = "bottomright", title.line = 1.5)
  # Hunan Soil (Soil)
  plotZC("MLL+19", "Bacteria", groupby = "Type", groups = c("Upland", "Paddy", "Sediment"), title.line = 1.5)

  if(pdf) dev.off()

}

# Figure 5: Local and global analysis of ZC-Eh7 correlations 20210828
orp16S5 <- function(pdf = FALSE) {

  if(pdf) pdf("Figure5.pdf", width = 10, height = 8)
  mat <- matrix(c(
    0, rep(1:2, each = 9), 0, 0, rep(3:4, each = 9), 0,
    rep(0, 40),
    18, 18, rep(5:10, each = 6), 19, 19,
    18, 18, rep(11:16, each = 6), 20, 20,
    rep(17, 40),
    rep(0, 40),
    rep(0, 5), rep(21, 13), rep(22, 13), rep(23, 9)
  ), nrow = 7, byrow = TRUE)
  layout(mat, heights = c(3, 0.5, 2, 2, 0.5, 0.5, 4))
  par(mar = c(4, 4, 1, 1))
  par(mgp = c(2.5, 1, 0))

  # Read linear fit coefficients
  dat <- read.csv(system.file("extdata/orp16S/EZlm.csv", package = "JMDplots"))
  # Use Bacteria only 20210913
  dat <- dat[dat$lineage == "Bacteria", ]

  ## Panel A: Linear fits for each datasets

  # River & seawater, lake & pond, hot spring, alkaline spring
  plot(c(-400, 750), c(-0.23, -0.12), type = "n", xlab = "Eh7 (mV)", ylab = cplab$ZC)
  i1 <- c(1, 2, 4, 5)
  for(i in i1) {
    studies <- envirotype[[i]]
    for(j in seq_along(studies)) {
      with(dat[dat$study == studies[j], ], {
        y <- function(x) intercept + slope * x
        Eh7 <- c(Eh7min, Eh7max)
        ZC <- y(Eh7)
        # Use dashed line for fits with low slope 20211002
        if(abs(slope * 1000) < 0.01) lty <- 2 else lty <- 1
        lines(Eh7, ZC, col = orp16Scol[i], lty = lty)
      })
    }
  }
  legend("bottomright", names(envirotype)[i1], lty = 1, col = orp16Scol[i1], cex = 0.8, bty = "n")
  label.figure("A", font = 2, cex = 2, xfrac = 0, yfrac = 0.9)

  # Groundwater, sediment, soil
  plot(c(-400, 750), c(-0.23, -0.12), type = "n", xlab = "Eh7 (mV)", ylab = cplab$ZC)
  i2 <- c(3, 6, 7)
  for(i in i2) {
    studies <- envirotype[[i]]
    for(j in seq_along(studies)) {
      with(dat[dat$study == studies[j], ], {
        y <- function(x) intercept + slope * x
        Eh7 <- c(Eh7min, Eh7max)
        ZC <- y(Eh7)
        if(abs(slope * 1000) < 0.01) lty <- 2 else lty <- 1
        lines(Eh7, ZC, col = orp16Scol[i], lty = lty)
      })
    }
  }
  legend("bottomright", names(envirotype)[i2], lty = 1, col = orp16Scol[i2], cex = 0.8, bty = "n")

  ## Panel B: Slope vs log10(number of samples)

  # Get color according to environment group
  env <- envirodat[match(dat$study, envirodat$study), ]
  # Get range for included samples 20210905
  j1 <- env$groupnum %in% i1
  j2 <- env$groupnum %in% i2
  xlim <- range(log10(dat$nsamp[j1 | j2]))
  # Multiply by 1000 to use V instead of mV 20210913
  ylim <- range(1000 * dat$slope[j1 | j2])
  # River & seawater, lake & pond, hot spring, alkaline spring
  plot(xlim, ylim, type = "n", xlab = "log10(Number of samples)", ylab = quote("Slope of linear fit"~(V^-1)))
  abline(h = 0, lty = 4, col = "gray")
  points(log10(dat$nsamp[j1]), 1000 * dat$slope[j1], pch = 19, col = orp16Scol[env$groupnum[j1]])
  label.figure("B", font = 2, cex = 2, xfrac = -0.05, yfrac = 0.9)
  # Groundwater, sediment, soil
  plot(xlim, ylim, type = "n", xlab = "log10(Number of samples)", ylab = quote("Slope of linear fit"~(V^-1)))
  abline(h = 0, lty = 4, col = "gray")
  points(log10(dat$nsamp[j2]), 1000 * dat$slope[j2], pch = 19, col = orp16Scol[env$groupnum[j2]])

  ## Panel C: Scatterplots and fits for each environment

  par(mar = c(0, 0, 1, 0))
  eachenv("Bacteria")
  par(mar = c(1, 0, 0, 0))
  eachenv("Archaea")
  # Add labels
  plot.new()
  text(0.5, -0.5, "Eh7 (mV)", cex = 1.2, xpd = NA)
  plot.new()
  text(0.2, 0.5, cplab$ZC, cex = 1.2, srt = 90, xpd = NA)
  label.figure("C", font = 2, cex = 2, xfrac = 0.2, yfrac = 1)
  plot.new()
  text(0.3, 0.5, "Bacteria", srt = 90, xpd = NA)
  plot.new()
  text(0.3, 0.5, "Archaea", srt = 90, xpd = NA)

  ## Panel D: Scatterplots and fits for Bacteria and Archaea in all environments 20210914
  # Read Eh7 - ZC data
  dat <- read.csv(system.file("extdata/orp16S/EZdat.csv", package = "JMDplots"))
  # Omit data for Hot Spring environment
  dat <- dat[dat$envirotype != "Hot Spring", ]

  # Start plot for Bacteria
  par(mar = c(4, 4, 1, 1))
  plot(c(-500, 650), c(-0.22, -0.12), type = "n", xlab = "Eh7 (mV)", ylab = cplab$ZC)
  # Use Bacteria only
  thisdat <- dat[dat$lineage == "Bacteria", ]
  # Add linear fit; include number of studies in legend 20210925
  nstudy <- length(unique(thisdat$study))
  add.linear(thisdat$Eh7, thisdat$ZC, nstudy)
  # Add points
  eachenv("Bacteria", add = TRUE, do.linear = FALSE)
  title("Bacteria", font.main = 1)
  label.figure("D", font = 2, cex = 2, xfrac = 0.02, yfrac = 1)

  # Now do Archaea
  plot(c(-500, 650), c(-0.22, -0.12), type = "n", xlab = "Eh7 (mV)", ylab = cplab$ZC)
  thisdat <- dat[dat$lineage == "Archaea", ]
  nstudy <- length(unique(thisdat$study))
  add.linear(thisdat$Eh7, thisdat$ZC, nstudy)
  eachenv("Archaea", add = TRUE, do.linear = FALSE)
  title("Archaea", font.main = 1)

  # Add legend
  par(mar = c(4, 1, 1, 1))
  plot.new()
  ienv = c(1, 2, 5, 3, 6, 7)
  legend("left", names(envirotype)[ienv], pch = 19, col = orp16Scol[ienv], bty = "n")

  if(pdf) dev.off()
}

# Figure 6: ZC-Eh7 fits for Bacteria and Archaea in hot springs 20210930
orp16S6 <- function(pdf = FALSE) {

  if(pdf) pdf("Figure6.pdf", width = 8, height = 3.5)
  par(mfrow = c(1, 3))
  par(mar = c(4, 4, 2, 1))
  par(mgp = c(2.5, 1, 0))

  # Read linear fit coefficients
  alldat <- read.csv(system.file("extdata/orp16S/EZlm.csv", package = "JMDplots"))
  # Use Hot Spring datasets
  i <- 4
  studies <- envirotype[[i]]
  # Assign colors
  # Use blue for neutral/alkaline
  col <- rep(orp16Scol[1], length(studies))
  # Use red for acidic
  iacid <- studies %in% c("LMG+20")
  col[iacid] <- orp16Scol[4]
  # Use turquoise for hypersaline
  ihyper <- studies %in% c("BMJ+19")
  col[ihyper] <- "turquoise3"
  # Use gray for sediments
  ised <- studies %in% c("PBU+20", "MWY+21", "OFY+19")
  col[ised] <- "gray"

  # Loop over Bacteria and Archaea
  lineages <- c("Bacteria", "Archaea")
  for(k in 1:2) {
    dat <- alldat[alldat$lineage == lineages[k], ]
    plot(c(-450, 400), c(-0.23, -0.08), type = "n", xlab = "Eh7 (mV)", ylab = cplab$ZC)
    for(j in seq_along(studies)) {
      with(dat[dat$study == studies[j], ], {
        y <- function(x) intercept + slope * x
        Eh7 <- c(Eh7min, Eh7max)
        ZC <- y(Eh7)
        if(length(Eh7) > 0) {
          if(abs(slope * 1000) < 0.01) lty <- 2 else lty <- 1
          lines(Eh7, ZC, col = col[j], lwd = 2, lty = lty)
          # Add number to identify dataset
          text(tail(Eh7, 1), tail(ZC, 1), j)
        }
      })
    }
    title(lineages[k], font.main = 1)
    label.figure(LETTERS[k], font = 2, cex = 1.5, yfrac = 0.96)
    if(k==1) {
      # Add legend
      lhyper <- paste0("Hypersaline (", paste(which(ihyper), collapse = ", "), ")")
      lsed <- paste0("Sediment (", paste(which(ised), collapse = ", "), ")")
      lacid <- paste0("Acidic (", paste(which(iacid), collapse = ", "), ")")
      lalk <- paste0("Neutral/Alkaline (", paste(which(col == orp16Scol[1]), collapse = ", "), ")")
      legend("topleft", c(lhyper, lsed, lacid, lalk), col = c("turquoise3", "gray", orp16Scol[4], orp16Scol[1]), lwd = 2, bty = "n")
    }
  }

  # Make a plot for both domains combined 20211002
  plot(c(-450, 400), c(-0.23, -0.08), type = "n", xlab = "Eh7 (mV)", ylab = cplab$ZC)
  for(j in seq_along(studies)) {
    # Skip study with no archaeal sequences
    if(studies[j] %in% c("PBU+20")) next
    pout <- plotZC(studies[j], show = "lm", lwd = 2, col.line = col[j], add = TRUE)
    text(pout$Eh7lim[2], pout$ZC[2], j)
  }
  title("Bacteria and Archaea", font.main = 1)
  label.figure("C", font = 2, cex = 1.5, yfrac = 0.96)

  if(pdf) dev.off()

}

add.linear <- function(Eh7, ZC, nstudy = NA) {
  # Plot linear fit and confidence interval 20210920
  # Use lighter colors for environments with few datasets 20210925
  text.col <- "black"
  polygon.col <- "gray80"
  line.col <- "gray62"
  if(!is.na(nstudy)) if(nstudy < 5) {
    text.col <- "gray60"
    polygon.col <- "gray90"
    line.col <- "gray70"
  }
  # Show number of studies and exit if there are zero 20210926
  if(!is.na(nstudy)) {
    legend("topleft", legend = nstudy, bty = "n", text.font = 2, text.col = text.col, inset = c(-0.05, 0))
    if(nstudy == 0) return()
  }
  # https://stackoverflow.com/questions/22717930/how-to-get-the-confidence-intervals-for-lowess-fit-using-r
  Ehvals <- seq(range(Eh7)[1], range(Eh7)[2], length.out = 100)
  EZlm <- lm(ZC ~ Eh7)
  plx <- predict(EZlm, newdata = data.frame(Eh7 = Ehvals), se = T)
  upper <- plx$fit + qt(0.975, plx$df) * plx$se
  lower <- plx$fit - qt(0.975, plx$df) * plx$se
  polygon(c(Ehvals, rev(Ehvals)), c(lower, rev(upper)), col = polygon.col, border = NA)
  lines(Ehvals, plx$fit, col = line.col)
  # Get the slope
  slope <- EZlm$coefficients[2]
  # Multiply by 1000 to convert from mV to V
  slope <- slope * 1000
  # Round to fixed number of decimal places
  slope <- formatC(slope, digits = 3, format = "f")
  stext <- bquote(.(slope)~V^-1)
  # Also show number of samples
  Ntext <- bquote(italic(N) == .(length(ZC)))
  # Add text to plot
  legend <- as.expression(c(Ntext, stext))
  legend("topright", legend = legend, bty = "n", text.col = text.col)

#  # Get P-value
#  Pval <- summary(EZlm)$coefficients[, "Pr(>|t|)"]["Eh7"]
#  # Divide by 2 to get one-sided p-value
#  Pval <- Pval / 2
#  Pval <- formatC(Pval, digits = 1, format = "e")
#  Ptext <- bquote(italic(P) == .(Pval))
#  legend("topleft", legend = Ptext, bty = "n")
}

# Scatterplots for all samples in each environment type 20210913
# NOTE: default ienv omits Hot Spring (4)
eachenv <- function(lineage = "Bacteria", add = FALSE, do.linear = TRUE, ienv = c(1, 2, 5, 3, 6, 7)) {
  # Read Eh7 - ZC data
  dat <- read.csv(system.file("extdata/orp16S/EZdat.csv", package = "JMDplots"))
  # Get overall x and y limits
  xlim <- range(dat$Eh7)
  ylim <- range(dat$ZC)
  # Use Bacteria or Archaea only
  dat <- dat[dat$lineage == lineage, ]
  # Get names of environment types
  envirotypes <- names(envirotype)
  # Loop over environment types
  for(i in ienv) {
    # Start plot
    if(!add) plot(xlim, ylim, type = "n", xlab = "", ylab = "", axes = FALSE)
    # Get Eh7 and ZC values
    thisdat <- dat[dat$envirotype == envirotypes[i], ]
    Eh7 <- thisdat$Eh7
    ZC <- thisdat$ZC
    if(do.linear) {
      # Include number of studies in legend 20210925
      nstudy <- length(unique(thisdat$study))
      add.linear(Eh7, ZC, nstudy)
    }
    if(!add) {
      # Add plot axes and title
      if(lineage == "Archaea") {
        axis(1, labels = NA)
        # Make rotated labels (modified from https://www.r-bloggers.com/rotated-axis-labels-in-r-plots/)
        x <- c(-400, 0, 600)
        text(x = x, y = par()$usr[3] - 1.5 * strheight("A"), labels = x, srt = 45, adj = 1, xpd = NA)
      }
      if(i == 1) axis(2)
      box()
      if(lineage == "Bacteria") title(envirotypes[i], font.main = 1, cex.main = 1)
    }
    # Plot points
    points(Eh7, ZC, pch = 19, cex = 0.2, col = orp16Scol[i])
  }
}

# Figure S1: ZC-Eh scatterplots for all studies 20210827
# This also creates files EZdat (Eh and ZC values) and
# EZlm (linear fits) for use by other plotting functions
orp16S_S1 <- function(pdf = FALSE) {

  # Setup figure
  if(pdf) pdf("Figure_S1.pdf", width = 12, height = 9)
  par(mfrow = c(3, 4))

  results <- c(

    message("\nRiver & Seawater"),
    plotZC("MLL+18", "Bacteria", groupby = "Season", groups = c("Summer", "Winter"), legend.x = "bottomright"),
    plotZC("SVH+19", "two", groupby = "Type", groups = c("Oxic", "Suboxic", "Euxinic")),
    plotZC("HXZ+20", "Bacteria", groupby = "Station", groups = c("SYBL", "C4")),
    plotZC("KLY+20", "Bacteria", groupby = "Year", groups = c(2018, 2019), legend.x = "bottomright"),
    plotZC("WHL+21", "two", groupby = "Season", groups = c("Spring", "Summer", "Autumn", "Winter"), legend.x = "bottomleft"),
    plotZC("LXH+20", "Bacteria", groupby = "Season", groups = c("Summer", "Winter"), legend.x = "bottomright"),
    plotZC("JVW+20", "two", groupby = "isolation_source", groups = c("Ulva laetevirens", "lagoon water"), legend.x = "topright"),
    plotZC("ZZL+21", "Bacteria", groupby = "Location", groups = c("Main Stream", "Animal Farm", "Hospital", "WWTP", "Tributary"), legend.x = "bottomright"),
    plotZC("GZL21", "Bacteria", groupby = "Type", groups = c("Surface water", "Middle water", "Bottom water"), legend.x = "bottomleft"),

    message("\nLake & Pond"),
    plotZC("SAR+13", "two", groupby = "Zone", groups = c("Photic-oxic", "Transition", "Anoxic")),
    plotZC("ECS+18", "Bacteria", groupby = "Lake", groups = c("Laguna Negra", "Lo Encanado")),
    plotZC("LLC+19", "Bacteria", groupby = "Size", groups = c("Free-living", "Particle-associated")),
    plotZC("BCA+21", "Bacteria", groupby = "Month", groups = c("Jul", "Nov", "Feb", "Apr")),
    plotZC("HLZ+18", "Bacteria", groupby = "Type", groups = c("Reservoir", "Pond"), legend.x = "bottomright"),
    plotZC("CNA+20", "Bacteria", groupby = "Season", groups = c("Summer", "Autumn", "Winter", "Spring"), legend.x = "topright"),
    plotZC("BWD+19", "two", groupby = "Cover", groups = c("Ice", "Ice Free"), legend.x = "bottomright"),
    plotZC("RSJ+21", "two", groupby = "Lake", groups = c("Kuiva", "Lovo"), legend.x = "topright"),
    plotZC("BOEM21", "Bacteria", groupby = "Stratum", groups = c("Upper", "Chemocline", "Lower")),
    plotZC("FAV+21", "two", groupby = "Type", groups = c("Oxic Surface", "Anoxic Surface", "Bottom")),

    # Hot Spring
    message("\nHot Spring"),
    plotZC("SMS+12", "two", groupby = "Type", groups = c("Chemotrophic", "Transition", "Phototrophic")),
    plotZC("BMJ+19", "two", groupby = "Environment", groups = c("Hydrothermal Pond", "Yellow Lake", "Black Lake", "Assale Lake", "Cave Water")),
    plotZC("LMG+20", "two", groupby = "Setting", groups = c("River", "Hot Spring"), legend.x = "bottomleft"),
    plotZC("GWSS21", "two", groupby = "Location", groups = c("Rehai", "Banglazhang")),
    plotZC("GWS+20", "two", groupby = "Field", groups = c("Batang", "Litang", "Kangding")),
    plotZC("PBU+20", "Bacteria", groupby = "Type", groups = c("Cauldron", "Sampling Pit", "Spring", "Geyser Valley (Control)"), legend.x = "bottomright"),
    plotZC("MWY+21", "two", groupby = "Location", groups = c("Quseyongba", "Moluojiang", "Daggyai", "Quzhuomu"), legend.x = "bottomright"),
    plotZC("OFY+19", "two", groupby = "Type", groups = c("Hot Spring", "Drain", "Lake")),

    message("\nAlkaline Spring"),
    plotZC("SBP+20", "two", groupby = "pH Group", groups = c("< 10", "> 10")),
    plotZC("RMB+17", "two", groupby = "pH Group", groups = c("< 10", "> 10")),
    plotZC("CTS+17", "two", groupby = "Type", groups = c("River", "Well", "Spring"), legend.x = "bottomleft"),
    plotZC("KSR+21", "Bacteria", groupby = "Location", groups = c("Lerone", "Branega", "Branega Creek Water")),
    plotZC("NTB+21", "two", groupby = "Well", groups = c("BA1A", "BA1D")),

    message("\nGroundwater"),
    plotZC("KLM+16", "Bacteria", groupby = "Day", groups = c(-1, 246, 448, 671)),
    plotZC("SDH+19", "two", groupby = "Type", groups = c("Freshwater", "Brackish", "Saltwater")),
    plotZC("SRM+19", "Bacteria", groupby = "Land Use", groups = c("Agriculture", "Community", "Landfill", "Mine")),
    plotZC("SKP+21", "Bacteria", groupby = "Type", groups = c("Groundwater", "Surface water"), legend.x = "topright"),
    plotZC("YHK+20", "Bacteria", groupby = "Location", groups = c("Upper Hillslope", "Middle Slope", "Lower Footslope")),
    plotZC("JDP+20", "Bacteria", groupby = "Roll-Front Setting", groups = c("Oxidized", "Intermediate", "Reduced"), legend.x = "bottomright"),
    plotZC("GWS+19", "Bacteria", groupby = "Type", groups = c("Un-confined", "Confined"), legend.x = "bottomleft"),
    plotZC("SRM+21", "Bacteria", groupby = "Depth", groups = c("Surface", "Shallow", "Deep"), legend.x = "bottomleft"),
    plotZC("ZCZ+21", "Bacteria", groupby = "Location", groups = c("LO", "CR1", "MN", "VA", "BS", "CR2"), legend.x = "topright"),

    message("\nSediment"),
    plotZC("JHL+12", "two", groupby = "Core", groups = c("GC6", "GC12")),
    plotZC("GFE+16", "Bacteria", groupby = "Station", groups = c(1, 4, 7, 18)),
#    plotZC("ZML+17", "two", groupby = "Season", groups = c("Winter", "Summer"), legend.x = "bottomleft"),
    plotZC("ZML+17", "two", groupby = "Type", groups = c("Mangrove Forest", "Intertidal Mudflat"), legend.x = "bottomleft"),
    plotZC("BYB+17", "two", groupby = "Type", groups = c("Freshwater", "Brackish", "Saltmarsh", "Hypersaline"), legend.x = "bottomright"),
    plotZC("BSPD17", "Bacteria", groupby = "Type", groups = c("Water", "Sediment"), legend.x = "bottomleft"),
    plotZC("ABT+17", "two", groupby = "Type", groups = c("Coastal", "Estuarine")),
    plotZC("HDZ+19", "Bacteria", groupby = "Type", groups = c("Water", "Sediment")),
    plotZC("SCM+18", "two", groupby = "Site", groups = c("Shallow", "Deep"), legend.x = "bottomleft"),
    plotZC("CLS+19", "two", groupby = "Type", groups = c("Water", "Sediment"), legend.x = "bottomright"),
    plotZC("ZDA+20", "Bacteria", groupby = "Type", groups = c("Water", "Sediment"), legend.x = "bottomleft"),
    plotZC("VMB+19", "two", groupby = "Type", groups = c("Water", "Sediment"), legend.x = "topright"),
    plotZC("HSF+19", "Bacteria", groupby = "Type", groups = c("Water", "Sediment-Water Interface", "Sediment"), legend.x = "topright"),
    plotZC("MCS+21", "two", groupby = "Type", groups = c("Water", "Sediment"), legend.x = "left"),
    plotZC("LMBA21", "two", groupby = "Year", groups = c(2013, 2015, 2016, 2017, 2018, 2019), legend.x = "bottomright"),
    plotZC("ZZLL21", "Bacteria", groupby = "Type", groups = c("Main stream", "Animal farm", "Hospital", "WWTP", "Tributary"), legend.x = "bottomleft"),
    plotZC("WFB+21", "Bacteria", groupby = "Treatment", groups = c("C. volutator", "H. diversicolor", "Cv & Hd", "MPB", "Manual turbation")),

    message("\nSoil"),
    plotZC("SBW+17", "Bacteria", groupby = "Treatment", groups = c("Control", "BC400", "BC600")),
    plotZC("MLL+19", "two", groupby = "Type", groups = c("Upland", "Paddy", "Sediment")),
    plotZC("BMOB18", "two", groupby = "Treatment", groups = c("Acetate", "No amendment", "Pre-incubation")),
    plotZC("ZZZ+18", "two", groupby = "Treatment", groups = c("None", "AQDS", "Biochar"), legend.x = "bottomleft"),
    plotZC("PMM+20", "two", groupby = "Type", groups = c("Fluctuating", "Static")),
    plotZC("ZHZ+19", "two", groupby = "Treatment", groups = c("Original", "Nitrate-reducing", "Ferric-reducing", "Sulfate-reducing", "Methanogenic")),
    plotZC("CWC+20", "Bacteria", groupby = "Management", groups = c("Flooding", "Draining")),
    plotZC("PSG+20", "two", groupby = "Treatment", groups = c("Initial", "NCC", "RB", "RGP", "TP")),
    plotZC("XLD+20", "two", groupby = "Treatment", groups = c("CK", "9K-Control", "Htt-sys", "Att-sys", "Co-sys"), legend.x = "bottomright"),
    plotZC("LJC+20", "two", groupby = "meanT", groups = c("MAT >= 21.5 degC", "MAT < 21.5 degC")),
    plotZC("DTJ+20", "two", groupby = "Zone", groups = c("Bulk Soil", "Mature", "Elongation", "Tip")),
    plotZC("LLL+21", "Bacteria", groupby = "Treatment", groups = c("CK", "FL", "EA", "SB", "BD")),
    plotZC("DLS21", "two", groupby = "Treatment", groups = c("control", "char", "silicate", "husk"))

  )
  # Done plotting!
  if(pdf) dev.off()

  # Assemble all data
  EZdat <- do.call(rbind, results[names(results) == "EZdat"])
  # Assemble all linear model results
  coefficients <- sapply(results[names(results) == "EZlm"], "[", "coefficients")
  intercept <- sapply(coefficients, "[[", "(Intercept)")
  slope <- sapply(coefficients, "[[", "Eh7")
  model <- lapply(results[names(results) == "EZlm"], "[[", "model")
  Eh7lim <- sapply(lapply(model, "[", "Eh7"), "range")
  Eh7min <- Eh7lim[1, ]
  Eh7max <- Eh7lim[2, ]
  study <- sapply(results[names(results) == "study"], "[")
  envirotype <- sapply(results[names(results) == "envirotype"], "[")
  lineage <- sapply(results[names(results) == "lineage"], "[")
  nsamp <- sapply(model, nrow)
  EZlm <- data.frame(study, envirotype, lineage, nsamp, Eh7min, Eh7max, slope, intercept)
  # Save data and results to files
  write.csv(EZdat, "EZdat.csv", row.names = FALSE, quote = FALSE)
  write.csv(EZlm, "EZlm.csv", row.names = FALSE, quote = FALSE)

}


# Plot ZC values vs Eh7 for a single study 20210827
# Use 'groupby' (name of column with sample groups) and 'groups' (names of sample groups) to apply the pch and col to individual samples
plotZC <- function(study, lineage = NULL, mincount = 100, pch = NULL, col = NULL, add = FALSE, type = "p", groupby = NULL, groups = NULL,
                   legend.x = "topleft", show = c("lm", "points"), col.line = "gray62", lwd = 1, cex = 1, mdat = NULL, title.line = NA) {

  if(identical(lineage, "two")) {
    # Make two plots for studies that have Bacteria and Archaea 20210913
    out1 <- plotZC(study, "Bacteria", mincount, pch, col, add, type, groupby, groups, legend.x, show, col.line, lwd, cex)
    # Don't show legend on second (Archaea) plot 20210914
    out2 <- plotZC(study, "Archaea", mincount, pch, col, add, type, groupby, groups, legend.x = NA, show, col.line, lwd, cex, mdat = out1$mdat)
    out <- c(out1, out2)
    return(invisible(out))
  }

  # Get metadata; use suppressMessages() to suppress messages
  if(is.null(mdat)) mdat <- suppressMessages(getmdat(study))
  mdat.orig <- mdat
  # For Bacteria or Archaea, use only runs that are labeled as such 20210920
  if("Domain" %in% colnames(mdat)) {
    if(identical(lineage, "Bacteria")) mdat <- mdat[mdat$Domain == "Bacteria", ]
    if(identical(lineage, "Archaea")) mdat <- mdat[mdat$Domain == "Archaea", ]
  }

  # Use capture.output to hide printed output
  null <- capture.output(
    # Use try() to capture errors (with no mapped sequences for lineage = "Archaea")
    met <- try(
      suppressMessages(
        getmetrics(study, mdat = mdat, lineage = lineage, mincount = mincount)
      ), silent = TRUE
    )
  )
  # Print message and skip dataset with no mapped sequences
  if(inherits(met, "try-error")) {
    print(paste0(study, ": no mapped sequences for ", lineage))
    return()
  }

  # Keep metadata only for samples with >= mincount counts 20201006
  mdat <- mdat[mdat$Run %in% met$Run, ]
  nsamp <- nrow(mdat)
  # Remove samples with NA Eh7 or ZC 20210822
  mdat <- mdat[!(is.na(mdat$Eh7) | is.na(met$ZC)), ]
  met <- met[met$Run %in% mdat$Run, ]
  stopifnot(all(mdat$Run == met$Run))
  # Print message about number of samples and Eh7 and ZC range
  ZCtext <- paste(range(round(met$ZC, 3)), collapse = " to ")
  if(!is.null(lineage)) ltext <- paste0(lineage, ": ") else ltext <- ""
  print(paste0(study, ": ", ltext, nrow(mdat), "/", nsamp, " samples, ZC ", ZCtext))

  # Assign pch and col to sample groups
  if(!is.null(groupby) & !is.null(groups)) {

    # Get default point symbols
    pchavail <- 21:25
    if(is.null(pch)) pch <- rep(pchavail, length(groups))

    # The pch and col for each sample type
    pchtype <- rep(pch, length.out = length(groups))
    coltype <- 1:length(groups)

    # The pch and col for individual samples
    pch <- col <- rep(NA, nrow(mdat))
    # The column with sample groups
    icol <- match(groupby, colnames(mdat))
    if(is.na(icol)) stop(paste(groupby, "is not a column name in metadata for", study))
    # Loop over sample groups
    for(i in seq_along(groups)) {
      # Find matching samples and set the pch and col
      itype <- mdat[, icol] == groups[i]
      pch[itype] <- pchtype[i]
      col[itype] <- coltype[i]
    }

  }

  # Defaults for pch and col if sample groups are not specified
  if(is.null(pch)) pch <- 19
  if(is.null(col)) col <- "#40404080"

  # Make data frame with Eh7 and ZC values
  EZdat <- data.frame(Eh = mdat$Ehorig, Eh7 = mdat$Eh7, ZC = round(met$ZC, 6))
  # Create subtitle for environment type 20210904
  sub <- envirotype <- envirodat$group[envirodat$study == study]
  if(!add) {
    # Start new plot
    plot(EZdat$Eh7, EZdat$ZC, xlab = "Eh7 (mV)", ylab = cplab$ZC, type = "n")
    # Take off suffix after underscore 20210914
    root <- strsplit(study, "_")[[1]][1]
    suffix <- strsplit(study, "_")[[1]][2]
    main <- paste0(na.omit(mdat$name)[1], " (", root, ")")
    title(main = main, font.main = 1, line = title.line)
    # Include suffix in subtite 20210914
    if(!is.na(suffix)) sub <- paste(sub, "-", suffix)
    # Add lineage 20210913
    if(!is.null(lineage)) sub <- paste(sub, "-", lineage)
    title(main = sub, line = 0.5, cex.main = 1)
  }
  # Add linear fit
  if("lm" %in% show) {
    EZlm <- lm(ZC ~ Eh7, EZdat)
    Eh7lim <- range(EZlm$model$Eh7)
    ZCpred <- predict.lm(EZlm, data.frame(Eh7 = Eh7lim))
    # Use solid or dashed line to indicate large or small slope 20210926
    slope <- EZlm$coefficients[2] * 1000
    if(is.na(slope)) lty <- 3 else if(abs(slope) < 0.01) lty <- 2 else lty <- 1
    lines(Eh7lim, ZCpred, col = col.line, lwd = lwd, lty = lty)
  }
  # Add points
  if("points" %in% show) {
    points(EZdat$Eh7, EZdat$ZC, pch = pch, col = orp16Scol[col], bg = orp16Scol[col], type = type, cex = cex)
  }

  if(!is.null(groupby) & !is.null(groups)) {
    # Add legend
    legend <- as.character(groups)
    legend(legend.x, legend, pch = pchtype, col = orp16Scol[coltype], pt.bg = coltype, title = groupby, cex = 0.9)
    # Add sample type (group) to output
    EZdat <- cbind(groupby = groupby, group = mdat[, icol], EZdat)
  } else {
    EZdat <- cbind(groupby = NA, group = NA, EZdat)
  }


  # Return values
  # Use first column name starting with "sample" or "Sample" 20210818
  sampcol <- grep("^sample", colnames(mdat), ignore.case = TRUE)[1]
  if(is.null(lineage)) lineage <- ""
  EZdat <- cbind(study = study, envirotype = envirotype, lineage = lineage, sample = mdat[, sampcol], Run = mdat$Run, EZdat)
  out <- list(study = study, envirotype = envirotype, lineage = lineage, mdat = mdat.orig, EZdat = EZdat)
  if("lm" %in% show) out <- c(out, list(EZlm = EZlm, Eh7lim = Eh7lim, ZCpred = ZCpred))
  invisible(out)

}

