# JMDplots/orp16S.R
# Plots for paper on ZC-ORP correlations 20210827

# Group studies by environment types 20210828
envirotype <- list(
  "River & Seawater" = c("MLL+18", "SVH+19", "HXZ+20", "KLY+20", "GSBT20_Prefilter", "GSBT20_Postfilter", "WHL+21", "LXH+20", "JVW+20", "ZZL+21",
                         "GZL21"),
  "Lake & Pond" = c("SAR+13", "LZR+17", "ECS+18", "LLC+19", "SCH+16", "BCA+21", "HLZ+18", "GRG+20", "CNA+20", "BWD+19",
                    "RSJ+21", "LRL+22", "BOEM21", "IBK+22", "GSY+20", "NLE+21", "SPA+21", "FAV+21", "PSV+22"),
  "Groundwater" = c("KLM+16", "YHK+19", "SDH+19", "SRM+19", "APV+20", "SKP+21", "YHK+20", "JDP+20", "GWS+19",
                    "SRM+21", "ZCZ+21", "CSW+22", "GXS+22"),
  # NOTE: Keep Hot Spring in 4th location to get red color 20210904
  "Hot Spring" = c("SMS+12", "PCL+18_Acidic", "PCL+18_Alkaline", "BMJ+19", "LMG+20", "GWSS21", "GWS+20", "PBU+20", "MWY+21"),
  "Hyperalkaline" = c("SBP+20", "RMB+17", "CTS+17", "SPH+21", "KSR+21", "PSB+21", "NTB+21"),
  "Sediment" = c("JHL+12", "GFE+16", "ZML+17", "BYB+17", "BSPD17", "HDZ+19", "TCN+17", "WHLH21", "SCM+18", "RSS+18",
                 "CLS+19", "ZDA+20", "VMB+19", "WHC+19", "HSF+19", "RBM+21", "ZHZ+19", "MCS+21", "LMBA21_2017", "HSF+22",
                 "ZZLL21", "BKR+22", "WFB+21", "HCW+22"),
  "Soil" = c("SBW+17", "MLL+19", "ZLH+22", "BMOB18", "ZZZ+18", "PMM+20", "WHLH21a", "CWC+20", "PSG+20", "XLD+20",
             "LJC+20", "DTJ+20", "ZWH+22", "LLL+21", "RKSK22", "DLS21_Bulk", "WKP+22", "CYG+22", "CKB+22")
)
# Turn the list into a data frame for easier lookup 20210904
envirodat <- do.call(rbind, lapply(seq_along(envirotype), function(i) data.frame(study = envirotype[[i]], groupnum = i)))
envirodat <- cbind(envirodat, group = names(envirotype)[envirodat$groupnum])

# Set the palette for R colors (numeric index of 'col') 20210914
orp16Scol <- palette.colors(n = length(envirotype), palette = "Classic Tableau", alpha = 0.75)

# Read Eh7 - ZC data
EZdat <- read.csv(system.file("extdata/orp16S/EZdat.csv", package = "JMDplots"))
# Read linear fit coefficients
EZlm <- read.csv(system.file("extdata/orp16S/EZlm.csv", package = "JMDplots"))

# Figure 1: Thermodynamic predictive framework 20210830
orp16S1 <- function(pdf = FALSE) {

  if(pdf) pdf("Figure1.pdf", width = 6, height = 4)

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
  par(mar = c(0, 0, 0, 0))
  plot(c(5, 95), c(1, 80), type = "n", axes = FALSE, xlab = "", ylab = "")
  # Uncomment this as a guide for making the 'grid' one 20210927
  #box(lwd = 2)
  grid.roundrect(0.5, 0.5, 0.99, 0.99, gp = gpar(fill = "azure"))

  dy <- -2.5

  # Plot lines to go behind shapes
  lines(c(20, 50), c(20, 60)+dy)
  lines(c(20, 50), c(60, 60)+dy)
  lines(c(50, 80), c(60, 60)+dy)
  # Add arrows along lines 20210927
  arrows(20, 20+dy, 20*0.35 + 50*0.65, (20*0.35 + 60*0.65)+dy, length = 0.1)
  arrows(20, 60+dy, 37, 60+dy, length = 0.1)
  arrows(50, 60+dy, 66, 60+dy, length = 0.1)

  # Plot shapes and text for biological methods
  text(20, 79+dy, "Biological Methods", col = RedText, font = 2)
  for(bg in c("white", Red80)) points(20, 60+dy, pch = 21, cex = 18, bg = bg)
  text(20, 66+dy, "RefSeq", font = 2, col = RedText)
  text(20, 57+dy, "Reference\nproteomes\nof taxa")
  for(bg in c("white", Red80)) points(20, 20+dy, pch = 21, cex = 18, bg = bg)
  text(20, 25+dy, "16S + RDP", font = 2, col = RedText)
  text(20, 18+dy, "Taxonomic\nabundances")

  # Plot shapes and text for chemical methods
  text(80, 79+dy, "Chemical Methods", col = BlueText, font = 2)
  for(bg in c("white", Blue80)) points(80, 60+dy, pch = 22, cex = 18, bg = bg)
  text(80, 68+dy, quote(bolditalic(Z)[bold(C)]), col = BlueText, cex = 1.2)
  text(80, 58.5+dy, "Carbon\noxidation\nstate")
  # Show multiple physicochemical variables 20210927
  # Function to draw rectangle at x,y with width and height w,h
  myrect <- function(x, y, w, h, ...) rect(x - w/2, y - h/2, x + w/2, y + h/2, ...)
  # T, Eh, pH, O2
  for(col in c("white", Blue80)) myrect(73, 16+dy, 5, 6, col = col)
  text(73, 16+dy, "T", col = BlueText)
  for(col in c("white", Blue80)) myrect(80, 20+dy, 7, 7, col = col)
  text(80, 20+dy, "Eh", font = 2, cex = 1.2, col = BlueText)
  for(col in c("white", Blue80)) myrect(87.5, 16+dy, 6, 6, col = col)
  text(87.5, 16+dy, "pH", col = BlueText)
  for(col in c("white", Blue80)) myrect(80, 11+dy, 5, 6, col = col)
  text(80, 11+dy, quote(O[2]), col = BlueText)

  # Plot inner rectangle and text
  third <- 100/3
  # Uncomment this to make the original rectangle (as a guide for the 'grid' one)
  #for(bg in c("white", Orange80)) rect(1.2*third, 50+dy, 1.8*third, 70+dy, col = bg)
  # Use this to get rounded corners 20210927
  for(fill in c("white", Orange80)) grid.roundrect(0.5, 0.7, 0.205, 0.23, gp = gpar(fill = fill))
  text(50, 67+dy, quote(bold(C[bolditalic(c)]*H[bolditalic(h)]*N[bolditalic(n)]*O[bolditalic(o)]*S[bolditalic(s)])), col = OrangeText)
  text(50, 58+dy, "Estimated\nCommunity\nProteomes")

  # Plot arrows and text labels
  arrows(1*third - 1, 20+dy, 1*third + 6, 20+dy, code = 1, lty = 1, length = 0.1)
  arrows(2*third - 4, 16+dy, 2*third + 3, 16+dy, code = 2, lty = 1, length = 0.1)
  text(50, 18+dy, "Community and\nenvironmental data")
  arrows(80, 25+dy, 80, 45+dy, code = 3, lwd = 1.5, length = 0.1, col = BlueText)
  text(46, 42+dy, "Thermodynamic prediction", font = 2, cex = 0.9, adj = c(0, 1))
  text(49, 42+dy, "\nCarbon oxidation state\nis positively correlated\nwith redox potential", font = 3, cex = 0.9, adj = c(0, 1))

  if(pdf) dev.off()

}

# Figure 2: Chemical and geobiochemical depth profiles in Winogradsky columns 20210829
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
  metrics.in <- getmetrics_orp16S("RBW+14")
  mdat <- getmdat_orp16S("RBW+14", metrics.in)
  metrics <- mdat$metrics
  # Get ZC values for each layer
  layers <- c("12 cm", "8 cm", "4 cm", "SWI", "Top")
  ZC <- lapply(layers, function(layer) metrics$ZC[mdat$metadata$layer == layer])
  # Make boxplots
  boxplot(ZC, horizontal = TRUE, show.names = FALSE, xlab = axis.label("ZC"), ylim = c(-0.18, -0.145), yaxs = "i", col = "azure3")
  axis(2, 1:5, labels = layers, las = 1)
  # Add sample sizes
  par(xpd = NA)
  for(i in 1:5) {
    label <- bquote(italic(N) == .(length(ZC[[i]])))
    text(-0.182, i - 0.25, label, adj = 1)
  }
  par(xpd = FALSE)
  label.figure("B", font = 2, cex = 1.5)

  # Reset layout to make orp16S3 in the examples run nicely 20211011
  if(pdf) dev.off() else layout(1)

}

# Figure 3: Sample locations on world map
orp16S3 <- function(pdf = FALSE) {

  if(pdf) pdf("Figure3.pdf", width = 26, height = 15)

  # Coordinates for orp16S datasets
  file <- tempfile()
  # Write spaces here (but don't save them in the file) to make this easier to read
  writeLines(con = file, text = gsub(" ", "", c(
    # Column names
    "study, latitude, longitude",
    ## River & Seawater (comments indicate source of coordinates from paper or SRA metadata)
    "MLL+18, 22.20, 113.09", # Fig. 1
    "SVH+19, 42.90, 30.68", # Materials and methods
    "HXZ+20, 16.52, 111.77", # Table 1
    "KLY+20, 37.60189, 126.8128", # Materials and methods
    # GSBT20 see below
    "WHL+21, 30.76, 115.36", # SAMN13242327
    "LXH+20, 25.42, 99.34", # SAMN15090995
    "JVW+20, 45.432025, 12.260878", # SAMN15796698
    "ZZL+21, 22.77, 113.79", # SAMN16964962
    # GZL21 see below
    ## Lake & Pond
    "SAR+13, 1.96, -157.33", # Materials and methods
    "LZR+17, 30.587, 104.310", # Table 1
    "ECS+18, -33.65, -70.117", # Materials and methods
    "LLC+19, 24.82, 118.15", # SAMN04549101
    "SCH+16, 67.081, -50.355", # Materials and methods
    "BCA+21, 46.3615, 25.0509", # SAMN07409474
    "HLZ+18, 24.795, 118.138", # SAMN07638080
    "GRG+20, 37.726, -6.555", # Materials and methods
    "CNA+20, 39.441, -77.371", # Web search for geographic center of Maryland --> https://sos.maryland.gov/mdkids/Pages/Geography.aspx
    "BWD+19, 47.120571, -88.545425", # SAMN09980099
    "RSJ+21, 61.833, 24.283", # Materials and methods
    "LRL+22, 20.8538, -87.1253", # SAMN16910034 
    "BOEM21, 43.051389, -75.965", # Materials and methods
    "IBK+22, 53.1516, 13.0262", # SAMN15366194
    "GSY+20, 37.59, -7.124", # Materials and methods
    "NLE+21, 32.833, 35.583", # SAMEA7280991
    "SPA+21, 45.8126, 8.7401", # SAMN17524543
    "FAV+21, 0.757, 36.372", # SAMN19267646
    "PSV+22, 50.178, 12.596", # Study site
    ## Hot Spring
    "SMS+12, 44.6, -110.9", # JGI IMG/M sample name 1_050719N
    "PCL+18_Acidic, -38.5, 176.0", # Fig. 1
    "BMJ+19, 14.089567, 40.348583", # SAMN11581539
    "LMG+20, -37.855, -71.158", # Table 1
    "GWSS21, 24.86, 98.33", # SAMN16802401
    "GWS+20, 30.12, 101.94", # SAMN13430433
    "PBU+20, 54.4395, 160.144194", # SAMN14538724
    # MWY+21 see below
    ## Hyperalkaline
    "SBP+20, 38.862, -122.414", # SAMN03850954
    "RMB+17, 22.9052, 58.6606", # SAMN05981641
    "CTS+17, 10.94323, -85.63485", # SAMN06226041
    "SPH+21, 38.8257, -122.35152", # SAMN16578990
    "KSR+21, 44.264340, 8.46442", # SAMN17101425
    "PSB+21, 38.8621, -122.4304", # SAMN17252996
    "NTB+21, 22.881, 58.701", # SAMN19998441
    ## Soil - put this group before Groundwater and Sediment for clearer visualization in GBA 20210927
    "SBW+17, 28.25, 116.92", # Materials and methods  ### Laboratory
    "MLL+19, 26.1, 112.5", # Materials and methods
    "ZLH+22, 30.50, 118.36", # SAMN07816355   ### Laboratory
    "BMOB18, 40.60842, -74.19258", # SAMN07828017  ### Laboratory
    "ZZZ+18, 21.816, 112.464", # Materials and methods  ### Laboratory
    "PMM+20, 43.397, -80.311",  # Web search for Cambridge, ON, Canada  ### Laboratory
    "WHLH21a, 37.53, 105.03", # Materials and methods
    "CWC+20, 28.226, 116.898", # Materials and methods  ### Laboratory
    "PSG+20, 36.61, -119.53", # Web search for Parlier, CA   ### Mesocosm
    "XLD+20, 27.35, 112.05", # Materials and methods   ### Laboratory
    # LJC+20 see below
    "DTJ+20, 26.45, 111.52", # SAMN14332759   ### Laboratory
    "ZWH+22, 29.18, 119.28", # SAMN16191137   ### Laboratory
    "LLL+21, 27.78, 113.13", # Materials and methods   ### Laboratory
    "RKSK22, 31.97283, -81.0304", # SAMN16678415
    "DLS21_Bulk, 39.39, -75.44", # SAMN17245435  ### Mesocosm
    "WKP+22, 53.29, 17.79", # SAMN23457780
    "CYG+22, 31.035, 118.8367", # Wikipedia Nanjing Agricultural University   ### Laboratory
    "CKB+22, 37.8635, 138.9426", # Wikipedia Niigata University   ### Laboratory
    ## Groundwater
    "KLM+16, 42.99, -82.30", # SAMN04423023
    # YHK+19 see below
    "SDH+19, 23.03, 113.38", # SAMN07692244
    "SRM+19, 12.67417, 101.3889", # Materials and methods
    "APV+20, 20.12, -99.23", # Materials and methods
    "SKP+21, 16.263306, 100.647778", # SAMN11191517
    "YHK+20, 51.209467, 10.791968", # SAMEA5714424
    "JDP+20, 44.8883, 110.1353", # SAMN12236980
    "GWS+19, 36.31, 94.81", # SAMN07765433
    "SRM+21, 14.83, 99.35", # SAMN14829351
    "ZCZ+21, 45.21, 9.57", # Table 1
    "CSW+22, 36.69, 109.86", # SAMN21040055
    "GXS+22, 36, 114", # Fig. 1
    ## Sediment
    "JHL+12, 73.566167, 8.1585", # SAMN00744894
    "GFE+16, -36.69, -73.07", # Table 1
    "ZML+17, 22.494, 114.029", # Table 1
    "BYB+17, 41.423, -112.087", # Table 1
    "BSPD17, 57.89297, 16.5855", # SAMN05163191   ### Laboratory
    "HDZ+19, 29.901, 113.52435", # SAMN05990289
    "TCN+17, 48.553, -4.536", # Materials and methods   ### Laboratory
    "WHLH21, 23.52, 113.495", # Materials and methods
    "SCM+18, 25.25, -97.23", # Table 1
    "RSS+18, 82, -71", # Introduction
    "CLS+19, 32.22, 118.83", # SAMN08683376
    "ZDA+20, -21.42996, -70.05874", # SAMEA4858706
    "VMB+19, 52.11, 79.17", # SAMN08987150
    "WHC+19, 30.12, 122.14", # Materials and methods
    "HSF+19, 47.803, 16.709", # methods
    "RBM+21, -23.25, -44.62", # SAMN10935837
    "ZHZ+19, 23.130, 113.671", # Materials and methods   ### Laboratory
    "MCS+21, -32.15, -71.1", # Materials and methods
    "LMBA17_2017, 43.42, -2.7", # Fig. 1
    "HSF+22, -9.42979, 46.49524", # SAMN14343437
    "ZZLL21, 22.68, 113.97", # SAMN16964887
    # BKR+22 see below
    "WFB+21, 56.440893, -2.863194", # methods   ### Mesocosm
    "HCW+22, 31.3, 119.98", # Materials and methods   ### Laboratory
    "NA, NA, NA"
  )))

  # This reads the data
  coords <- read.csv(file, as.is = TRUE)

  # Start map
  # https://cran.r-project.org/web/packages/oce/vignettes/map_projections.html
  par(mar = c(2, 0.5, 0, 0.5))
  # We don't need data(coastlineWorld) ... it's the default map 20211003
  # A color between azure2 and azure3 20220516
  azure23 <- rgb(t((col2rgb("azure2") + col2rgb("azure3")) / 2), maxColorValue = 255)
  mapPlot(col = azure23, projection = "+proj=wintri", border = "white", drawBox = FALSE)

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

  # Get colors for studies
  icol <- envirodat$groupnum[match(coords$study, envirodat$study)]
  # Identify studies that use samples from laboratory or mesocosm experiments
  lab <- c(
    # Sediment
    "BSPD17", "TCN+17", "ZHZ+19", "WFB+21", "HCW+22",
    # Soil
    "SBW+17", "ZLH+22", "BMOB18", "ZZZ+18", "PMM+20", "CWC+20", "PSG+20", "XLD+20", "DTJ+20",
    "ZWH+22", "LLL+21", "DLS21_Bulk", "CYG+22", "CKB+22"
  )
  pch <- ifelse(coords$study %in% lab, 15, 19)
  # Use smaller points for high-density regions 20210915
  cex <- ifelse(coords$study %in% c(
    "MLL+19", "XLD+20", "LLL+21", "DTJ+20", # Hunan
    "ZZL+21", "MLL+18", "SDH+19", "ZML+17", "ZZLL21", "ZZZ+18", "ZHZ+19", "WHLH21", # GD-HK-MO GBA
    "ZLH+22", "CYG+22", "CLS+19", "ZWH+22", "HCW+22", "WHC+19", "SBW+17" # Jiangsu-Gansu-Zhejiang
  ), 1.5, 2.5)
  # Plot sample locations
  mapPoints(coords$longitude, coords$latitude, pch = pch, col = orp16Scol[icol], cex = cex)

  # Plot transects 20210929
  # Coordinates for East Asia Paddy Soil dataset are from
  # Sourcedata.xlsx from https://doi.org/10.6084/m9.figshare.12622829
  dat <- getmdat_orp16S("LJC+20")
  mapPoints(dat$longitude, dat$latitude, col = orp16Scol[7], lwd = 2)
  # Coordinates for Southern Tibetan Plateau dataset are from Table 1 of MWY+21
  dat <- getmdat_orp16S("MWY+21")
  latlon <- paste(dat$latitude, dat$longitude)
  isuniq <- !duplicated(latlon)
  dat <- dat[isuniq, ]
  mapPoints(dat$longitude, dat$latitude, col = orp16Scol[4], lwd = 2)
  # Coordinates for Three Gorges Reservoir are from Table S1 of GZL21
  dat <- getmdat_orp16S("GZL21")
  dat <- dat[!is.na(dat$"ORP (mV)"), ]
  latlon <- paste(dat$Latitude, dat$Longitude)
  isuniq <- !duplicated(latlon)
  dat <- dat[isuniq, ]
  mapPoints(dat$Longitude, dat$Latitude, col = orp16Scol[1], lwd = 1)
  # Coordinates for Mahomet Aquifer are from Table S1 of YHK+19
  dat <- getmdat_orp16S("YHK+19")
  latlon <- paste(dat$Latitude, dat$Longitude)
  mapPoints(dat$Longitude, dat$Latitude, col = orp16Scol[3], lwd = 1)
  # Coordinates for Port Microbes are from https://github.com/rghannam/portmicrobes/data/metadata/pm_metadata.csv
  dat <- getmdat_orp16S("GSBT20")
  # Use first sample for each port
  dat <- dat[!duplicated(dat$Port), ]
  mapPoints(dat$Longitude, dat$Latitude, col = orp16Scol[1], pch = 17)
  # Coordinates for Barents Sea are from BioSamples of BioProject PRJNA624280
  dat <- getmdat_orp16S("BKR+22")
  dat <- dat[!duplicated(dat$Station), ]
  mapPoints(dat$Longitude, dat$Latitude, col = orp16Scol[6], pch = 19, cex = 1.5)

  # Add legend
  ienv = c(1, 2, 4, 5, 3, 6, 7)
  par(xpd = NA)
  legend("bottomleft", names(envirotype)[ienv], pch = 19, col = orp16Scol[ienv], bty = "n", cex = 2, inset = c(0, -0.03))
  ltext <- c("Field sites", "Laboratory or mesocosm", "Smaller symbols for", "densely sampled areas", "Port sites", "Open symbols for transects")
  legend("bottomright", ltext, pch = c(19, 15, 20, NA, 17, 1), col = c(1, 1, 1, NA, orp16Scol[1], 1),
         bty = "n", cex = 2, pt.cex = c(2, 2, 2, 2, 1, 1), inset = c(0, -0.03))
  par(xpd = FALSE)

  if(pdf) dev.off()

}

# Linear regressions between ZC and Eh7 at local scales 20220517
orp16S4 <- function(pdf = FALSE) {
  if(pdf) pdf("Figure4.pdf", width = 8, height = 6)
  mat <- matrix(c(1,1,1,1, 2,2,2,2, 3,3,3,3,
                  0,0,0,0,0,0,0,0,0,0,0,0,
                  4,4,4,4,4, 0, 6,6,6, 7,7,7,
                  5,5,5,5,5, 0, 6,6,6, 7,7,7),
                nrow = 4, byrow = TRUE)
  layout(mat, heights = c(3, 0.75, 3, 3))
  par(mar = c(3, 4, 3, 1))
  par(mgp = c(2.5, 1, 0))

  ## Panel A: Analysis of selected datasets 20211003
  # Daya Bay (Sediment)
  plotEZ("WHLH21", "Bacteria", groupby = "Position", groups = c("Surface", "Middle", "Bottom"),
    legend.x = "bottomleft", title.line = NULL, dxlim = c(-20, 0))
  title("Daya Bay\n(Sediment bacteria)", font.main = 1)
  label.figure("A", font = 2, cex = 1.5, yfrac = 0.9)
  # Bay of Biscay (Sediment)
  plotEZ("LMBA21_2017", "Bacteria", groupby = "Season", groups = c("Summer", "Winter"),
    legend.x = "bottomright", title.line = NULL, dxlim = c(0, 150))
  title("Bay of Biscay\n(Sediment bacteria)", font.main = 1)
  # Hunan Soil (Soil)
  plotEZ("MLL+19", "Bacteria", groupby = "Type", groups = c("Upland", "Paddy", "Sediment"),
    title.line = NULL, dylim = c(0, 0.005))
  title("Hunan Province\n(Soil and sediment bacteria)", font.main = 1)

  ## Panel B: Slope vs log10(number of samples) for all datasets
  par(mar = c(3.6, 4, 1, 1))
  # Use Bacteria only 20210913
  lmbac <- EZlm[EZlm$lineage == "Bacteria", ]
  # Get color according to environment group
  env <- envirodat[match(lmbac$study, envirodat$study), ]
  # Get range for included samples 20210905
  i1 <- c(1, 2, 4, 5)
  j1 <- env$groupnum %in% i1
  i2 <- c(3, 6, 7)
  j2 <- env$groupnum %in% i2
  xlim <- range(log10(lmbac$nsamp[j1 | j2]))
  # Multiply by 1e3 to use V-1 instead of mV-1 20210913
  # NOTE: conversion to V-1 is moved to orp16S_S2() 20220520
  ymaxabs <- max(abs(lmbac$slope[j1 | j2]))
  ylim <- c(-ymaxabs*1.2, ymaxabs)
  # River & seawater, lake & pond, hot spring, alkaline spring
  plot(xlim, ylim, type = "n", xlab = "log10(Number of samples)", ylab = quote("Slope of linear fit"~(V^-1)))
  abline(h = 0, lty = 2, lwd = 1.5, col = "gray50")
  points(log10(lmbac$nsamp[j1]), lmbac$slope[j1], pch = 19, col = orp16Scol[env$groupnum[j1]])
  # Add legend
  ltext <- names(envirotype)[i1[1:2]]
  legend("bottomleft", ltext, pch = 19, col = orp16Scol[i1[1:2]])
  ltext <- names(envirotype)[i1[3:4]]
  legend("bottomright", ltext, pch = 19, col = orp16Scol[i1[3:4]])
  title("Linear regressions for bacterial communities\n(all datasets)", font.main = 1, xpd = NA, line = 0.7)
  label.figure("B", font = 2, cex = 1.5, yfrac = 1.05)
  # Groundwater, sediment, soil
  plot(xlim, ylim, type = "n", xlab = "log10(Number of samples)", ylab = quote("Slope of linear fit"~(V^-1)))
  abline(h = 0, lty = 2, lwd = 1.5, col = "gray50")
  points(log10(lmbac$nsamp[j2]), lmbac$slope[j2], pch = 19, col = orp16Scol[env$groupnum[j2]])
  ltext <- names(envirotype)[i2]
  legend("bottomright", ltext, pch = 19, col = orp16Scol[i2])

  ## Panel C: Distinctions in carbon oxidation state estimated for different hot springs 20210930
  # Use Hot Spring datasets
  i <- 4
  studies <- envirotype[[i]]
  # Assign colors
  # Use blue for neutral/alkaline
  col <- rep(orp16Scol[1], length(studies))
  # Use red for acidic
  iacid <- studies %in% c("PCL+18_Acidic", "LMG+20")
  col[iacid] <- orp16Scol[4]
  # Use turquoise for hypersaline
  ihyper <- studies %in% c("BMJ+19")
  col[ihyper] <- "turquoise3"
  # Use gray for sediments
  ised <- studies %in% c("PBU+20", "MWY+21")
  col[ised] <- "gray"

  # Offset for labels 20211012
  dx <- list(
    c(-80, 40, -190, 40, -150, 40, -220, 40, -180),
    c(40, -380, -130, 20, 40, 40, -120, NA, -200)
  )
  dy <- list(
    c(-0.022, 0, -0.015, 0, 0.004, 0, -0.03, 0, -0.01),
    c(0, -0.002, 0.03, 0, 0, 0, -0.007, NA, 0.007)
  )

  # Loop over Bacteria and Archaea
  lineages <- c("Bacteria", "Archaea")
  for(k in 1:2) {
    dat <- EZlm[EZlm$lineage == lineages[k], ]
    plot(c(-450, 400), c(-0.23, -0.08), type = "n", xlab = "Eh7 (mV)", ylab = cplab$ZC)
    for(j in seq_along(studies)) {
      with(dat[dat$study == studies[j], ], {
        # Divide by 1e3 to convert slope from V-1 to mV-1 20220520
        y <- function(x) intercept + slope * x / 1e3
        Eh7 <- c(Eh7min, Eh7max)
        ZC <- y(Eh7)
        if(length(Eh7) > 0) {
          if(abs(slope) < 0.01) lty <- 2 else lty <- 1
          lines(Eh7, ZC, col = col[j], lwd = 2, lty = lty)
          # Add number to identify dataset
          text(tail(Eh7, 1) + dx[[k]][j], tail(ZC, 1) + dy[[k]][j], j)
        }
      })
    }
    title(paste0("Hot spring\n", tolower(lineages[k])), font.main = 1, xpd = NA, line = 0.7)
    if(k==1) {
      label.figure("C", font = 2, cex = 1.5, yfrac = 1.025)
      # Add legend
      lhyper <- paste0("Hypersaline (", paste(which(ihyper), collapse = ", "), ")")
      lsed <- paste0("Sediment (", paste(which(ised), collapse = ", "), ")")
      lacid <- paste0("Acidic (", paste(which(iacid), collapse = ", "), ")")
      lneut <- "Circumneutral to"
      lalk <- paste0("Alkaline (", paste(which(col == orp16Scol[1]), collapse = ", "), ")")
      legend("topleft", c(lhyper, lsed, lacid, lneut, lalk), col = c("turquoise3", "gray", orp16Scol[4], orp16Scol[1], NA), lwd = 2, bty = "n")
    }
  }

  if(pdf) dev.off()

}

# Figure 5: Linear regressions between ZC and Eh7 at a global scale 20210828
orp16S5 <- function(pdf = FALSE) {

  if(pdf) pdf("Figure5.pdf", width = 10, height = 7)
  mat <- matrix(c(
    16, 16, rep(1:7, each = 6), 17, 17,
    16, 16, rep(8:14, each = 6), 18, 18,
    rep(15, 46),
    rep(0, 46),
    rep(0, 5), rep(19, 16), rep(20, 16), rep(21, 9)
  ), nrow = 5, byrow = TRUE)
  layout(mat, heights = c(2, 2, 0.5, 0.5, 4))
  par(mar = c(4, 4, 1, 1))
  par(mgp = c(2.5, 1, 0))

  ## Panel A: Scatterplots and fits for Bacteria and Archaea in each environment
  par(mar = c(0, 0, 1, 0))
  xlim <- range(EZdat$Eh7)
  ylim <- range(EZdat$ZC)
  eedat <- EZdat[EZdat$lineage == "Bacteria", ]
  eachenv(eedat, xlim = xlim, ylim = ylim, lineage = "Bacteria")
  par(mar = c(1, 0, 0, 0))
  eedat <- EZdat[EZdat$lineage == "Archaea", ]
  eachenv(eedat, xlim = xlim, ylim = ylim, lineage = "Archaea")
  # Add labels
  plot.new()
  text(0.5, -0.5, "Eh7 (mV)", cex = 1.2, xpd = NA)
  plot.new()
  text(0.2, 0.5, cplab$ZC, cex = 1.2, srt = 90, xpd = NA)
  label.figure("A", font = 2, cex = 2, xfrac = 0.2, yfrac = 0.97)
  plot.new()
  text(0.3, 0.5, "Bacteria", srt = 90, xpd = NA)
  plot.new()
  text(0.3, 0.5, "Archaea", srt = 90, xpd = NA)

  ## Panel B: Scatterplots and fits for Bacteria and Archaea in all environments 20210914
  # Start plot for Bacteria
  par(mar = c(4, 4, 1, 1))
  plot(c(-500, 650), range(EZdat$ZC), type = "n", xlab = "Eh7 (mV)", ylab = cplab$ZC)
  # Use Bacteria only
  thisdat <- EZdat[EZdat$lineage == "Bacteria", ]
  # Add linear fit; include number of studies in legend 20210925
  nstudy <- length(unique(thisdat$study))
  add.linear(thisdat$Eh7, thisdat$ZC, nstudy)
  # Add points
  eachenv(thisdat, add = TRUE, do.linear = FALSE)
  title("Bacteria", font.main = 1)
  label.figure("B", font = 2, cex = 2, xfrac = 0.02, yfrac = 1)

  # Now do Archaea
  plot(c(-500, 650), range(EZdat$ZC), type = "n", xlab = "Eh7 (mV)", ylab = cplab$ZC)
  thisdat <- EZdat[EZdat$lineage == "Archaea", ]
  nstudy <- length(unique(thisdat$study))
  add.linear(thisdat$Eh7, thisdat$ZC, nstudy)
  eachenv(thisdat, add = TRUE, do.linear = FALSE)
  title("Archaea", font.main = 1)

  # Add legend
  par(mar = c(4, 1, 1, 1))
  plot.new()
  ienv = c(1, 2, 4, 5, 3, 6, 7)
  ltext <- names(envirotype)[ienv]
  legend("left", ltext, pch = 19, col = orp16Scol[ienv], bty = "n")

  if(pdf) dev.off()

}

# Figure S2: ZC-Eh scatterplots for all studies 20210827
# This also creates files EZdat (Eh and ZC values) and
# EZlm (linear fits) for use by other plotting functions
orp16S_S2 <- function(pdf = FALSE) {

  # Setup figure
  if(pdf) pdf("Figure_S2.pdf", width = 12, height = 9)
  par(mfrow = c(3, 4))

  # To reproduce subsampling with size= argument 20220509
  set.seed(1234)
  results <- c(

    message("\nRiver & Seawater"),
    plotEZ("MLL+18", "Bacteria", groupby = "Season", groups = c("Summer", "Winter"), legend.x = "bottomright"),
    plotEZ("SVH+19", "two", groupby = "Type", groups = c("Oxic", "Suboxic", "Euxinic")),
    plotEZ("HXZ+20", "Bacteria", groupby = "Station", groups = c("SYBL", "C4")),
    plotEZ("KLY+20", "Bacteria", groupby = "Year", groups = c(2018, 2019), legend.x = "bottomright"),
    plotEZ("GSBT20_Prefilter", "two", groupby = "Region", groups = c("West Coast U.S.", "Great Lakes", "East Coast U.S.", "Europe", "Asia"), legend.x = "bottomright"),
    plotEZ("GSBT20_Postfilter", "two", groupby = "Region", groups = c("West Coast U.S.", "Great Lakes", "East Coast U.S.", "Europe", "Asia"), legend.x = "bottomright"),
    plotEZ("WHL+21", "Bacteria", groupby = "Season", groups = c("Spring", "Summer", "Autumn", "Winter"), legend.x = "bottomleft"),
    plotEZ("LXH+20", "Bacteria", groupby = "Season", groups = c("Summer", "Winter"), legend.x = "bottomright"),
    plotEZ("JVW+20", "Bacteria", groupby = "isolation_source", groups = c("Ulva laetevirens", "lagoon water"), legend.x = "topright"),
    plotEZ("ZZL+21", "Bacteria", groupby = "Location", groups = c("Main Stream", "Animal Farm", "Hospital", "WWTP", "Tributary"), legend.x = "bottomright"),
    plotEZ("GZL21", "Bacteria", groupby = "Type", groups = c("Surface water", "Middle water", "Bottom water"), legend.x = "bottomleft"),

    message("\nLake & Pond"),
    plotEZ("SAR+13", "two", groupby = "Zone", groups = c("Photic-oxic", "Transition", "Anoxic")),
    plotEZ("LZR+17", "Bacteria", groupby = "Elevation", groups = c("< 1000 m", "1000 - 4000 m", "> 4000 m"), legend.x = "bottomleft"),
    plotEZ("ECS+18", "Bacteria", groupby = "Lake", groups = c("Laguna Negra", "Lo Enca\u00F1ado")),
    plotEZ("LLC+19", "Bacteria", groupby = "Size", groups = c("Free-living", "Particle-associated")),
    plotEZ("SCH+16", "two", groupby = "Type", groups = c("Oxic", "Oxycline", "Anoxic")),
    plotEZ("BCA+21", "Bacteria", groupby = "Month", groups = c("Jul", "Nov", "Feb", "Apr")),
    plotEZ("HLZ+18", "Bacteria", groupby = "Type", groups = c("Reservoir", "Pond"), legend.x = "bottomright"),
    plotEZ("GRG+20", "two", groupby = "Type", groups = c("Oxic", "Anoxic"), legend.x = "bottomleft"),
    plotEZ("CNA+20", "Bacteria", groupby = "Season", groups = c("Summer", "Autumn", "Winter", "Spring"), legend.x = "topright"),
    plotEZ("BWD+19", "Bacteria", groupby = "Cover", groups = c("Ice", "Ice Free"), legend.x = "bottomright"),
    plotEZ("RSJ+21", "two", groupby = "Lake", groups = c("Kuiva", "Lovo"), legend.x = "topright"),
    plotEZ("LRL+22", "two", groupby = "Zone", groups = c("Freshwater", "Redoxcline", "Halocline")),
    plotEZ("BOEM21", "Bacteria", groupby = "Stratum", groups = c("Upper", "Chemocline", "Lower")),
    plotEZ("IBK+22", "two", groupby = "Land Use", groups = c("Arable", "Forest", "Grassland")),
    plotEZ("GSY+20", "two", groupby = "Lake", groups = c("La Zarza", "Filon Centro"), legend.x = "bottomright"),
    plotEZ("NLE+21", "Bacteria", groupby = "Year", groups = c("2017", "2018"), legend.x = "bottomleft"),
    plotEZ("SPA+21", "Bacteria", groupby = "Depth", groups = c("Epi", "Secchi", "Meso"), legend.x = "bottomleft"),
    plotEZ("FAV+21", "two", groupby = "Type", groups = c("Oxic Surface", "Anoxic Surface", "Bottom")),
    plotEZ("PSV+22", "Bacteria", groupby = "Location", groups = c("Center", "West"), legend.x = "bottomright"),

    # Hot Spring
    message("\nHot Spring"),
    plotEZ("SMS+12", "two", groupby = "Type", groups = c("Chemotrophic", "Transition", "Phototrophic")),
    plotEZ("PCL+18_Acidic", "two", legend.x = "bottomright"),
    plotEZ("PCL+18_Alkaline", "two"),
    plotEZ("BMJ+19", "two", groupby = "Environment", groups = c("Hydrothermal Pond", "Yellow Lake", "Black Lake", "Assale Lake", "Cave Water")),
    plotEZ("LMG+20", "two", groupby = "Setting", groups = c("River", "Hot Spring"), legend.x = "bottomleft"),
    plotEZ("GWSS21", "two", groupby = "Location", groups = c("Rehai", "Banglazhang")),
    plotEZ("GWS+20", "two", groupby = "Hydrothermal Field", groups = c("Batang", "Litang", "Kangding")),
    plotEZ("PBU+20", "Bacteria", groupby = "Type", groups = c("Cauldron", "Sampling Pit", "Spring", "Geyser Valley (Control)"), legend.x = "bottomright"),
    plotEZ("MWY+21", "two", groupby = "Location", groups = c("Quseyongba", "Moluojiang", "Daggyai", "Quzhuomu"), legend.x = "bottomright"),

    message("\nHyperalkaline"),
    plotEZ("SBP+20", "Bacteria", groupby = "pH Group", groups = c("< 10", "> 10")),
    plotEZ("RMB+17", "two", groupby = "pH Group", groups = c("< 10", "> 10")),
    plotEZ("CTS+17", "two", groupby = "Type", groups = c("River", "Well", "Spring"), legend.x = "bottomleft"),
    plotEZ("SPH+21", "Bacteria"),
    plotEZ("KSR+21", "Bacteria", groupby = "Location", groups = c("Lerone", "Branega", "Branega Creek Water")),
    plotEZ("PSB+21", "Bacteria", groupby = "O2 range", groups = c("> 0.5 mg/L", "0.2-0.5 mg/L", "< 0.2 mg/L")),
    plotEZ("NTB+21", "two", groupby = "Well", groups = c("BA1A", "BA1D")),

    message("\nGroundwater"),
    plotEZ("KLM+16", "Bacteria", groupby = "Day", groups = c(-1, 246, 448, 671)),
    plotEZ("YHK+19", "two"),
    plotEZ("SDH+19", "Bacteria", groupby = "Type", groups = c("Freshwater", "Brackish", "Saltwater")),
    plotEZ("SRM+19", "Bacteria", groupby = "Land Use", groups = c("Agriculture", "Community", "Landfill", "Mine")),
    plotEZ("APV+20", "two", groupby = "Type", groups = c("Canal", "Piezometer", "Well", "Spring")),
    plotEZ("SKP+21", "Bacteria", groupby = "Type", groups = c("Groundwater", "Surface water"), legend.x = "topright"),
    plotEZ("YHK+20", "Bacteria", groupby = "Location", groups = c("Upper Hillslope", "Middle Slope", "Lower Footslope")),
    plotEZ("JDP+20", "Bacteria", groupby = "Roll-Front Setting", groups = c("Oxidized", "Intermediate", "Reduced"), legend.x = "bottomright"),
    plotEZ("GWS+19", "Bacteria", groupby = "Type", groups = c("Un-confined", "Confined"), legend.x = "bottomleft"),
    plotEZ("SRM+21", "Bacteria", groupby = "Depth", groups = c("Surface", "Shallow", "Deep"), legend.x = "bottomleft"),
    plotEZ("ZCZ+21", "Bacteria", groupby = "Location", groups = c("LO", "CR1", "MN", "VA", "BS", "CR2"), legend.x = "topright"),
    plotEZ("CSW+22", "two", groupby = "BTEX", groups = c("High", "Low", "No"), legend.x = "topright"),
    plotEZ("GXS+22", "Bacteria", groupby = "Subarea", groups = c("A", "B", "C")),

    message("\nSediment"),
    plotEZ("JHL+12", "two", groupby = "Core", groups = c("GC6", "GC12"), legend.x = "bottomright"),
    plotEZ("GFE+16", "Bacteria", groupby = "Station", groups = c(1, 4, 7, 18)),
    plotEZ("ZML+17", "two", groupby = "Type", groups = c("Mangrove Forest", "Intertidal Mudflat"), legend.x = "bottomleft"),
    plotEZ("BYB+17", "two", groupby = "Type", groups = c("Freshwater", "Brackish", "Saltmarsh", "Hypersaline"), legend.x = "bottomright"),
    plotEZ("BSPD17", "Bacteria", groupby = "Type", groups = c("Water", "Sediment")),
    plotEZ("HDZ+19", "Bacteria", groupby = "Type", groups = c("Water", "Sediment")),
    plotEZ("TCN+17", "Bacteria", groupby = "Treatment", groups = c("Oxic", "Anoxic/Oxic", "Anoxic"), legend.x = "bottomright"),
    plotEZ("WHLH21", "Bacteria", groupby = "Position", groups = c("Surface", "Middle", "Bottom"), legend.x = "bottomleft"),
    plotEZ("SCM+18", "two", groupby = "Site", groups = c("Shallow", "Deep"), legend.x = "bottomleft"),
    plotEZ("RSS+18", "Bacteria", groupby = "Site", groups = c("Deep Hole", "Snowgoose Bay", "John's Island", "Skeleton Lake"), dylim = c(0, 0.001)),
    plotEZ("CLS+19", "two", groupby = "Type", groups = c("Water", "Sediment"), legend.x = "bottomright"),
    plotEZ("ZDA+20", "Bacteria", groupby = "Type", groups = c("Water", "Sediment"), legend.x = "bottomleft"),
    plotEZ("VMB+19", "two", groupby = "Type", groups = c("Water", "Sediment"), legend.x = "topright"),
    plotEZ("WHC+19", "Archaea", groupby = "Type", groups = c("Mid-tide", "Low-tide", "Subtidal"), legend.x = "topright"),
    plotEZ("HSF+19", "Bacteria", groupby = "Type", groups = c("Water", "Sediment-Water Interface", "Sediment"), legend.x = "topright"),
    plotEZ("RBM+21", "Bacteria", groupby = "Site", groups = c("3 (external)", "5 (middle)", "7 (inner)"), legend.x = "topright"),
    plotEZ("ZHZ+19", "two", groupby = "Treatment", groups = c("Original", "Nitrate-reducing", "Ferric-reducing", "Sulfate-reducing", "Methanogenic")),
    plotEZ("MCS+21", "two", groupby = "Type", groups = c("Water", "Sediment"), legend.x = "topright"),
    plotEZ("LMBA21_2017", "two", groupby = "Season", groups = c("Summer", "Winter"), legend.x = "bottomright"),
    plotEZ("HSF+22", "Bacteria", groupby = "Location", groups = c("West Lagoon", "North Lagoon", "South Lagoon", "Cinq Cases"), legend.x = "bottomright"),
    plotEZ("ZZLL21", "Bacteria", groupby = "Type", groups = c("Main stream", "Animal farm", "Hospital", "WWTP", "Tributary"), legend.x = "bottomleft"),
    plotEZ("BKR+22", "two", groupby = "Station", groups = c(5407, 5412, 5424, 5441, 5454)),
    plotEZ("WFB+21", "Bacteria", groupby = "Treatment", groups = c("C. volutator", "H. diversicolor", "Cv & Hd", "MPB", "Manual turbation")),
    plotEZ("HCW+22", "Bacteria", groupby = "Condition", groups = c("Static", "Weak", "Strong")),

    message("\nSoil"),
    plotEZ("SBW+17", "Bacteria", groupby = "Treatment", groups = c("Control", "BC400", "BC600")),
    plotEZ("MLL+19", "two", groupby = "Type", groups = c("Upland", "Paddy", "Sediment")),
    plotEZ("ZLH+22", "Bacteria", groupby = "Treatment", groups = c("Control", "Low FeCl2", "High FeCl2")),
    plotEZ("BMOB18", "two", groupby = "Treatment", groups = c("Acetate", "No amendment", "Pre-incubation")),
    plotEZ("ZZZ+18", "two", groupby = "Treatment", groups = c("None", "AQDS", "Biochar"), legend.x = "bottomleft"),
    plotEZ("PMM+20", "Bacteria", groupby = "Type", groups = c("Fluctuating", "Static")),
    plotEZ("WHLH21a", "Bacteria", groupby = "Stage", groups = c("Algae", "Cyanolichen", "Chlorolichen", "Moss"), legend.x = "bottomright"),
    plotEZ("CWC+20", "Bacteria", groupby = "Management", groups = c("Flooding", "Draining")),
    plotEZ("PSG+20", "two", groupby = "Treatment", groups = c("Initial", "NCC", "RB", "RGP", "TP")),
    plotEZ("XLD+20", "two", groupby = "Treatment", groups = c("Original Soil", "CK", "9K-Control", "Htt-sys", "Att-sys", "Co-sys"), legend.x = "bottomright"),
    plotEZ("LJC+20", "two", groupby = "meanT", groups = c("MAT >= 21.5 degC", "MAT < 21.5 degC")),
    plotEZ("DTJ+20", "two", groupby = "Zone", groups = c("Bulk Soil", "Mature", "Elongation", "Tip")),
    plotEZ("ZWH+22", "Bacteria", groupby = "Treatment", groups = c("Control", "Flooded", "Flooded + 0.5% silkworm excrement", "Flooded + 1% silkworm excrement"), legend.x = "bottomright"),
    plotEZ("LLL+21", "Bacteria", groupby = "Treatment", groups = c("CK", "FL", "EA", "SB", "BD")),
    plotEZ("RKSK22", "two", groupby = "Compartment", groups = c("Bulk sediment", "Rhizosphere", "Root"), legend.x = "bottomright"),
    plotEZ("DLS21_Bulk", "Bacteria", groupby = "Treatment", groups = c("control", "char", "silicate", "husk")),
    plotEZ("WKP+22", "Bacteria", groupby = "Type", groups = c("Intercropping", "Monoculture"), legend.x = "bottomleft"),
    plotEZ("CYG+22", "two", groupby = "Soil", groups = c("Fengyang", "Tancheng", "Shuyang")),
    plotEZ("CKB+22", "two", groupby = "Treatment", groups = c("FM 5 Mg/ha", "FM 10 Mg/ha", "RS 5 Mg/ha", "RS 10 Mg/ha"), legend.x = "bottomright")

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
  pearson.r <- unlist(lapply(results[names(results) == "pearson"], "[[", "estimate"))
  P.value <- unlist(lapply(results[names(results) == "pearson"], "[[", "p.value"))
  # Note: slope is mutiplied by 1e3 to convert from mV-1 to V-1
  EZlm <- data.frame(study, envirotype, lineage, nsamp, Eh7min, Eh7max,
    slope = signif(slope * 1e3, 6), intercept = signif(intercept, 6), pearson.r = signif(pearson.r, 6), P.value = signif(P.value, 6))
  # Save data and results to files
  write.csv(EZdat, "EZdat.csv", row.names = FALSE, quote = FALSE)
  write.csv(EZlm, "EZlm.csv", row.names = FALSE, quote = FALSE)

}


############################
### UNEXPORTED FUNCTIONS ###
############################

add.linear <- function(Eh7, ZC, nstudy = NA, O2 = FALSE) {
  # Plot linear fit and confidence interval 20210920
  text.col <- "black"
  line.col <- "gray62"
  # Use lighter colors for environments with few datasets 20210925
  if(!is.na(nstudy)) if(nstudy < 5) {
    text.col <- "gray60"
    line.col <- "gray70"
  }
  # Show number of studies and exit if there are zero 20210926
  if(!is.na(nstudy)) {
    legend("topleft", legend = nstudy, bty = "n", text.font = 2, text.col = text.col, inset = c(-0.05, 0))
    if(nstudy == 0) return()
  }
  # Fit data with linear model
  EZlm <- lm(ZC ~ Eh7)
  # Draw linear fit
  Ehvals <- range(Eh7)
  plx <- predict(EZlm, newdata = data.frame(Eh7 = Ehvals), se = T)
  lines(Ehvals, plx$fit, col = line.col)
  # Get the slope
  slope <- EZlm$coefficients[2]
  # Multiply by 1e3 to convert from mV-1 to V-1
  slope <- slope * 1e3
  # Round to fixed number of decimal places
  slope <- formatC(slope, digits = 3, format = "f")
  # Units for Eh7
  stext <- bquote(.(slope)~V^-1)
  # Units for O2
  if(O2) stext <- bquote(.(slope)~L~mu*mol^-1)
  # Also show number of samples
  Ntext <- bquote(italic(N) == .(length(ZC)))
  # Add text to plot
  legend <- as.expression(c(Ntext, stext))
  legend("topright", legend = legend, bty = "n", text.col = text.col)

  # Calculate Pearson correlation 20211009
  pearson <- cor.test(Eh7, ZC, method = "pearson")
  # Get P-value
  pval <- pearson$p.value
  # Format correlation coefficient
  rtext <- formatC(pearson$estimate, digits = 2, format = "f")
  rtext <- bquote(italic(r) == .(rtext))
  # Format P-value
  ptext <- formatC(signif(pval, 1))
  ptext <- bquote(italic(P) == .(ptext))
  ltext <- c(rtext, ptext)
  legend("bottomright", legend = ltext, bty = "n", text.col = text.col)
}

# Scatterplots for all samples in each environment type 20210913
eachenv <- function(eedat, add = FALSE, do.linear = TRUE, ienv = c(1, 2, 4, 5, 3, 6, 7), cols = orp16Scol,
  lineage = NULL, xlim = NULL, ylim = NULL, O2 = FALSE) {
  # Decide whether x variable is Eh7 or O2
  if(O2) Xvar <- eedat$O2_umol_L else Xvar <- eedat$Eh7
  eedat <- cbind(eedat, Xvar)
  # Get overall x and y limits
  if(is.null(xlim)) xlim <- range(eedat$Xvar)
  if(is.null(ylim)) ylim <- range(eedat$ZC)
  # Get names of environment types
  envirotypes <- names(envirotype)
  # Loop over environment types
  for(i in ienv) {
    # Start plot
    if(!add) plot(xlim, ylim, type = "n", xlab = "", ylab = "", axes = FALSE)
    # Get Eh7/O2 and ZC values
    thisdat <- eedat[eedat$envirotype == envirotypes[i], ]
    Xvar <- thisdat$Xvar
    ZC <- thisdat$ZC
    if(do.linear) {
      # Include number of studies in legend 20210925
      nstudy <- length(unique(thisdat$study))
      add.linear(Xvar, ZC, nstudy)
    }
    if(!isTRUE(add)) {
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
    points(Xvar, ZC, pch = 19, cex = 0.2, col = cols[i])
  }
}

# Get metadata for a study, appending columns for pch and col 20200914
getmdat_orp16S <- function(study, metrics = NULL, dropNA = TRUE, size = NULL, quiet = TRUE) {
  # Read metadata file
  # Remove suffix after underscore 20200929
  studyfile <- gsub("_.*", "", study)
  datadir <- system.file("extdata/orp16S", package = "JMDplots")
  file <- file.path(datadir, "metadata", paste0(studyfile, ".csv"))
  metadata <- read.csv(file, as.is = TRUE, check.names = FALSE)

  if(dropNA) {
    # Exclude samples with NA name 20200916
    iname <- match("name", tolower(colnames(metadata)))
    noname <- is.na(metadata[, iname])
    if(any(noname)) {
      if(!quiet) print(paste0("getmdat_orp16S [", study, "]: dropping ", sum(noname), " samples with NA name"))
      metadata <- metadata[!is.na(metadata[, iname]), ]
    }
  }
  # Use NULL pch as flag for unavailable dataset 20210820
  pch <- NULL
  infotext <- NULL

  ## Datasets for orp16S paper 20211003
  shortstudy <- study
  if(study == "RBW+14") {
    type <- rep("reducing", nrow(metadata))
    type[metadata$layer == "Top"] <- "oxidizing"
    type[metadata$layer == "SWI"] <- "transition"
    pch <- sapply(type, switch, oxidizing = 24, transition = 20, reducing = 25)
    col <- sapply(type, switch, oxidizing = 4, transition = 1, reducing = 2)
  }
  if(grepl("PCL\\+18", study)) {
    # PCL+18, PCL+18_Acidic, PCL+18_Alkaline
    Type <- sapply(strsplit(study, "_"), "[", 2)
    if(!is.na(Type)) metadata <- metadata[metadata$Type == Type, ]
    shortstudy <- "PCL+18"
  }
  if(grepl("LMBA21", study)) {
    # LMBA21, LMBA21_2013, LMBA21_2015, LMBA21_2016, LMBA21_2017, LMBA21_2018, LMBA21_2019
    Year <- sapply(strsplit(study, "_"), "[", 2)
    if(!is.na(Year)) metadata <- metadata[metadata$Year == Year, ]
    shortstudy <- "LMBA21"
  }
  if(grepl("DLS21", study)) {
    # DLS21, DLS21_Bulk, DLS21_Rhizosphere, DLS21_Plaque
    source <- sapply(strsplit(study, "_"), "[", 2)
    if(!is.na(source)) {
      if(source == "Bulk") metadata <- metadata[metadata$Source == "bulk soil", ]
      if(source == "Rhizosphere") metadata <- metadata[metadata$Source == "rhizosphere soil", ]
      if(source == "Plaque") metadata <- metadata[metadata$Source == "iron plaque", ]
    }
    shortstudy <- "DLS21"
  }
  if(grepl("GSBT20", study)) {
    # GSBT20, GSBT20_Prefilter, GSBT20_Postfilter
    filter <- sapply(strsplit(study, "_"), "[", 2)
    if(!is.na(filter)) {
      if(filter == "Postfilter") metadata <- metadata[grepl("post", metadata$Sample), ]
      if(filter == "Prefilter") metadata <- metadata[!grepl("post", metadata$Sample), ]
    }
    shortstudy <- "GSBT20"
  }
  if(shortstudy %in% c(
    "MLL+19", "HXZ+20", "BCA+21", "RSJ+21", "RMB+17", "SBP+20", "NTB+21", "MWY+21", "SAR+13", "CTS+17",
    "SCM+18", "HDZ+19", "BOEM21", "ZHZ+19", "YHK+20", "CNA+20", "BMJ+19", "SRM+19", "HLZ+18", "XLD+20",
    "JHL+12", "PSG+20", "KSR+21", "ZCZ+21", "SKP+21", "ZZL+21", "PBU+20", "GWS+19", "KLY+20", "SRM+21",
    "MLL+18", "JDP+20", "BWD+19", "LXH+20", "LMG+20", "WHL+21", "LLL+21", "SDH+19", "GWSS21", "HSF+19",
    "ZML+17", "DTJ+20", "WFB+21", "SBW+17", "KLM+16", "LMBA21", "ZDA+20", "ZZZ+18", "BSPD17", "CWC+20",
    "BMOB18", "JVW+20", "LJC+20", "GFE+16", "ECS+18", "FAV+21", "VMB+19", "DLS21", "ZZLL21", "GWS+20",
    "CLS+19", "SMS+12", "BYB+17", "MCS+21", "SVH+19", "PMM+20", "GZL21", "LLC+19", "NLE+21", "GSY+20",
    "SCH+16", "LZR+17", "GRG+20", "APV+20", "YHK+19", "WHC+19", "WHLH21", "PCL+18", "GSBT20", "SPA+21",
    "PSV+22", "CSW+22", "GXS+22", "IBK+22", "HSF+22", "RBM+21", "HCW+22", "BKR+22", "WKP+22", "ZLH+22",
    "CKB+22", "ZWH+22", "LRL+22", "WHLH21a", "CYG+22", "RSS+18", "SPH+21", "PSB+21", "RKSK22", "TCN+17"
  )) {
    # General processing of metadata for orp16S datasets 20210820
    # Get Eh or ORP values (uses partial name matching, can match a column named "Eh (mV)")
    Eh <- metadata$Eh
    Ehname <- "Eh"
    if(is.null(Eh)) {
      Eh <- metadata$ORP
      Ehname <- "ORP"
    }
    # Also match a column starting with "redox" (case-insensitive)
    iEh <- grep("^redox", colnames(metadata), ignore.case = TRUE)
    if(length(iEh) > 0) {
      Eh <- metadata[, iEh[1]]
      Ehname <- strsplit(colnames(metadata)[iEh[1]], " ")[[1]][1]
    }
    if(is.null(Eh)) stop("can't find Eh or ORP column")

    # Append Ehorig column (normalized column name used for exporting data) 20210821
    metadata <- cbind(metadata, Ehorig = Eh)
    
    # Get temperature for dEh/dpH calculation 20210829
    iT <- match("T", colnames(metadata))  # matches "T"
    if(is.na(iT)) iT <- grep("^T\\ ", colnames(metadata))[1]  # matches "T (°C)" but not e.g. "Treatment"
    if(is.na(iT)) iT <- grep("^Temp", colnames(metadata))[1]  # matches "Temperature (°C)"
    if(!is.na(iT)) {
      T <- metadata[, iT]
      Ttext <- paste(round(range(na.omit(T)), 1), collapse = " to ")
    } else {
      T <- rep(25, nrow(metadata))
      Ttext <- "assumed 25"
    }

    # Adjust Eh to pH 7 20210828
    # Find pH column
    ipH <- match("pH", colnames(metadata))
    if(!is.na(ipH)) {
      # pH values
      pH <- metadata[, ipH]
      pHtext <- paste(round(range(na.omit(pH)), 1), collapse = " to ")
      ## Eh(mV)-pH slope at 25 °C
      #dEhdpH <- -59.2
      # Find dEh/dpH as a function of T 20210829
      TK <- T + 273.15
      R <- 0.0083147
      F <- 96.4935
      dEhdpH <- - (log(10) * R * TK) / F * 1000
      # Difference to pH 7
      pHdiff <- 7 - pH
      # Adjust to pH 7
      Eh7 <- round(Eh + pHdiff * dEhdpH, 2)
    } else {
      Eh7 <- Eh
      pHtext <- "assumed 7"
    }
    metadata <- cbind(metadata, Eh7 = Eh7)
    # Print message about T, pH and Eh ranges
    Ehtext <- paste(round(range(na.omit(Eh))), collapse = " to ")
    Eh7text <- paste(round(range(na.omit(Eh7))), collapse = " to ")
    infotext <- paste0(study, ": T ", Ttext, ", pH ", pHtext, ", Eh ", Ehtext, ", Eh7 ", Eh7text)
    if(!quiet) print(infotext)

    # Get O2 values 20220517
    O2_umol_L <- NA
    iO2 <- grep("^O2", colnames(metadata))[1]
    if(is.na(iO2)) iO2 <- grep("^DO", colnames(metadata))[1]
    if(is.na(iO2)) iO2 <- grep("^oxygen", tolower(colnames(metadata)))[1]
    if(!is.na(iO2)) {
      # Get the units
      O2name <- colnames(metadata)[iO2]
      # Remove text up to and including left parenthesis, and right parenthesis
      units <- gsub(".*\\(|\\)", "", O2name)
      if(units %in% c("\u03BCmol/L", "\u03BCmol L-1", "umol/L", "umol L-1", "umol kg-1", "\u03BCM")) O2_umol_L <- metadata[, iO2]
      else if(units %in% c("mM", "mmol L-1")) O2_umol_L <- metadata[, iO2] * 1000
      else if(units %in% c("mg/L", "mg L-1")) O2_umol_L <- metadata[, iO2] * 1000 / 32
      else if(units %in% c("mL/L", "mL L-1")) O2_umol_L <- metadata[, iO2] * 1.42905 * 1000 / 32
      else if(units == "%") {
        # Calculate logK for O2(gas) = O2(aq)
        logK <- rep(NA, nrow(metadata))
        not_NA <- !is.na(T)
        logK[not_NA] <- subcrt(c("O2", "O2"), c("gas", "aq"), c(-1, 1), T = T[not_NA])$out$logK
        # Calculate saturated oxygen concentration from logK = logaO2(aq) - logfO2(gas)
        # Assume fugacity is partial pressure of oxygen in sea-level atmosphere
        logaO2 <- logK + log10(0.21 * 1.01325)
        # Remove logarithm and convert mol to umol
        O2_umol_L_sat <- 10^logaO2 * 1e6
        # Multiply by percent to get concentration
        O2_umol_L <- metadata[, iO2] / 100 * O2_umol_L_sat
      }
      else if(identical(units, O2name)) warning(paste0(study, ": no units given for ", O2name))
      else stop(paste0(study, ": unrecognized units for oxygen concentration: ", units))
    }
    metadata <- cbind(metadata, O2_umol_L = O2_umol_L)

    # Cluster samples by Eh7 value
    if(nrow(metadata) > 3) {
      # Keep NA values out of clusters
      ina <- is.na(Eh7)
      # Divide data into 3 clusters
      cl <- kmeans(Eh7[!ina], 3, nstart = 100)
      # Find the clusters with the lowest and highest means
      imin <- which.min(cl$centers[, 1])
      imax <- which.max(cl$centers[, 1])
      imid <- (1:3)[-c(imin, imax)]
      # Name the clusters
      cluster <- rep(NA, nrow(metadata))
      cluster[!ina][cl$cluster == imin] <- "reducing"
      cluster[!ina][cl$cluster == imid] <- "transition"
      cluster[!ina][cl$cluster == imax] <- "oxidizing"
      # Get pch and col for plot
      # Use is.null to let previous assignments take precedence (for SVH+19 in geo16S) 20211008
      if(is.null(pch)) {
        pch <- sapply(cluster, switch, oxidizing = 24, transition = 21, reducing = 25, NA)
        col <- sapply(cluster, switch, oxidizing = 4, transition = 1, reducing = 2, NA)
      }
      if(!quiet) {
        # Print message about Eh7 or ORP ranges of clusters
        ro <- range(na.omit(Eh7[cluster=="oxidizing"]))
        no <- sum(na.omit(cluster=="oxidizing"))
        message(paste0("getmdat_orp16S: oxidizing Eh7 ", paste(ro, collapse = " to "), " mV (n = ", no, ")"))
        rt <- range(na.omit(Eh7[cluster=="transition"]))
        nt <- sum(na.omit(cluster=="transition"))
        message(paste0("getmdat_orp16S: transition Eh7 ", paste(rt, collapse = " to "), " mV (n = ", nt, ")"))
        rr <- range(na.omit(Eh7[cluster=="reducing"]))
        nr <- sum(na.omit(cluster=="reducing"))
        message(paste0("getmdat_orp16S: reducing Eh7 ", paste(rr, collapse = " to "), " mV (n = ", nr, ")"))
      }
      # Append cluster names to data frame
      metadata <- cbind(metadata, cluster)
    } else {
      # For very small datasets, don't create clusters but still assign pch and col
      pch <- rep(21, length(Eh7))
      pch[which.max(Eh7)] <- 24
      pch[which.min(Eh7)] <- 25
      col <- rep(1, length(Eh7))
      col[which.max(Eh7)] <- 4
      col[which.min(Eh7)] <- 2
    }
  }

  if(is.null(pch)) stop(paste(study, "metadata file exists, but not set up for processing"))

  metadata <- cbind(metadata, pch, col)
  # Use the infotext as an attribute for printing by orp16S_info 20220513
  attr(metadata, "infotext") <- infotext
  # Return both metadata and metrics, if provided 20220506
  if(is.null(metrics)) metadata else {
    # Keep metadata only for samples with metrics 20201006
    metadata <- metadata[metadata$Run %in% metrics$Run, ]
    # Keep specified number of random samples 20220509
    if(!is.null(size)) if(nrow(metadata) > size) {
      isample <- sample(1:nrow(metadata), size = size)
      metadata <- metadata[isample, ]
    }
    # Put metrics in same order as metadata 20220505
    imet <- match(metadata$Run, metrics$Run)
    metrics <- metrics[imet, ]
    # Insert sample column in metrics
    # Use first column name starting with "sample" or "Sample" 20210818
    sampcol <- grep("^sample", colnames(metadata), ignore.case = TRUE)[1]
    metrics <- data.frame(Run = metrics$Run, sample = metadata[, sampcol], nH2O = metrics$nH2O, ZC = metrics$ZC)
    list(metadata = metadata, metrics = metrics)
  }
}

# Function to calculate metrics for a given study 20220506
getmetrics_orp16S <- function(study, mincount = 100, quiet = TRUE, ...) {
  # Remove suffix after underscore 20200929
  studyfile <- gsub("_.*", "", study)
  datadir <- system.file("extdata/orp16S/RDP", package = "JMDplots")
  RDPfile <- file.path(datadir, paste0(studyfile, ".tab.xz"))
  # If there is no .xz file, look for a .tab file 20210607
  if(!file.exists(RDPfile)) RDPfile <- file.path(datadir, paste0(studyfile, ".tab"))
  RDP <- readRDP(RDPfile, mincount = mincount, quiet = quiet, ...)
  map <- mapRDP(RDP, quiet = quiet)
  getmetrics(RDP, map)
}

# Function to gather info (study name, number of samples with Bacteria and Archaea,
# T, pH, Eh, and Eh7 ranges) for filling in Table S1 20220513
orp16S_info <- function(study) {
  # Get all samples with mincount = 100
  metrics <- getmetrics_orp16S(study)
  # Get metadata for these samples
  mdat <- getmdat_orp16S(study, metrics)
  metadata <- mdat$metadata
  iname <- match("name", tolower(colnames(metadata)))
  name <- na.omit(metadata[, iname])[1]
  print(name)
  # Number of samples with Bacteria and Archaea (note mincount = 100 for each one)
  metrics.Bac <- try(getmetrics_orp16S(study, lineage = "Bacteria"), TRUE)
  if(inherits(metrics.Bac, "try-error")) nBac <- 0 else {
    mdat.Bac <- getmdat_orp16S(study, metrics.Bac)
    nBac <- nrow(mdat.Bac$metrics)
  }
  metrics.Arc <- try(getmetrics_orp16S(study, lineage = "Archaea"), TRUE)
  if(inherits(metrics.Arc, "try-error")) nArc <- 0 else {
    mdat.Arc <- getmdat_orp16S(study, metrics.Arc)
    nArc <- nrow(mdat.Arc$metrics)
  }
  print(paste0("nBac: ", nBac, "; nArc: ", nArc))
  # Now print the T, pH, Eh, Eh7 ranges
  print(attributes(metadata)$infotext)
  # Read linear fit coefficients
  dat <- EZlm[EZlm$study == study, ]
  for(lineage in c("Bacteria", "Archaea")) {
    idat <- dat$lineage == lineage
    if(any(idat)) {
      slope <- dat$slope[idat]
      slopetxt <- "-- (close to zero)"
      if(slope > 0.01) slopetxt <- "positive (> 0.01 V-1)"
      if(slope < -0.01) slopetxt <- "negative (< -0.01 V-1)"
      slopetxt <- paste("Slope of ZC-Eh7 correlation for", lineage, "is", slopetxt)
      print(slopetxt)
    }
  }
}

# Summary table of regression slopes 20220516
orp16S_T2 <- function() {
  # All environment types
  envirotype <- unique(EZlm$envirotype)
  # Initialize output
  out <- matrix(nrow = length(envirotype), ncol = 6)
  # Loop over Bacteria and Archaea
  lineage <- c("Bacteria", "Archaea")
  for(ilin in 1:2) {
    # Column offset
    if(lineage[ilin] == "Archaea") dcol <- 3 else dcol <- 0
    # Loop over environment type
    for(ienv in 1:length(envirotype)) {
      # Get fitting results
      thisdat <- EZlm[EZlm$lineage == lineage[ilin] & EZlm$envirotype == envirotype[ienv], ]
      # Get number of datasets and count those with positive and negative slopes
      ndat <- nrow(thisdat)
      npos <- sum(thisdat$slope > 0.01)
      nneg <- sum(thisdat$slope < -0.01)
      # Calcualate percentages
      ppos <- round(npos / ndat * 100)
      pneg <- round(nneg / ndat * 100)
      # Enter values into table
      out[ienv, 1 + dcol] <- ndat
      out[ienv, 2 + dcol] <- ppos
      out[ienv, 3 + dcol] <- pneg
    }
  }
  # Put in names
  rownames(out) <- envirotype
  colnames(out) <- paste(rep(c("N", "Pos", "Neg"), 2), rep(c("Bac", "Arc"), each = 3), sep = "_")
  out
}

# Global Eh-pH diagram 20220516
orp16S_S1 <- function(pdf = FALSE) {
  if(pdf) pdf("Figure_S1.pdf", width = 8, height = 6)
  layout(t(matrix(1:2)), widths = c(3, 1))
  # Get data for unique samples
  ssr <- paste(EZdat$study, EZdat$sample, EZdat$Run, sep = "_")
  idup <- duplicated(ssr)
  thisdat <- EZdat[!idup, ]
  # Remove NA pH
  thisdat <- thisdat[!is.na(thisdat$pH), ]
  # Get colors
  ienv <- match(thisdat$envirotype, names(envirotype))
  col <- orp16Scol[ienv]
  # Make plot
  par(mar = c(4, 4, 1, 1))
  plot(thisdat$pH, thisdat$Eh, col = col, pch = 19, cex = 0.5, xlab = "pH", ylab = "Eh (mV)",
       xlim = c(-2, 14), ylim = c(-600, 1000), xaxs = "i", yaxs = "i", axes = FALSE)
  box()
  axis(1, at = seq(-2, 14, 2))
  axis(2, at = seq(-600, 1000, 200))
  # Add outline from Bass Becking et al. (1960)
  BKM60 <- read.csv(system.file("extdata/orp16S/BKM60.csv", package = "JMDplots"))
  segment <- NULL
  top <- subset(BKM60, segment == "top")
  lines(top$pH, top$Eh)
  bottom <- subset(BKM60, segment == "bottom")
  lines(bottom$pH, bottom$Eh)
  left <- subset(BKM60, segment == "left")
  for(i in 1:(nrow(left)/2)) {
    ii <- c(i*2-1, i*2)
    lines(left$pH[ii], left$Eh[ii])
  }
  right <- subset(BKM60, segment == "right")
  for(i in 1:(nrow(right)/2)) {
    ii <- c(i*2-1, i*2)
    lines(right$pH[ii], right$Eh[ii])
  }
  O2 <- subset(BKM60, segment == "O2")
  lines(O2$pH, O2$Eh, lwd = 2)
  H2 <- subset(BKM60, segment == "H2")
  lines(H2$pH, H2$Eh, lwd = 2)
  text(0, 55, quote(H^"+"), cex = 1.2, srt = -30)
  text(-0.5, -35, quote(H[2]), cex = 1.2, srt = -30)
  text(5.5, 840, quote(H[2]*O), cex = 1.2, srt = -30)
  text(6, 930, quote(O[2]), cex = 1.2, srt = -30)
  # Add legend
  par(mar = c(4, 1, 1, 1))
  plot.new()
  ienv = c(1, 2, 4, 5, 3, 6, 7)
  ltext <- names(envirotype)[ienv]
  legend("left", ltext, pch = 19, col = orp16Scol[ienv], bty = "n")
  if(pdf) dev.off()
}

# Compare regressions with Eh and O2 20220517
orp16S_S3 <- function(pdf = FALSE) {

  if(pdf) pdf("Figure_S3.pdf", width = 7, height = 6)
  mat <- matrix(1:6, nrow = 2, byrow = TRUE)
  layout(mat, widths = c(2, 2, 1))

  # Use only samples with non-NA O2
  thisdat <- EZdat[!is.na(EZdat$O2_umol_L), ] 

  ## Bacteria only
  par(mar = c(4, 4, 2.5, 1))
  bacdat <- thisdat[thisdat$lineage == "Bacteria", ]
  # Start ZC-Eh7 plot
  plot(c(-500, 650), range(bacdat$ZC), type = "n", xlab = "Eh7 (mV)", ylab = cplab$ZC)
  # Add linear fit; include number of studies in legend 20210925
  nstudy <- length(unique(bacdat$study))
  add.linear(bacdat$Eh7, bacdat$ZC, nstudy)
  # Add points
  eachenv(bacdat, add = TRUE, do.linear = FALSE)
#  label.figure("A", cex = 1.5, font = 2)
  # Start ZC-O2 plot
  plot(c(0, 750), range(bacdat$ZC), type = "n", xlab = quote(O[2]~"("*mu*"mol L"^{-1}*")"), ylab = cplab$ZC)
  add.linear(bacdat$O2_umol_L, bacdat$ZC, nstudy, O2 = TRUE)
  eachenv(bacdat, add = TRUE, do.linear = FALSE, O2 = TRUE)
  mtext("Bacteria", adj = -0.42, line = 1)
  # Add legend
  par(mar = c(0, 0, 0, 0))
  plot.new()
  ienv = c(1, 2, 4, 5, 3, 6, 7)
  ltext <- names(envirotype)[ienv]
  # Add number of samples in each environment 20220518
  nsamp <- table(bacdat$envirotype)[ltext]
  ltext <- paste0(ltext, " (", nsamp, ")")
  legend("left", ltext, pch = 19, col = orp16Scol[ienv], bty = "n")

  ## Archaea only
  par(mar = c(4, 4, 2.5, 1))
  arcdat <- thisdat[thisdat$lineage == "Archaea", ]
  # Start ZC-Eh7 plot
  plot(c(-500, 650), range(arcdat$ZC), type = "n", xlab = "Eh7 (mV)", ylab = cplab$ZC)
  # Add linear fit; include number of studies in legend 20210925
  nstudy <- length(unique(arcdat$study))
  add.linear(arcdat$Eh7, arcdat$ZC, nstudy)
  # Add points
  eachenv(arcdat, add = TRUE, do.linear = FALSE)
#  label.figure("B", cex = 1.5, font = 2)
  # Start ZC-O2 plot
  plot(c(0, 750), range(arcdat$ZC), type = "n", xlab = quote(O[2]~"("*mu*"mol L"^{-1}*")"), ylab = cplab$ZC)
  add.linear(arcdat$O2_umol_L, arcdat$ZC, nstudy, O2 = TRUE)
  eachenv(arcdat, add = TRUE, do.linear = FALSE, O2 = TRUE)
  mtext("Archaea", adj = -0.42, line = 1)
  # Add legend
  par(mar = c(0, 0, 0, 0))
  plot.new()
  ienv = c(1, 2, 4, 5, 3, 6, 7)
  ltext <- names(envirotype)[ienv]
  # Add number of samples in each environment 20220518
  nsamp <- table(arcdat$envirotype)[ltext]
  ltext <- paste0(ltext, " (", nsamp, ")")
  legend("left", ltext, pch = 19, col = orp16Scol[ienv], bty = "n")

  if(pdf) dev.off()

}
