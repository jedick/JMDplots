# JMDplots/orp16S.R
# Plots for paper on ZC-ORP correlations 20210827

## Uncomment to source and run these functions interactively (developer mode)
#source("chem16S.R")
#options(chem16Sdir = system.file("extdata/orp16S", package = "JMDplots"))

# Group studies by environment types 20210828
envirotype <- list(
  "River & Seawater" = c("MLL+18", "SVH+19", "HXZ+20", "KLY+20", "WHL+21", "LXH+20", "JVW+20", "ZZL+21", "GZL21"),
  "Lake & Pond" = c("SAR+13", "LZR+17", "ECS+18", "LLC+19", "SCH+16", "BCA+21", "HLZ+18", "CNA+20",
                    "BWD+19", "RSJ+21", "BOEM21", "GSY+20", "NLE+21", "FAV+21", "GRG+20"),
  "Groundwater" = c("KLM+16", "YHK+19", "SDH+19", "SRM+19", "APV+20", "SKP+21", "YHK+20", "JDP+20", "GWS+19", "SRM+21", "ZCZ+21"),
  # NOTE: Keep Hot Spring at #4 to get red color 20210904
  "Hot Spring" = c("SMS+12", "PCL+18_Acidic", "PCL+18_Alkaline", "BMJ+19", "LMG+20", "GWSS21", "GWS+20", "PBU+20", "MWY+21", "OFY+19"),
  "Alkaline Spring" = c("SBP+20", "RMB+17", "CTS+17", "KSR+21", "NTB+21"),
  "Sediment" = c("JHL+12", "GFE+16", "ZML+17", "BYB+17", "BSPD17", "HDZ+19", "WHLH21", "SCM+18",
                 "CLS+19", "ZDA+20", "VMB+19", "WHC+19", "HSF+19", "MCS+21", "LMBA21_2017", "ZZLL21", "WFB+21"),
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

# Figure 2: Chemical and geobiochemical depth profiles in Winogradsky columns 20210829
orp16S2 <- function(pdf = FALSE) {

  if(!grepl("orp16S", options("chem16Sdir")[[1]])) 
    stop('Please run this first: options(chem16Sdir = system.file("extdata/orp16S", package = "JMDplots"))')

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

  if(!grepl("orp16S", options("chem16Sdir")[[1]])) 
    stop('Please run this first: options(chem16Sdir = system.file("extdata/orp16S", package = "JMDplots"))')

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
    "LZR+17, 30.587, 104.310", # Table 1
    "ECS+18, -33.65, -70.117", # Materials and methods
    "LLC+19, 24.82, 118.15", # SAMN04549101
    "SCH+16, 67.081, -50.355", # Materials and methods
    "BCA+21, 46.3615, 25.0509", # SAMN07409474
    "HLZ+18, 24.795, 118.138", # SAMN07638080
    "CNA+20, 39.441, -77.371", # Web search for geographic center of Maryland --> https://sos.maryland.gov/mdkids/Pages/Geography.aspx
    "BWD+19, 47.120571, -88.545425", # SAMN09980099
    "RSJ+21, 61.833, 24.283", # Materials and methods
    "BOEM21, 43.051389, -75.965", # Materials and methods
    "GSY+20, 37.59, -7.124", # Materials and methods
    "NLE+21, 32.833, 35.583", # SAMEA7280991
    "FAV+21, 0.757, 36.372", # SAMN19267646
    "GRG+20, 37.726, -6.555", # Materials and methods
    # Hot Spring
    "SMS+12, 44.6, -110.9", # JGI IMG/M sample name 1_050719N
    "PCL+18_Acidic, -38.5, 176.0", # Fig. 1
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
#    "YHK+19, 40.460, -87.764", # SAMD00089561
    "SDH+19, 23.03, 113.38", # SAMN07692244
    "SRM+19, 12.67417, 101.3889", # Materials and methods
    "APV+20, 20.12, -99.23", # Materials and methods
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
    "HDZ+19, 29.901, 113.52435", # SAMN05990289
    "WHLH21, 23.52, 113.495", # Materials and methods
    "SCM+18, 25.25, -97.23", # Table 1
    "CLS+19, 32.22, 118.83", # SAMN08683376
    "ZDA+20, -21.42996, -70.05874", # SAMEA4858706
    "VMB+19, 52.11, 79.17", # SAMN08987150
    "WHC+19, 30.12, 122.14", # Materials and methods
    "HSF+19, 47.803, 16.709", # methods
    "MCS+21, -32.15, -71.1", # Materials and methods
    "LMBA17_2017, 43.42, -2.7", # Fig. 1
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
  mapPlot(col = "slategray2", projection = "+proj=wintri", border = "white", drawBox = FALSE)

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
  latlon <- paste(dat$latitude, dat$longitude)
  isuniq <- !duplicated(latlon)
  dat <- dat[isuniq, ]
  mapPoints(dat$longitude, dat$latitude, col = orp16Scol[4], lwd = 2)
  # Coordinates for Three Gorges Reservoir are from Table S1 of GZL21
  dat <- getmdat("GZL21")
  dat <- dat[!is.na(dat$"ORP (mV)"), ]
  latlon <- paste(dat$Latitude, dat$Longitude)
  isuniq <- !duplicated(latlon)
  dat <- dat[isuniq, ]
  mapPoints(dat$Longitude, dat$Latitude, col = orp16Scol[1], lwd = 1)
  # Coordinates for Mahomet Aquifer are from Table S1 of YHK+19
  dat <- getmdat("YHK+19")
  latlon <- paste(dat$Latitude, dat$Longitude)
  mapPoints(dat$Longitude, dat$Latitude, col = orp16Scol[3], lwd = 1)

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
    "ZZL+21", "MLL+18", "SDH+19", "ZML+17", "ZZLL21", "ZZZ+18", "ZHZ+19", "WHLH21" # GD-HK-MO GBA
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

# Figure 4: Analysis of selected datasets for each environment type 20211003
orp16S4 <- function(pdf = FALSE) {

  if(!grepl("orp16S", options("chem16Sdir")[[1]])) 
    stop('Please run this first: options(chem16Sdir = system.file("extdata/orp16S", package = "JMDplots"))')

  if(pdf) pdf("Figure4.pdf", width = 12, height = 6)
  par(mfrow = c(2, 4))
  par(mar = c(4, 4, 3, 1))
  par(cex.lab = 1.1)
  par(mgp = c(2.5, 1, 0))

  # Black Sea (River & Seawater)
  plotEZ("SVH+19", "Bacteria", groupby = "Type", groups = c("Oxic", "Suboxic", "Euxinic"), title.line = 1.5)
  # Fayetteville Green Lake (Lake & Pond)
  plotEZ("BOEM21", "Bacteria", groupby = "Stratum", groups = c("Upper", "Chemocline", "Lower"), title.line = 1.5)
  # Samail Ophiolite (Alkaline Spring)
  plotEZ("RMB+17", "Bacteria", groupby = "pH Group", groups = c("< 10", "> 10"), title.line = 1.5)
  # Eastern Tibetan Plateau (Hot Spring)
  plotEZ("GWS+20", "Bacteria", groupby = "Hydrothermal Field", groups = c("Batang", "Litang", "Kangding"), title.line = 1.5)
  # Hainich Critical Zone (Groundwater)
  plotEZ("YHK+20", "Bacteria", groupby = "Location", groups = c("Upper Hillslope", "Middle Slope", "Lower Footslope"), title.line = 1.5)
  # Phetchabun-Pichit Gold Mine (Groundwater)
  plotEZ("SKP+21", "Bacteria", groupby = "Type", groups = c("Groundwater", "Surface water"), legend.x = "topright", title.line = 1.5)
  # Bay of Biscay (Sediment)
  plotEZ("LMBA21_2017", "Bacteria", groupby = "Season", groups = c("Summer", "Winter"), legend.x = "bottomright", title.line = 1.5)
  # Hunan Soil (Soil)
  plotEZ("MLL+19", "Bacteria", groupby = "Type", groups = c("Upland", "Paddy", "Sediment"), title.line = 1.5)

  if(pdf) dev.off()

}

# Figure 5: Local and global analysis of ZC and Eh7 data 20210828
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
        if(length(slope) > 0) {
          y <- function(x) intercept + slope * x
          Eh7 <- c(Eh7min, Eh7max)
          ZC <- y(Eh7)
          # Use dashed line for fits with low slope 20211002
          if(abs(slope * 1000) < 0.01) lty <- 2 else lty <- 1
          lines(Eh7, ZC, col = orp16Scol[i], lty = lty)
        }
      })
    }
  }
  l1names <- names(envirotype)[i1][c(1, NA, 2, 3, 4)]
  l1names[1] <- "River &"
  l1names[2] <- "Seawater"
  l1names[5] <- "Alkaline Spr."
  l1col <- orp16Scol[i1][c(1, NA, 2, 3, 4)]
  legend("bottomright", l1names, lty = 1, col = l1col, cex = 0.8, bty = "n")
  label.figure("A", font = 2, cex = 2, xfrac = 0, yfrac = 0.9)

  # Groundwater, sediment, soil
  plot(c(-400, 750), c(-0.23, -0.12), type = "n", xlab = "Eh7 (mV)", ylab = cplab$ZC)
  i2 <- c(3, 6, 7)
  for(i in i2) {
    studies <- envirotype[[i]]
    for(j in seq_along(studies)) {
      with(dat[dat$study == studies[j], ], {
        # Need this 'if' to catch study without Bacteria (WHC+19) 20211005
        if(length(slope) > 0) {
          y <- function(x) intercept + slope * x
          Eh7 <- c(Eh7min, Eh7max)
          ZC <- y(Eh7)
          if(abs(slope * 1000) < 0.01) lty <- 2 else lty <- 1
          lines(Eh7, ZC, col = orp16Scol[i], lty = lty)
        }
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

  ## Panel D: Scatterplots and fits for Bacteria and Archaea in all environments except hot springs 20210914
  # Read Eh7 - ZC data
  dat <- read.csv(system.file("extdata/orp16S/EZdat.csv", package = "JMDplots"))
  # Omit data for Hot Spring environment
  dat <- dat[dat$envirotype != "Hot Spring", ]

  # Start plot for Bacteria
  par(mar = c(4, 4, 1, 1))
  plot(c(-500, 650), range(dat$ZC), type = "n", xlab = "Eh7 (mV)", ylab = cplab$ZC)
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
  plot(c(-500, 650), range(dat$ZC), type = "n", xlab = "Eh7 (mV)", ylab = cplab$ZC)
  thisdat <- dat[dat$lineage == "Archaea", ]
  nstudy <- length(unique(thisdat$study))
  add.linear(thisdat$Eh7, thisdat$ZC, nstudy)
  eachenv("Archaea", add = TRUE, do.linear = FALSE)
  title("Archaea", font.main = 1)

  # Add legend
  par(mar = c(4, 1, 1, 1))
  plot.new()
  ienv = c(1, 2, 5, 3, 6, 7, NA)
  ltext <- names(envirotype)[ienv]
  ltext[7] <- "(No hot springs)"
  legend("left", ltext, pch = 19, col = orp16Scol[ienv], bty = "n")

  if(pdf) dev.off()
}

# Figure 6: Distinctions in carbon oxidation state estimated for different hot springs, and global fits for all environments. 20210930
orp16S6 <- function(pdf = FALSE) {

  if(!grepl("orp16S", options("chem16Sdir")[[1]])) 
    stop('Please run this first: options(chem16Sdir = system.file("extdata/orp16S", package = "JMDplots"))')

  if(pdf) pdf("Figure6.pdf", width = 8, height = 6)
  mat <- matrix(c(1,1,1,1, 2,2,2,2, 3,3,3,3, 0,4,4,4,4,4, 5,5,5,5,5,0), nrow = 2, byrow = TRUE)
  layout(mat, heights = c(3, 2))
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
  iacid <- studies %in% c("PCL+18_Acidic", "LMG+20")
  col[iacid] <- orp16Scol[4]
  # Use turquoise for hypersaline
  ihyper <- studies %in% c("BMJ+19")
  col[ihyper] <- "turquoise3"
  # Use gray for sediments
  ised <- studies %in% c("PBU+20", "MWY+21", "OFY+19")
  col[ised] <- "gray"

  # Offset for labels 20211012
  dx <- list(
    c(-60, 20, -190, 20, 20, 20, 20, 20, 20, 35),
    c(20, 20, 20, 20, 20, 20, 20, 20, -90, 20),
    c(-80, 20, -150, 20, 20, 20, 20, 20, -90, 35)
  )
  dy <- list(
    c(-0.022, 0, -0.015, 0, 0.003, 0, 0, 0, 0, 0),
    c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    c(-0.025, 0, -0.01, 0, 0, 0, -0.003, 0, -0.0025, 0)
  )

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
          text(tail(Eh7, 1) + dx[[k]][j], tail(ZC, 1) + dy[[k]][j], j)
        }
      })
    }
    title(lineages[k], font.main = 1)
    if(k==1) {
      label.figure("A", font = 2, cex = 1.5, yfrac = 0.96)
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
    pout <- plotEZ(studies[j], show = "lm", lwd = 2, col.line = col[j], add = TRUE)
    text(pout$Eh7lim[2] + dx[[3]][j], pout$ZC[2] + dy[[3]][j], j)
  }
  title("Bacteria and Archaea", font.main = 1)

  # Scatterplots and fits for Bacteria and Archaea in all environments incluing hot springs 20211009
  # Read Eh7 - ZC data
  dat <- read.csv(system.file("extdata/orp16S/EZdat.csv", package = "JMDplots"))
  # Loop over Bacteria and Archaea
  for(k in 1:2) {
    par(mar = c(4, 4, 1, 1))
    plot(c(-500, 650), range(dat$ZC), type = "n", xlab = "Eh7 (mV)", ylab = cplab$ZC)
    # Use Bacteria only
    thisdat <- dat[dat$lineage == lineages[k], ]
    # Add linear fit; include number of studies in legend 20210925
    nstudy <- length(unique(thisdat$study))
    add.linear(thisdat$Eh7, thisdat$ZC, nstudy, pvalue_upper_right = TRUE)
    # Add points - highlight hot springs in red
    ienv <- c(1, 2, 5, 3, 6, 7, 4)
    # Use red for hot spring, gray for other environments
    cols <- rep("slategray2", length(ienv))
    cols[4] <- orp16Scol[4]
    eachenv(lineages[k], add = TRUE, do.linear = FALSE, ienv = ienv, cols = cols)
    title(lineages[k], font.main = 1)
    if(k == 1) {
      label.figure("B", font = 2, cex = 1.5, yfrac = 0.96)
      legend("bottomright", c("Hot Spring", "Other environments"), col = c(orp16Scol[4], "slategray2"), pch = 19, bty = "n", cex = 0.9)
    }
  }

  if(pdf) dev.off()

}

# Figure S1: ZC-Eh scatterplots for all studies 20210827
# This also creates files EZdat (Eh and ZC values) and
# EZlm (linear fits) for use by other plotting functions
orp16S_S1 <- function(pdf = FALSE) {

  if(!grepl("orp16S", options("chem16Sdir")[[1]])) 
    stop('Please run this first: options(chem16Sdir = system.file("extdata/orp16S", package = "JMDplots"))')

  # Setup figure
  if(pdf) pdf("Figure_S1.pdf", width = 12, height = 9)
  par(mfrow = c(3, 4))

  results <- c(

    message("\nRiver & Seawater"),
    plotEZ("MLL+18", "Bacteria", groupby = "Season", groups = c("Summer", "Winter"), legend.x = "bottomright"),
    plotEZ("SVH+19", "two", groupby = "Type", groups = c("Oxic", "Suboxic", "Euxinic")),
    plotEZ("HXZ+20", "Bacteria", groupby = "Station", groups = c("SYBL", "C4")),
    plotEZ("KLY+20", "Bacteria", groupby = "Year", groups = c(2018, 2019), legend.x = "bottomright"),
    plotEZ("WHL+21", "Bacteria", groupby = "Season", groups = c("Spring", "Summer", "Autumn", "Winter"), legend.x = "bottomleft"),
    plotEZ("LXH+20", "Bacteria", groupby = "Season", groups = c("Summer", "Winter"), legend.x = "bottomright"),
    plotEZ("JVW+20", "Bacteria", groupby = "isolation_source", groups = c("Ulva laetevirens", "lagoon water"), legend.x = "topright"),
    plotEZ("ZZL+21", "Bacteria", groupby = "Location", groups = c("Main Stream", "Animal Farm", "Hospital", "WWTP", "Tributary"), legend.x = "bottomright"),
    plotEZ("GZL21", "Bacteria", groupby = "Type", groups = c("Surface water", "Middle water", "Bottom water"), legend.x = "bottomleft"),

    message("\nLake & Pond"),
    plotEZ("SAR+13", "two", groupby = "Zone", groups = c("Photic-oxic", "Transition", "Anoxic")),
    plotEZ("LZR+17", "Bacteria", groupby = "Elevation", groups = c("< 1000 m", "1000 - 4000 m", "> 4000 m"), legend.x = "bottomleft"),
    plotEZ("ECS+18", "Bacteria", groupby = "Lake", groups = c("Laguna Negra", "Lo Encanado")),
    plotEZ("LLC+19", "Bacteria", groupby = "Size", groups = c("Free-living", "Particle-associated")),
    plotEZ("SCH+16", "two", groupby = "Type", groups = c("Oxic", "Oxycline", "Anoxic")),
    plotEZ("BCA+21", "Bacteria", groupby = "Month", groups = c("Jul", "Nov", "Feb", "Apr")),
    plotEZ("HLZ+18", "Bacteria", groupby = "Type", groups = c("Reservoir", "Pond"), legend.x = "bottomright"),
    plotEZ("CNA+20", "Bacteria", groupby = "Season", groups = c("Summer", "Autumn", "Winter", "Spring"), legend.x = "topright"),
    plotEZ("BWD+19", "Bacteria", groupby = "Cover", groups = c("Ice", "Ice Free"), legend.x = "bottomright"),
    plotEZ("RSJ+21", "two", groupby = "Lake", groups = c("Kuiva", "Lovo"), legend.x = "topright"),
    plotEZ("BOEM21", "Bacteria", groupby = "Stratum", groups = c("Upper", "Chemocline", "Lower")),
    plotEZ("GSY+20", "two", groupby = "Lake", groups = c("La Zarza", "Filon Centro"), legend.x = "bottomright"),
    plotEZ("NLE+21", "Bacteria", groupby = "Year", groups = c("2017", "2018"), legend.x = "bottomleft"),
    plotEZ("FAV+21", "two", groupby = "Type", groups = c("Oxic Surface", "Anoxic Surface", "Bottom")),
    plotEZ("GRG+20", "two", groupby = "Type", groups = c("Oxic", "Anoxic"), legend.x = "bottomleft"),

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
    plotEZ("OFY+19", "two", groupby = "Type", groups = c("Hot Spring", "Drain", "Lake")),

    message("\nAlkaline Spring"),
    plotEZ("SBP+20", "Bacteria", groupby = "pH Group", groups = c("< 10", "> 10")),
    plotEZ("RMB+17", "two", groupby = "pH Group", groups = c("< 10", "> 10")),
    plotEZ("CTS+17", "two", groupby = "Type", groups = c("River", "Well", "Spring"), legend.x = "bottomleft"),
    plotEZ("KSR+21", "Bacteria", groupby = "Location", groups = c("Lerone", "Branega", "Branega Creek Water")),
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

    message("\nSediment"),
    plotEZ("JHL+12", "two", groupby = "Core", groups = c("GC6", "GC12"), legend.x = "bottomright"),
    plotEZ("GFE+16", "Bacteria", groupby = "Station", groups = c(1, 4, 7, 18)),
    plotEZ("ZML+17", "two", groupby = "Type", groups = c("Mangrove Forest", "Intertidal Mudflat"), legend.x = "bottomleft"),
    plotEZ("BYB+17", "two", groupby = "Type", groups = c("Freshwater", "Brackish", "Saltmarsh", "Hypersaline"), legend.x = "bottomright"),
    plotEZ("BSPD17", "Bacteria", groupby = "Type", groups = c("Water", "Sediment")),
    plotEZ("HDZ+19", "Bacteria", groupby = "Type", groups = c("Water", "Sediment")),
    plotEZ("WHLH21", "Bacteria", groupby = "Position", groups = c("Surface", "Middle", "Bottom"), legend.x = "bottomleft"),
    plotEZ("SCM+18", "two", groupby = "Site", groups = c("Shallow", "Deep"), legend.x = "bottomleft"),
    plotEZ("CLS+19", "two", groupby = "Type", groups = c("Water", "Sediment"), legend.x = "bottomright"),
    plotEZ("ZDA+20", "Bacteria", groupby = "Type", groups = c("Water", "Sediment"), legend.x = "bottomleft"),
    plotEZ("VMB+19", "two", groupby = "Type", groups = c("Water", "Sediment"), legend.x = "topright"),
    plotEZ("WHC+19", "Archaea", groupby = "Type", groups = c("Mid-tide", "Low-tide", "Subtidal"), legend.x = "topright"),
    plotEZ("HSF+19", "Bacteria", groupby = "Type", groups = c("Water", "Sediment-Water Interface", "Sediment"), legend.x = "topright"),
    plotEZ("MCS+21", "two", groupby = "Type", groups = c("Water", "Sediment"), legend.x = "topright"),
    plotEZ("LMBA21_2017", "two", groupby = "Season", groups = c("Summer", "Winter"), legend.x = "bottomright"),
    plotEZ("ZZLL21", "Bacteria", groupby = "Type", groups = c("Main stream", "Animal farm", "Hospital", "WWTP", "Tributary"), legend.x = "bottomleft"),
    plotEZ("WFB+21", "Bacteria", groupby = "Treatment", groups = c("C. volutator", "H. diversicolor", "Cv & Hd", "MPB", "Manual turbation")),

    message("\nSoil"),
    plotEZ("SBW+17", "Bacteria", groupby = "Treatment", groups = c("Control", "BC400", "BC600")),
    plotEZ("MLL+19", "two", groupby = "Type", groups = c("Upland", "Paddy", "Sediment")),
    plotEZ("BMOB18", "two", groupby = "Treatment", groups = c("Acetate", "No amendment", "Pre-incubation")),
    plotEZ("ZZZ+18", "two", groupby = "Treatment", groups = c("None", "AQDS", "Biochar"), legend.x = "bottomleft"),
    plotEZ("PMM+20", "Bacteria", groupby = "Type", groups = c("Fluctuating", "Static")),
    plotEZ("ZHZ+19", "two", groupby = "Treatment", groups = c("Original", "Nitrate-reducing", "Ferric-reducing", "Sulfate-reducing", "Methanogenic")),
    plotEZ("CWC+20", "Bacteria", groupby = "Management", groups = c("Flooding", "Draining")),
    plotEZ("PSG+20", "two", groupby = "Treatment", groups = c("Initial", "NCC", "RB", "RGP", "TP")),
    plotEZ("XLD+20", "two", groupby = "Treatment", groups = c("Original Soil", "CK", "9K-Control", "Htt-sys", "Att-sys", "Co-sys"), legend.x = "bottomright"),
    plotEZ("LJC+20", "two", groupby = "meanT", groups = c("MAT >= 21.5 degC", "MAT < 21.5 degC")),
    plotEZ("DTJ+20", "two", groupby = "Zone", groups = c("Bulk Soil", "Mature", "Elongation", "Tip")),
    plotEZ("LLL+21", "Bacteria", groupby = "Treatment", groups = c("CK", "FL", "EA", "SB", "BD")),
    plotEZ("DLS21", "Bacteria", groupby = "Treatment", groups = c("control", "char", "silicate", "husk"))

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


############################
### UNEXPORTED FUNCTIONS ###
############################

add.linear <- function(Eh7, ZC, nstudy = NA, pvalue_upper_right = FALSE) {
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

  # Calculate Pearson correlation 20211009
  pearson <- cor.test(Eh7, ZC)
  # Get p-value stars
  pval <- pearson$p.value
  stars <- ""
  if(pval < 1e-2) stars <- "*"
  if(pval < 1e-5) stars <- "**"
  if(pval < 1e-20) stars <- "***"
  if(pval < 1e-30) stars <- "****"
  # Format correlation coefficient and paste stars
  rtext <- formatC(pearson$estimate, digits = 2, format = "f")
  rtext <- paste0(rtext, stars)
  rtext <- bquote(italic(r) == .(rtext))
  if(pvalue_upper_right) {
    legend("topright", legend = rtext, bty = "n", text.col = text.col, inset = c(0.02, 0.2))
  } else {
    legend("bottomleft", legend = rtext, bty = "n", text.col = text.col)
  }
}

# Scatterplots for all samples in each environment type 20210913
# NOTE: default ienv omits Hot Spring (4)
eachenv <- function(lineage = "Bacteria", add = FALSE, do.linear = TRUE, ienv = c(1, 2, 5, 3, 6, 7), cols = orp16Scol) {
  # Read Eh7 - ZC data
  dat <- read.csv(system.file("extdata/orp16S/EZdat.csv", package = "JMDplots"))
  # Get overall x and y limits
  xlim <- range(dat$Eh7)
  ylim <- range(dat$ZC)
  # Use Bacteria or Archaea only
  if(!is.null(lineage)) dat <- dat[dat$lineage == lineage, ]
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
    points(Eh7, ZC, pch = 19, cex = 0.2, col = cols[i])
  }
}

