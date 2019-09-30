# gradox/MG.R
# prepare metagenomic DNA data files and make plots of ZC
# 20180214 first version (jmd in chiang mai, thailand)
# 20180831 updated for paper submission
# 20181215 updated for paper revision

## requires CHNOSZ 1.1.3 and R 3.6.0 (for 'gap.axis' argument in axis())
#library(CHNOSZ)
#data(thermo)

# usage:
# mprep("BalticSea_Sediment", "SRA_MTD")  # sample DNA compositions
# mprep("BalticSea_Sediment", "SRA_MTR")  # sample RNA compositions
# mplot("BalticSea_Sediment", "SRA_MT")   # plot DNA and RNA ZC
# mprep("BalticSea_Sediment", "SRA_MTP")  # sample protein compositions
# mplot("BalticSea_Sediment", "SRA_MTP")  # plot protein ZC
# mprep(...)         # as above, for all other datasets
# mout <- mpage()    # plot DNA and RNA ZC for all datasets
# pout <- ppage()    # plot protein ZC for all datasets
# mcomp(mout)        # compare DNA ZC with RNA ZC or GC content
# mcomp(mout, pout)  # compare DNA ZC with protein ZC

# usage for sequences extracted for species classified by Kraken
# mprep("Menez_Gwen", "SRA_MGD", taxid=1427364)
# mprep("Menez_Gwen", "SRA_MGD", taxid=1420917)
# mplot("Menez_Gwen", "SRA_MG", taxid=c(1427364, 1420917))
# run mprep() for all of the species shown in Fig. 5
# mprep.species()

# naming conventions:
# dataset = study_seqtype
# study = basestudy[|-substudy]
# basestudy = name1_name2
# seqtype = [SRA|IMG|MGRAST|GenBank]_[MG|MT]
# append [D|R|P] to seqtype to indicate processed data (DNA, RNA, protein)

# identify the datasets used for the figure in the paper
# (all datasets are shown in the Supporting Information)
#usedinpaper <- c("BalticSea_Sediment", "Diffuse_Vents", "ETNP_OMZ", "Bison_Pool",
#  "Shimokita_Peninsula", "Menez_Gwen", "ETSP_OMZ", "Organic_Lake")
usedinpaper <- c("BalticSea_Sediment", "Shimokita_Peninsula",  # sediment
                 "Diffuse_Vents", "Menez_Gwen",                # hydrothermal
                 "ETNP_OMZ", "ETSP_OMZ",                       # ocean
                 "Mono_Lake", "Organic_Lake",                  # hypersaline
                 "Bison_Pool", "Guerrero_Negro")               # microbial mat
# mfrow setting for mpage() and ppage() (plot all ZC comparisons on a single page)
mfrow <- c(6, 3)
mfrow.usedinpaper <- c(2, 4)

# groups: hydrothermal, rock, mat, marine, hypersaline 20180315
# added OMZ, plume 20180321
studies <- list(
  BalticSea_Sediment = list( # 20180324 [ZMB+17, TFLS16], 20180608 [MKNJ18]
    `depth (mbsf)` = c("10cm", "12m", "15m", "42m", "10cm", "11m", "41m", "47m", "67m", "81m"),
    xlabels = c(0.1, 12, 15, 42, 0.1, 11, 41, 47, 67, 81),
    SRA_MG = c(NA, NA, NA, NA, "ERX509195", "SRX1537899", "SRX1537877", "SRX1537900", "SRX1537878", "SRX1537882"),
    SRA_MT = c("ERX509194", "SRX2869655", "SRX2869654", "SRX2869653", NA, NA, NA, NA, NA, NA),
    MT_range = c(0.615, 0.595),
    MTP_range = c(-0.16, -0.2),
    MG_range = c(0.635, 0.605),
    MGP_range = c(-0.12, -0.18),
    MG_srange = c(0.65, 0.6),
    abbrev = "BS",
    group = c("sediment0", "sediment", "sediment", "sediment", "sediment0", "sediment", "sediment", "sediment", "sediment", "sediment"),
    dx = list(MG = 0.0185, MT = 0.003),
    dy = list(MG = -0.003, MT = -0.003),
    plot_real_x = TRUE
  ),
  Bison_Pool = list( # 20180215 [DS11, HRM+11]
    `distance from source (m)` = c("N", "S", "R", "Q", "P"),
    xlabels = c(0, 6, 11, 14, 22),
    IMG_MG = c(2009439003, 2009439000, 2010170001, 2010170002, 2010170003),
    MG_range = c(0.61, 0.645),
    MGP_range = c(-0.22, -0.14),
    abbrev = "BP",
    group = c("yellowstone", "yellowstone", "yellowstone", "yellowstone1", "yellowstone1"),
    dx = list(MG = 0.0023),
    dy = list(MG = 0.0019),
    plot_real_x = TRUE
  ),
  Diffuse_Vents = list( # 20180320 [RRM+16, FLBH18]
    sample = c("Ginger_Castle", "Shrimp_Gulley_2", "Marker_33", "Marker_113", "Anemone", "Plume", "Background"),
    xlabels = c("GC", "SG", "33", "113", "An", "Pl", "Bd"),
    SRA_MG = c("ERX947644", "ERX947703", "ERX2080779", "ERX2080783", "ERX2080781", "ERX2080787", "ERX2080789"),
    SRA_MT = c(NA, NA, "ERX2080780", "ERX2080784", "ERX2080782", "ERX2080788", "ERX2080790"),
    MG_range = c(0.575, 0.62),
    MG_srange = c(0.57, 0.64),
    MT_range = c(0.60, 0.62),
    MGP_range = c(-0.2, -0.14),
    MTP_range = c(-0.18, -0.13),
    abbrev = "DV",
    group = c("vent", "vent", "vent", "vent", "vent", "plume", "plume0"),
    dx = list(MG = -0.030, MT = 0.005),
    dy = list(MG = -0.049, MT = -0.001),
    plot_real_x = FALSE
  ),
  ETNP_OMZ = list( # 20180301 [GKG+15] (MG), [GBL+15] (size-fractionated MT)
    `depth (m)` = c("30m", "85m", "100m", "125m", "300m"),
    xlabels = c(30, 85, 100, 125, 300),
    SRA_MG = c("SRX648440", "SRX648442", "SRX648443", "SRX648444", "SRX648445"),
    # sequences for small size fraction
    SRA_MT = c("SRX854074", "SRX854078", "SRX854081", "SRX854083", "SRX854084"),
    ## sequences for large size fraction
    #SRA_MT = c("SRX895490", "SRX895491", "SRX895493", "SRX895494", "SRX895496"),
    MG_range = c(0.615, 0.595),
    MG_srange = c(0.61, 0.57),
    MT_range = c(0.615, 0.60),
    MGP_range = c(-0.14, -0.165),
    MTP_range = c(-0.15, -0.165),
    abbrev = "EN",
    group = c("OMZ0", "OMZ", "OMZ", "OMZ", "OMZ"),
    dx = list(MG = -0.0025, MT = -0.001),
    dy = list(MG = 0, MT = 0.0035),
    plot_real_x = TRUE
  ),
  ETSP_OMZ = list( # 20180218 [SUD12]
    `depth (m)` = c("15m", "50m", "65m", "85m", "110m", "200m", "500m", "800m"),
    xlabels = c(15, 50, 65, 85, 110, 200, 500, 800),
    SRA_MG = c("SRX080962", "SRX025906", "SRX080938", "SRX025908", "SRX025910", "SRX025912", "SRX080950", "SRX080951"),
    MG_range = c(0.615, 0.58),
    MG_srange = c(0.61, 0.57),
    MGP_range = c(-0.16, -0.18),
    abbrev = "ES",
    group = c("OMZ0", "OMZ0", "OMZ", "OMZ", "OMZ", "OMZ", "OMZ", "OMZ"),
    plot_real_x = TRUE
  ),
  Guerrero_Negro = list( # 20180214 [KRH+08]
    layer = c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10"),
    xlabels = 1:10,
    IMG_MG = c(2004247000, 2004247001, 2004247002, 2004247003, 2004247004, 2004247005, 2004247006, 2004247007, 2004247008, 2004247009),
    MG_range = c(0.64, 0.62),
    MGP_range = c(-0.12, -0.15),
    abbrev = "GN",
    group = c("mat1", "mat1", "mat1", "mat", "mat", "mat", "mat", "mat", "mat", "mat"),
    dx = list(MG = -0.0038),
    dy = list(MG = 0.002),
    plot_real_x = FALSE
  ),
  HOT_ALOHA = list( # 20180312 [STED11]
    `depth (m)` = c("25m", "75m", "125m", "500m"),
    xlabels = c(25, 75, 125, 500),
    SRA_MG = c("SRX002155", "SRX000174", "SRX002157", "SRX002159"),
    SRA_MT = c("SRX002156", "SRX000175", "SRX002158", "SRX002160"),
    MG_range = c(0.62, 0.58),
    MT_range = c(0.62, 0.58),
    MGP_range = c(-0.15, -0.21),
    MTP_range = c(-0.17, -0.22),
    abbrev = "HA",
    group = c("oxic0", "oxic", "oxic", "oxic"),
    dx = list(MG = -0.019, MT = 0.002),
    dy = list(MG = -0.021, MT = -0.003),
    plot_real_x = TRUE
  ),
  `HOT_ALOHA-2010` = list( # 20181117 [MBA+17]
    `depth (m)` = c("25m", "45m", "75m", "125m", "200m", "500m", "770m", "1000m"),
    xlabels = c(25, 45, 75, 125, 200, 500, 770, 1000),
    SRA_MG = c("SRX2334639", "SRX2334707", "SRX2334640", "SRX2334668", "SRX2334697", "SRX2334652", "SRX2334669", "SRX2334641"),
    MG_srange = c(0.61, 0.57),
    abbrev = "HA",
    group = c("oxic0", "oxic", "oxic", "oxic", "oxic", "oxic", "oxic", "oxic"),
    plot_real_x = TRUE
  ),
  Menez_Gwen = list( # 20180407 [MBG+16]
    "distance from crack (cm)" = c("WdCr-f2", "Wd-10L" ,"Wd-40UP"),
    xlabels = c(0, 10, 40),
    SRA_MG = c("ERX1158314", "ERX1158315", "ERX1158316"),
    MG_range = c(0.6025, 0.6225),
    MG_srange = c(0.57, 0.64),
    MGP_range = c(-0.16, -0.13),
    abbrev = "MZ",
    group = c("vent", "vent", "plume"),
    dx = list(MG = -0.0022),
    dy = list(MG = 0.0025),
    plot_real_x = TRUE
  ),
  Mono_Lake = list( # 20180216 [EH17]
    `depth (m)` = c("10m", "15m", "18m", "25m", "31m"),
    xlabels = c(10, 15, 18, 25, 31),
    # SRA experiments for "A" samples
    SRA_MT = c("SRX1527166", "SRX1527805", "SRX1527807", "SRX1527809", "SRX1527811"),
    ## SRA experiments for "B" samples
    #SRA_MT = c("SRX1527804", "SRX1527806", "SRX1527808", "SRX1527810", "SRX1527812"),
    MT_range = c(0.625, 0.605),
    MTP_range = c(-0.15, -0.175),
    abbrev = "ML",
    group = c("hypersaline0", "hypersaline", "hypersaline", "hypersaline", "hypersaline"),
    dx = list(MT = 0.0035),
    dy = list(MT = 0.001),
    plot_real_x = TRUE
  ),
  Mud_Volcano = list( # 20180320 [CCT+12]
    `depth (cm)` = c("3cm", "13cm", "23cm", "31cm"),
    xlabels = c(3, 13, 23, 31),
    # SRA experiments corresponding to run IDs given in the paper (SRR389129, SRR389128, SRR389127, SRR389126)
    SRA_MG = c("SRX111280", "SRX111279", "SRX111278", "SRX111277"),
    MG_range = c(0.62, 0.60),
    MGP_range = c(-0.14, -0.19),
    abbrev = "MV",
    group = c("rock0", "rock", "rock", "rock"),
    dy = list(MG = -0.017),
    dx = list(MG = -0.0057),
    plot_real_x = TRUE
  ),
  Organic_Lake = list( # 20180215 [YLW+13]
    `depth (m)` = c("1.7m", "4.2m", "5.7m", "6.5m", "6.7m"),
    xlabels = c(1.7, 4.2, 5.7, 6.5, 6.7),
    # SRA experiments for samples GS374, GS375, GS376, GS377, GS378 (0.8 um size fraction)
    SRA_MG = c("SRX024781", "SRX024782", "SRX024783", "SRX024784", "SRX024843"),
    MG_range = c(0.625, 0.61),
    MG_srange = c(0.64, 0.57),
    MGP_range = c(-0.125, -0.155),
    abbrev = "OL",
    group = c("hypersaline0", "hypersaline", "hypersaline", "hypersaline", "hypersaline"),
    dx = list(MG = 0.0023),
    dy = list(MG = 0.0018),
    plot_real_x = TRUE
  ),
  Serpentinite_Springs = list( # 20180318 [BB09, BNS12]
    sample = c("LostCity", "WHC2B", "TLE"),
    MGRAST_MG = c("4461585.3", "4460690.3", "4460689.3"),
    MG_range = c(0.58, 0.65),
    MGP_range = c(-0.2, -0.08),
    abbrev = "SS",
    group = c("vent", "rock", "rock0"),
    dx = list(MG=-0.0025),
    dy = list(MG=0),
    plot_real_x = FALSE
  ),
  Shimokita_Peninsula = list( # 20180613 [KFT+14]
    `depth (mbsf)` = c("0.8", "5.1", "18.6", "48.5", "107.0"),
    xlabels = c(0.8, 5.1, 18.6, 48.5, 107.0),
    GenBank_MG = c("BARS", "BART", "BARU", "BARV", "BARW"),
    MG_range = c(0.62, 0.58),
    MGP_range = c(-0.15, -0.21),
    abbrev = "SP",
    group = c("sediment0", "sediment", "sediment", "sediment", "sediment"),
    dx = list(MG=-0.013),
    dy = list(MG=-0.032),
    plot_real_x = TRUE
  ),
  Yellowstone_Park = list( # 20180310 [IJT+13]
    samples = c("Ar-red", "Aq-red", "Aq", "Ph-red", "Ph"),
    xlabels = c("Ar(red)  ", " Aq(red)", "Aq", "Ph(red)", "Ph"),
    # Archaeal-reducing: Sites 01, 02, 03, 19, 04, 18
    # IMG IDs (Celera): 2022920009, 2022920014, 2014031003, 2022920017, 2022920008, 2022920019
    #   NOTE: 2022920002 for Site 03 is indexed but not present on IMG server, so we use the ID for PGA assembly given by IJT+13
    # Aquificales-reducing: Sites 09, 10, 12; IMG IDs: 2022920010, 2022920015, 2022920011
    # Aquificales: Sites 14, 11, 13; IMG IDs: 2022920007, 2022920012, 2022920006
    # Phototroph-reducing: Sites 05, 20; IMG IDs: 2022920003, 2022920020
    # Phototroph: Sites 06, 07, 15, 16; IMG IDs: 2022920004, 2022920013, 2022920016, 2022920018
    IMG_MG = c("Archaeal-reducing", "Aquificales-reducing", "Aquificales", "Phototroph-reducing", "Phototroph"),
    MG_range = c(0.59, 0.63),
    MGP_range = c(-0.23, -0.16),
    abbrev = "YP",
    group = c("yellowstone", "yellowstone", "yellowstone", "yellowstone1", "yellowstone1"),
    dy = c(MG = -0.004),
    plot_real_x = FALSE
  )
)

# get metadata (location names and sequencing IDs) for a study
# extracted from mprep/mplot 20180312
mdata <- function(study, seqtype) {
  samples <- studies[[study]][[1]]
  if(is.null(samples)) stop("metadata for ", study, " study not available")
  xlabels <- studies[[study]]$xlabels
  if(is.null(xlabels)) xlabels <- samples
  group <- studies[[study]][["group"]]
  seqtype.for.ID <- seqtype
  # change e.g. SRA_MGP to SRA_MG
  seqtype.for.ID <- gsub("_MG.$", "_MG", seqtype.for.ID)
  seqtype.for.ID <- gsub("_MT.$", "_MT", seqtype.for.ID)
  IDs <- studies[[study]][[seqtype.for.ID]]
  if(is.null(IDs)) stop(seqtype.for.ID, " IDs not available for ", study, " study")
  # remove NA IDs and corresponding samples, xlabels, and groups
  samples[is.na(IDs)] <- NA
  samples <- na.omit(samples)
  xlabels[is.na(IDs)] <- NA
  xlabels <- na.omit(xlabels)
  if(length(group) > 1) {
    group[is.na(IDs)] <- NA
    group <- na.omit(group)
  }
  IDs <- na.omit(IDs)
  abbrev <- studies[[study]][["abbrev"]]
  return(list(samples=samples, xlabels=xlabels, IDs=IDs, group=group, abbrev=abbrev))
}

# make page of plots for MG/MT 20180225
mpage <- function(plottype="bars", dsDNA=TRUE, forSI=FALSE, set.par=TRUE) {
  if(!forSI) mfrow <- mfrow.usedinpaper
  if(set.par) par(mfrow=mfrow, mar=c(4, 3.5, 2, 1), mgp=c(2.5, 1, 0))
  # initialize output of ZC values
  mout <- list()
  # extract only the datasets used for the paper
  if(!forSI) studies <- studies[match(usedinpaper, names(studies))]
  iletter <- 1
  for(i in 1:length(studies)) {
    # name of the study without seqtype (e.g. Guerrero_Negro)
    study <- names(studies[i])
    # skip HOT_ALOHA-2010
    if(grepl("-", study)) next
    for(j in 1:length(studies[[i]])) {
      # the seqtype (e.g. IMG_MG)
      seqtype <- names(studies[[i]][j])
      # the seqtype isn't the samples (first position) or the xlabels, ZC range, abbreviation, group, dx, or dx
      if(j==1 | grepl("xlabels", seqtype) | grepl("range", seqtype) |
         seqtype=="abbrev" | seqtype=="group" | seqtype=="dx" | seqtype=="dy" | seqtype=="plot_real_x") next
      # for the figure in the paper, take only one dataset for each study (i.e. metagenome, except for Mono Lake)
      if(!forSI & (grepl("_MT", seqtype) & study!="Mono_Lake")) next
      ZC <- list(mplot(study, seqtype, plottype, dsDNA=dsDNA, add.label=forSI))
      names(ZC) <- paste0(study, "_", seqtype)
      mout <- c(mout, ZC)
    }
    # add figure label 20181210
    if(!forSI) {
      label.figure(LETTERS[iletter], cex=1.6, font=2, yfrac=0.936)
      iletter <- iletter + 2
    }
  }
  mout
}

# scatterplot of DNA and RNA ZC 20180307
# mout <- mpage(); mcomp(mout)
# yvar can be RNA or GC
mcomp <- function(mout, yvar="RNA") {
  # set up plot
  if(yvar=="RNA") plot(c(0.57, 0.65), c(0.22, 0.32), xlab="ZC(DNA)", ylab="DZC(RNA - DNA)", type="n")
  if(yvar=="GC") plot(c(0.57, 0.65), c(0.3, 0.65), xlab="ZC(DNA)", ylab="GC", type="n")
  for(i in 1:length(mout)) {
    DNA <- mout[[i]]$DNA
    if(yvar=="RNA") Y <- mout[[i]]$RNA
    if(yvar=="GC") Y <- mout[[i]]$GC
    # order points by increasing DNA value
    ord <- order(DNA)
    DNA <- DNA[ord]
    Y <- Y[ord]
    # for RNA, plot difference from DNA
    if(yvar=="RNA") Y <- Y - DNA
    # color: red for MG, blue for MT
    if(grepl("_MG", names(mout[i]))) col <- "red"
    if(grepl("_MT", names(mout[i]))) col <- "blue"
    # plot lines
    lines(DNA, Y, col=col)
  }
}

# make page of plots for MGP/MTP 20180225
ppage <- function(forSI=FALSE, set.par=TRUE) {
  if(!forSI) mfrow <- mfrow.usedinpaper
  if(set.par) par(mfrow=mfrow, mar=c(4, 3.5, 2, 1), mgp=c(2.5, 1, 0))
  # initialize output of ZC values
  pout <- list()
  # extract only the datasets used for the paper
  if(!forSI) studies <- studies[match(usedinpaper, names(studies))]
  iletter <- 2
  for(i in 1:length(studies)) {
    # name of the study without seqtype (e.g. Columbia_River)
    study <- names(studies[i])
    # we only handle studies with proteins (MGP or MTP)
    if(!any(c("MGP_range", "MTP_range") %in% names(studies[[i]]))) next
    for(j in 1:length(studies[[i]])) {
      # the seqtype
      seqtype <- names(studies[[i]][j])
      # the seqtype isn't the samples (first position) or the xlabels, ZC range, abbreviation, group, dx, or dx
      if(j==1 | grepl("xlabels", seqtype) | grepl("range", seqtype) |
         seqtype=="abbrev" | seqtype=="group" | seqtype=="dx" | seqtype=="dy" | seqtype=="plot_real_x") next
      # for the figure in the paper, take only one dataset for each study (i.e. metagenome, except for Mono Lake)
      if(!forSI & (grepl("_MT", seqtype) & study!="Mono_Lake")) next
      # add "P" for proteins
      seqtype <- paste0(seqtype, "P")
      ZC <- list(mplot(study, seqtype, add.label=forSI))
      names(ZC) <- paste0(study, "_", seqtype)
      pout <- c(pout, ZC)
    }
    # add figure label 20181210
    if(!forSI) {
      label.figure(LETTERS[iletter], cex=1.6, font=2, yfrac=0.936)
      iletter <- iletter + 2
    }
  }
  pout
}

# scatterplot of DNA and protein ZC 20180307
# mout <- mpage(); pout <- ppage(); pcomp(mout, pout)
pcomp <- function(mout, pout, seqtype="MG", xaxis="DNA", yaxis="AA", parts=c("plot", "legend")) {
  if("plot" %in% parts) {
    # set up plot
    if(xaxis=="DNA") {
      if(seqtype=="MG") xlim <- c(0.58, 0.65)
      if(seqtype=="MT") xlim <- c(0.585, 0.63)
      xlab <- quote(italic(Z)[C]~of~DNA)
      dx <- -0.004
    }
    if(xaxis=="RNA") {xlim <- c(1.85, 1.92); xlab <- "ZC(RNA)"; dx <- -0.004}
    if(xaxis=="CM") {xlim <- c(0.024, 0.037); xlab <- "Cys+Met fraction"; dx <- -0.0003}
    if(yaxis=="AA") {
      if(seqtype=="MG") ylim <- c(-0.22, -0.098)
      if(seqtype=="MT") ylim <- c(-0.21, -0.14)
      ylab <- quote(italic(Z)[C]~of~protein)
    }
    if(yaxis=="CM") {ylim <- c(0.024, 0.037); ylab <- "Cys+Met fraction"}
    plot(xlim, ylim, xlab=xlab, ylab=NA, type="n")
    mtext(ylab, side=2, line=3.5, las=0, cex=par("cex"))
    ## show MT range on MG plot
    #if(seqtype=="MG") {
    #  xs <- extendrange(c(0.595, 0.64))
    #  ys <- extendrange(c(-0.193, -0.117))
    #  rect(xs[1], ys[1], xs[2], ys[2], lty=3)
    #}
    for(i in 1:length(pout)) {
      # get study name and change MGP or MTP to MG or MT
      study <- names(pout[i])
      study <- gsub("MGP", "MG", study)
      study <- gsub("MTP", "MT", study)
      # select MG or MT
      if(seqtype=="MG" & !grepl("_MG$", study)) next
      if(seqtype=="MT" & !grepl("_MT$", study)) next
      # find this study in mout
      imout <- match(study, names(mout))
      # assemble ZC values and get environment group
      if(yaxis=="AA") yvals <- pout[[i]]$AA
      if(yaxis=="CM") yvals <- pout[[i]]$CM
      if(xaxis=="DNA") xvals <- mout[[imout]]$DNA
      if(xaxis=="RNA") xvals <- mout[[imout]]$RNA
      if(xaxis=="CM") xvals <- pout[[i]]$CM
      group <- mout[[imout]]$group
      # color: by group
      col <- ""
      col[group %in% c("yellowstone", "yellowstone1")] <- "orange"
      col[group %in% c("rock", "rock0")] <- "brown"
      col[group %in% c("vent", "vent0")] <- "red"
      col[group %in% c("mat", "mat1")] <- "green3"
      col[group %in% c("oxic", "oxic0")] <- "blue"
      col[group %in% c("hypersaline", "hypersaline0")] <- "turquoise3"
      col[group %in% c("OMZ", "OMZ0")] <- "black"
      col[group %in% c("plume", "plume0")] <- "purple1"
      col[group %in% c("sediment", "sediment0")] <- "slategrey"
      # order points by increasing DNA/RNA ZC value
      ord <- order(xvals)
      xvals <- xvals[ord]
      yvals <- yvals[ord]
      group <- group[ord]
      # plot lines and points
      # for lines, use the color of most of the points 20180501
      # change this to gray 20181114
      lines(xvals, yvals, col="dimgray", lwd=0.8, lty=2)
      # filled circles for marine, filled squares for terrestrial 20181114
      pch <- rep(19, length(xvals))
      pch[group %in% c("yellowstone", "yellowstone1", "rock", "rock0", "mat", "mat1", "hypersaline", "hypersaline0")] <- 15
      points(xvals, yvals, col=col, pch=pch)
      # outline circle or square for groups with "0" (upper part of OMZ etc)
      pch0 <- rep(1, length(xvals))
      pch0[group %in% c("yellowstone1", "rock0", "mat1", "hypersaline0")] <- 0
      cex0 <- rep(1.8, length(xvals))
      cex0[group %in% c("yellowstone1", "rock0", "mat1", "hypersaline0")] <- 1.6
      i0 <- grepl("0", group)
      if(any(i0)) points(xvals[i0], yvals[i0], pch=pch0[i0], col=col[i0], cex=cex0[i0])
      # outline green circle or square for groups with "1" (phototrophic mat in Yellowstone / upper part of Guerrero Negro mat)
      i1 <- grepl("1", group)
      if(any(i1)) points(xvals[i1], yvals[i1], pch=pch0[i1], col="green3", cex=cex0[i1])
      # get abbreviation for this study
      studyname <- paste(strsplit(study, "_")[[1]][1:2], collapse="_")
      istudy <- match(studyname, names(studies))
      abbrev <- studies[[istudy]][["abbrev"]]
      dx <- studies[[istudy]][["dx"]][[seqtype]]
      dy <- studies[[istudy]][["dy"]][[seqtype]]
      if(is.null(dx)) dx <- 0
      if(is.null(dy)) dy <- 0.003
      # add text label
      #text(xvals[1] + dx, yvals[1], abbrev, cex=0.7)
      imax <- which.max(yvals)
      text(xvals[imax] + dx, yvals[imax] + dy, abbrev, cex=0.7)
    }
  }
  # add legend
  if("legend" %in% parts) {
    if("plot" %in% parts) {
      marine.x <- "topleft"
      terrestrial.x <- "bottomright"
      bty <- "o"
    } else {
      marine.x <- "left"
      terrestrial.x <- "right"
      bty <- "n"
    }
    # marine
    legend(-0.063, 1.014, c("", "", "", "", "", "", "", "", ""), lty=2, bty="n", col=c(NA, rep("dimgray", 8)))
    legend(marine.x, lty=2, lwd=NA, bty=bty, pch=19,
           legend=c("MARINE", "OMZ", "ocean surface", "oxic", "vent", "plume", "background SW", "sediment", "sediment surface"),
           col=c(NA, "black", "black", "blue", "red", "purple1", "purple1", "slategrey", "slategrey"))
    # overlay ocean surface, seawater, and sediment surface points (blue circles)
    legend(marine.x, lty=0, lwd=0, bty=bty, pch=1, pt.cex=1.8, pt.lwd=1,
           legend=c("", "", "", "", "", "", "", "", ""),
           col=c(NA, NA, "black", NA, NA, NA, "purple1", NA, "slategrey"))
    # terrestrial
    legend(0.47, 1.014, c("", "", "", "", "", "", "", "", ""), lty=2, bty="n", col=c(NA, rep("dimgray", 8)))
    legend(terrestrial.x, lty=2, lwd=NA, pch=15,
           legend=c("TERRESTRIAL", "hypersaline", "lake surface", "mat", "upper mat", "rock-derived", "surface mixing", "Yellowstone", "phototrophic mat"),
           col=c(NA, "turquoise3", "turquoise3", "green3", "green3", "brown", "brown", "orange", "orange"), bty=bty)
    # overlay lake surface and phototrophic mat points (blue and green)
    legend(0.492, 1.014, lty=0, lwd=0, bty=bty, pch=0, pt.cex=1.6, pt.lwd=1,
           legend=c("", "", "", "", "", "", "", "", ""),
           col=c(NA, NA, "turquoise3", NA, "green3", NA, "brown", NA, "green3"))
  }
}

# function to plot ZC of metagenomic DNA 20180215
# or ZC of metagenomic proteins 20180228
# NOTE: dataset should end in "_MG" or "_MT" (plot DNA and RNA compositions)
# or "_MGP" or "_MTP" (plot protein compositions)
plotMG <- function(dataset="Guerrero_Negro_IMG_MG", plottype="bars",
  samples=formatC(10:1, width=2, flag="0"), labels=formatC(10:1, width=2, flag="0"),
  group="mat", xlab="layer", ylim=c(1.67, 1.77), abbrev=NULL, dsDNA=TRUE, plot.RNA=TRUE,
  taxid=NULL, lwd=1, lty=2, lwd.bars=2, col=NULL, extendrange=FALSE, add.label=TRUE,
  plot_real_x=FALSE, maxdepth=NULL) {
  # samples: (used for suffixes on file names)
  # labels: (used for labeling x-axis ticks)
  # xlab: "layer", ...
  isprotein <- grepl("_MGP$", dataset) | grepl("_MTP$", dataset)
  # where to keep mean and high/lo (+/- SD) of ZC at each site
  CM <- GC <- ZClo <- ZCmean <- ZChi <- numeric()
  # gradox data location in JMDplots package 20190928
  datadir <- system.file("extdata/gradox", package = "JMDplots")
  # set up for proteins or DNA
  if(isprotein) {
    filestart <- paste0(datadir, "/MGP/", dataset)
    ZCfun <- ZCAA
    col <- "darkgreen"
  } else {
    # data directory for DNA
    filestart <- paste0(datadir, "/MGD/", dataset, "D")
    ZCfun <- ZCnuc
    if(is.null(col)) col <- "red"
  }
  # keep sample values for violin plot 20180515
  DNA <- data.frame(ZC=numeric(length(samples)*100), sample=character(length(samples)*100), stringsAsFactors=FALSE)
  meancomp <- NULL
  for(i in 1:length(samples)) {
    sample <- samples[i]
    if(!is.null(taxid)) file <- paste0(filestart, "_", sample, "-", taxid, ".csv")
    else file <- paste0(filestart, "_", sample, ".csv")
    # use NA if the file is missing 20180529
    if(!file.exists(file)) {
      ZCmean <- c(ZCmean, NA)
      ZClo <- c(ZClo, NA)
      ZChi <- c(ZChi, NA)
      message(paste("missing file", file))
    } else {
      mycomp <- read.csv(file, as.is=TRUE)
      if(isprotein) {
        # calculate Cys+Met fraction 20180324
        CM <- c(CM, CMAA(mycomp))
      } else {
        # use base-paired (double-stranded) DNA
        if(dsDNA) mycomp <- make_dsDNA(mycomp)
        # calculate GC ratio 20180309
        if(!isprotein) GC <- c(GC, GCnuc(mycomp))
      }
      myZC <- ZCfun(mycomp, "deoxyribose")
      ZCmean <- c(ZCmean, mean(myZC))
      ZClo <- c(ZClo, mean(myZC) - sd(myZC))
      ZChi <- c(ZChi, mean(myZC) + sd(myZC))
      # initialize data frame for mean compositions
      if(is.null(meancomp)) {
        meancomp <- mycomp[1:length(samples), ]
        meancomp[] <- NA
        rownames(meancomp) <- samples
      }
      # get mean compositions 20180505
      meancomp[i, ] <- colMeans(mycomp)
      # keep sample values for violin plot 20180515
      istart <- (i-1) * 100 + 1
      iend <- istart + 99
      DNA$ZC[istart:iend] <- myZC
      DNA$sample[istart:iend] <- paste0("X", letters[i])
    }
  }
  # make plot
  isdeep <- logical(length(labels))
  if(plot_real_x) {
    if(!is.null(maxdepth)) isdeep <- labels > maxdepth
    xlim <- range(labels[!isdeep])
    if(extendrange) xlim <- extendrange(xlim)
    if(!grepl("Bison_Pool", dataset) & !grepl("Menez_Gwen", dataset)) xlim <- rev(xlim)
    at <- labels
  } else {
    nsamp <- length(samples)
    xlim <- c(1, nsamp)
    if(extendrange) xlim <- extendrange(xlim)
    if(plottype=="violin") xlim[1] <- xlim[1] - 0.5
    if(plottype=="violin") xlim[length(xlim)] <- xlim[length(xlim)] + 0.5
    at <- 1:nsamp
  }
  if(is.null(taxid) | identical(taxid, 0)) {
    plot(0, 0, xlim=xlim, ylim=ylim, xlab=xlab, ylab=NA, xaxt="n")
    # loop over non-numeric labels so that R doesn't omit any (because they're crowded)
    # and lower gap.axis to show more numeric labels 20181215 (requires R 3.6.0)
    if(!is.numeric(labels)) for(i in 1:nsamp) axis(1, at=i, labels=labels[i])
    else axis(1, at=at, labels=labels, gap.axis=0.02)
    mtext(quote(italic(Z)[C]), side=2, line=2, las=0)
  }
  if(plottype=="lines") {
    lines(at[!isdeep], ZClo[!isdeep], col=col, lty=3)
    lines(at[!isdeep], ZCmean[!isdeep], col=col)
    lines(at[!isdeep], ZChi[!isdeep], col=col, lty=3)
  }
  if(plottype=="bars") {
    # error bar plots 20180515
    # apply small offset to x-position to separate DNA and RNA
    if(!is.null(taxid)) dx <- 0
    else dx <- abs(diff(par("usr")[1:2])) / 120
    # bars (whiskers) at one SD from mean (arrows trick from https://stackoverflow.com/questions/13032777/scatter-plot-with-error-bars)
    arrows(at-dx, ZClo, at-dx, ZChi, length = 0.03, angle = 90, code = 3, col=col, lwd=lwd.bars)
    # remove NAs so we can draw lines between all sites 20180529
    iNA <- is.na(ZCmean)
    lines((at)[!iNA & !isdeep]-dx, ZCmean[!iNA & !isdeep], col=col, lty=lty, lwd=lwd)
  }
  # add points to show > 1% species abundance 20181118
  if(!is.null(taxid) & !identical(taxid, 0)) {
    # gradox data location in JMDplots package 20190928
    datadir <- system.file("extdata/gradox", package = "JMDplots")
    file <- paste0(datadir, "/one_percent/", dataset, ".csv")
    dat <- read.csv(file)
    # which samples have this taxid with at least 1% abundance?
    dat <- dat[dat$taxid==taxid, ]
    dat <- dat[dat$percentage >= 1, ]
    ioneperc <- samples %in% dat$sample
    points(at[ioneperc], ZCmean[ioneperc], pch=19, cex=0.5, col=col)
  }
  # now do it for RNA
  if(!isprotein) {
    RNA_ZClo <- RNA_ZCmean <- RNA_ZChi <- numeric()
    # keep sample values for violin plot 20180515
    RNA <- data.frame(ZC=numeric(length(samples)*100), sample=character(length(samples)*100), stringsAsFactors=FALSE)
    for(i in 1:length(samples)) {
      sample <- samples[i]
      file <- paste0("data/MGR/", dataset, "R_", sample, ".csv")
      if(file.exists(file)) {
        myRNA <- read.csv(file, as.is=TRUE)
        myZC <- ZCfun(myRNA, "ribose")
        RNA_ZCmean <- c(RNA_ZCmean, mean(myZC))
        RNA_ZClo <- c(RNA_ZClo, mean(myZC) - sd(myZC))
        RNA_ZChi <- c(RNA_ZChi, mean(myZC) + sd(myZC))
      }
      # keep sample values for violin plot 20180515
      istart <- (i-1) * 100 + 1
      iend <- istart + 99
      RNA$ZC[istart:iend] <- myZC
      RNA$sample[istart:iend] <- paste0("X", letters[i])
    }
    if(plot.RNA & length(RNA_ZCmean) > 0) {
      # offset to apply to ZC of RNA
      dZC <- -0.28
      if(plottype=="lines") {
        lines(at, RNA_ZClo + dZC, col="blue", lty=3)
        lines(at, RNA_ZCmean + dZC, col="blue")
        lines(at, RNA_ZChi + dZC, col="blue", lty=3)
      }
      if(plottype=="bars") {
        # apply small offset to x-position to separate DNA and RNA
        dx <- abs(diff(par("usr")[1:2])) / 200
        arrows(at+dx, RNA_ZClo + dZC, at+dx, RNA_ZChi + dZC, length = 0.03, angle = 90, code = 3, col="blue", lwd=lwd.bars)
        lines(at+dx, RNA_ZCmean + dZC, col="blue", lty=2, lwd=lwd)
      }
    }
    # return ZC values
    outval <- list(DNA=ZCmean, RNA=RNA_ZCmean, GC=GC, group=group, meancomp=meancomp, abbrev=abbrev)
  } else {
    outval <- list(AA=ZCmean, CM=CM, meancomp=meancomp)
  }
  # add title 20180225
  if(is.null(taxid) | identical(taxid, 0)) {
    main <- dataset2main(dataset, abbrev)
    title(main=main, font.main=1)
    # add in-plot label: MG or MT 20180829
    if(add.label) {
      label <- "MG"
      if(grepl("_MT.*", dataset)) label <- "MT"
      # change position for some datasets
      xfrac <- 0.1
      if(grepl("OMZ", dataset)) xfrac <- 0.9
      if(grepl("Negro", dataset)) xfrac <- 0.9
      if(grepl("ALOHA", dataset)) xfrac <- 0.9
      CHNOSZ::label.plot(label, xfrac=xfrac, yfrac=0.9)
    }
  }
  # done!
  return(outval)
}

# parse dataset name to plot title 20180505
dataset2main <- function(dataset, abbrev=NULL) {
  main <- dataset
  main <- gsub("_SRA", "", main)
  main <- gsub("_IMG", "", main)
  main <- gsub("_MGRAST", "", main)
  main <- gsub("_GenBank", "", main)
  # strip _MG or _MGP
  main <- gsub("_MG.*", "", main)
  main <- gsub("_MT.*", "", main)
  # replace underscore with space, and put space in Baltic Sea
  main <- gsub("_", " ", main)
  main <- gsub("BalticSea", "Baltic Sea", main)
  main <- gsub("Mud Volcano", "SYNH Mud Volcano", main)
  # add dataset abbreviation at end 20180829
  if(!is.null(abbrev)) main <- paste0(main, " (", abbrev, ")")
  main
}

### utility functions for ZC and GC calculations ###

# calculate GC ratio for given nucleobase compositions 20180309
GCnuc <- function(nuccomp) {
  # find columns with names for the nucleobases
  isbase <- colnames(nuccomp) %in% c("A", "C", "G", "T", "U")
  # keep only these columns
  nuccomp <- nuccomp[, isbase]
  # find GC columns and not-GC columns
  isGC <- colnames(nuccomp) %in% c("G", "C")
  isnotGC <- colnames(nuccomp) %in% c("A", "T", "U")
  # calculate GC ratio
  sum(nuccomp[, isGC]) / sum(nuccomp)
}

# calculate ZC for given nucleobase compositions 20171223
# changed to nucleosides 20180512
ZCnuc <- function(nuccomp, sugar="deoxyribose") {
  ## ZC and nC of nucleobases
  #ZC_nuc <- c(A=2, C=1.5, G=2.4, T=0.8, U=1.5)
  #nC_nuc <- c(A=5, C=4, G=5, T=5, U=4)
  # the number of carbons of the nucleosides (same in DNA and RNA)
  nC_nuc <- c(A=10, C=9, G=10, T=10, U=9)
  if(sugar=="deoxyribose") {
    # ZC of the nucleosides in DNA
    # DNAsides <- c("deoxyadenosine", "deoxycytidine", "deoxyguanosine", "deoxythymidine", "deoxyuridine")
    # ZC_nuc <- ZC(info(info(DNAsides))$formula)
    ZC_nuc <- c(A=0.8, C=4/9, G=1, T=0.2, U=4/9)
  }
  if(sugar=="ribose") {
    # ZC of the nucleosides in RNA
    # RNAsides <- c("adenosine", "cytidine", "guanosine", "thymidine", "uridine")
    # ZC_nuc <- ZC(info(info(RNAsides))$formula)
    ZC_nuc <- c(A=1, C=2/3, G=1.2, T=0.4, U=2/3)
  }
  # find columns with names for the nucleosides
  isbase <- colnames(nuccomp) %in% c("A", "C", "G", "T", "U")
  ibase <- match(colnames(nuccomp)[isbase], names(ZC_nuc))
  # calculate the nC from the counts
  multC <- t(t(nuccomp[, isbase]) * nC_nuc[ibase])
  # multiply nC by ZC to get total formal charge on carbon (Z)
  multZ <- t(t(multC) * ZC_nuc[ibase])
  # calculate the Z and nC in each composition (row)
  Ztot <- rowSums(multZ)
  nCtot <- rowSums(multC)
  # the ZC in each composition (row)
  Ztot / nCtot
}

# calculate ZC for amino acid compositions 20180228
ZCAA <- function(AAcomp, nothing=NULL) {
  # a dummy second argument is needed because of how this function is used in plotMG()
  # the number of carbons of the amino acids
  nC_AA <- c(Ala = 3, Cys = 3, Asp = 4, Glu = 5, Phe = 9, Gly = 2, His = 6, 
    Ile = 6, Lys = 6, Leu = 6, Met = 5, Asn = 4, Pro = 5, Gln = 5, 
    Arg = 6, Ser = 3, Thr = 4, Val = 5, Trp = 11, Tyr = 9)
  # the Ztot of the amino acids == CHNOSZ::ZC(info(info(aminoacids("")))$formula) * nC_AA
  Ztot_AA <- c(Ala = 0, Cys = 2, Asp = 4, Glu = 2, Phe = -4, Gly = 2, His = 4, 
    Ile = -6, Lys = -4, Leu = -6, Met = -2, Asn = 4, Pro = -2, Gln = 2, 
    Arg = 2, Ser = 2, Thr = 0, Val = -4, Trp = -2, Tyr = -2)
  # the ZC of the amino acids == CHNOSZ::ZC(info(info(aminoacids("")))$formula)
  ZC_AA <- Ztot_AA / nC_AA
  # find columns with names for the amino acids
  isAA <- colnames(AAcomp) %in% c("Ala", "Cys", "Asp", "Glu", "Phe", "Gly", "His", "Ile", "Lys", 
    "Leu", "Met", "Asn", "Pro", "Gln", "Arg", "Ser", "Thr", "Val", "Trp", "Tyr")
  iAA <- match(colnames(AAcomp)[isAA], names(ZC_AA))
  # calculate the nC for all occurrences of each amino acid
  multC <- t(t(AAcomp[, isAA]) * nC_AA[iAA])
  # multiply nC by ZC
  multZC <- t(t(multC) * ZC_AA[iAA])
  # calculate the total ZC and nC, then the overall ZC
  ZCtot <- rowSums(multZC)
  nCtot <- rowSums(multC)
  ZCtot / nCtot
}

# function to convert ssDNA base counts to dsDNA 20180318
make_dsDNA <- function(nuccomp) {
  nuccomp[, c("A", "C", "G", "T")] <- nuccomp[, c("A", "C", "G", "T")] + nuccomp[, c("T", "G", "C", "A")]
  nuccomp
}

# calculate Cys+Met fraction of amino acid compositions 20180324
CMAA <- function(AAcomp) {
  # find columns with names for the amino acids
  AA <- c("Ala", "Cys", "Asp", "Glu", "Phe", "Gly", "His", "Ile", "Lys",
    "Leu", "Met", "Asn", "Pro", "Gln", "Arg", "Ser", "Thr", "Val", "Trp", "Tyr")
  isAA <- colnames(AAcomp) %in% AA
  # columns for Cys and Met
  isCM <- colnames(AAcomp) %in% c("Cys", "Met")
  # calculate (Cys+Met) / (total AA)
  sum(AAcomp[, isCM]) / sum(AAcomp[, isAA])
}

