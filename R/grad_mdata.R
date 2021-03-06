# JMDplots/grad_mdata.R
# metadata for studies used in gradox and gradH2O papers 20181230
# added to JMDplots package 20190930

gradH2O <- list(
  `Baltic_Sea-0.1s` = list( # 20181228 [DLY+14], 20190104 [ZMI+17]
    # surface 0.1 um samples
    station = c("GS667", "GS665", "GS669", "GS673", "GS675", "GS659", "GS679", "GS681", "GS683", "GS685", "GS687", "GS694"),
    iMicrobe_MG = c("GS667_0.1", "GS665_0.1", NA, "GS673_0.1", NA, "GS659_0.1", "GS679_0.1", "GS681_0.1", "GS683_0.1", "GS685_0.1", "GS687_0.1", "GS694_0.1"),
    SRA_MT = c("SRR3747412", NA, "SRR3747762", NA, "SRR3747289", NA, "SRR3747362", NA, "SRR3747402", NA, NA, "SRR3747415"),
    MG_range = c(0.59, 0.63),
    MGP_range = c(-0.17, -0.13),
    MT_range = c(0.59, 0.63),
    MTP_range = c(-0.17, -0.13),
    abbrev = "GS",
    group = c("lake", "lake", "lake", "lake", "lake", "lake", "oxic", "oxic", "oxic", "oxic", "oxic", "oxic"),
    plot_real_x = FALSE,
    techtype = "454" # note: MT is Illumina
  ),
  `Baltic_Sea-0.8s` = list( # 20181228 [DLY+14], 20190104 [ZMI+17]
    # surface 0.8 um samples
    station = c("GS667", "GS665", "GS669", "GS673", "GS675", "GS659", "GS679", "GS681", "GS683", "GS685", "GS687", "GS694"),
    iMicrobe_MG = c("GS667_0.8", "GS665_0.8", NA, "GS673_0.8", NA, "GS659_0.8", "GS679_0.8", "GS681_0.8", "GS683_0.8", "GS685_0.8", "GS687_0.8", "GS694_0.8"),
    SRA_MT = c("SRR3747359", NA, "SRR3747763", NA, "SRR3747315", NA, "SRR3747363", NA, "SRR3747403", NA, NA, "SRR3747416"),
    MG_range = c(0.59, 0.63),
    MGP_range = c(-0.17, -0.13),
    MT_range = c(0.59, 0.63),
    MTP_range = c(-0.17, -0.13),
    abbrev = "GS",
    group = c("lake", "lake", "lake", "lake", "lake", "lake", "oxic", "oxic", "oxic", "oxic", "oxic", "oxic"),
    plot_real_x = FALSE,
    techtype = "454" # note: MT is Illumina
  ),
  `Baltic_Sea-3.0s` = list( # 20181228 [DLY+14], 20190104 [ZMI+17]
    # surface 3.0 um samples
    station = c("GS667", "GS665", "GS669", "GS673", "GS675", "GS659", "GS679", "GS681", "GS683", "GS685", "GS687", "GS694"),
    iMicrobe_MG = c("GS667_3.0", "GS665_3.0", NA, "GS673_3.0", NA, "GS659_3.0", "GS679_3.0", "GS681_3.0", "GS683_3.0", "GS685_3.0", "GS687_3.0", "GS694_3.0"),
    SRA_MT = c("SRR3747331", NA, "SRR3747764", NA, "SRR3746048", NA, "SRR3747364", NA, "SRR3747404", NA, NA, "SRR3747417"),
    MG_range = c(0.59, 0.63),
    MGP_range = c(-0.17, -0.13),
    MT_range = c(0.59, 0.63),
    MTP_range = c(-0.17, -0.13),
    abbrev = "GS",
    group = c("lake", "lake", "lake", "lake", "lake", "lake", "oxic", "oxic", "oxic", "oxic", "oxic", "oxic"),
    plot_real_x = FALSE,
    techtype = "454" # note: MT is Illumina
  ),
  `Baltic_Sea-0.1d` = list( # 20181228 [DLY+14], 20190104 [ZMI+17]
    # deep 0.1 um samples
    station = c("GS666", "GS670", "GS674", "GS660", "GS676", "GS677", "GS680", "GS682", "GS684", "GS686", "GS688", "GS695"),
    iMicrobe_MG = c("GS666_0.1", NA, "GS674_0.1", "GS660_0.1", NA, "GS677_0.1", NA, "GS682_0.1", "GS684_0.1", "GS686_0.1", "GS688_0.1", "GS695_0.1"),
    SRA_MT = c(NA, "SRR3747765", NA, NA, "SRR3747324", "SRR3747327", "SRR3747368", NA, "SRR3747413", NA, NA, "SRR3747418"),
    MG_range = c(0.59, 0.63),
    MGP_range = c(-0.17, -0.13),
    MT_range = c(0.59, 0.63),
    MTP_range = c(-0.17, -0.13),
    abbrev = "GS",
    group = c("lake", "lake", "lake", "lake", "lake", "lake", "oxic", "oxic", "oxic", "oxic", "oxic", "oxic"),
    plot_real_x = FALSE,
    techtype = "454" # note: MT is Illumina
  ),
  `Baltic_Sea-0.8d` = list( # 20181228 [DLY+14], 20190104 [ZMI+17]
    # deep 0.8 um samples
    station = c("GS666", "GS670", "GS674", "GS660", "GS676", "GS677", "GS680", "GS682", "GS684", "GS686", "GS688", "GS695"),
    iMicrobe_MG = c("GS666_0.8", NA, "GS674_0.8", "GS660_0.8", NA, "GS677_0.8", NA, "GS682_0.8", "GS684_0.8", "GS686_0.8", "GS688_0.8", "GS695_0.8"),
    SRA_MT = c(NA, "SRR3747775", NA, NA, "SRR3747325", "SRR3747328", "SRR3747400", NA, "SRR3747414", NA, NA, "SRR3747419"),
    MG_range = c(0.59, 0.63),
    MGP_range = c(-0.17, -0.13),
    MT_range = c(0.59, 0.63),
    MTP_range = c(-0.17, -0.13),
    abbrev = "GS",
    group = c("lake", "lake", "lake", "lake", "lake", "lake", "oxic", "oxic", "oxic", "oxic", "oxic", "oxic"),
    plot_real_x = FALSE,
    techtype = "454" # note: MT is Illumina
  ),
  `Baltic_Sea-3.0d` = list( # 20181228 [DLY+14], 20190104 [ZMI+17]
    # deep 3.0 um samples
    station = c("GS666", "GS670", "GS674", "GS660", "GS676", "GS677", "GS680", "GS682", "GS684", "GS686", "GS688", "GS695"),
    iMicrobe_MG = c("GS666_3.0", NA, "GS674_3.0", "GS660_3.0", NA, "GS677_3.0", NA, "GS682_3.0", "GS684_3.0", "GS686_3.0", "GS688_3.0", "GS695_3.0"),
    SRA_MT = c(NA, "SRR3747777", NA, NA, "SRR3747326", "SRR3747329", "SRR3747401", NA, "SRR3747405", NA, NA, "SRR3747420"),
    MG_range = c(0.59, 0.63),
    MGP_range = c(-0.17, -0.13),
    MT_range = c(0.59, 0.63),
    MTP_range = c(-0.17, -0.13),
    abbrev = "GS",
    group = c("lake", "lake", "lake", "lake", "lake", "lake", "oxic", "oxic", "oxic", "oxic", "oxic", "oxic"),
    plot_real_x = FALSE,
    techtype = "454" # note: MT is Illumina
  ),
  Eiler_Freshwater = list( # 20190719 [EZM+14]
    sample = c("DamariscottaSP", "DamariscottaSU", "Ekoln", "Erken", "Lanier", "MendotaSP", "MendotaSU",
               "Spark", "Trout", "Vattern", "Yellowstone1", "Yellowstone2"),
    SRA_MG = c("ERR358545", "ERR358546", "ERR358543", "ERR358542", "SRR063691", "ERR358550", "ERR358549",
               "ERR358548", "ERR358547", "ERR358544", "SRR077348", "SRR078855"),
    MG_range = c(0.57, 0.64),
    MGP_range = c(-0.20, -0.13),
    abbrev = "EF",
    group = c("lake"),
    plot_real_x = FALSE,
    techtype = "454"
  ),
  Eiler_Marine = list( # 20190719 [EZM+14]
    sample = c("BATS0", "BATS200", "BATS250", "BATS40", "EqDP35155", "NPTG35179", "PNEq35163", "PNEqCc35171",
               "SPSG35131", "SPSG35139", "SPSG35147", "WChannelApr", "WChannelJan"),
    iMicrobe_MG = c("BATS_SMPL_BATS-167-0", "BATS_SMPL_BATS-167-200", "BATS_SMPL_BATS-167-250", "BATS_SMPL_BATS-167-40",
                    "BACTERIOPLANKTON_SMPL_S_35155", "BACTERIOPLANKTON_SMPL_S_35179", "BACTERIOPLANKTON_SMPL_S_35163", "BACTERIOPLANKTON_SMPL_S_35171",
                    "BACTERIOPLANKTON_SMPL_S_35131", "BACTERIOPLANKTON_SMPL_S_35139", "BACTERIOPLANKTON_SMPL_S_35147",
                    "WESTERNCHANNELOMM_SMPL_APR_DAY", "WESTERNCHANNELOMM_SMPL_JAN_DAY"),
    MG_range = c(0.57, 0.64),
    MGP_range = c(-0.20, -0.13),
    abbrev = "EM",
    group = c("oxic"),
    plot_real_x = FALSE,
    techtype = "454"
  ),
  Santa_Pola = list( # 20181224 [GPF+11, FGM+13]
    sample = c("SS13", "SS19", "SS33", "SS37"),
    SRA_MG = c("SRR944625", "SRR328982", "SRR979792", "SRR328983"),
    MG_range = c(0.61, 0.64),
    MGP_range = c(-0.13, -0.08),
    abbrev = "SA",
    group = c("hypersaline_low", "hypersaline_low", "hypersaline", "hypersaline"),
    plot_real_x = FALSE,
    techtype = "454"
  ),
  SouthBay_Water = list( # 20181225 [KBW+18]
    pond = c("2C", "1C", "A23"),
    SRA_MG = c("SRR4030050", "SRR4030051", "SRR4030044"),
    MG_range = c(0.61, 0.65),
    MGP_range = c(-0.15, -0.09),
    abbrev = "SB",
    group = c("hypersaline_low", "hypersaline_low", "hypersaline"),
    plot_real_x = FALSE,
    techtype = "Illumina"
  ),
  Kulunda_Steppe = list( # 20190104 [VGR+16]
    Lake = c("T5", "PL", "Tc", "B1"),
    SRA_MG = c("SRR490135", "SRR490126", "SRR490146", "SRR490137"),
    MG_range = c(0.61, 0.65),
    MGP_range = c(-0.15, -0.09),
    abbrev = "KS",
    group = c("hypersaline_low", "hypersaline_low", "hypersaline", "hypersaline"),
    plot_real_x = FALSE,
    techtype = "Illumina"
  ),
  `Amazon_River-FL` = list( # 20190101 [SFD+15, SZD+14]
    site = c("OB", "TAPS", "TAPD", "MCPN", "MCPS", "BLM", "S10", "S3", "S2", "S23", "S25", "S27"),
    SRA_MG = c("SRR1790676", "SRR1796116", "SRR1792674", "SRR1786279", "SRR1788318", "SRR1790489",
               "SRR1199271", "SRR1205253", "SRR1199270", "SRR1186214", "SRR1204581", "SRR1202091"),
    SRA_MT = c("SRR1781945", "SRR1777513", "SRR1781714", "SRR1785351", "SRR1784305", "SRR1781804",
               "SRR1186930", "SRR1199283", "SRR1193205", "SRR1193632", "SRR1205247", "SRR1193627"),
    MG_range = c(0.56, 0.62),
    MGP_range = c(-0.20, -0.14),
    MT_range = c(0.60, 0.63),
    MTP_range = c(-0.16, -0.12),
    abbrev = "ARF",
    group = c("riverFL", "riverFL", "riverFL", "riverFL", "riverFL", "riverFL",
              "plumeFL", "plumeFL", "plumeFL", "plumeFL", "plumeFL", "plumeFL"),
    plot_real_x = FALSE,
    techtype = "Illumina"
  ),
  `Amazon_River-PA` = list( # 20190101 [SFD+15, SZD+14]
    site = c("OB", "TAPS", "TAPD", "MCPN", "MCPS", "BLM", "S10", "S3", "S2", "S23", "S25", "S27"),
    SRA_MG = c("SRR1790678", "SRR1796118", "SRR1792852", "SRR1786281", "SRR1790487", "SRR1790644",
               "SRR1209978", "SRR1199272", "SRR1202081", "SRR1202089", "SRR1209976", "SRR1202095"),
    SRA_MT = c("SRR1782579", "SRR1778024", "SRR1781802", "SRR1785352", "SRR1785207", "SRR1781811",
               "SRR1193190", "SRR1199284", "SRR1193177", "SRR1193237", "SRR1205800", "SRR1193629"),
    MG_range = c(0.58, 0.63),
    MGP_range = c(-0.18, -0.13),
    MT_range = c(0.60, 0.63),
    MTP_range = c(-0.16, -0.12),
    abbrev = "ARP",
    group = c("riverPA", "riverPA", "riverPA", "riverPA", "riverPA", "riverPA",
              "plumePA", "plumePA", "plumePA", "plumePA", "plumePA", "plumePA"),
    plot_real_x = FALSE,
    techtype = "Illumina"
  ),
  `Rodriguez_Brito-mic` = list( # 20190723 [RLW+10]
    sample = c("A-M", "B-M", "C-M", "B2-M",
               "D-M", "E-M", "F-M",
               "G-M", "H-M", "I-M", "K-M",
               "P-M", "M-M"),
    MGRAST_MG = c("mgm4440440.3", "mgm4440413.3", "mgm4440422.3", "mgm4440411.3",
                  "mgm4440437.3", "mgm4440324.3", "mgm4440426.3",
                  "mgm4440435.3", "mgm4440434.3", "mgm4440425.3", "mgm4440416.3",
                  "mgm4440438.3", "mgm4440419.3"
    ),
    MG_range = c(0.57, 0.64),
    MGP_range = c(-0.20, -0.13),
    abbrev = "RM",
    group = c("lake", "lake", "lake", "lake",
              "oxic", "oxic", "oxic",
              "hypersaline", "hypersaline", "hypersaline", "hypersaline",
              "hypersaline", "hypersaline"),
    plot_real_x = FALSE,
    techtype = "454"
  ),
  `Rodriguez_Brito-vir` = list( # 20190723 [RLW+10]
    sample = c("A-V", "B-V", "C-V", "B2-V",
               "D-V", "E-V", "F-V",
               "G-V", "H-V", "I-V", "J-V", "K-V",
               "L-V", "M-V", "N-V"),
    MGRAST_MG = c("mgm4440439.3", "mgm4440412.3", "mgm4440424.3", "mgm4440414.3",
                  "mgm4440436.3", "mgm4440432.3", "mgm4440420.3",
                  "mgm4440431.3", "mgm4440325.3", "mgm4440428.3", "mgm4440417.3", "mgm4440427.3",
                  "mgm4440421.3", "mgm4440144.4", "mgm4440145.4"
    ),
    MG_range = c(0.57, 0.64),
    MGP_range = c(-0.20, -0.13),
    abbrev = "RV",
    group = c("lake", "lake", "lake", "lake",
              "oxic", "oxic", "oxic",
              "hypersaline", "hypersaline", "hypersaline", "hypersaline", "hypersaline",
              "hypersaline", "hypersaline", "hypersaline"),
    plot_real_x = FALSE,
    techtype = "454"
  )
)

# groups: hydrothermal, rock, mat, marine, hypersaline 20180315
# added OMZ, plume 20180321
gradox <- list(
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
    plot_real_x = TRUE,
    techtype = "Illumina"
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
    plot_real_x = TRUE,
    techtype = "Sanger"
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
    plot_real_x = FALSE,
    techtype = "Illumina"
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
    plot_real_x = TRUE,
    techtype = "Illumina"
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
    plot_real_x = TRUE,
    techtype = "454"
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
    plot_real_x = FALSE,
    techtype = "Sanger"
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
    plot_real_x = TRUE,
    techtype = "454"
  ),
  `HOT_ALOHA-2010` = list( # 20181117 [MBA+17]
    `depth (m)` = c("25m", "45m", "75m", "125m", "200m", "500m", "770m", "1000m"),
    xlabels = c(25, 45, 75, 125, 200, 500, 770, 1000),
    SRA_MG = c("SRX2334639", "SRX2334707", "SRX2334640", "SRX2334668", "SRX2334697", "SRX2334652", "SRX2334669", "SRX2334641"),
    MG_srange = c(0.61, 0.57),
    abbrev = "HA",
    group = c("oxic0", "oxic", "oxic", "oxic", "oxic", "oxic", "oxic", "oxic"),
    plot_real_x = TRUE,
    techtype = "Illumina"
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
    plot_real_x = TRUE,
    techtype = "Illumina"
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
    plot_real_x = TRUE,
    techtype = "Illumina"
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
    plot_real_x = TRUE,
    techtype = "454"
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
    plot_real_x = TRUE,
    techtype = "454"
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
    plot_real_x = FALSE,
    techtype = "454"  # note: Lost City is Sanger sequencing [BNS12]
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
    plot_real_x = TRUE,
    techtype = "Sanger"
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
    plot_real_x = FALSE,
    techtype = "Sanger"
  )
)

studies <- c(gradox, gradH2O)

