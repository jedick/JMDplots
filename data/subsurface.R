# JMDplots/subsurface.R
# metadata for a couple of subsurface metagenomes
# used to exemplify user-provided data sets 20191004

subsurface <- list(
  Peru_Margin = list( # 20180313 [OECB13]
    depth = c("5m", "30m", "50m", "70m", "91m", "159m"),
    # FIXME: read.fasta() has spurious unrecognized letters in SRX187180, SRX330851, SRX331958
    SRA_MT = c("SRX187180", "SRX330851", "SRX331883", "SRX331884", "SRX331956", "SRX331958"),
    MT_range = c(1.64, 1.72),
    MTP_range = c(-0.18, -0.13),
    abbrev = "PM",
    group = "rock",
    plot_real_x = FALSE
  ),
  North_Pond = list( # 20180320 [TWGH18]
    # time points of samples from U1832A 
    sample = c("TP1", "TP2", "TP3", "TP4", "TP5", "TP6", "TP7", "TP8", "TP9", "CTD"),
    SRA_MG = c("SRX3143893", "SRX3143894", "SRX3143895", "SRX3143896", "SRX3143897",
               "SRX3143898", "SRX3143899", "SRX3143900", "SRX3143886", "SRX3143890"),
    MG_range = c(1.63, 1.71),
    MGP_range = c(-0.18, -0.14),
    abbrev = "NP",
    group = c("oxygenated", "oxygenated", "oxygenated", "oxygenated", "oxygenated",
              "oxygenated", "oxygenated", "oxygenated", "oxygenated", "oxygenated0"),
    plot_real_x = FALSE
  )
)
