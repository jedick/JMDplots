# JMDplots/pdat_osmotic_halo.R
# retrieve protein IDs for hyper- and hypoosmotic stress in halophiles
# new datasets added 20191104
# renamed to pdat_osmotic_halo 20200418
# 20240426 moved from canprot package

pdat_osmotic_halo <- function(dataset = 2020) {
  if(identical(dataset, 2020)) {
    return(c("LRB+09_2.6=hypoosmotic", "LRB+09_5.1",
             "ZLZ+16_10=hypoosmotic", "ZLZ+16_17.5",
             "LLYL17_0=hypoosmotic", "LLYL17_3.5",
             "JSP+19_LoS=hypoosmotic", "JSP+19_HiS" #, "JSP+19_LoT", "JSP+19_HiT"
             ))
  }
  # remove tags
  dataset <- strsplit(dataset, "=")[[1]][1]
  # get study and stage/condition
  study <- strsplit(dataset, "_")[[1]][1]
  stage <- paste(strsplit(dataset, "_")[[1]][-1], collapse="_")
  extdatadir <- system.file("extdata", package="JMDplots")
  datadir <- paste0(extdatadir, "/diffexpr/halophiles/")
  if(study=="ZLZ+16") {
    # 20191103 Nocardiopsis xinjiangensis, Zhang et al., 2016
    # ZLZ+16_10, ZLZ+16_17.5
    dat <- read.csv(file.path(datadir, "ZLZ+16.csv.xz"), as.is = TRUE)
    if(stage=="10") {
      description <- "_Nocardiopsis xinjiangensis_ 6% / 10% NaCl"
      pval <- dat$p.Value..6..vs..10..
      icolratio <- grep("ratio_6_over_10", colnames(dat))
    }
    if(stage=="17.5") {
      description <- "_Nocardiopsis xinjiangensis_ 17.5% / 10% NaCl"
      pval <- dat$p.Value..10..vs..17.5..
      icolratio <- grep("ratio_17.5_over_10", colnames(dat))
    }
    # use selected condition
    idiff <- (dat[, icolratio] > 1.3 | dat[, icolratio] < 1/1.3) & pval < 0.05
    dat <- dat[idiff, ]
    up2 <- dat[, icolratio] > 1.3
    pcomp <- protcomp(dat$Entry, aa_file = file.path(extdatadir, "aa/bacteria/ZLZ+16_aa.csv.xz"))
  } else if(study=="LRB+09") {
    # 20191101 Halobacterium salinarum NaCl adjustment, Leuko et al., 2009
    # LRB+09_2.6, LRB+09_5.1
    dat <- read.csv(file.path(datadir, "LRB+09.csv.xz"), as.is=TRUE)
    if(stage=="2.6") description <- "_Halobacterium salinarium_ 2.6 M / 4.3 M NaCl"
    if(stage=="5.1") description <- "_Halobacterium salinarium_ 5.1 M / 4.3 M NaCl"
    # use selected dataset
    icol <- grep(stage, colnames(dat))
    dat <- dat[!is.na(dat[, icol]), ]
    # up-expressed proteins in low or high
    up2 <- dat[, icol] > 0
    # drop missing proteins
    dat <- cleanup(dat, "Entry", up2)
    pcomp <- protcomp(dat$Entry, aa_file=file.path(extdatadir, "aa/archaea/LRB+09_aa.csv.xz"))
  } else if(study=="LLYL17") {
    # 20191102 Tetragenococcus halophilus NaCl adjustment, Lin et al., 2017
    # LLYL17_0, LLYL17_3.5
    dat <- read.csv(file.path(datadir, "LLYL17.csv.xz"), as.is=TRUE)
    if(stage=="0") description <- "_Tetragenococcus halophilus_ 0 M / 1 M NaCl"
    if(stage=="3.5") description <- "_Tetragenococcus halophilus_ 3.5 M / 1 M NaCl"
    # use selected dataset
    icol <- grep(stage, colnames(dat))
    dat <- dat[!is.na(dat[, icol]), ]
    # up-expressed proteins in low or high salinity (0 M / 1 M or 3.5 M / 1 M)
    if(stage=="0") up2 <- dat[, icol] > 0
    if(stage=="3.5") up2 <- dat[, icol] < 0
    pcomp <- protcomp(dat$UniProtKB.Entry, aa_file=file.path(extdatadir, "aa/bacteria/LLYL17_aa.csv.xz"))
  } else if(study=="JSP+19") {
    # 20191102 Haloferax volcanii salt and temperature, Jevtić et al., 2019
    # JSP+19_LoS, JSP+19_HiS, JSP+19_LoT, JSP+19_HiT
    dat <- read.csv(file.path(datadir, "JSP+19.csv.xz"), as.is=TRUE)
    if(stage=="LoS") description <- "_Haloferax volcanii_ 10.8% / 15% NaCl"
    if(stage=="HiS") description <- "_Haloferax volcanii_ 19.2% / 15% NaCl"
    # use selected condition
    icol <- grep(stage, colnames(dat))
    # at least two-fold, significant difference
    idiff <- abs(dat[, icol[1]]) >= 1 & dat[, icol[2]]==1
    dat <- dat[idiff, ]
    up2 <- dat[, icol[1]] > 0
    # remove NA accessions
    dat <- cleanup(dat, "UniProt.Accession", up2)
    pcomp <- protcomp(dat$UniProt.Accession, aa_file=file.path(extdatadir, "aa/archaea/JSP+19_aa.csv.xz"))
  } else stop(paste("osmotic_halo dataset", dataset, "not available"))
  print(paste0("pdat_osmotic_halo: ", description, " [", dataset, "]"))
  # use the up2 from the cleaned-up data, if it exists 20190407
  if("up2" %in% colnames(dat)) up2 <- dat$up2
  return(list(dataset=dataset, pcomp=pcomp, up2=up2, names=names, description=description))
}
