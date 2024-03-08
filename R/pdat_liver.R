# JMDplots/pdat_liver.R
# retrieve protein IDs for liver cancer
# 20200330 first version jmd
# 20240426 moved from canprot package

pdat_liver <- function(dataset = 2020) {
  if(identical(dataset, 2020)) { 
    return(c(
             "LHT+04",
             "BLP+05", "LTZ+05",
             "DTS+07", "SXS+07",
             "CHN+08",
             "RLA+10=rat",
             "LMG+11_nuclear", "LMG+11_cytoskeletal",
             "LRL+12",
             "KOK+13", "MBK+13",
             "XWS+14",
             "BSG15=mouse", "RPM+15",
             "NBM+16_G1", "NBM+16_G2", "NBM+16_G3", "NMB+16", "QXC+16_T1", "QXC+16_T2", "QXC+16_T3",
             "GJZ+17", "GWS+17", "QPP+17", "WLL+17_small", "WLL+17_medium", "WLL+17_large", "WLL+17_huge",
             "BOK+18", "YXZ+18",
             "BEM+20=mouse", "GZD+19_protein", "GZD+19_phosphoprotein", "JSZ+19", "ZZL+19",
             "GZL+20", "SCL+20_differential", "SCL+20_unique"
             ))
  }
  # remove tags
  dataset <- strsplit(dataset, "=")[[1]][1]
  # get study and stage/condition
  study <- strsplit(dataset, "_")[[1]][1]
  stage <- paste(strsplit(dataset, "_")[[1]][-1], collapse="_")
  extdatadir <- system.file("extdata", package="JMDplots")
  datadir <- paste0(extdatadir, "/diffexpr/liver/")
  if(study=="MBK+13") {
    # 20160417 Megger et al., 2013
    dat <- read.csv(paste0(datadir, "MBK+13.csv.xz"), as.is=TRUE)
    description <- "T / N"
    up2 <- dat$Highest.mean.condition == "HCC"
    dat <- cleanup(dat, "Entry", up2)
    pcomp <- protcomp(dat$Entry)
  } else if(study=="CHN+08") {
    # 20160419 Chaerkady et al., 2008
    dat <- read.csv(paste0(datadir, "CHN+08.csv.xz"), as.is=TRUE)
    description <- "T / N"
    up2 <- dat$Regulated == "up"
    dat <- cleanup(dat, "Entry", up2)
    pcomp <- protcomp(dat$Entry)
  } else if(study=="NMB+16") {
    # 20160717 Naboulsi et al., 2016
    dat <- read.csv(paste0(datadir, "NMB+16.csv.xz"), as.is=TRUE)
    description <- "T / N"
    dat <- check_IDs(dat, "Accession")
    up2 <- dat$Highest.mean.condition == "HCC"
    pcomp <- protcomp(dat$Accession)
  } else if(study=="KOK+13") {
    # 20170114 Kimura et al., 2013
    dat <- read.csv(paste0(datadir, "KOK+13.csv.xz"), as.is=TRUE)
    description <- "T / N"
    dat <- check_IDs(dat, "Accession.no.")
    up2 <- dat$Fold.difference..Ratio.of.means. > 1
    dat <- cleanup(dat, "Accession.no.", up2)
    pcomp <- protcomp(dat$Accession.no.)
  } else if(study=="XWS+14") {
    # 20170115 Xu et al., 2014
    dat <- read.csv(paste0(datadir, "XWS+14.csv.xz"), as.is=TRUE)
    description <- "T / N"
    dat <- check_IDs(dat, "Uniprot")
    up2 <- dat$Fold.change > 0
    dat <- cleanup(dat, "Uniprot", up2)
    pcomp <- protcomp(dat$Uniprot)
  } else if(study=="LRL+12") {
    # 20170617 Li et al., 2012
    dat <- read.csv(paste0(datadir, "LRL+12.csv.xz"), as.is=TRUE)
    description <- "T / N"
    up2 <- dat$HCC.non.HCC.ratio > 1
    dat <- cleanup(dat, "Entry", up2)
    pcomp <- protcomp(dat$Entry)
  } else if(study=="GZD+19") {
    # 20200330 Gao et al., 2019
    # GZD+19_protein, GZD+19_phosphoprotein
    dat <- read.csv(paste0(datadir, "GZD+19.csv.xz"), as.is=TRUE)
    description <- paste("T / N", stage)
    icol <- grep(paste0("^", stage), colnames(dat))
    dat <- dat[!is.na(dat[, icol]), ]
    up2 <- dat[, icol] > 1
    pcomp <- protcomp(dat$Entry)
  } else if(study=="LHT+04") {
    # 20200330 Li et al., 2004
    dat <- read.csv(paste0(datadir, "LHT+04.csv.xz"), as.is=TRUE)
    description <- "T / N"
    up2 <- dat$Protein.ratio > 1
    dat <- cleanup(dat, "Entry", up2)
    pcomp <- protcomp(dat$Entry)
  } else if(study=="SXS+07") {
    # 20200330 Sun et al., 2007
    dat <- read.csv(paste0(datadir, "SXS+07.csv.xz"), as.is=TRUE)
    description <- "T / N"
    up2 <- dat$Regulation == "Up"
    dat <- check_IDs(dat, "Uniprot")
    dat <- cleanup(dat, "Uniprot", up2)
    pcomp <- protcomp(dat$Uniprot)
  } else if(study=="RPM+15") {
    # 20200330 Reis et al., 2015
    dat <- read.csv(paste0(datadir, "RPM+15.csv.xz"), as.is=TRUE)
    description <- "T / N"
    dat <- check_IDs(dat, "Accession")
    up2 <- dat$Highest.mean.condition == "HCC"
    pcomp <- protcomp(dat$Accession)
  } else if(study=="GWS+17") {
    # 20200330 Gao et al., 2017
    dat <- read.csv(paste0(datadir, "GWS+17.csv.xz"), as.is=TRUE)
    description <- "T / N"
    dat <- check_IDs(dat, "Accession")
    up2 <- dat$average.FC > 1
    dat <- cleanup(dat, "Accession", up2)
    pcomp <- protcomp(dat$Accession)
  } else if(study=="WLL+17") {
    # 20200330 Wang et al., 2017
    # WLL+17_small, WLL+17_medium, WLL+17_large, WLL+17_huge
    dat <- read.csv(paste0(datadir, "WLL+17.csv.xz"), as.is=TRUE)
    description <- paste("T / N", stage)
    # keep proteins that are significantly different in this tumor size
    icol <- grep(stage, colnames(dat))
    ifold <- dat[, icol[1]] > 2 | dat[, icol[1]] < 0.5
    ip <- dat[, icol[2]] < 0.05
    dat <- dat[ifold & ip, ]
    up2 <- dat[, icol[1]] > 2
    dat <- cleanup(dat, "Entry", up2)
    pcomp <- protcomp(dat$Entry)
  } else if(study=="BOK+18") {
    # 20200331 Buczak et al., 2018
    dat <- read.csv(paste0(datadir, "BOK+18.csv.xz"), as.is=TRUE)
    description <- "T / N"
    dat <- check_IDs(dat, "ID")
    up2 <- dat$log2FC > 1
    dat <- cleanup(dat, "ID", up2)
    pcomp <- protcomp(dat$ID)
  } else if(study=="NBM+16") {
    # 20200331 Naboulsi et al., 2016
    # NBM+16_G1, NBM+16_G2, NBM+16_G3
    dat <- read.csv(paste0(datadir, "NBM+16.csv.xz"), as.is=TRUE)
    description <- paste("T / N", stage)
    dat <- check_IDs(dat, "Accession")
    # keep proteins that are differentially expressed in this grade
    icol <- grep(stage, colnames(dat))
    dat <- dat[!is.na(dat[, icol]), ]
    up2 <- dat[, icol[1]] > 1
    dat <- cleanup(dat, "Accession", up2)
    pcomp <- protcomp(dat$Accession)
  } else if(study=="JSZ+19") {
    # 20200331 Jiang et al., 2019
    dat <- read.csv(paste0(datadir, "JSZ+19.csv.xz"), as.is=TRUE)
    description <- "T / N"
    dat <- check_IDs(dat, "Uniprot_ID")
    up2 <- dat$logFC > 0
    pcomp <- protcomp(dat$Uniprot_ID)
  } else if(study=="ZZL+19") {
    # 20200331 Zhu et al., 2019
    dat <- read.csv(paste0(datadir, "ZZL+19.csv.xz"), as.is=TRUE)
    description <- "T / N"
    up2 <- dat$fold.change > 2
    pcomp <- protcomp(dat$UniProt)
  } else if(study=="GJZ+17") {
    # 20200401 Guo et al., 2017
    dat <- read.csv(paste0(datadir, "GJZ+17.csv.xz"), as.is=TRUE)
    description <- "T / N"
    up2 <- dat$iTRAQ.ratio..HCC.HC. > 1
    pcomp <- protcomp(dat$Entry)
  } else if(study=="BLP+05") {
    # 20200401 Blanc et al., 2005
    dat <- read.csv(paste0(datadir, "BLP+05.csv.xz"), as.is=TRUE)
    description <- "T / N"
    up2 <- dat$Regulation == "Up"
    dat <- cleanup(dat, "Entry", up2)
    pcomp <- protcomp(dat$Entry)
  } else if(study=="LTZ+05") {
    # 20200401 Li et al., 2005
    dat <- read.csv(paste0(datadir, "LTZ+05.csv.xz"), as.is=TRUE)
    description <- "T / N"
    dat <- check_IDs(dat, "Entry")
    up2 <- dat$Regulation == "Up"
    dat <- cleanup(dat, "Entry", up2)
    pcomp <- protcomp(dat$Entry)
  } else if(study=="LMG+11") {
    # 20200401 Lee et al., 2011
    # LMG+11_nuclear, LMG+11_cytoskeletal
    dat <- read.csv(paste0(datadir, "LMG+11.csv.xz"), as.is=TRUE)
    description <- paste("T / N", stage)
    # keep proteins that are differentially expressed in this location
    icol <- grep(stage, colnames(dat))
    dat <- dat[!is.na(dat[, icol]), ]
    dat <- check_IDs(dat, "UniProt")
    up2 <- dat[, icol] > 2
    dat <- cleanup(dat, "UniProt", up2)
    pcomp <- protcomp(dat$UniProt)
  } else if(study=="QPP+17") {
    # 20200401 Qiao et al., 2017
    dat <- read.csv(paste0(datadir, "QPP+17.csv.xz"), as.is=TRUE)
    description <- "T / N"
    dat <- check_IDs(dat, "UniProt")
    up2 <- dat$FC > 1.5
    dat <- cleanup(dat, "UniProt", up2)
    pcomp <- protcomp(dat$UniProt)
  } else if(study=="QXC+16") {
    # 20200403 Qi et al., 2016
    # QXC+16_T1, QXC+16_T2, QXC+16_T3
    dat <- read.csv(paste0(datadir, "QXC+16.csv.xz"), as.is=TRUE)
    description <- paste("T / N", stage)
    dat <- check_IDs(dat, "Protein.IDs")
    icol <- grep(stage, colnames(dat))
    dat <- dat[!is.na(dat[, icol]), ]
    up2 <- dat[, icol] > 1
    pcomp <- protcomp(dat$Protein.IDs)
  } else if(study=="BEM+20") {
    # 20200405 mouse liver cancer, Berndt et al., 2020
    dat <- read.csv(paste0(datadir, "BEM+20.csv.xz"), as.is=TRUE)
    description <- paste("T / N mouse")
    dat <- check_IDs(dat, "majority.protein.ids", aa_file = paste0(extdatadir, "/aa/mouse/BEM+20_aa.csv.xz"))
    up2 <- dat$fold.change > 2
    pcomp <- protcomp(dat$majority.protein.ids, aa_file = paste0(extdatadir, "/aa/mouse/BEM+20_aa.csv.xz"))
  } else if(study=="BSG15") {
    # 20200416 EGF transgenic mice, Borlak et al., 2015
    dat <- read.csv(paste0(datadir, "BSG15.csv.xz"), as.is=TRUE)
    description <- "T / N EGF transgenic mice"
    up2 <- dat$Ratio.T.C > 1
    dat <- cleanup(dat, "Entry", up2)
    pcomp <- protcomp(dat$Entry, aa_file = paste0(extdatadir, "/aa/mouse/BSG15_aa.csv.xz"))
  } else if(study=="DTS+07") {
    # 20200416 Dos Santos et al., 2007
    dat <- read.csv(paste0(datadir, "DTS+07.csv.xz"), as.is=TRUE)
    description <- "T / N"
    dat <- check_IDs(dat, "SwissProt.accession.no.")
    # include values for both tumor homogenates and LM samples
    up2 <- sapply(apply(sign(dat[, 2:3]), 1, unique), na.omit) == 1
    dat <- cleanup(dat, "SwissProt.accession.no.", up2)
    pcomp <- protcomp(dat$SwissProt.accession.no.)
  } else if(study=="RLA+10") {
    # 20200417 transitional endoplasmic reticulum, Roy et al., 2010
    dat <- read.csv(paste0(datadir, "RLA+10.csv.xz"), as.is=TRUE)
    description <- "T / N rat transitional endoplasmic reticulum"
    up2 <- dat$Ratio > 1
    pcomp <- protcomp(dat$Entry, aa_file = paste0(extdatadir, "/aa/rat/RLA+10_aa.csv.xz"))
  } else if(study=="SCL+20") {
    # 20200417 Shin et al., 2020
    # SCL+20_differential, SCL+20_unique
    dat <- read.csv(paste0(datadir, "SCL+20.csv.xz"), as.is=TRUE)
    description <- paste("T / N mitochondrial", stage)
    if(stage == "differential") dat <- dat[!is.infinite(dat$Ratio) & !dat$Ratio==0, ]
    if(stage == "unique") dat <- dat[is.infinite(dat$Ratio) | dat$Ratio==0, ]
    dat <- check_IDs(dat, "ID", aa_file = paste0(extdatadir, "/aa/human/SCL+20_aa.csv.xz"))
    up2 <- dat$Ratio > 1
    dat <- cleanup(dat, "ID", up2)
    pcomp <- protcomp(dat$ID, aa_file = paste0(extdatadir, "/aa/human/SCL+20_aa.csv.xz"))
  } else if(study=="GZL+20") {
    # 20200418 Gao et al., 2020
    dat <- read.csv(paste0(datadir, "GZL+20.csv.xz"), as.is=TRUE)
    description <- "T / N"
    up2 <- dat$Expr.Fold.Change > 0
    pcomp <- protcomp(dat$Accession)
  } else if(study=="YXZ+18") {
    # 20200420 liver cancer, Yang et al., 2018
    dat <- read.csv(paste0(datadir, "YXZ+18.csv.xz"), as.is=TRUE)
    description <- "T / N"
    up2 <- dat$Ratio..tumor.adjacent. > 1
    dat <- cleanup(dat, "Entry", up2)
    pcomp <- protcomp(dat$Entry)
  } else stop(paste("liver dataset", dataset, "not available"))
  print(paste0("pdat_liver: ", description, " [", dataset, "]"))
  # use the up2 from the cleaned-up data, if it exists 20190407
  if("up2" %in% colnames(dat)) up2 <- dat$up2
  return(list(dataset = dataset, pcomp = pcomp, up2 = up2, description = description))
}
