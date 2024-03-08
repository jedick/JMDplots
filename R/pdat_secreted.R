# JMDplots/pdat_secreted.R
# retrieve IDs for proteins secreted in hypoxia
# 20190325 extracted from pdat_hypoxia.R
# 20240426 moved from canprot package

pdat_secreted <- function(dataset = 2020) {
  if(identical(dataset, 2020)) {
    return(c("BRA+10", "PTD+10_Hx48=cancer", "PTD+10_Hx72=cancer",
             "JVC+12",
             "SKA+13", "SRS+13a_3", "SRS+13a_8",
             "LRS+14_Hy",
             "YKK+14_soluble=cancer", "YKK+14_exosome=cancer",
             "CRS+15_wt=cancer", "CRS+15_BT=cancer", "RTA+15=cancer",
             "RSE+16",
             "CGH+17_exosomes", "CGH+17_secretome",
             "CLY+18_secretome=cancer", "DWW+18=cancer", "FPR+18", "ODS+18",
             "CWG+19=cancer", "KAN+19_secretome=cancer", "NJVS19_CAM=cancer", "NJVS19_NTM", "PDT+19=cancer"))
  }
  # remove tags
  dataset <- strsplit(dataset, "=")[[1]][1]
  # get study and stage/condition
  study <- strsplit(dataset, "_")[[1]][1]
  stage <- paste(strsplit(dataset, "_")[[1]][-1], collapse="_")
  extdatadir <- system.file("extdata", package="JMDplots")
  datadir <- paste0(extdatadir, "/diffexpr/secreted/")
  if(study=="LRS+14") {
    # 20160717 rat heart myoblast secretome, Li et al., 2014
    # LRS+14_Hy, LRS+14_Re
    dat <- read.csv(paste0(datadir, "LRS+14.csv.xz"), as.is=TRUE)
    if(stage=="Hy") description <- "myoblast secretome"
    if(stage=="Re") description <- "myoblast secretome reoxygenation / normoxia"
    # select proteins with differential expression in Hy or Re
    icol <- grep(paste0(stage, ".Ctrl_iTRAQ"), colnames(dat))
    up2 <- dat[, icol] > 1
    dat <- cleanup(dat, "Entry", up2)
    pcomp <- protcomp(dat$Entry, aa_file=paste0(extdatadir, "/aa/rat/LRS+14_aa.csv.xz"))
  } else if(study=="RSE+16") {
    # 20160729 adipose-derived stem cells, Riis et al., 2016
    dat <- read.csv(paste0(datadir, "RSE+16.csv.xz"), as.is=TRUE)
    description <- "adipose-derived stem cells"
    pcomp <- protcomp(dat$Entry)
    up2 <- dat$Regulated == "up"
  } else if(study=="PTD+10") {
    # 20160801 A431 hypoxic / reoxygenated, Park et al., 2010
    # PTD+10_Hx48, PTD+10_Hx72, PTD+10_ReOx
    dat <- read.csv(paste0(datadir, "PTD+10.csv.xz"), as.is=TRUE)
    description <- paste("A431 squamous carcinoma cells", stage)
    if(stage=="Hx48") icol <- grep("115", colnames(dat))
    if(stage=="Hx72") icol <- grep("117", colnames(dat))
    if(stage=="ReOx") icol <- grep("116", colnames(dat))
    # filter ratio, p-value, EF value
    irat <- dat[, icol[1]] > 1.3 | dat[, icol[1]] < 1/1.3
    ipval <- dat[, icol[2]] < 0.05
    ief <- dat[, icol[3]] < 2
    dat <- dat[irat & ipval & ief, ]
    # drop missing entries
    up2 <- dat[, icol[1]] > 1.3
    dat <- cleanup(dat, "Entry", up2)
    pcomp <- protcomp(dat$Entry)
  } else if(study=="BRA+10") {
    # 20160805 placental tissue secretome, Blankley et al., 2010
    dat <- read.csv(paste0(datadir, "BRA+10.csv.xz"), as.is=TRUE)
    description <- "placental secretome"
    dat$UniProt.accession <- sapply(strsplit(dat$UniProt.accession, "|", fixed=TRUE), "[", 2)
    dat <- check_IDs(dat, "UniProt.accession")
    pcomp <- protcomp(dat$UniProt.accession)
    up2 <- dat$Fold.change > 0
  } else if(study=="DWW+18") {
    # 20190322 hypoxia-induced exosomes, Dorayappan et al., 2018
    dat <- read.csv(paste0(datadir, "DWW+18.csv.xz"), as.is=TRUE)
    description <- "ovarian cancer cell exosomes"
    dat <- check_IDs(dat, "Genes.symbol")
    pcomp <- protcomp(dat$Genes.symbol)
    up2 <- dat$FC > 1
  } else if(study=="CGH+17") {
    # 20190324 mouse cardiac fibroblast exosomes, secretome, Cosme et al., 2017
    # CGH+17_exosomes, CGH+17_secretome
    return(.pdat_multi(dataset))
  } else if(study=="CLY+18") {
    # 20190324 HCT116 cells, Chen et al., 2018
    # CLY+18_secretome
    return(.pdat_multi(dataset))
  } else if(study=="PDT+19") {
    # 20190326 tumor exosomes, Park et al., 2019
    dat <- read.csv(paste0(datadir, "PDT+19.csv.xz"), as.is=TRUE)
    description <- "mouse melanoma B16-F0 exosomes"
    # drop isoform suffixes
    dat$Accession <- sapply(strsplit(dat$Accession, "-"), "[", 1)
    dat <- check_IDs(dat, "Accession", aa_file=paste0(extdatadir, "/aa/mouse/PDT+19_aa.csv.xz"))
    up2 <- dat$Log2.127.126. > 0
    dat <- cleanup(dat, "Accession", up2)
    pcomp <- protcomp(dat$Accession, aa_file=paste0(extdatadir, "/aa/mouse/PDT+19_aa.csv.xz"))
  } else if(study=="SRS+13a") {
    # 20190327 placental mesenchymal stem cells 3% and 8% vs 1% O2, Salomon et al., 2013
    # SRS+13a_3, SRS+13a_8
    dat <- read.csv(paste0(datadir, "SRS+13a.csv.xz"), as.is=TRUE)
    description <- paste("pMSC", stage, "/ 1 % O2")
    # use selected dataset
    dat <- dat[dat$O2.percent %in% c(1, stage), ]
    pcomp <- protcomp(dat$Entry)
    up2 <- dat$O2.percent == stage
  } else if(study=="JVC+12") {
    # 20191204 endothelial cell-derived exosomes, de Jong et al., 2012
    dat <- read.csv(paste0(datadir, "JVC+12.csv.xz"), as.is=TRUE)
    description <- "endothelial cell-derived exosomes"
    # keep highly differential proteins
    dat <- dat[abs(dat$Hypoxia.median) > 0.2, ]
    up2 <- dat$Hypoxia.median > 0
    pcomp <- protcomp(dat$Entry)
  } else if(study=="SKA+13") {
    # 20191207 cytotrophoblast-derived exosomes, Salomon et al., 2013
    dat <- read.csv(paste0(datadir, "SKA+13.csv.xz"), as.is=TRUE)
    description <- "cytotrophoblast-derived exosomes"
    # compare 8% / 1% conditions
    dat <- dat[dat$O2.1_percent | dat$O2.8_percent, ]
    dat <- dat[xor(dat$O2.1_percent, dat$O2.8_percent), ]
    up2 <- dat$O2.8_percent
    pcomp <- protcomp(dat$Entry)
  } else if(study=="YKK+14") {
    # 20191207 U373MG glioma cells, Yoon et al., 2014
    # YKK+14_soluble, YKK+14_exosome
    dat <- read.csv(paste0(datadir, "YKK+14.csv.xz"), as.is=TRUE)
    description <- paste("U373MG glioma cells", stage)
    # get differential proteins for specified condition
    icol <- grep(stage, colnames(dat))
    dat <- dat[abs(dat[, icol]) > 0.5, ]
    up2 <- dat[, icol] > 0.5
    dat <- cleanup(dat, "Uniprot.Acc", up2)
    pcomp <- protcomp(dat$Uniprot.Acc)
  } else if(study=="NJVS19") {
    # 20191226 cancer-associated and normal tissue myofibroblasts, Najgebauer et al., 2019
    # NJVS19_CAM, NJVS19_NTM
    dat <- read.csv(paste0(datadir, "NJVS19.csv.xz"), as.is=TRUE)
    if(stage=="CAM") description <- "cancer-associated myofibroblasts"
    if(stage=="NTM") description <- "normal tissue myofibroblasts"
    # use selected dataset
    dat <- dat[!is.na(dat[, stage]), ]
    dat <- check_IDs(dat, "Majority.protein.IDs")
    pcomp <- protcomp(dat$Majority.protein.IDs)
    up2 <- dat[, stage] > 1
  } else if(study=="FPR+18") {
    # 20191226 endothelial progenitor cells, Felice et al., 2018
    dat <- read.csv(paste0(datadir, "FPR+18.csv.xz"), as.is=TRUE)
    description <- "endothelial progenitor cells"
    up2 <- dat$Modulation == "UP"
    pcomp <- protcomp(dat$Accession.number)
  } else if(study=="KAN+19") {
    # 20191226 human umbilical vein ECs, Kugeratski et al., 2019
    # KAN+19_secretome
    return(.pdat_multi(dataset))
  } else if(study=="CRS+15") {
    # 20200116 breast cancer MDA-MB-231 breast cancer parental and bone tropic cells, Cox et al., 2015
    # CRS+15_wt, CRS+15_BT
    dat <- read.csv(paste0(datadir, "CRS+15.csv.xz"), as.is=TRUE)
    if(stage=="wt") description <- "MDA-MB-231 breast cancer parental cells"
    if(stage=="BT") description <- "MDA-BT breast cancer bone tropic cells"
    icol <- grep(stage, colnames(dat))
    # calculate fold-change and keep highly differential proteins
    log2FC <- dat[, icol[2]] - dat[, icol[1]]
    dat <- cbind(dat, log2FC = log2FC)
    dat <- dat[abs(dat$log2FC) > 0.2, ]
    up2 <- dat$log2FC > 0.2
    pcomp <- protcomp(dat$Entry)
  } else if(study=="RTA+15") {
    # 20200117 LNCaP and PC3 cells, Ramteke et al., 2015
    dat <- read.csv(paste0(datadir, "RTA+15.csv.xz"), as.is=TRUE)
    description <- "LNCaP and PC3 prostate cancer cell exosomes"
    # remove proteins identified in both normoxic and hypoxic conditions
    dat <- dat[!(dat$Normoxic & dat$Hypoxic), ]
    up2 <- dat$Hypoxic & !dat$Normoxic
    pcomp <- protcomp(dat$Entry)
  } else if(study=="ODS+18") {
    # 20200117 AC10 ventricular cardiomyocyte extracellular vesicles, Ontoria-Oviedo et al., 2018
    dat <- read.csv(paste0(datadir, "ODS+18.csv.xz"), as.is=TRUE)
    description <- "AC10 ventricular cardiomyocyte extracellular vesicles"
    dat$Accession <- sapply(strsplit(dat$Accession, "\\|"), "[", 2)
    # remove proteins identified in both hypoxia and normoxia
    hyp <- dat$Accession[dat$Identified.in == "hypoxia"]
    nor <- dat$Accession[dat$Identified.in == "normoxia"]
    both <- intersect(hyp, nor)
    dat <- dat[!dat$Accession %in% both, ]
    dat <- check_IDs(dat, "Accession")
    up2 <- dat$Identified.in == "hypoxia"
    dat <- cleanup(dat, "Accession", up2)
    pcomp <- protcomp(dat$Accession)
  } else if(study=="CWG+19") {
    # 20200405 U87-MG glioma cell secretome
    dat <- read.csv(paste0(datadir, "CWG+19.csv.xz"), as.is=TRUE)
    description <- "U87-MG glioma cell extracellular vesicles"
    # drop proteins that are identified in both hypoxia and normoxia
    dat <- dat[!(dat$hypoxia & dat$normoxia), ]
    dat <- check_IDs(dat, "Accession", aa_file=paste0(extdatadir, "/aa/human/CWG+19_aa.csv.xz"))
    up2 <- dat$hypoxia
    pcomp <- protcomp(dat$Accession, aa_file = paste0(extdatadir, "/aa/human/CWG+19_aa.csv.xz"))
  } else stop(paste("secreted dataset", dataset, "not available"))
  print(paste0("pdat_secreted: ", description, " [", dataset, "]"))
  # use the up2 from the cleaned-up data, if it exists 20191120
  if("up2" %in% colnames(dat)) up2 <- dat$up2
  return(list(dataset=dataset, pcomp=pcomp, up2=up2, description=description))
}
