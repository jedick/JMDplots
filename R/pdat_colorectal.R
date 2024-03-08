# JMDplots/pdat_colorectal.R
# get protein data for colorectal cancer
# 20160703 jmd
# 20161011 updated with new data [LXM+16]; add =AD tag (adenoma as n2)
# 20170904 add =NT tag (normal tissue as n1)
# 20190318-20200117 updates for 2020 compilation
# 20240426 moved from canprot package

pdat_colorectal <- function(dataset = 2020) {
  # list available datasets in 2020 compilation
  if(identical(dataset, 2020)) { 
    return(c(
             "WTK+08",
             "XZC+10_I", "XZC+10_II", "ZYS+10",
             "BPV+11_stage.I", "BPV+11_stage.II", "BPV+11_stage.III", "BPV+11_stage.IV",
             "JCF+11", "MRK+11_AC.NC", "SHHS11",
             "FGW+12", "KYK+12", "WOD+12",
             "CZD+14",
             "STK+15", "WDO+15_C.N",
             "LXM+16", "PHL+16_CIS", "PHL+16_ICC",
             "CTW+17", "HZW+17", "LLL+17", "NKG+17", "QMB+17", "TMS+17", "ZLY+17",
             "AKG+18",
             "STA+19_CC.NM", "STA+19_CC.M", "VHW+19", "WYL+19"
             ))
  }
  # list available datasets in 2017 compilation
  if(identical(dataset, 2017)) { 
    return(c("WTK+08=NT",
             "AKP+10_CRC", "AKP+10_CIN", "AKP+10_MIN", "JKMF10", "XZC+10_I=NT", "XZC+10_II=NT", "ZYS+10=NT",
             "BPV+11_adenoma=AD=NT", "BPV+11_stage.I=NT", "BPV+11_stage.II=NT", "BPV+11_stage.III=NT", "BPV+11_stage.IV=NT",
             "JCF+11=NT", "MRK+11_AD.NC=AD=NT", "MRK+11_AC.AD", "MRK+11_AC.NC=NT",
             "KKL+12", "KYK+12=NT", "WOD+12=NT", "YLZ+12",
             "MCZ+13=NT",
             "KWA+14", "UNS+14=AD=NT", "WKP+14",
             "STK+15=NT", "WDO+15_A.N=AD=NT", "WDO+15_C.A", "WDO+15_C.N=NT",
             "LPL+16_ACP=AD=NT", "LPL+16_CIS=NT", "LPL+16_ICC=NT", "LXM+16=NT",
             "PHL+16_AD=AD=NT", "PHL+16_CIS=NT", "PHL+16_ICC=NT"))
  }
  # remove tags
  dataset <- strsplit(dataset, "=")[[1]][1]
  # get study and stage/condition
  study <- strsplit(dataset, "_")[[1]][1]
  stage <- paste(strsplit(dataset, "_")[[1]][-1], collapse="_")
  extdatadir <- system.file("extdata", package="JMDplots")
  datadir <- paste0(extdatadir, "/diffexpr/colorectal/")
  if(study=="JKMF10") {
    # 20150520 up- and down-regulated CRC-associated proteins reported in 4 or more studies, from Jimenez et al., 2010
    description <- "serum biomarkers up / down"
    dat <- read.csv(paste0(datadir, "JKMF10.csv.xz"), as.is=TRUE)
    # get amino acid compositions
    pcomp <- protcomp(dat$Uniprot.ID)
    up2 <- dat$Change=="UP"
  } else if(study=="KWA+14") {
    # 20150908 chromatin-binding fraction, Knol et al., 2014
    dat <- read.csv(paste0(datadir, "KWA+14.csv.xz"), as.is=TRUE)
    description <- "chromatin-binding C / A"
    pcomp <- protcomp(dat$Entry)
    up2 <- dat$Fold.change > 0 | dat$Only.in == "CRC"
    up2[is.na(up2)] <- FALSE
  } else if(study=="STK+15") {
    # 20151004 CRC membrane-enriched proteome, Sethi et al., 2015
    dat <- read.csv(paste0(datadir, "STK+15.csv.xz"), as.is=TRUE)
    description <- "membrane enriched T / N"
    up2 <- dat$invratio > 1
    dat <- cleanup(dat, "uniprot", up2)
    pcomp <- protcomp(dat$uniprot)
  } else if(study=="UNS+14") {
    # 20151005 epithelial cell signature, Uzozie et al., 2014
    dat <- read.csv(paste0(datadir, "UNS+14.csv.xz"), as.is=TRUE)
    description <- "epithelial adenoma / normal"
    dat <- check_IDs(dat, "uniprot")
    up2 <- dat$log2_fold > 0
    dat <- cleanup(dat, "uniprot", up2)
    pcomp <- protcomp(dat$uniprot)
  } else if(study=="BPV+11") {
    # 20160414 CRC Besson et al., 2015
    # BPV+11_adenoma, BPV+11_stage.I, BPV+11_stage.II, BPV+11_stage.III, BPV+11_stage.IV
    dat <- read.csv(paste0(datadir, "BPV+11.csv.xz"), as.is=TRUE)
    description <- paste0(gsub("\\.", " ", stage), " / normal")
    # keep signifcantly changed proteins for the cancer stage
    istage <- match(tolower(stage), tolower(colnames(dat)))
    dat <- dat[dat[, istage] != 1, ]
    pcomp <- protcomp(dat$Entry)
    up2 <- dat[, istage] > 1
  } else if(study=="WDO+15") {
    # 20160414 Wisniewski et al., 2015
    # WDO+15_A.N, WDO+15_C.A, WDO+15_C.N
    dat <- read.csv(paste0(datadir, "WDO+15.csv.xz"), as.is=TRUE)
    if(stage=="A.N") description <- "LCM FFPE adenoma / normal"
    if(stage=="C.A") description <- "LCM FFPE carcinoma / adenoma"
    if(stage=="C.N") description <- "LCM FFPE T / adjacent N"
    # columns with the signficance and ratio
    isig <- grep(paste("Significant", stage, sep="."), colnames(dat), fixed=TRUE)
    irat <- grep(paste("ratio", stage, sep="."), colnames(dat), fixed=TRUE)
    # keep only the proteins marked with a significant change
    dat <- dat[dat[, isig]=="+", ]
    # remove the CON__ prefix
    dat$Majority.protein.IDs <- gsub("CON__", "", dat$Majority.protein.IDs)
    dat <- check_IDs(dat, "Majority.protein.IDs")
    up2 <- dat[, irat] > 0
    dat <- cleanup(dat, "Majority.protein.IDs", up2)
    pcomp <- protcomp(dat$Majority.protein.IDs)
  } else if(study=="WOD+12") {
    # 20160418 CRC tumor tissue, Wisniewski et al., 2012
    dat <- read.csv(paste0(datadir, "WOD+12.csv.xz"), as.is=TRUE)
    description <- "LCM FFPE T / adjacent N"
    dat <- check_IDs(dat, "Uniprot")
    up2 <- dat$Median.Ratio.C.N > 1
    dat <- cleanup(dat, "Uniprot", up2)
    pcomp <- protcomp(dat$Uniprot)
  } else if(study=="JCF+11") {
    # 20160422 tumor vs normal, Jankova et al., 2011
    dat <- read.csv(paste0(datadir, "JCF+11.csv.xz"), as.is=TRUE)
    description <- "T / N"
    pcomp <- protcomp(dat$Accession)
    up2 <- dat$Av..Fold.Change > 0
  } else if(study=="XZC+10") {
    # 20160426 stage I and II vs normal, Xie et al., 2010
    # XZC+10_I, XZC+10_II
    dat <- read.csv(paste0(datadir, "XZC+10.csv.xz"), as.is=TRUE)
    description <- paste("stage", stage, "/ normal")
    # use data for the specified stage
    icol <- grep(paste0("Log2.", stage, ".N"), colnames(dat))
    dat <- dat[!is.na(dat[, icol]), ]
    dat <- check_IDs(dat, "UniProt")
    up2 <- dat[, icol] > 0
    dat <- cleanup(dat, "UniProt", up2)
    pcomp <- protcomp(dat$UniProt)
  } else if(study=="AKP+10") {
    # 20160427 adenoma ADE vs CRC, CIN, MIN, Albrethsen et al., 2010
    # AKP+10_CRC, AKP+10_CIN, AKP+10_MIN
    dat <- read.csv(paste0(datadir, "AKP+10.csv.xz"), as.is=TRUE)
    description <- paste(stage, "nuclear matrix C / A")
    # use the specified data set
    icol <- grep(paste0("Fold.Change.ADE.", stage), colnames(dat))
    dat <- dat[!is.na(dat[, icol]), ]
    up2 <- dat[, icol] > 0
    dat <- cleanup(dat, "Entry", up2)
    pcomp <- protcomp(dat$Entry)
  } else if(study=="KKL+12") {
    # 20160428 poor / good prognosis, Kim et al., 2012
    dat <- read.csv(paste0(datadir, "KKL+12.csv.xz"), as.is=TRUE)
    description <- "poor / good prognosis"
    dat <- check_IDs(dat, "UniProt")
    up2 <- dat$protein.ratio..G.P. < 1
    dat <- cleanup(dat, "UniProt", up2)
    pcomp <- protcomp(dat$UniProt)
  } else if(study=="WKP+14") {
    # 20160428 tissue secretome, de Wit et al., 2014
    dat <- read.csv(paste0(datadir, "WKP+14.csv.xz"), as.is=TRUE)
    description <- "tissue secretome T / N"
    pcomp <- protcomp(dat$Entry)
    up2 <- dat$Fold.change > 0
  } else if(study=="KYK+12") {
    # 20160428 MSS-type CRC, Kang et al., 2012
    dat <- read.csv(paste0(datadir, "KYK+12.csv.xz"), as.is=TRUE)
    description <- "MSS-type T / N"
    pcomp <- protcomp(dat$Entry)
    up2 <- dat$mTRAQ.ratio..N.C.a < 0.5
  } else if(study=="ZYS+10") {
    # 20160430 microdissected T / N, Zhang et al., 2010
    dat <- read.csv(paste0(datadir, "ZYS+10.csv.xz"), as.is=TRUE)
    description <- "microdissected T / N"
    dat <- check_IDs(dat, "UniProt")
    up2 <- dat$Ratio..cancer.normal. > 1
    dat <- cleanup(dat, "UniProt", up2)
    pcomp <- protcomp(dat$UniProt)
  } else if(study=="MRK+11") {
    # 20160509 T / N, Mikula et al., 2011
    # MRK+11_AD.NC, MRK+11_AC.AD, MRK+11_AC.NC
    dat <- read.csv(paste0(datadir, "MRK+11.csv.xz"), as.is=TRUE)
    if(stage=="AD.NC") description <- "adenoma / normal"
    if(stage=="AC.AD") description <- "adenocarcinoma / adenoma"
    if(stage=="AC.NC") description <- "adenocarcinoma / normal"
    # the fold change and false discovery rate columns
    iFC <- grep(paste0(stage, ".FC"), colnames(dat))
    iFDR <- grep(paste0(stage, ".FDR"), colnames(dat))
    # apply the cutoffs
    isFC <- dat[, iFC] >= 3/2 | dat[, iFC] <= 2/3
    isFDR <- dat[, iFDR] <= 0.01
    dat <- dat[isFC & isFDR, ]
    dat <- check_IDs(dat, "Swiss.ID")
    up2 <- dat[, iFC] > 1
    dat <- cleanup(dat, "Swiss.ID", up2)
    pcomp <- protcomp(dat$Swiss.ID)
  } else if(study=="YLZ+12") {
    # 20160511 conditioned media T / N, Yao et al., 2012
    dat <- read.csv(paste0(datadir, "YLZ+12.csv.xz"), as.is=TRUE)
    description <- "CM T / N"
    dat <- check_IDs(dat, "UniProt")
    pcomp <- protcomp(dat$UniProt)
    up2 <- dat$Rsca > 0
  } else if(study=="WTK+08") {
    # 20160511 T / N, Watanabe et al., 2008
    dat <- read.csv(paste0(datadir, "WTK+08.csv.xz"), as.is=TRUE)
    description <- "T / N"
    dat <- check_IDs(dat, "Accession.No.")
    up2 <- dat$Average.T.N.ratio > 1
    dat <- cleanup(dat, "Accession.No.", up2)
    pcomp <- protcomp(dat$Accession.No.)
  } else if(study=="PHL+16") {
    # 20160602 AD/NC, CIS/NC, ICC/NC, Peng et al., 2016
    # PHL+16_AD, PHL+16_CIS, PHL+16_ICC
    dat <- read.csv(paste0(datadir, "PHL+16.csv.xz"), as.is=TRUE)
    if(stage=="AD") description <- "AD / N"
    if(stage=="CIS") description <- "CIS / N"
    if(stage=="ICC") description <- "ICC / N"
    # AD/NC, CIS/NC, ICC/NC
    if(stage=="AD") ratio <- 2^dat$log.of.113.114 
    if(stage=="CIS") ratio <- 2^dat$log.of.115.114
    if(stage=="ICC") ratio <- 2^dat$log.of.116.114
    idiff <- ratio < 0.67 | ratio > 1.5
    dat <- dat[idiff, ]
    pcomp <- protcomp(dat$Entry)
    up2 <- ratio[idiff] > 1.5
  } else if(study=="LPL+16") {
    # 20160602 stromal ACP/NNCM, CIS/NNCM, ICC/NNCM, Li et al., 2016
    # LPL+16_ACP, LPL+16_CIS, LPL+16_ICC
    dat <- read.csv(paste0(datadir, "LPL+16.csv.xz"), as.is=TRUE)
    if(stage=="ACP") description <- "stromal AD / NC"
    if(stage=="CIS") description <- "stromal CIS / NC"
    if(stage=="ICC") description <- "stromal ICC / NC"
    icol <- grep(stage, colnames(dat))
    # keep only significantly changed proteins
    dat <- dat[dat[, icol] < 0.67 | dat[, icol] > 1.5, ]
    pcomp <- protcomp(dat$Entry)
    up2 <- dat[, icol] > 1.5
  } else if(study=="MCZ+13") {
    # 20160602 stromal T / N, Mu et al., 2013
    dat <- read.csv(paste0(datadir, "MCZ+13.csv.xz"), as.is=TRUE)
    description <- "stromal T / N"
    up2 <- dat$CS.vs..NS > 0
    dat <- cleanup(dat, "Entry", up2)
    pcomp <- protcomp(dat$Entry)
  } else if(study=="LXM+16") {
    # 20160728 CRC, Liu et al., 2016
    dat <- read.csv(paste0(datadir, "LXM+16.csv.xz"), as.is=TRUE)
    description <- "biopsy T / N"
    dat <- check_IDs(dat, "UniProt")
    up2 <- dat$Ratio..C.N. > 1
    dat <- cleanup(dat, "UniProt", up2)
    pcomp <- protcomp(dat$UniProt)
  } else if(study=="CTW+17") {
    # 20190318 colorectal organoid, Cristobal et al., 2017
    dat <- read.csv(paste0(datadir, "CTW+17.csv.xz"), as.is=TRUE)
    description <- "organoid T / N"
    dat <- check_IDs(dat, "Accession")
    pcomp <- protcomp(dat$Accession)
    up2 <- dat$Tumor.Healthy > 0
  } else if(study=="SHHS11") {
    # 20160422 tumor vs normal, Shi et al., 2011
    # re-added 20190326
    dat <- read.csv(paste0(datadir, "SHHS11.csv.xz"), as.is=TRUE)
    description <- "LCM T / N"
    pcomp <- protcomp(dat$Accession.no.)
    up2 <- dat$Fold.of.change > 0
  } else if(study=="CZD+14") {
    # 20190328 tumor vs normal, Chen et al., 2014
    dat <- read.csv(paste0(datadir, "CZD+14.csv.xz"), as.is=TRUE)
    description <- "T / N"
    pcomp <- protcomp(dat$Entry)
    up2 <- dat$CRC.vs..Normal > 1
  } else if(study=="NKG+17") {
    # 20190328 tumor vs normal, Nishio et al., 2017
    dat <- read.csv(paste0(datadir, "NKG+17.csv.xz"), as.is=TRUE)
    description <- "T / N"
    up2 <- dat$Average.ratio.of.115.114..cancer.normal. > 1
    dat <- cleanup(dat, "Accession.no.", up2)
    pcomp <- protcomp(dat$Accession.no.)
  } else if(study=="QMB+17") {
    # 20190331 cancer vs control, Quesada-Calvo et al., 2017
    dat <- read.csv(paste0(datadir, "QMB+17.csv.xz"), as.is=TRUE)
    description <- "FFPE T / N"
    # use first Uniprot ID
    dat$Uniprot.Accession.number. <- sapply(strsplit(dat$Uniprot.Accession.number., ";"), "[", 1)
    pcomp <- protcomp(dat$Uniprot.Accession.number.)
    up2 <- grepl("ADK", dat$Highest.mean.condition)
  } else if(study=="ZLY+17") {
    # 20190403 tumor / normal, Zhang et al., 2017
    dat <- read.csv(paste0(datadir, "ZLY+17.csv.xz"), as.is=TRUE)
    description <- "T / N"
    dat <- check_IDs(dat, "UniProt")
    up2 <- dat$Ratio..cancer.normal. > 1
    dat <- cleanup(dat, "UniProt", up2)
    pcomp <- protcomp(dat$UniProt)
  } else if(study=="TMS+17") {
    # 20191203 cancer / normal, Tu et al., 2017
    dat <- read.csv(paste0(datadir, "TMS+17.csv.xz"), as.is=TRUE)
    description <- "epithelial T / N"
    up2 <- dat$Regulation=="up"
    pcomp <- protcomp(dat$Accession.Number)
  } else if(study=="AKG+18") {
    # 20191223 cancer / normal, Atak et al., 2018
    dat <- read.csv(paste0(datadir, "AKG+18.csv.xz"), as.is=TRUE)
    description <- "T / N"
    # keep highly differential proteins
    dat <- dat[dat$Average.Fold.change..n.11. > 1.5 | dat$Average.Fold.change..n.11. < 1/1.5, ]
    up2 <- dat$Average.Fold.change..n.11. > 1.5
    pcomp <- protcomp(dat$accession_number)
  } else if(study=="HZW+17") {
    # 20191224 cancer / adjacent normal, Hao et al., 2017
    dat <- read.csv(paste0(datadir, "HZW+17.csv.xz"), as.is=TRUE)
    description <- "T / adjacent N"
    dat <- check_IDs(dat, "SP")
    up2 <- dat$Ratio > 1
    pcomp <- protcomp(dat$SP)
  } else if(study=="LLL+17") {
    # 20191228 LCM cancer / non-neoplastic mucosa, Li et al., 2017
    dat <- read.csv(paste0(datadir, "LLL+17.csv.xz"), as.is=TRUE)
    description <- "LCM cancer / non-neoplastic mucosa"
    dat$Accession <- sapply(strsplit(dat$Accession, "\\|"), "[", 2)
    dat <- check_IDs(dat, "Accession")
    up2 <- dat$coloncancer.NNCM > 1
    dat <- cleanup(dat, "Accession", up2)
    pcomp <- protcomp(dat$Accession)
  } else if(study=="WYL+19") {
    # 20191228 tumor-associated / normal vascular endothelial cells, Wang et al., 2019
    dat <- read.csv(paste0(datadir, "WYL+19.csv.xz"), as.is=TRUE)
    description <- "tumor-associated / normal vascular endothelial cells"
    dat <- check_IDs(dat, "Accession")
    up2 <- dat$Ratio.T.N. > 1
    pcomp <- protcomp(dat$Accession)
  } else if(study=="STA+19") {
    # 20200117 colorectal cancer, Saleem et al., 2019
    # STA+19_NAP, STA+19_CC.NM, STA+19_CC.M
    dat <- read.csv(paste0(datadir, "STA+19.csv.xz"), as.is=TRUE)
    if(stage == "NAP") description <- "non-cancer non-adenomatous polyp"
    if(stage == "CC.NM") description <- "non-metastatic colon cancer, T / N"
    if(stage == "CC.M") description <- "metastatic colon cancer, T / N"
    icol <- match(stage, colnames(dat))
    dat <- dat[!is.na(dat[, icol]), ]
    up2 <- dat[, icol] > 0
    pcomp <- protcomp(dat$Accession)
  } else if(study=="VHW+19") {
    # 20200328 colorectal cancer, Vasaikar et al., 2019
    dat <- read.csv(paste0(datadir, "VHW+19.csv.xz"), as.is=TRUE)
    description <- "T / adjacent N"
    dat <- dat[abs(dat$median_log2FC) > 1, ]
    up2 <- dat$median_log2FC > 1
    pcomp <- protcomp(dat$Entry)
  } else if(study=="FGW+12") {
    # 20200422 colorectal cancer, Fan et al., 2012
    dat <- read.csv(paste0(datadir, "FGW+12.csv.xz"), as.is=TRUE)
    description <- "T / matched N"
    up2 <- dat$expression.level.in.cancer == "Increased"
    dat <- check_IDs(dat, "UniProt")
    dat <- cleanup(dat, "UniProt", up2)
    pcomp <- protcomp(dat$UniProt)
  } else stop(paste("colorectal dataset", dataset, "not available"))
  print(paste0("pdat_colorectal: ", description, " [", dataset, "]"))
  # use the up2 from the cleaned-up data, if it exists 20190429
  if("up2" %in% colnames(dat)) up2 <- dat$up2
  return(list(dataset=dataset, description=description, pcomp=pcomp, up2=up2))
}
