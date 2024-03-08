# JMDplots/pdat_lung.R
# retrieve protein IDs for lung cancer studies
# 20160720-20191230 assemble data for 2020 compilation
# 20240426 moved from canprot package

pdat_lung <- function(dataset = 2020) {
  if(identical(dataset, 2020)) {
    return(c(
             "LXC+06",
             "KHA+12_ADC", "KHA+12_SCC", "YLL+12", "ZZD+12",
             "ZZY+13",
             "LLY+14", "LWT+14", "ZLH+14", "ZLS+14",
             "KNT+15",
             "BLL+16_protein", "FGP+16", "JCP+16=mouse",
             "HHH+16_pN0", "HHH+16_pN1", "HHH+16_pN2.M1", "TLB+16",
             "FGW+17", "LZW+17", "SFS+17_LF", "WLC+17",
             "YCC+17_SqCC.Oncogene", "YCC+17_SqCC.TSG", "YCC+17_SqCC.Glycoprotein",
             "YCC+17_ADC.Oncogene", "YCC+17_ADC.TSG", "YCC+17_ADC.Glycoprotein",
             "KPS+20_early", "KPS+20_advanced", "XZW+20"
             ))
  }
  # remove tags
  dataset <- strsplit(dataset, "=")[[1]][1]
  # get study and stage/condition
  study <- strsplit(dataset, "_")[[1]][1]
  stage <- paste(strsplit(dataset, "_")[[1]][-1], collapse="_")
  extdatadir <- system.file("extdata", package="JMDplots")
  datadir <- paste0(extdatadir, "/diffexpr/lung/")
  if(study=="HHH+16") {
    # 20160720 lung stages, Hu et al., 2016
    # HHH+16_pN0, HHH+16_pN1, HHH+16_pN2.M1
    dat <- read.csv(paste0(datadir, "HHH+16.csv.xz"), as.is=TRUE)
    description <- paste(stage, "/ normal")
    # get data for the selected experiment
    icol <- match(paste0(stage, ".Normal"), colnames(dat))
    # include proteins with a large expression change
    dat <- dat[dat[, icol] <= 0.66 | dat[, icol] >= 1.5, ]
    dat <- check_IDs(dat, "Uniprot.accession.No.")
    pcomp <- protcomp(dat$Uniprot.accession.No.)
    up2 <- dat[, icol] > 1
  } else if(study=="YLL+12") {
    # 20170113 human lung squamous carcinoma, Yan et al., 2012
    dat <- read.csv(paste0(datadir, "YLL+12.csv.xz"), as.is=TRUE)
    description <- "LCM SCC / NBE"
    dat <- check_IDs(dat, "AC")
    up2 <- dat$X115.113 > 1
    dat <- cleanup(dat, "AC", up2)
    pcomp <- protcomp(dat$AC)
  } else if(study=="ZLH+14") {
    # 20170113 lung adenocarcinoma, Zhang et al., 2014
    dat <- read.csv(paste0(datadir, "ZLH+14.csv.xz"), as.is=TRUE)
    description <- "membrane microdissected ADC / ANT"
    dat <- check_IDs(dat, "Accession")
    up2 <- dat$X117.118 > 1
    dat <- cleanup(dat, "Accession", up2)
    pcomp <- protcomp(dat$Accession)
  } else if(study=="LXC+06") {
    # 20190317 lung squamous carcinoma, Li et al., 2006
    dat <- read.csv(paste0(datadir, "LXC+06.csv.xz"), as.is=TRUE)
    description <- "SCC / NBE"
    dat <- check_IDs(dat, "Accession.no.")
    up2 <- grepl("^T", dat$Spot.no.)
    dat <- cleanup(dat, "Accession.no.", up2)
    pcomp <- protcomp(dat$Accession.no.)
  } else if(study=="FGP+16") {
    # 20190318 lung adenocarcinoma, Fahrmann et al., 2016
    dat <- read.csv(paste0(datadir, "FGP+16.csv.xz"), as.is=TRUE)
    description <- "adenocarcinoma / ANT"
    up2 <- dat$Fold.Change > 1
    dat <- cleanup(dat, "Entry", up2)
    pcomp <- protcomp(dat$Entry)
  } else if(study=="JCP+16") {
    # 20190318 mouse tumor / normal endothelial cells, Jin et al., 2016
    dat <- read.csv(paste0(datadir, "JCP+16.csv.xz"), as.is=TRUE)
    description <- "mouse endothelial tumor / normal"
    up2 <- dat$Ratio.TEC.NEC > 1
    dat <- cleanup(dat, "Entry", up2)
    pcomp <- protcomp(dat$Entry, aa_file=paste0(extdatadir, "/aa/mouse/JCP+16_aa.csv.xz"))
  } else if(study=="FGW+17") {
    # 20190325 lung adenocarcinoma, Fahrmann et al., 2017
    dat <- read.csv(paste0(datadir, "FGW+17.csv.xz"), as.is=TRUE)
    description <- "adenocarcinoma / ANT"
    pcomp <- protcomp(dat$Entry)
    up2 <- dat$Fold.Change..Tumor.Control. > 1
  } else if(study=="KHA+12") {
    # 20190325 ADC or SCC versus normal, Kikuchi et al., 2012
    # KHA+12_ADC, KHA+12_SCC
    dat <- read.csv(paste0(datadir, "KHA+12.csv.xz"), as.is=TRUE)
    description <- paste(stage, "/ pooled normal")
    # use ADC or SCC dataset
    dat <- dat[dat$Higher.in == stage | dat$Lower.in == stage, ]
    up2 <- dat$Higher.in == stage
    dat <- cleanup(dat, "Entry", up2)
    pcomp <- protcomp(dat$Entry)
  } else if(study=="LLY+14") {
    # 20190325 SCC / normal, Lihong et al., 2014
    dat <- read.csv(paste0(datadir, "LLY+14.csv.xz"), as.is=TRUE)
    description <- "SCC / normal"
    # remove version suffixes
    dat$Protein.ID <- sapply(strsplit(dat$Protein.ID, "\\."), "[", 1)
    dat <- check_IDs(dat, "Protein.ID")
    up2 <- dat$Average.Ratio > 0
    dat <- cleanup(dat, "Protein.ID", up2)
    pcomp <- protcomp(dat$Protein.ID)
  } else if(study=="LWT+14") {
    # 20190325 NSCLC tumor / normal, Li et al., 2014
    dat <- read.csv(paste0(datadir, "LWT+14.csv.xz"), as.is=TRUE)
    description <- "NSCLC / ANT"
    dat <- check_IDs(dat, "Swissprot.ID")
    pcomp <- protcomp(dat$Swissprot.ID)
    up2 <- dat$logFC > 0
  } else if(study=="SFS+17") {
    # 20190326 lung tumor / normal (method comparison), Stewart et al., 2017
    # SFS+17_LF, SFS+17_TMT, SFS+17_DIA
    dat <- read.csv(paste0(datadir, "SFS+17.csv.xz"), as.is=TRUE)
    description <- paste("SCC / ANT", stage)
    # use data for specified method
    icol <- grep(stage, colnames(dat))
    up2 <- dat[, icol] > 0
    dat <- cleanup(dat, "UniProt", up2)
    pcomp <- protcomp(dat$UniProt)
  } else if(study=="TLB+16") {
    # 20190402 tumor/control, Tenzer et al., 2016
    dat <- read.csv(paste0(datadir, "TLB+16.csv.xz"), as.is=TRUE)
    description <- "NSCLC / ANT"
    # ratios in table are control/cancer: up in cancer is less than 1
    up2 <- dat$fold.change < 1
    dat <- cleanup(dat, "Entry", up2)
    pcomp <- protcomp(dat$Entry)
  } else if(study=="WLC+17") {
    # 20190403 NSCLC, Wang et al., 2017
    dat <- read.csv(paste0(datadir, "WLC+17.csv.xz"), as.is=TRUE)
    description <- "NSCLC / ANT"
    # use first protein ID
    dat$Entry <- sapply(strsplit(dat$Entry, ";"), "[", 1)
    pcomp <- protcomp(dat$Entry)
    up2 <- dat$Log2.ratio..N.T.. < 0
  } else if(study=="YCC+17") {
    # 20190403 SqCC and ADC, TSG/oncogenes and glycproteins, Yang et al., 2017
    # YCC+17_SqCC.Oncogene, YCC+17_SqCC.TSG, YCC+17_SqCC.Glycoprotein
    # YCC+17_ADC.Oncogene, YCC+17_ADC.TSG, YCC+17_ADC.Glycoprotein
    dat <- read.csv(paste0(datadir, "YCC+17.csv.xz"), as.is=TRUE)
    description <- paste(gsub("SqCC", "SCC", stage), "/ ANT")
    # use specified dataset
    icol <- grep(stage, colnames(dat))
    dat <- dat[!is.na(dat[, icol]), ]
    pcomp <- protcomp(dat$Entry)
    up2 <- dat[, icol] >= 1.5
  } else if(study=="ZZD+12") {
    # 20190408 SCC, Zeng et al., 2012
    dat <- read.csv(paste0(datadir, "ZZD+12.csv.xz"), as.is=TRUE)
    description <- "LCM LSCC / NBE"
    dat <- check_IDs(dat, "UniProt")
    up2 <- dat$NBE.vs..LSCC > 1
    dat <- cleanup(dat, "UniProt", up2)
    pcomp <- protcomp(dat$UniProt)
  } else if(study=="ZZY+13") {
    # 20190408 plasma membrane AdC vs normal, Zhang et al., 2013
    dat <- read.csv(paste0(datadir, "ZZY+13.csv.xz"), as.is=TRUE)
    description <- "plasma membrane ADC / ANT"
    pcomp <- protcomp(dat$Accession.no.)
    up2 <- dat$AdC.vs..PNLT > 1
  } else if(study=="BLL+16") {
    # 20191228 non small-cell lung cancer transcriptomics and proteomics, Backes et al., 2016
    # BLL+16_gene, BLL+16_protein
    dat <- read.csv(paste0(datadir, "BLL+16.csv.xz"), as.is=TRUE)
    if(stage=="gene") description <- "NSCLC transcriptome"
    if(stage=="protein") description <- "NSCLC proteome"
    # use selected dataset
    dat <- dat[!is.na(dat[, stage]), ]
    dat <- dat[abs(dat$log.2..fold.change) > 1, ]
    dat <- check_IDs(dat, "Entry")
    up2 <- dat$log.2..fold.change > 0
    dat <- cleanup(dat, "Entry", up2)
    pcomp <- protcomp(dat$Entry)
  } else if(study=="KNT+15") {
    # 20191228 FFPE lepidic predominant invasive adenocarcinoma / pseudo-Normal, Kato et al., 2015
    dat <- read.csv(paste0(datadir, "KNT+15.csv.xz"), as.is=TRUE)
    description <- "FFPE LPIA / pseudo-normal"
    # remove HIP000 accession nubmers
    dat$Accession.Number.Code[grepl("HIP000", dat$Accession.Number.Code)] <- NA
    # remove _1 suffixes
    dat$Accession.Number.Code <- sapply(strsplit(dat$Accession.Number.Code, "_"), "[", 1)
    up2 <- dat$Fold.changeLPIA.vs.pN > 0
    dat <- cleanup(dat, "Accession.Number.Code", up2)
    pcomp <- protcomp(dat$Accession.Number.Code)
  } else if(study=="LZW+17") {
    # 20191228 mitochondria-related proteins adenocarcinoma / normal, Li et al., 2017
    dat <- read.csv(paste0(datadir, "LZW+17.csv.xz"), as.is=TRUE)
    description <- "mitochondria-related proteins adenocarcinoma / normal"
    dat <- check_IDs(dat, "Accession.number")
    up2 <- dat$C.N > 1
    dat <- cleanup(dat, "Accession.number", up2)
    pcomp <- protcomp(dat$Accession.number)
  } else if(study=="ZLS+14") {
    # 20191230 lung SCC, Zhou et al., 2014
    dat <- read.csv(paste0(datadir, "ZLS+14.csv.xz"), as.is=TRUE)
    description <- "endothelial SCC / normal"
    up2 <- dat$Ratio.114.113 > 1.5
    dat <- cleanup(dat, "Entry", up2)
    pcomp <- protcomp(dat$Entry)
  } else if(study=="KPS+20") {
    # 20201016 lung cancer, Kelemen et al., 2020
    # KPS+20_early, KPS+20_advanced
    dat <- read.csv(paste0(datadir, "KPS+20.csv.xz"), as.is=TRUE)
    description <- paste("ADC / adjacent normal", stage)
    icol <- grep(stage, colnames(dat))
    dat <- dat[!is.na(dat[, icol]), ]
    up2 <- dat[, icol] >= 2
    dat <- cleanup(dat, "Accession", up2)
    pcomp <- protcomp(dat$Accession)
  } else if(study=="XZW+20") {
    # 20201016 lung cancer, Xu et al., 2020
    dat <- read.csv(paste0(datadir, "XZW+20.csv.xz"), as.is=TRUE)
    description <- "ADC / non-cancerous adjacent"
    dat <- check_IDs(dat, "Majority.protein.IDs")
    up2 <- dat$NAT.vs.tumor == "up"
    dat <- cleanup(dat, "Majority.protein.IDs", up2)
    pcomp <- protcomp(dat$Majority.protein.IDs)
  } else stop(paste("lung dataset", dataset, "not available"))
  print(paste0("pdat_lung: ", description, " [", dataset, "]"))
  # use the up2 from the cleaned-up data, if it exists 20190407
  if("up2" %in% colnames(dat)) up2 <- dat$up2
  return(list(dataset=dataset, pcomp=pcomp, up2=up2, description=description))
}
