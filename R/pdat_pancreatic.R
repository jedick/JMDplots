# JMDplots/pdat_pancreatic.R
# retrieve protein IDs for pancreatic cancer studies
# 20160827 jmd
# 20170904 add =NT tag (comparison between cancer and normal tissue)
# 20190239-20191206 updates for 2020 compilation
# 20240426 moved from canprot package

pdat_pancreatic <- function(dataset = 2020) {
  if(identical(dataset, 2020)) {
    return(c("LHE+04",
             "CYD+05", "CGB+05",
             "CTZ+09",
             "MLC+11", "PCS+11_PDAC", "TMW+11",
             "KBK+12",
             "KHO+13", "KPC+13_all", "WLL+13a_PC_NT", "WLL+13a_PC.DM_NT.DM", "YKK+13", "ZNWL13",
             "ISI+14", "BZQ+14", "MZH+14=mouse",
             "BHB+15_T=mouse",
             "KKC+16_T1=mouse",
             "CHO+18", "SWW+18",
             "ZAH+19"
             ))
  }
  if(identical(dataset, 2017)) {
    return(c("LHE+04=NT",
             "CYD+05=NT", "CGB+05=NT",
             "CBP+07=low",
             "CTZ+09=NT",
             "MLC+11=NT", "PCS+11_PDAC=NT", #"PCS+11_MCP", "PCS+11_SCP", 
             "TMW+11=NT",
             "KBK+12=NT",
             "KHO+13=NT", "KPC+13_all=NT", #"KPC+13_2-fold-signif=NT",
             "PKB+13_AIP", "PKB+13_CP",
             "WLL+13_low=low=NT", "WLL+13_high=NT", "WLL+13a_PC_NT=NT", "WLL+13a_PC.DM_NT.DM=NT",
             "ZNWL13=NT", "ISI+14=NT",
             "KKC+16_T4=mouse=NT", "KKC+16_T3=mouse=NT", "KKC+16_T2=mouse=NT", "KKC+16_T1=mouse=NT"))
  }
  # remove tags
  dataset <- strsplit(dataset, "=")[[1]][1]
  # get study and stage/condition
  study <- strsplit(dataset, "_")[[1]][1]
  stage <- paste(strsplit(dataset, "_")[[1]][-1], collapse="_")
  extdatadir <- system.file("extdata", package="JMDplots")
  datadir <- paste0(extdatadir, "/canH2O/pancreatic/")
  if(study=="KKC+16") {
    # 20160717 mouse PDAC, Kuo et al., 2016
    # KKC+16_T1, KKC+16_T2, KKC+16_T3, KKC+16_T4 (10, 5, 3.5, 2.5 weeks)
    dat <- read.csv(paste0(datadir, "KKC+16.csv.xz"), as.is=TRUE)
    if(stage=="T1") description <- "mouse 10 w T / N"
    if(stage=="T2") description <- "mouse 5 w T / N"
    if(stage=="T3") description <- "mouse 3.5 w T / N"
    if(stage=="T4") description <- "mouse 2.5 w T / N"
    # use only proteins differentially expressed at this stage
    icol <- grep(paste0(stage, ".N"), colnames(dat))
    dat <- dat[!is.na(dat[, icol]), ]
    pcomp <- protcomp(dat$Entry, aa_file=paste0(extdatadir, "/aa/mouse/KKC+16_aa.csv.xz"))
    up2 <- dat[, icol] > 1
  } else if(study=="ISI+14") {
    # 20160827 PDAC, Iuga et al., 2014
    dat <- read.csv(paste0(datadir, "ISI+14.csv.xz"), as.is=TRUE)
    description <- "T / adjacent N"
    dat <- check_IDs(dat, "Protein.accesion.number")
    pcomp <- protcomp(dat$Protein.accesion.number)
    up2 <- dat$Simple.ratio > 1
  } else if(study=="MLC+11") {
    # 20160827 PDAC, McKinney et al., 2011
    dat <- read.csv(paste0(datadir, "MLC+11.csv.xz"), as.is=TRUE)
    description <- "T / adjacent N"
    dat <- check_IDs(dat, "UniProt")
    up2 <- dat$Regulated == "up"
    dat <- cleanup(dat, "UniProt", up2)
    pcomp <- protcomp(dat$UniProt)
  } else if(study=="PCS+11") {
    # 20160828 PDAC, Pan et al., 2011
    # PCS+11_MCP, PCS+11_SCP, PCS+11_PDAC
    dat <- read.csv(paste0(datadir, "PCS+11.csv.xz"), as.is=TRUE)
    if(stage=="PDAC") stext <- "T" else stext <- stage
    description <- paste("FFPE", stext, " / N")
    # keep proteins with reported ratio
    icol <- grep(stage, colnames(dat))
    dat <- dat[!is.na(dat[, icol]), ]
    up2 <- dat[, icol] > 1
    dat <- cleanup(dat, "Entry", up2)
    pcomp <- protcomp(dat$Entry)
  } else if(study=="WLL+13") {
    # 20160829 PDAC, Wang et al., 2013
    # WLL+13_low, WLL+13_high
    dat <- read.csv(paste0(datadir, "WLL+13.csv.xz"), as.is=TRUE)
    description <- paste0(stage, "-grade T / N")
    # which columns hold the expression data
    if(stage=="low") icol <- 2:5
    if(stage=="high") icol <- 6:9
    # 2 of the 4 ratios must be >=1.5 or <=0.667
    is.signif <- dat[, icol] >= 1.5 | dat[, icol] <= 0.667
    is.signif <- rowSums(is.signif) >= 2
    # 4 of the 4 ratios must be in the same direction
    is.up <- dat[, icol] > 1
    is.all.up <- apply(is.up, 1, all)
    is.all.down <- apply(!is.up, 1, all)
    # the final list
    is.diff <- is.signif & (is.all.down | is.all.up)
    dat <- dat[is.diff, ]
    pcomp <- protcomp(dat$Entry)
    up2 <- apply(dat[, icol] > 1, 1, all)
  } else if(study=="CGB+05") {
    # 20160829 PDAC, Crnogorac-Jurcevic et al., 2005
    dat <- read.csv(paste0(datadir, "CGB+05.csv.xz"), as.is=TRUE)
    description <- "T / N"
    pcomp <- protcomp(dat$Entry)
    up2 <- dat$Regulated == "up"
  } else if(study=="CTZ+09") {
    # 20160829 PDAC, Cui et al., 2009
    dat <- read.csv(paste0(datadir, "CTZ+09.csv.xz"), as.is=TRUE)
    description <- "T / adjacent N"
    up2 <- dat$C.N > 1
    dat <- cleanup(dat, "Swissprot.ID", up2)
    pcomp <- protcomp(dat$Swissprot.ID)
  } else if(study=="KBK+12") {
    # 20160830 PDAC, Kojima et al., 2012
    dat <- read.csv(paste0(datadir, "KBK+12.csv.xz"), as.is=TRUE)
    description <- "FFPE T / N"
    dat <- check_IDs(dat, "Sequence.Id")
    pcomp <- protcomp(dat$Sequence.Id)
    up2 <- !(grepl("-", dat$Fold.Change..PDAC.Control.) | grepl("Adjacent", dat$Fold.Change..PDAC.Control.))
  } else if(study=="ZNWL13") {
    # 20160830 PDAC, Zhu et al., 2013
    dat <- read.csv(paste0(datadir, "ZNWL13.csv.xz"), as.is=TRUE)
    description <- "LCM T / adjacent N"
    pcomp <- protcomp(dat$Accession)
    up2 <- dat$Up.Down == "Up"
  } else if(study=="KPC+13") {
    # 20160831 Kosanam et al., 2013
    # KPC+13_all, KPC+13_2-fold, KPC+13_2-fold-signif
    dat <- read.csv(paste0(datadir, "KPC+13.csv.xz"), as.is=TRUE)
    if(stage=="all") stext <- "" else stext <- paste0(stage, " ")
    description <- paste0(stext, "T / adjacent N")
    # use all listed proteins or at least 2-fold differentially expressed ones
    if(grepl("2-fold", stage)) dat <- dat[dat$PDAC.Benign.fold.change. >= 2 | dat$PDAC.Benign.fold.change. <= 0.5, ]
    if(grepl("signif", stage)) dat <- dat[dat$t.test.pvalue < 0.1, ]
    up2 <- dat$PDAC.Benign.fold.change. > 1
    dat <- cleanup(dat, "Entry", up2)
    pcomp <- protcomp(dat$Entry)
  } else if(study=="CYD+05") {
    # 20160907 Chen et al., 2005
    dat <- read.csv(paste0(datadir, "CYD+05.csv.xz"), as.is=TRUE)
    description <- "T / N"
    up2 <- dat$Ratio..cancer.normal. > 1
    dat <- cleanup(dat, "Entry", up2)
    pcomp <- protcomp(dat$Entry)
  } else if(study=="CBP+07") {
    # 20160909 Chen et al., 2007
    dat <- read.csv(paste0(datadir, "CBP+07.csv.xz"), as.is=TRUE)
    description <- "CP / N"
    pcomp <- protcomp(dat$Entry)
    up2 <- dat$Ratio.CP.NL > 1
  } else if(study=="KHO+13") {
    # 20160910 Kawahara et al., 2013
    dat <- read.csv(paste0(datadir, "KHO+13.csv.xz"), as.is=TRUE)
    description <- "T / N"
    up2 <- rowMeans(dat[, 6:12]) > 1
    dat <- cleanup(dat, "Entry", up2)
    pcomp <- protcomp(dat$Entry)
  } else if(study=="LHE+04") {
    # 20160910 Lu et al., 2004
    dat <- read.csv(paste0(datadir, "LHE+04.csv.xz"), as.is=TRUE)
    description <- "T / adjacent N"
    dat <- check_IDs(dat, "Acc.no.")
    up2 <- dat$Higher.in == "cancer"
    dat <- cleanup(dat, "Acc.no.", up2)
    pcomp <- protcomp(dat$Acc.no.)
  } else if(study=="PKB+13") {
    # 20160910 Paulo et al., 2013
    # proteins exclusively identified in
    # autoimmune pancreatitis (AIP), chronic pancreatitis (CP), and pancreatic cancer (PC) cohorts
    # PKB+13_AIP, PKB+13_CP
    dat <- read.csv(paste0(datadir, "PKB+13.csv.xz"), as.is=TRUE)
    description <- paste("FFPE PC /", stage)
    # keep only proteins for the indicated comparison
    dat <- dat[dat$Cohort %in% c(stage, "PC"), ]
    dat <- check_IDs(dat, "UniProt.ID")
    pcomp <- protcomp(dat$UniProt.ID)
    up2 <- dat$Cohort == "PC"
  } else if(study=="TMW+11") {
    # 20160910 Turtoi et al., 2011
    dat <- read.csv(paste0(datadir, "TMW+11.csv.xz"), as.is=TRUE)
    description <- "accessible T / N"
    pcomp <- protcomp(dat$Entry)
    up2 <- dat$tumor - dat$normal > 1
  } else if(study=="WLL+13a") {
    # 20161110 PC +/- diabetes mellitus, Wang et al., 2013
    # WLL+13a_PC_NT, WLL+13a_PC.DM_NT.DM
    dat <- read.csv(paste0(datadir, "WLL+13a.csv.xz"), as.is=TRUE)
    description <- stage
    if(stage=="PC_NT") description <- "T / adjacent N with DM"
    if(stage=="PC.DM_NT.DM") description <- "T / adjacent N without DM"
    # which columns hold the expression data
    icol <- grep(stage, colnames(dat))
    # 3 of 4 experiments must have ratios >=1.5 or <=0.667
    nup <- rowSums(dat[, icol] > 3/2, na.rm=TRUE)
    ndn <- rowSums(dat[, icol] < 2/3, na.rm=TRUE)
    is.diff <- (nup >=3 & ndn==0) | (ndn >= 3 & nup==0)
    dat <- dat[is.diff, ]
    pcomp <- protcomp(dat$Entry)
    up2 <- apply(dat[, icol] > 1, 1, sum) >= 3
  } else if(study=="BHB+15") {
    # 20190329 mouse organoids, Boj et al., 2015
    # BHB+15_P, BHB+15_T
    # P - pancreatic epithelial neoplasm organoid
    # T - tumor organoid
    dat <- read.csv(paste0(datadir, "BHB+15.csv.xz"), as.is=TRUE)
    description <- paste("mouse organoids", stage, "/ N")
    # use selected dataset
    icol <- grep(paste0("logFC.", stage, "vsN"), colnames(dat))
    dat <- dat[!is.na(dat[, icol]), ]
    # remove isoform suffixes
    dat$ProteinID <- sapply(strsplit(dat$ProteinID, "-"), "[", 1)
    dat <- check_IDs(dat, "ProteinID", aa_file=paste0(extdatadir, "/aa/mouse/BHB+15_aa.csv.xz"))
    up2 <- dat[, icol] > 0
    dat <- cleanup(dat, "ProteinID", up2)
    pcomp <- protcomp(dat$ProteinID, aa_file=paste0(extdatadir, "/aa/mouse/BHB+15_aa.csv.xz"))
  } else if(study=="BZQ+14") {
    # 20190329 PDAC / normal, Britton et al., 2014
    dat <- read.csv(paste0(datadir, "BZQ+14.csv.xz"), as.is=TRUE)
    description <- "T / matched N"
    dat <- check_IDs(dat, "Uniprot.ID")
    pcomp <- protcomp(dat$Uniprot.ID)
    up2 <- dat$log2.T.NT.ratios > 0 
  } else if(study=="MZH+14") {
    # 20190407 tumor / healthy, Mirus et al., 2014
    dat <- read.csv(paste0(datadir, "MZH+14.csv.xz"), as.is=TRUE)
    description <- "mouse tumor / healthy"
    up2 <- dat$Coefficient > 0
    dat <- cleanup(dat, "Entry", up2)
    pcomp <- protcomp(dat$Entry, aa_file=paste0(extdatadir, "/aa/mouse/MZH+14_aa.csv.xz"))
  } else if(study=="YKK+13") {
    # 20190408 PDAC / normal, Yu et al., 2013
    dat <- read.csv(paste0(datadir, "YKK+13.csv.xz"), as.is=TRUE)
    description <- "T / adjacent N"
    # keep proteins identified in at least 2 patients
    dat <- dat[dat$No..patients > 1, ]
    dat <- check_IDs(dat, "Accession")
    up2 <- dat$Mean > 1
    dat <- cleanup(dat, "Accession", up2)
    pcomp <- protcomp(dat$Accession)
  } else if(study=="SWW+18") {
    # 20191204 PDAC tumor / adjacent, Song et al., 2018
    dat <- read.csv(paste0(datadir, "SWW+18.csv.xz"), as.is=TRUE)
    description <- "T / adjacent N"
    # calculate fold change in each sample
    FC1 <- log2(dat$t1 / dat$n1)
    FC2 <- log2(dat$t2 / dat$n2)
    FC3 <- log2(dat$t3 / dat$n3)
    dat <- cbind(dat, FC1 = FC1, FC2 = FC2, FC3 = FC3)
    # remove proteins with ambiguous (up- and down-) expression
    iup <- dat$FC1 > 0 | dat$FC2 > 0 | dat$FC3 > 0
    idn <- dat$FC1 < 0 | dat$FC2 < 0 | dat$FC3 < 0
    iambi <- iup & idn
    iambi[is.na(iambi)] <- FALSE
    dat <- dat[!iambi, ]
    dat <- check_IDs(dat, "Accession")
    up2 <- dat$FC1 > 0 | dat$FC2 > 0 | dat$FC3 > 0
    up2[is.na(up2)] <- FALSE
    dat <- cleanup(dat, "Accession", up2)
    pcomp <- protcomp(dat$Accession)
  } else if(study=="ZAH+19") {
    # 20191204 cancer / normal, Zhou et al., 2019
    dat <- read.csv(paste0(datadir, "ZAH+19.csv.xz"), as.is=TRUE)
    description <- "T / N"
    up2 <- dat$Regulation == "Up"
    pcomp <- protcomp(dat$Entry)
  } else if(study=="CHO+18") {
    # 20191206 tumor / adjacent normal, Coleman et al., 2018
    dat <- read.csv(paste0(datadir, "CHO+18.csv.xz"), as.is=TRUE)
    description <- "T / adjacent N"
    up2 <- dat$Highest.mean.condition == "Tumour"
    pcomp <- protcomp(dat$Accession)
  } else stop(paste("pancreatic dataset", dataset, "not available"))
  print(paste0("pdat_pancreatic: ", description, " [", dataset, "]"))
  # use the up2 from the cleaned-up data, if it exists 20190407
  if("up2" %in% colnames(dat)) up2 <- dat$up2
  return(list(dataset=dataset, pcomp=pcomp, up2=up2, description=description))
}

