# JMDplots/pdat_saltygenes.R
# retrieve protein IDs for differentially expressed genes in osmotic stress
# 20200418 jmd

pdat_saltygenes <- function(dataset=2020, basis="rQEC") {
  if(identical(dataset, 2020)) {
    return(c(
             "KSA+02_NaCl", "KSA+02_sorbitol", "WJ02",
             "HZP+05_HSS", "HZP+05_HOS", "LGW+05", "SLA+05",
             "GCP08_30", "GCP08_43",
             "SBB+09_NaCl", "SBB+09_Sucrose",
             "HMO+10_transcriptomics",
             # don't include yeast 20200424
             #"LTH+11_RNA_30", "LTH+11_RNA_60", "LTH+11_RNA_90", "LTH+11_RNA_120", "LTH+11_RNA_240",
             "BBWB12_37_2.5", "BBWB12_37_5", "BBWB12_37_10", "BBWB12_37_20", #"BBWB12_37_exp", 29 up, 6 down
             "LB12_NaCl",
             "QHT+13_Gene.24.h", "QHT+13_Gene.48.h", "QHT+13_Gene.72.h",
             "WGB+13_N", "WGB+13_U",
             "ADW+14_Gene", "KKG+14_Gene_30min", "KKG+14_Gene_80min", "KKG+14_Gene_310min",
             "KSM+14_NaCl", "KSM+14_GB",
             "MGM+14_3.5", "MGM+14_4.5", "MGM+14_5", "MGM+14_5.5",
             "SLM+14_5", "SLM+14_30", "SLM+14_60",
             "FRH+15_NaCl_1h", "FRH+15_NaCl_6h", "FRH+15_NaCl_24h",
             "FRH+15_KCl_1h", "FRH+15_KCl_6h", "FRH+15_KCl_24h",
             "FRH+15_glycerol_1h", "FRH+15_glycerol_6h", "FRH+15_glycerol_24h",
             "KLB+15_trans-NaCl", "KLB+15_trans-suc",
             "HLL17_45min", "HLL17_14h",
             "MWZ+18"
             ))
  }
  # remove tags
  dataset <- strsplit(dataset, "=")[[1]][1]
  # get study and stage/condition
  study <- strsplit(dataset, "_")[[1]][1]
  stage <- paste(strsplit(dataset, "_")[[1]][-1], collapse="_")
  extdatadir <- system.file("extdata", package="JMDplots")
  datadir <- paste0(extdatadir, "/expression/osmotic/")
  if(study %in% c("KLB+15", "KKG+14", "ADW+14", "QHT+13")) {
    return(pdat_osmotic_bact(dataset, basis))
  #} else if(study =="LTH+11") {
  #  return(pdat_osmotic_euk(dataset, basis))
  } else if(study=="HMO+10") {
    # 20191102 Bacillus subtilis, Hahne et al., 2010
    # HMO+10_prot-cytosol, HMO+10_prot-membrane, HMO+10_transcriptomics
    dat <- read.csv(paste0(datadir, "HMO+10.csv.xz"), as.is=TRUE)
    #description <- paste("Bacillus subtilis", stage)
    description <- "Bacillus subtilis in 6% w/v NaCl vs control"
    # use selected dataset
    dat <- dat[dat$Experiment==stage, ]
    # count up- or down-expression at each time point
    diff <- apply(sign(dat[, 4:7]), 1, sum)
    dat <- cbind(dat, diff = diff)
    # keep proteins that have same expression in at least 3 out of 4 time points
    dat <- dat[abs(dat$diff) >= 2, ]
    up2 <- dat$diff > 0
    # update IDs / remove duplicated IDs
    aa_file <- file.path(extdatadir, "aa/bacteria/HMO+10_aa.csv.xz")
    dat <- check_IDs(dat, "UniProtKB", aa_file)
    dat <- cleanup(dat, "UniProtKB", up2)
    pcomp <- protcomp(dat$UniProtKB, basis = basis, aa_file = aa_file)
  } else if(study=="WJ02") {
    # 20200419 E. coli, Weber and Jung, 2002
    dat <- read.csv(paste0(datadir, "WJ02.csv.xz"), as.is=TRUE)
    description <- "E. coli in 0.4 M NaCl vs control"
    up2 <- dat$Regulation == "up"
    pcomp <- protcomp(dat$Entry, basis = basis, aa_file = file.path(extdatadir, "aa/bacteria/WJ02_aa.csv.xz"))
  } else if(study=="KSA+02") {
    # 20200419 Synechocystis, Kanesaki et al., 2002
    # KSA+02_NaCl, KSA+02_sorbitol
    dat <- read.csv(paste0(datadir, "KSA+02.csv.xz"), as.is=TRUE)
    description <- paste("Synechocystis sp. PCC 6803 in 0.5 M", stage, "vs control")
    # use highly differential proteins for the indicated experiment
    icol <- grep(stage, colnames(dat))
    dat <- dat[abs(dat[, icol]) >= 3, ]
    up2 <- dat[, icol] >= 3
    pcomp <- protcomp(dat$Entry, basis = basis, aa_file = file.path(extdatadir, "aa/bacteria/KSA+02_aa.csv.xz"))
  } else if(study=="GCP08") {
    # 20200419 E. coli, Gunasekera et al., 2008
    # GCP08_30, GCP08_43
    dat <- read.csv(paste0(datadir, "GCP08.csv.xz"), as.is=TRUE)
    description <- paste("E. coli in 0.5 M vs 0 M NaCl at", stage, "deg C")
    if(stage == "30") icol <- grep("C2.C1", colnames(dat))
    if(stage == "43") icol <- grep("C5.C4", colnames(dat))
    dat <- dat[abs(dat[, icol]) > 1, ]
    up2 <- dat[, icol] > 1
    dat <- cleanup(dat, "Entry", up2)
    pcomp <- protcomp(dat$Entry, basis = basis, aa_file = file.path(extdatadir, "aa/bacteria/GCP08_aa.csv.xz"))
  } else if(study=="LGW+05") {
    # 20200419 Shewanella oneidensis, Liu et al., 2005
    dat <- read.csv(paste0(datadir, "LGW+05.csv.xz"), as.is=TRUE)
    description <- "Shewanella oneidensis MR-1 in 0.5 vs 0.1 M NaCl"
    up2 <- dat$AVG > 1
    dat <- cleanup(dat, "Entry", up2)
    pcomp <- protcomp(dat$Entry, basis = basis, aa_file = file.path(extdatadir, "aa/bacteria/LGW+05_aa.csv.xz"))
  } else if(study=="HLL17") {
    # 20200419 Methylocystis, Han et al., 2017
    # HLL17_45min, HLL17_14h
    dat <- read.csv(paste0(datadir, "HLL17.csv.xz"), as.is=TRUE)
    description <- paste("Methylocystis sp. strain SC2 in 0.75 % NaCl vs control at", stage)
    icol <- grep(stage, colnames(dat))
    dat <- dat[!is.na(dat[, icol]), ]
    up2 <- dat[, icol] > 0
    pcomp <- protcomp(dat$Entry, basis = basis, aa_file = file.path(extdatadir, "aa/bacteria/HLL17_aa.csv.xz"))
  } else if(study=="BBWB12") {
    # 20200419 Listeria monocytogenes, Bergholz et al., 2012
    # BBWB12_37_2.5, 37_5, 37_10, 37_20, 37_exp, 7_2.5, 7_5, 7_10, 7_20, 7_exp
    dat <- read.csv(paste0(datadir, "BBWB12.csv.xz"), as.is=TRUE)
    description <- paste("L. monocytogenes strain H7858 in 6 % NaCl vs control at", stage)
    icol <- grep(paste0("X", stage), colnames(dat))
    dat <- dat[abs(dat[, icol]) > 1, ]
    up2 <- dat[, icol] > 1
    dat <- cleanup(dat, "Entry", up2)
    pcomp <- protcomp(dat$Entry, basis = basis, aa_file = file.path(extdatadir, "aa/bacteria/BBWB12_aa.csv.xz"))
  } else if(study=="KSM+14") {
    # 20200419 Bacillus subtilis, Kohlstedt et al., 2014
    # KSM+14_NaCl, KSM+14_GB
    dat <- read.csv(paste0(datadir, "KSM+14.csv.xz"), as.is=TRUE)
    if(stage=="NaCl") description <- "Bacillus subtilis in 1.2 M NaCl vs control"
    if(stage=="GB") description <- "Bacillus subtilis with vs without glycine betaine in 1.2 M NaCl"
    icol <- grep(stage, colnames(dat))
    dat <- dat[!is.na(dat[, icol]), ]
    up2 <- dat[, icol] > 0
    dat <- cleanup(dat, "Entry", up2)
    pcomp <- protcomp(dat$Entry, basis = basis, aa_file = file.path(extdatadir, "aa/bacteria/KSM+14_aa.csv.xz"))
  } else if(study=="SLM+14") {
    # 20200419 Enterococcus faecalis, Solheim et al., 2014
    # SLM+14_5, SLM+14_30, SLM+14_60
    dat <- read.csv(paste0(datadir, "SLM+14.csv.xz"), as.is=TRUE)
    description <- paste("Enterococcus faecalis in 6.5 % NaCl vs control at", stage, "min")
    icol <- grep(stage, colnames(dat))
    dat <- dat[!is.na(dat[, icol]), ]
    up2 <- dat[, icol] > 0
    dat <- cleanup(dat, "Entry", up2)
    pcomp <- protcomp(dat$Entry, basis = basis, aa_file = file.path(extdatadir, "aa/bacteria/SLM+14_aa.csv.xz"))
  } else if(study=="MWZ+18") {
    # 20200419 Lactobacillus paracasei, Ma et al., 2018
    dat <- read.csv(paste0(datadir, "MWZ+18.csv.xz"), as.is=TRUE)
    description <- "Lactobacillus paracasei L9 with vs without 0.13% ox bile"
    up2 <- dat$Log2Fold_change > 0
    dat <- cleanup(dat, "Entry", up2)
    pcomp <- protcomp(dat$Entry, basis = basis, aa_file = file.path(extdatadir, "aa/bacteria/MWZ+18_aa.csv.xz"))
  } else if(study=="HZP+05") {
    # 20200419 Yersinia pestis, Han et al., 2005
    # HZP+05_HOS, HZP+05_HSS
    dat <- read.csv(paste0(datadir, "HZP+05.csv.xz"), as.is=TRUE)
    if(stage=="HOS") description <- paste("Yersinia pestis in 0.5 M sorbitol vs control")
    if(stage=="HSS") description <- paste("Yersinia pestis in 0.5 M NaCl vs control")
    icol <- grep(stage, colnames(dat))
    dat <- dat[!is.na(dat[, icol]), ]
    dat <- dat[dat[, icol] >= 2 | dat[, icol] <= 0.5, ]
    up2 <- dat[, icol] >= 2
    dat <- cleanup(dat, "Entry", up2)
    pcomp <- protcomp(dat$Entry, basis = basis, aa_file = file.path(extdatadir, "aa/bacteria/HZP+05_aa.csv.xz"))
  } else if(study=="SLA+05") {
    # 20200420 Synechocystis, Shapiguzov et al., 2005
    dat <- read.csv(paste0(datadir, "SLA+05.csv.xz"), as.is=TRUE)
    description <- "Synechocystis sp. PCC 6803 in 0.5 M sorbitol vs control"
    up2 <- dat$IF > 0
    pcomp <- protcomp(dat$Entry, basis = basis, aa_file = file.path(extdatadir, "aa/bacteria/SLA+05_aa.csv.xz"))
  } else if(study=="LB12") {
    # 20200420 Synechococcus, Ludwig and Bryant, 2012
    # LB12_NaCl, LB12_glycerol
    dat <- read.csv(paste0(datadir, "LB12.csv.xz"), as.is=TRUE)
    if(stage=="NaCl") description <- "Synechococcus sp. strain PCC 7002 in 1.5 M NaCl vs control"
    if(stage=="glycerol") description <- "Synechococcus sp. strain PCC 7002 in 10 mM glycerol vs control"
    icol <- grep(stage, colnames(dat))
    dat <- dat[!is.na(dat[, icol]), ]
    dat <- dat[dat[, icol] > 2 | dat[, icol] < 0.5, ]
    up2 <- dat[, icol] > 2
    dat <- cleanup(dat, "Entry", up2)
    pcomp <- protcomp(dat$Entry, basis = basis, aa_file = file.path(extdatadir, "aa/bacteria/LB12_aa.csv.xz"))
  } else if(study=="FRH+15") {
    # 20200420 Salmonella enterica, Finn et al., 2015
    # FRH+15_NaCl_1h, FRH+15_NaCl_6h, FRH+15_NaCl_24h
    # FRH+15_KCl_1h, FRH+15_KCl_6h, FRH+15_KCl_24h
    # FRH+15_glycerol_1h, FRH+15_glycerol_6h, FRH+15_glycerol_24h
    dat <- read.csv(paste0(datadir, "FRH+15.csv.xz"), as.is=TRUE)
    description <- paste("Salmonella enterica", stage)
    icol <- grep(stage, colnames(dat))
    dat <- dat[!is.na(dat[, icol]), ]
    up2 <- dat[, icol] > 1
    dat <- cleanup(dat, "Entry", up2)
    pcomp <- protcomp(dat$Entry, basis = basis, aa_file = file.path(extdatadir, "aa/bacteria/FRH+15_aa.csv.xz"))
  } else if(study=="MGM+14") {
    # 20200423 E. coli, Metris et al., 2014
    # MGM+14_3.5, MGM+14_4.5, MGM+14_5, MGM+14_5.5 
    dat <- read.csv(file.path(datadir, "MGM+14.csv.xz"), as.is=TRUE)
    description <- paste("E. coli in", stage, "vs 2 % NaCl (with glycine betaine)")
    icol <- grep(paste0("_", stage, "_"), colnames(dat))
    dat <- dat[abs(dat[, icol]) > 1, ]
    up2 <- dat[, icol] > 1
    dat <- cleanup(dat, "Entry", up2)
    pcomp <- protcomp(dat$Entry, basis, aa_file = file.path(extdatadir, "aa/bacteria/MGM+14_aa.csv.xz"))
  } else if(study=="SBB+09") {
    # 20200423 E. coli, Shabala et al., 2009
    # SBB+09_NaCl, SBB+09_Sucrose
    dat <- read.csv(file.path(datadir, "SBB+09.csv.xz"), as.is=TRUE)
    description <- paste("Escherichia coli in ~2.7 Os/kg", gsub("Sucrose", "sucrose", stage), "vs control")
    icol <- grep(stage, colnames(dat))
    dat <- dat[!is.na(dat[, icol]), ]
    up2 <- dat[, icol] == "up"
    dat <- cleanup(dat, "Entry", up2)
    pcomp <- protcomp(dat$Entry, basis, aa_file = file.path(extdatadir, "aa/bacteria/SBB+09_aa.csv.xz"))
  } else if(study=="WGB+13") {
    # 20200423 E. coli in NaCl or urea, Withman et al., 2013
    # WGB+13_N, WGB+13_U
    dat <- read.csv(file.path(datadir, "WGB+13.csv.xz"), as.is=TRUE)
    if(stage=="N") description <- "Escherichia coli in 0.3 M NaCl vs control"
    if(stage=="U") description <- "Escherichia coli in 0.6 M urea vs control"
    icol <- grep(paste0(stage, ".vs.K"), colnames(dat))
    dat <- dat[!is.na(dat[, icol]), ]
    up2 <- dat[, icol] > 0
    dat <- cleanup(dat, "Entry", up2)
    pcomp <- protcomp(dat$Entry, basis, aa_file = file.path(extdatadir, "aa/bacteria/WGB+13_aa.csv.xz"))
  } else stop(paste("saltygenes dataset", dataset, "not available"))
  print(paste0("pdat_saltygenes: ", description, " [", dataset, "]"))
  # use the up2 from the cleaned-up data, if it exists 20190407
  if("up2" %in% colnames(dat)) up2 <- dat$up2
  return(list(dataset=dataset, basis=basis, pcomp=pcomp, up2=up2, description=description))
}
