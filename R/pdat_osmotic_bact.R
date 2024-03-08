# JMDplots/pdat_osmotic_bact.R
# retrieve protein IDs for hyperosmotic experiments
# 20160926 jmd
# 20200418 extract data for bacteria from pdat_osmotic.R
# 20240426 moved from canprot package

pdat_osmotic_bact <- function(dataset = 2020) {
  if(identical(dataset, 2020)) {
    return(c(
             "PNWB09",
             "FTR+10",
             "LPK+13", "QHT+13_Protein.24.h", "QHT+13_Protein.48.h",
             "ADW+14_Protein",
             "KKG+14_Protein_immediate", "KKG+14_Protein_30min", "KKG+14_Protein_80min", "KKG+14_Protein_310min",
             "PBP+14_4.C", "PBP+14_37.C",
             "KLB+15_prot-suc", "KLB+15_prot-NaCl",
             "SKV+16_Glucose_LB", "SKV+16_Osmotic.stress.glucose_LB",
             "KAK+17", "LYS+17",
             "HGC+18", "KSK+18", "LJC+18_wt", "LJC+18_mutant", "TSC18_WT", "TSC18_GsrN",
             "LWS+19", "MGF+19_10", "MGF+19_20",
             "AST+20", "GBR+20_CIRM129", "GBR+20_CIRM1025"
             ))
  }
  # remove tags
  dataset <- strsplit(dataset, "=")[[1]][1]
  # get study and stage/condition
  study <- strsplit(dataset, "_")[[1]][1]
  stage <- paste(strsplit(dataset, "_")[[1]][-1], collapse="_")
  extdatadir <- system.file("extdata", package="JMDplots")
  datadir <- paste0(extdatadir, "/diffexpr/osmotic/")
  if(study=="KKG+12") {
    # 20160918 Escherichia coli, Kocharunchitt et al., 2012
    # KKG+12_25C_aw0.985, KKG+12_14C_aw0.985, KKG+12_25C_aw0.967, KKG+12_14C_aw0.967
    dat <- read.csv(paste0(datadir, "KKG+12.csv.xz"), as.is=TRUE)
    description <- paste("_Escherichia coli_ in NaCl", stage)
    # use specified temperature and subcellular fraction
    icol <- grep(stage, colnames(dat))
    dat <- dat[!is.na(dat[, icol]), ]
    up2 <- dat[, icol] > 0
    dat <- cleanup(dat, "Entry", up2)
    pcomp <- protcomp(dat$Entry, aa_file=paste0(extdatadir, "/aa/bacteria/KKG+12_aa.csv.xz"))
  } else if(study=="KLB+15") {
    # 20160926 Caulobacter crescentus, Kohler et al., 2015
    # KLB+15_trans-suc, KLB+15_trans-NaCl, KLB+15_prot-suc, KLB+15_prot-NaCl
    dat <- read.csv(paste0(datadir, "KLB+15.csv.xz"), as.is=TRUE)
    if(grepl("suc", stage)) osmoticum <- "200 mM sucrose vs M2 minimal salts medium"
    if(grepl("NaCl", stage)) osmoticum <- "40/50 mM NaCl vs M2 minimal salts medium"
    if(grepl("trans", stage)) molecule <- "Gene"
    if(grepl("prot", stage)) molecule <- "Protein"
    description <- paste("_Caulobacter crescentus_", molecule, "in", osmoticum)
    # use protein identified in given experiment
    icol <- grep(gsub("-", ".*", stage), colnames(dat))
    dat <- dat[!is.na(dat[, icol]), ]
    pcomp <- protcomp(dat$Entry, aa_file=paste0(extdatadir, "/aa/bacteria/KLB+15_aa.csv.xz"))
    up2 <- dat[, icol] > 0
  } else if(study=="LJC+18") {
    # 20191102 Listeria monocytogenes membrane vesicles, Lee et al., 2018
    # LJC+18_wt, LJC+18_mutant
    dat <- read.csv(file.path(datadir, "LJC+18.csv.xz"), as.is=TRUE)
    description <- paste("_Listeria monocytogenes_", gsub("wt", "WT", stage), "in 0.5 M NaCl vs control medium")
    # use selected dataset
    icol <- grep(stage, colnames(dat))
    dat <- dat[dat[, icol], ]
    pcomp <- protcomp(dat$Entry, aa_file=file.path(extdatadir, "aa/bacteria/LJC+18_aa.csv.xz"))
    up2 <- dat$condition=="salt"
  } else if(study=="KAK+17") {
    # 20191102 Lactobacillus fermentum NCDC 400 bile salt exposure, Kaur et al., 2017
    dat <- read.csv(paste0(datadir, "KAK+17.csv.xz"), as.is=TRUE)
    description <- "_Lactobacillus fermentum_ with vs without 1.2% w/v bile salts"
    up2 <- dat$Fold.Change > 1
    dat <- cleanup(dat, "Protein.IDs", up2)
    pcomp <- protcomp(dat$Protein.IDs, aa_file=paste0(extdatadir, "/aa/bacteria/KAK+17_aa.csv.xz"))
  } else if(study=="FTR+10") {
    # 20161112 Corynebacterium glutamicum, Fränzel et al., 2010
    dat <- read.csv(file.path(datadir, "FTR+10.csv.xz"), as.is = TRUE)
    description <- "_Corynebacterium glutamicum_ in 750 mM NaCl vs control medium"
    # exclude entries with any NA protein expression data
    dat <- dat[!rowSums(is.na(dat[, 4:6])) > 0, ]
    # use only proteins with consistent expression at 3 time points
    dat <- dat[abs(rowSums(sign(dat[, 4:6]))) == 3, ]
    pcomp <- protcomp(dat$Entry, aa_file = file.path(extdatadir, "aa/bacteria/FTR+10_aa.csv.xz"))
    up2 <- rowSums(sign(dat[, 4:6])) == 3
  } else if(study=="MGF+19") {
    # 20200216 Staphylococcus aureus 10 and 20% NaCl, Ming et al., 2019
    # MGF+19_10, MGF+19_20
    dat <- read.csv(paste0(datadir, "MGF+19.csv.xz"), as.is=TRUE)
    description <- paste0("_Staphylococcus aureus_ in ", stage, "% vs 0% NaCl")
    icol <- grep(stage, colnames(dat))
    # keep proteins with differential expression in the selected experiment
    idiff <- dat[, icol] > 2 | dat[, icol] < 0.5
    idiff[is.na(idiff)] <- FALSE
    dat <- dat[idiff, ]
    up2 <- dat[, icol] > 2
    pcomp <- protcomp(dat$Entry, aa_file = paste0(extdatadir, "/aa/bacteria/MGF+19_aa.csv.xz"))
  } else if(study=="AST+20") {
    # 20200407 Lactobacillus fermentum, Ali et al., 2020
    dat <- read.csv(paste0(datadir, "AST+20.csv.xz"), as.is=TRUE)
    description <- "_Lactobacillus fermentum_ with vs without 0.3% to 1.5% w/v bile salts"
    up2 <- dat$Folds.change > 1
    pcomp <- protcomp(dat$Protein.IDs, aa_file = paste0(extdatadir, "/aa/bacteria/AST+20_aa.csv.xz"))
  } else if(study=="QHT+13") {
    # 20200408 Synechocystis sp. PCC 6803, Qiao et al., 2013
    # QHT+13_Protein.24.h, QHT+13_Protein.48.h
    # QHT+13_Gene.24.h, QHT+13_Gene.48.h, QHT+13_Gene.72.h
    dat <- read.csv(paste0(datadir, "QHT+13.csv.xz"), as.is=TRUE)
    hours <- strsplit(stage, "\\.")[[1]][2]
    molecule <- strsplit(stage, "\\.")[[1]][1]
    description <- paste("_Synechocystis_ sp. PCC 6803", molecule, "in 4% w/v vs 0% added NaCl for", hours, "h")
    icol <- grep(stage, colnames(dat))
    dat <- dat[!is.na(dat[, icol]), ]
    up2 <- dat[, icol] > 1
    dat <- cleanup(dat, "Entry", up2)
    pcomp <- protcomp(dat$Entry, aa_file = paste0(extdatadir, "/aa/bacteria/QHT+13_aa.csv.xz"))
  } else if(study=="LPK+13") {
    # 20200411 Lactobacillus johnsonii, Lee et al., 2013
    dat <- read.csv(paste0(datadir, "LPK+13.csv.xz"), as.is=TRUE)
    description <- "_Lactobacillus johnsonii_ with vs without 0.1-0.3% bile salt"
    up2 <- dat$Regulation == "up"
    dat <- cleanup(dat, "UniProt", up2)
    pcomp <- protcomp(dat$UniProt, aa_file = paste0(extdatadir, "/aa/bacteria/LPK+13_aa.csv.xz"))
  } else if(study=="LYS+17") {
    # 20200412 Lactobacillus salivarius LI01, Lv et al., 2017
    dat <- read.csv(paste0(datadir, "LYS+17.csv.xz"), as.is=TRUE)
    description <- "_Lactobacillus salivarius_ LI01 with vs without 0.15% ox bile"
    up2 <- dat$log2fold_change > 0
    dat <- cleanup(dat, "Entry", up2)
    pcomp <- protcomp(dat$Entry, aa_file = paste0(extdatadir, "/aa/bacteria/LYS+17_aa.csv.xz"))
  } else if(study=="LWS+19") {
    # 20200412 Lactobacillus plantarum FS5-5, Li et al., 2019
    dat <- read.csv(paste0(datadir, "LWS+19.csv.xz"), as.is=TRUE)
    description <- "_Lactobacillus plantarum_ FS5-5 in 6-8% w/v vs 0% NaCl"
    up2 <- dat$X118.113 > 1
    dat <- cleanup(dat, "No.", up2)
    pcomp <- protcomp(dat$No., aa_file = paste0(extdatadir, "/aa/bacteria/LWS+19_aa.csv.xz"))
  } else if(study=="KSK+18") {
    # 20200413 Acidihalobacter prosperus DSM 14174, Khaleque et al., 2018
    dat <- read.csv(file.path(datadir, "KSK+18.csv.xz"), as.is=TRUE)
    description <- "_Acidihalobacter prosperus_ DSM 14174 30 g/L / 5 g/L NaCl"
    dat <- check_IDs(dat, "Protein", aa_file = file.path(extdatadir, "aa/bacteria/KSK+18_aa.csv.xz"))
    up2 <- dat$FCProteins > 1
    pcomp <- protcomp(dat$Protein, aa_file = file.path(extdatadir, "aa/bacteria/KSK+18_aa.csv.xz"))
  } else if(study=="PNWB09") {
    # 20200416 Synechocystis sp. PCC6803, Pandhal et al., 2009
    dat <- read.csv(file.path(datadir, "PNWB09.csv.xz"), as.is=TRUE)
    description <- "_Synechocystis_ sp. PCC6803 in 6% w/v NaCl vs no added salt"
    up2 <- dat$Expression == "increased"
    pcomp <- protcomp(dat$Entry, aa_file = file.path(extdatadir, "aa/bacteria/PNWB09_aa.csv.xz"))
  } else if(study=="GBR+20") {
    # 20200416 Propionibacterium freudenreichii, Gaucher et al., 2020
    # GBR+20_CIRM129, GBR+20_CIRM1025
    dat <- read.csv(file.path(datadir, "GBR+20.csv.xz"), as.is=TRUE)
    description <- paste("_Propionibacterium freudenreichii_", stage, "in NaCl vs MMO")
    icol <- match(stage, colnames(dat))
    dat <- dat[!is.na(dat[, icol]), ]
    up2 <- dat[, icol] > 1
    pcomp <- protcomp(dat$Entry, aa_file = file.path(extdatadir, "aa/bacteria/GBR+20_aa.csv.xz"))
  } else if(study=="HGC+18") {
    # 20200418 Lactobacillus casei, Huang et al., 2018
    dat <- read.csv(paste0(datadir, "HGC+18.csv.xz"), as.is=TRUE)
    description <- "_Lactobacillus casei_ BL23 in hyper-concentrated vs isotonic sweet whey"
    up2 <- dat$Regulation == "up"
    pcomp <- protcomp(dat$Entry, aa_file = paste0(extdatadir, "/aa/bacteria/HGC+18_aa.csv.xz"))
  } else if(study=="SKV+16") {
    # 20200420 E. coli, Schmidt et al., 2016
    # SKV+16_Glucose_LB, SKV+16_Osmotic.stress.glucose_LB
    dat <- read.csv(file.path(datadir, "SKV+16.csv.xz"), as.is=TRUE)
    description <- paste("_Escherichia coli_ in", gsub("_", " vs ", stage))
    icol <- grep(stage, colnames(dat))
    dat <- dat[!is.na(dat[, icol]), ]
    dat <- dat[dat[, icol] > 2 | dat[, icol] < 0.5, ]
    up2 <- dat[, icol] > 2
    dat <- check_IDs(dat, "Uniprot.Accession", aa_file = file.path(extdatadir, "aa/bacteria/SKV+16_aa.csv.xz"))
    dat <- cleanup(dat, "Uniprot.Accession", up2)
    pcomp <- protcomp(dat$Uniprot.Accession, aa_file = file.path(extdatadir, "aa/bacteria/SKV+16_aa.csv.xz"))
  } else if(study=="KKG+14") {
    # 20200420 E. coli, Kocharunchitt et al., 2014
    # KKG+14_Gene_immediate, KKG+14_Gene_30min, KKG+14_Gene_80min, KKG+14_Gene_310min
    # KKG+14_Protein_immediate, KKG+14_Protein_30min, KKG+14_Protein_80min, KKG+14_Protein_310min
    dat <- read.csv(file.path(datadir, "KKG+14.csv.xz"), as.is=TRUE)
    molecule <- strsplit(stage, "_")[[1]][[1]]
    time <- strsplit(stage, "_")[[1]][[2]]
    description <- paste("_Escherichia coli_", molecule, "in NaCl (0.967 aw) vs control for", gsub("min", " min", time))
    icol <- grep(stage, colnames(dat))
    dat <- dat[!is.na(dat[, icol]), ]
    up2 <- dat[, icol] > 0
    dat <- cleanup(dat, "Entry", up2)
    pcomp <- protcomp(dat$Entry, aa_file = file.path(extdatadir, "aa/bacteria/KKG+14_aa.csv.xz"))
  } else if(study=="ADW+14") {
    # 20200420 Bifidobacterium longum, An et al., 2014
    # ADW+14_Gene, ADW+14_Protein
    dat <- read.csv(file.path(datadir, "ADW+14.csv.xz"), as.is=TRUE)
    description <- paste("_Bifidobacterium longum_ BBMN68", stage, "with vs without 0.75 g/l ox bile")
    icol <- grep(paste0(stage, "_"), colnames(dat))
    dat <- dat[!is.na(dat[, icol]), ]
    up2 <- dat[, icol] > 0
    pcomp <- protcomp(dat$Entry, aa_file = file.path(extdatadir, "aa/bacteria/ADW+14_aa.csv.xz"))
  } else if(study=="TSC18") {
    # 20200423 Caulobacter crescentus, Tien et al., 2018
    # TSC18_GsrN, TSC18_WT
    dat <- read.csv(file.path(datadir, "TSC18.csv.xz"), as.is=TRUE)
    description <- paste("_Caulobacter crescentus_", stage, "in 300 mM sucrose vs control")
    icol <- grep(stage, colnames(dat))
    dat <- dat[!is.na(dat[, icol]), ]
    dat <- dat[dat[, icol] > 1.5 | dat[, icol] < 2/3, ]
    # exclude 0 or infinite ratios (i.e. ones that involve unquantified proteins)
    dat <- dat[dat[, icol]!=0 & !is.infinite(dat[, icol]), ]
    up2 <- dat[, icol] > 1.5
    pcomp <- protcomp(dat$Majority.protein.IDs, aa_file = file.path(extdatadir, "aa/bacteria/TSC18_aa.csv.xz"))
  } else if(study=="PBP+14") {
    # 20200426 Listeria monocytogenes, Pittman et al., 2014
    # PBP+14_4.C, PBP+14_37.C
    dat <- read.csv(file.path(datadir, "PBP+14.csv.xz"), as.is=TRUE)
    description <- paste("_Listeria monocytogenes_ in 3% NaCl vs control at", stage)
    icol <- grep(stage, colnames(dat))
    dat <- dat[!is.na(dat[, icol]), ]
    up2 <- dat[, icol] > 1
    pcomp <- protcomp(dat$Entry, aa_file = file.path(extdatadir, "aa/bacteria/PBP+14_aa.csv.xz"))
  } else stop(paste("osmotic_bact dataset", dataset, "not available"))
  print(paste0("pdat_osmotic_bact: ", description, " [", dataset, "]"))
  # use the up2 from the cleaned-up data, if it exists 20191120
  if("up2" %in% colnames(dat)) up2 <- dat$up2
  return(list(dataset=dataset, pcomp=pcomp, up2=up2, description=description))
}
