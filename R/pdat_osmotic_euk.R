# JMDplots/pdat_osmotic_euk.R
# retrieve protein IDs for osmotic stress (NaCl only) experiments in eukaryotes
# 20160926 jmd
# 20200418 eukaryotes extracted from pdat_osmotic.R

pdat_osmotic_euk <- function(dataset = 2020) {
  if(identical(dataset, 2020)) {
    return(c(
             "DAA+05",
             "MHN+08_Hyper",
             "LTH+11_Protein_30=yeast", "LTH+11_Protein_60=yeast", "LTH+11_Protein_90=yeast", "LTH+11_Protein_120=yeast", "LTH+11_Protein_240=yeast",
             "OBBH11",
             "LFY+12_C1h", "LFY+12_C8h", "LFY+12_C2p", "LFY+12_N1h", "LFY+12_N8h", "LFY+12_N2p",
             "CLG+15", "SCG+15_nodelay", "SCG+15_delayed", "YDZ+15=yeast",
             "GAM+16_HTS", "GAM+16_HTS.Cmx", "RBP+16=yeast",
             "JBG+18=yeast"
             ))
  }
  # remove tags
  dataset <- strsplit(dataset, "=")[[1]][1]
  # get study and stage/condition
  study <- strsplit(dataset, "_")[[1]][1]
  stage <- paste(strsplit(dataset, "_")[[1]][-1], collapse="_")
  extdatadir <- system.file("extdata", package="JMDplots")
  datadir <- paste0(extdatadir, "/canH2O/osmotic/")
  if(study=="CLG+15") {
    # 20160925 conjunctival epithelial cells, Chen et al., 2015
    dat <- read.csv(paste0(datadir, "CLG+15.csv.xz"), as.is=TRUE)
    description <- paste("human conjunctival epithelial cells in 380 or 480 mOsm vs 280 mOsm NaCl")
    # use proteins that have same direction of change in both conditions
    dat <- dat[(dat$T1 > 1 & dat$T2 > 1) | (dat$T1 < 1 & dat$T2 < 1), ]
    pcomp <- protcomp(dat$accession..UniProtKB.Swiss.Prot.)
    up2 <- dat$T1 > 1
  } else if(study=="OBBH11") {
    # 20160925 adipose-derived stem cells, Oswald et al., 2011
    dat <- read.csv(paste0(datadir, "OBBH11.csv.xz"), as.is=TRUE)
    description <- "adipose-derived stem cells in 400 mOsm vs 300 mOsm NaCl"
    dat <- check_IDs(dat, "Uniprot.Protein.Code")
    pcomp <- protcomp(dat$Uniprot.Protein.Code)
    up2 <- dat$Elucidator.Expression.Ratio..Treated.Control. > 1
  } else if(study=="YDZ+15") {
    # 20160926 Yarrowia lipolytica, Yang et al., 2015
    dat <- read.csv(paste0(datadir, "YDZ+15.csv.xz"), as.is=TRUE)
    description <- paste("_Yarrowia lipolytica_ in 4.21 osmol/kg vs 3.17 osmol/kg NaCl")
    up2 <- dat$Av..ratio..high.low. > 0
    dat <- cleanup(dat, "Accession.No.", up2)
    pcomp <- protcomp(substr(dat$Accession.No., 4, 12), aa_file=paste0(extdatadir, "/aa/yeast/YDZ+15_aa.csv.xz"))
  } else if(study=="RBP+16") {
    # 20161112 Paracoccidioides lutzii, da Silva Rodrigues et al., 2016
    dat <- read.csv(paste0(datadir, "RBP+16.csv.xz"), as.is=TRUE)
    description <- "_Paracoccidioides lutzii_ in 0.1 M KCl vs medium with no added KCl"
    up2 <- dat$Fold.change > 1
    dat <- cleanup(dat, "Entry", up2)
    pcomp <- protcomp(dat$Entry, aa_file=paste0(extdatadir, "/aa/yeast/RBP+16_aa.csv.xz"))
  } else if(study=="JBG+18") {
    # 20200406 Candida albicans 1 M NaCl, Jacobsen et al., 2018
    dat <- read.csv(paste0(datadir, "JBG+18.csv.xz"), as.is=TRUE)
    description <- "_Candida albicans_ in 1 M NaCl vs medium with no added NaCl"
    up2 <- dat$log2.ratio > 0
    dat <- cleanup(dat, "Entry", up2)
    pcomp <- protcomp(dat$Entry, aa_file = paste0(extdatadir, "/aa/yeast/JBG+18_aa.csv.xz"))
  } else if(study=="LFY+12") {
    # 20200417 HEK293 cells, Li et al., 2012
    # LFY+12_C1h, LFY+12_C8h, LFY+12_C2p, LFY+12_N1h, LFY+12_N8h, LFY+12_N2p
    dat <- read.csv(file.path(datadir, "LFY+12.csv.xz"), as.is=TRUE)
    compartment <- ifelse(grepl("C", stage), "cytoplasm", "nucleus")
    time <- substr(stage, 2, 2)
    units <- ifelse(grepl("h", stage), "h", "passages")
    description <- paste(compartment, "of HEK293 cells in 500 (NaCl added) vs 300 mosmol/kg medium for", time, units)
    icol <- grep(stage, colnames(dat))
    dat <- dat[!is.na(dat[, icol]), ]
    dat <- check_IDs(dat, "UniProt.ID")
    up2 <- dat[, icol] > 0
    pcomp <- protcomp(dat$UniProt.ID)
  } else if(study=="GAM+16") {
    # 20200418 human small airway epithelial cells, Gamboni et al., 2016
    # GAM+16_HTS, GAM+16_HTS.Cmx
    dat <- read.csv(file.path(datadir, "GAM+16.csv.xz"), as.is=TRUE)
    control <- paste0("isotonic", substr(stage, 4, 7))
    description <- paste("human small airway epithelial cells in", stage, "vs", control)
    icol <- grep(paste0(stage, ".Iso"), colnames(dat))
    dat <- dat[!is.na(dat[, icol]), ]
    dat <- dat[dat[, icol] > 2 | dat[, icol] < 0.5, ]
    up2 <- dat[, icol] > 2
    pcomp <- protcomp(dat$Uniprot.ID)
  } else if(study=="DAA+05") {
    # 20200418 thick ascending limb of Henle's loop cells, Dihazi et a., 2005
    dat <- read.csv(file.path(datadir, "DAA+05.csv.xz"), as.is=TRUE)
    description <- "thick ascending limb of Henle's loop cells in 600 (NaCl added) vs 300 mosmol/kg medium"
    up2 <- dat$Regulation == "up"
    dat <- cleanup(dat, "Entry", up2)
    pcomp <- protcomp(dat$Entry, aa_file = paste0(extdatadir, "/aa/mouse/DAA+05_aa.csv.xz"))
  } else if(study=="LTH+11") {
    # 20200419 S. cerevisiae, Lee et al., 2011
    # LTH+11_RNA_30, 60, 90, 120, 240
    # LTH+11_Protein_30, 60, 90, 120, 240
    dat <- read.csv(file.path(datadir, "LTH+11.csv.xz"), as.is=TRUE)
    molecule <- strsplit(stage, "_")[[1]][1]
    time <- strsplit(stage, "_")[[1]][2]
    description <- paste("_Saccharomyces cerevisiae_", molecule, "in 0.7 M NaCl vs control medium for", time, "min")
    # remove rows that have high q-value
    iq <- grep(paste0(molecule, ".Change"), colnames(dat))
    dat[is.na(dat[, iq]), iq] <- 1
    dat <- dat[dat[, iq] < 0.05, ]
    # keep rows that have high fold change
    icol <- grep(stage, colnames(dat))
    if(molecule=="protein") dat <- dat[abs(dat[, icol]) > 0.2, ]
    if(molecule=="RNA") dat <- dat[abs(dat[, icol]) > 0.5, ]
    up2 <- dat[, icol] > 0
    dat <- cleanup(dat, "Entry", up2)
    pcomp <- protcomp(dat$Entry, aa_file = paste0(extdatadir, "/aa/yeast/LTH+11_aa.csv.xz"))
  } else if(study=="SCG+15") {
    # 20200419 S. cerevisiae, Selevsek et al., 2015
    # SCG+15_nodelay, SCG+15_delayed
    dat <- read.csv(file.path(datadir, "SCG+15.csv.xz"), as.is=TRUE)
    description <- paste("_Saccharomyces cervisiae_ in 0.4 M NaCl vs control -", stage)
    if(stage=="nodelay") {
      dat <- dat[dat$Cluster %in% c(1, 4), ]
      up2 <- dat$Cluster == 1
    }
    if(stage=="delayed") {
      dat <- dat[dat$Cluster %in% c(2, 3), ]
      up2 <- dat$Cluster == 2
    }
    pcomp <- protcomp(dat$Entry, aa_file = paste0(extdatadir, "/aa/yeast/SCG+15_aa.csv.xz"))
  } else if(study=="MHN+08") {
    # 20200422 mouse, Mao et al., 2008
    # MHN+08_Hyper, MHN+08_Hypo
    dat <- read.csv(file.path(datadir, "MHN+08.csv.xz"), as.is=TRUE)
    if(stage=="Hyper") description <- "mouse CGR8 embryonic stem cells in 500 (NaCl added) vs 340 mOsM medium"
    if(stage=="Hypo") description <- "mouse CGR8 embryonic stem cells in 200 (lower NaCl) vs 340 mOsM medium"
    icol <- grep(stage, colnames(dat))
    dat <- dat[!is.na(dat[, icol]), ]
    dat <- check_IDs(dat, "Entry", aa_file = paste0(extdatadir, "/aa/mouse/MHN+08_aa.csv.xz"))
    up2 <- dat[, icol] > 1
    dat <- cleanup(dat, "Entry", up2)
    pcomp <- protcomp(dat$Entry, aa_file = paste0(extdatadir, "/aa/mouse/MHN+08_aa.csv.xz"))
  } else stop(paste("osmotic_euk dataset", dataset, "not available"))
  print(paste0("pdat_osmotic_euk: ", description, " [", dataset, "]"))
  # use the up2 from the cleaned-up data, if it exists 20191120
  if("up2" %in% colnames(dat)) up2 <- dat$up2
  return(list(dataset=dataset, pcomp=pcomp, up2=up2, description=description))
}
