# canprot/R/pdat_hypoxia.R
# retrieve protein IDs for hypoxia data sets
# 20160414 jmd first version
# 20190322-20200118 assemble data for 2020 compilation
# - exclude datasets with < 30 either up- or down-regulated proteins
# - exclude ReOx datasets
# - move secretome data to pdat_secreted.R
# - move 3D data to pdat_3D.R
# - add new datasets
# 20240426 moved from canprot package

pdat_hypoxia <- function(dataset = 2020) {
  if(identical(dataset, 2020)) {
    return(c(
             "SBB+06=cancer",
             "FWH+13", "RHD+13_Hx48=cancer", "RHD+13_Hx72=cancer", "VTMF13",
             "DCH+14=cancer", "DYL+14_Hx48-S=cancer", "DYL+14_Hx72-S=cancer", "DYL+14_Hx48-P=cancer", "DYL+14_Hx72-P=cancer",
             "BSA+15=cancer",
             "HWA+16=cancer",
             "LCS16_translation=cancer",
             "CGH+17_whole", "ZXS+17=cancer",
             "CLY+18_proteome", "GBH+18=cancer", "LKK+18", "WTG+18",
             "CSK+19=cancer", "GPT+19_Light.S=cancer", "GPT+19_Heavy.S=cancer", "KAN+19_proteome=cancer", "LLL+19",
             "BCMS20=cancer", "RVN+20_DMSO=cancer", "RVN+20_NO.sul", "RVN+20_sul", "RVN+20_DMSO.4Gy", "RVN+20_NO.sul.4Gy", "RVN+20_sul.4Gy",
             "SPJ+20_POS=cancer", "SPJ+20_HMPOS=cancer"
             ))
  }
  if(identical(dataset, 2017)) {
    return(c("HXS+06",
             "BRA+10", "DPL+10",
             "BMJ+11", "CBW+11",
             "LAR+12", "MHG+12_P5=SPH", "MHG+12_P2=SPH",
             "MVC+12_perinecrotic=SPH", "MVC+12_necrotic=SPH",
             "FWH+13", "RHD+13_Hx48", "RHD+13_Hx72", "RHD+13_ReOx=ReOx", "VTMF13",
             "DYL+14_Hx48-S", "DYL+14_Hx72-S", "DYL+14_ReOx-S=ReOx",
             "DYL+14_Hx48-P", "DYL+14_Hx72-P", "DYL+14_ReOx-P=ReOx",
             "RKP+14=SPH", "WRK+14=SPH",
             "BSA+15",
             "HWA+16",
             "LCS16_transcription", "LCS16_translation",
             "RSE+16=ASC", "XCJ+16_CoCl2", "XCJ+16_SAL=ReOx", "YLW+16=SPH"))
  }
  # remove tags
  dataset <- strsplit(dataset, "=")[[1]][1]
  # get study and stage/condition
  study <- strsplit(dataset, "_")[[1]][1]
  stage <- paste(strsplit(dataset, "_")[[1]][-1], collapse="_")
  extdatadir <- system.file("extdata", package="JMDplots")
  datadir <- paste0(extdatadir, "/diffexpr/hypoxia/")
  if(study=="BSA+15") {
    # 20160412 HeLa hypoxia, Bousquet et al., 2005
    dat <- read.csv(paste0(datadir, "BSA+15.csv.xz"), as.is=TRUE)
    description <- "HeLa cervical cancer cells"
    # use updated UniProt IDs
    inew <- which(dat$Uniprot.new != "")
    dat$Uniprot.accession[inew] <- dat$Uniprot.new[inew]
    # remove Q8IWE2, which is duplicated with opposite expression ratio
    dat <- dat[dat$Uniprot.accession!="Q8IWE2", ]
    pcomp <- protcomp(dat$Uniprot.accession)
    up2 <- dat$Ratio..H.L. > 1
  } else if(study=="MVC+12") {
    # 20160413 spheriod hypoxia, McMahon et al., 2012
    # MVC+12_perinecrotic, MVC+12_necrotic
    return(pdat_3D(dataset))
  } else if(study=="MHG+12") {
    # 20160415 MCF-7 tumourspheres, Morrison et al., 2012
    # MHG+12_P5, MHG+12_P2
    return(pdat_3D(dataset))
  } else if(study=="HXS+06") {
    # 20160415 leukemic U937 cells, Han et al., 2006
    dat <- read.csv(paste0(datadir, "HXS+06.csv.xz"), as.is=TRUE)
    description <- "U937"
    # update IDs with new ones
    inew <- dat$UniProt.new != ""
    dat$Swiss.Prot.accession.no.[inew] <- dat$UniProt.new[inew]
    up2 <- dat$Mean.fold.H.N > 1
    dat <- cleanup(dat, "Swiss.Prot.accession.no.", up2)
    pcomp <- protcomp(dat$Swiss.Prot.accession.no.)
  } else if(study=="RHD+13") {
    # 20160419 A431 cells, Ren et al., 2013
    # RHD+13_Hx48, RHD+13_Hx72, RHD+13_ReOx
    dat <- read.csv(paste0(datadir, "RHD+13.csv.xz"), as.is=TRUE)
    description <- paste("A431 epithelial carcinoma cells", stage)
    # columns with the ratios
    if(stage=="Hx48") icol <- grep("115", colnames(dat))
    if(stage=="Hx72") icol <- grep("116", colnames(dat))
    if(stage=="ReOx") icol <- grep("117", colnames(dat))
    # remove NA values (indicates p-value >= 0.05)
    dat <- dat[!is.na(dat[, icol]), ]
    # get highly differential proteins
    dat <- dat[dat[, icol] > sqrt(2) | dat[, icol] < 1/sqrt(2), ]
    # get known UniProt IDs
    dat <- check_IDs(dat, "Accession")
    up2 <- dat[, icol[1]] > 1
    dat <- cleanup(dat, "Accession", up2)
    pcomp <- protcomp(dat$Accession)
  } else if(study=="BMJ+11") {
    # 20160713 DU145 cells prolonged hypoxia, van den Beucken et al., 2011
    dat <- read.csv(paste0(datadir, "BMJ+11.csv.xz"), as.is=TRUE)
    description <- "DU145"
    # keep proteins detected in prolonged hypoxia
    dat <- dat[dat$induced_prolonged | dat$repressed_prolonged, ]
    up2 <- dat$induced_prolonged
    dat <- cleanup(dat, "uniprot", up2)
    pcomp <- protcomp(dat$uniprot)
  } else if(study=="FWH+13") {
    # 20160716 THP-1 macrophages CV (control virus) hypoxia, Fuhrmann et al., 2013
    dat <- read.csv(paste0(datadir, "FWH+13.csv.xz"), as.is=TRUE)
    description <- "THP-1 macrophages"
    up2 <- dat$Norm.CH > 0
    pcomp <- protcomp(dat$UniProt)
  } else if(study=="HWA+16") {
    # 20160716 U87MG and 786-O translatome, Ho et al., 2016
    dat <- read.csv(paste0(datadir, "HWA+16.csv.xz"), as.is=TRUE)
    # keep those with fold change < 0.5 or > 2 (log2 < -1 or > 1)
    dat <- dat[ dat$Hypoxia.Heavy - dat$Normoxia.Heavy > 1 |
                dat$Hypoxia.Heavy - dat$Normoxia.Heavy < -1, ]
    description <- "U87MG and 786-O cancer cells"
    up2 <- dat$Hypoxia.Heavy - dat$Normoxia.Heavy > 0
    dat <- check_IDs(dat, "Uniprot.Accession")
    pcomp <- protcomp(dat$Uniprot.Accession)
  } else if(study=="RKP+14") {
    # 20160718 organotypic spheroids, Rajcevic et al., 2014
    return(pdat_3D(dataset))
  } else if(study=="CBW+11") {
    # 20160720 neuroblastoma cells, Cifani et al., 2011
    dat <- read.csv(paste0(datadir, "CBW+11.csv.xz"), as.is=TRUE)
    description <- "SK-N-BE(2)c; IMR-32"
    # remove "-1" suffix for isoform 1
    dat$UniProt <- gsub("-1", "", dat$UniProt)
    dat <- check_IDs(dat, "UniProt")
    up2 <- dat$Be2c > 1
    dat <- cleanup(dat, "UniProt", up2)
    pcomp <- protcomp(dat$UniProt)
  } else if(study=="WRK+14") {
    # 20160721 3D spheroids / 2D culture, Wrzesinski et al., 2014
    return(pdat_3D(dataset))
  } else if(study=="DPL+10") {
    # 20160722 B104 rat neuroblastoma cells, Datta et al., 2010
    dat <- read.csv(paste0(datadir, "DPL+10.csv.xz"), as.is=TRUE)
    description <- "B104"
    # select highly changed proteins
    dat <- dat[dat$HYP.LSC > 1.2 | dat$HYP.LSC < 0.83, ]
    dat <- check_IDs(dat, "UniProt", aa_file=paste0(extdatadir, "/aa/rat/DPL+10_aa.csv.xz"))
    up2 <- dat$HYP.LSC > 1
    dat <- cleanup(dat, "UniProt", up2)
    pcomp <- protcomp(dat$UniProt, aa_file=paste0(extdatadir, "/aa/rat/DPL+10_aa.csv.xz"))
  } else if(study=="LCS16") {
    # 20160728 HCT116 transcription and translation, Lai et al., 2016
    # LCS16_transcription, LCS16_translation
    dat <- read.csv(paste0(datadir, "LCS16.csv.xz"), as.is=TRUE)
    description <- paste("HCT116 colon cancer", stage)
    # select the experiment
    icol <- grep(stage, tolower(colnames(dat)))
    idiff <- sapply(dat[, icol[1]] | dat[, icol[2]], isTRUE)
    dat <- dat[idiff, ]
    # which proteins (genes) are up-regulated
    iicol <- icol[grep("up", colnames(dat)[icol])]
    up2 <- sapply(dat[, iicol], isTRUE)
    dat <- cleanup(dat, "UniProt", up2)
    pcomp <- protcomp(dat$UniProt)
  } else if(study=="DYL+14") {
    # 20160729 A431 cells, Dutta et al., 2014
    # DYL+14_Hx48-S, DYL+14_Hx72-S, DYL+14_ReOx-S,
    # DYL+14_Hx48-P, DYL+14_Hx72-P, DYL+14_ReOx-P
    dat <- read.csv(paste0(datadir, "DYL+14.csv.xz"), as.is=TRUE)
    description <- paste("A431 epithelial carcinoma cells", stage)
    # -S (supernatant) and -P (pellet) datasets
    if(stage=="Hx48-S") icol <- grep("114", colnames(dat))
    if(stage=="Hx72-S") icol <- grep("115", colnames(dat))
    if(stage=="ReOx-S") icol <- grep("116", colnames(dat))
    if(stage=="Hx48-P") icol <- grep("118", colnames(dat))
    if(stage=="Hx72-P") icol <- grep("119", colnames(dat))
    if(stage=="ReOx-P") icol <- grep("121", colnames(dat))
    # now deal with -P (pellet) datasets
    # divide (expt / Nx-S) by (Nx-P / Nx-S)
    if(grepl("-P", stage)) dat[, icol[1]] <- dat[, icol[1]] / dat$X117.113
    # keep significantly changed proteins based on p-value and ratio
    if(!grepl("-P", stage)) dat <- dat[dat[, icol[2]] < 0.05, ]
    dat <- dat[dat[, icol[1]] > sqrt(2) | dat[, icol[1]] < 1/sqrt(2), ]
    dat <- check_IDs(dat, "Accession")
    pcomp <- protcomp(dat$Accession)
    up2 <- dat[, icol[1]] > 1
  } else if(study=="RSE+16") {
    # 20160729 adipose-derived stem cells, Riis et al., 2016
    return(pdat_secreted(dataset))
  } else if(study=="VTMF13") {
    # 20160804 neuroblastoma cell line, Villeneuve et al., 2013
    dat <- read.csv(paste0(datadir, "VTMF13.csv.xz"), as.is=TRUE)
    description <- "SH-SY5Y neuroblastoma cells"
    # keep proteins with large expression ratio
    dat <- dat[dat$Ratio.H.L.Normalized > 1.2 | dat$Ratio.H.L.Normalized < 0.83, ]
    # find known UniProt IDs
    dat <- check_IDs(dat, "Uniprot")
    up2 <- dat$Ratio.H.L.Normalized > 1.2
    dat <- cleanup(dat, "Uniprot", up2)
    pcomp <- protcomp(dat$Uniprot)
  } else if(study=="BRA+10") {
    # 20160805 placental tissue secretome, Blankley et al., 2010
    return(pdat_secreted(dataset))
  } else if(study=="LAR+12") {
    # 20160826 rat heart ischemia, Li et al., 2012
    dat <- read.csv(paste0(datadir, "LAR+12.csv.xz"), as.is=TRUE)
    description <- "H9C2"
    up2 <- dat$Isch.Ctrl > 1
    dat <- cleanup(dat, "UniProt", up2)
    pcomp <- protcomp(dat$UniProt, aa_file=paste0(extdatadir, "/aa/mouse/LAR+12_aa.csv.xz"))
  } else if(study=="YLW+16") {
    # 20161109 HT29 colon cancer cell 3D/2D, Yue et al., 2011
    return(pdat_3D(dataset))
  } else if(study=="XCJ+16") {
    # 20161119 cardiomyocytes CoCl2 (hypoxia mimetic) or SAL (anti-hypoxic), Xu et al., 2016
    # XCJ+16_CoCl2, XCJ+16_SAL
    dat <- read.csv(paste0(datadir, "XCJ+16.csv.xz"), as.is=TRUE)
    description <- paste("cardiomyocytes", stage)
    # use selected dataset
    icol <- grep(paste0("Log2.", stage), colnames(dat))
    dat <- dat[!is.na(dat[, icol]), ]
    pcomp <- protcomp(dat$Entry, aa_file=paste0(extdatadir, "/aa/rat/XCJ+16_aa.csv.xz"))
    up2 <- dat[, icol] > 0
  } else if(study=="CGH+17") {
    # 20190324 mouse cardiac whole-cell lysate, Cosme et al., 2017
    # CGH+17_whole
    return(.pdat_multi(dataset))
  } else if(study=="CLY+18") {
    # 20190324 HCT116 cells, Chen et al., 2018
    # CLY+18_proteome
    return(.pdat_multi(dataset))
  } else if(study=="DCH+14") {
    # 20160109 hypoxia-regulated proteins in mouse breast cancer, Djidja et al., 2014
    # 20191204 use all proteins in Table S1
    dat <- read.csv(paste0(datadir, "DCH+14.csv.xz"), as.is=TRUE)
    description <- "mouse 4T1 cells"
    # remove isoform suffixes
    dat$Accession <- substr(dat$Accession, 1, 6)
    dat <- check_IDs(dat, "Accession", aa_file=paste0(extdatadir, "/aa/mouse/DCH+14_aa.csv.xz"))
    up2 <- dat$Regulation == "Up"
    pcomp <- protcomp(dat$Accession, aa_file=paste0(extdatadir, "/aa/mouse/DCH+14_aa.csv.xz"))
  } else if(study=="ZXS+17") {
    # 20190407 glioblastoma cells, Zhang et al., 2017
    dat <- read.csv(paste0(datadir, "ZXS+17.csv.xz"), as.is=TRUE)
    description <- "U87 and U251 glioblastoma cells"
    pcomp <- protcomp(dat$Accession)
    up2 <- dat$Fold.Change > 0
  } else if(study=="GBH+18") {
    # 20190409 SW620 cells 1% / 21% O2, Greenhough et al., 2018
    dat <- read.csv(paste0(datadir, "GBH+18.csv.xz"), as.is=TRUE)
    description <- "SW620 colorectal cancer cells"
    up2 <- dat$Av..Fold.change..hyp.norm. > 1
    dat <- cleanup(dat, "Entry", up2)
    pcomp <- protcomp(dat$Entry)
  } else if(study=="WTG+18") {
    # 20190413 mesenchymal stem cells, Wobma et al., 2018
    dat <- read.csv(paste0(datadir, "WTG+18.csv.xz"), as.is=TRUE)
    description <- "adipose-derived mesehnchymal stem cells"
    dat <- check_IDs(dat, "UniProt.Accession")
    pcomp <- protcomp(dat$UniProt.Accession)
    up2 <- dat$Normalized.Ratio..Hypoxia.MSC..Control.MSC. > 1
  } else if(study=="LKK+18") {
    # 20191127 mesenchymal stem cells, Lee et al., 2018
    dat <- read.csv(paste0(datadir, "LKK+18.csv.xz"), as.is=TRUE)
    description <- "hUCB mesehnchymal stem cells"
    dat <- check_IDs(dat, "Entry")
    up2 <- dat$Regulation == "up"
    dat <- cleanup(dat, "Entry", up2)
    pcomp <- protcomp(dat$Entry)
  } else if(study=="CSK+19") {
    # 20191204 HeLa cells, Chachami et al., 2019
    dat <- read.csv(paste0(datadir, "CSK+19.csv.xz"), as.is=TRUE)
    description <- "HeLa cervical cancer cells"
    up2 <- dat$input_Log2ratioHL_firstIP > 0
    pcomp <- protcomp(dat$Entry)
  } else if(study=="KAN+19") {
    # 20191226 human umbilical vein ECs, Kugeratski et al., 2019
    # KAN+19_proteome
    return(.pdat_multi(dataset))
  } else if(study=="SBB+06") {
    # 20161109 mouse hypoxia-adapated malignant melanoma, Stockwin et al., 2006
    # 20200115 added to canprot
    dat <- read.csv(paste0(datadir, "SBB+06.csv.xz"), as.is = TRUE)
    description <- "mouse malignant melanoma plasma membrane"
    dat <- check_IDs(dat, "PAN", aa_file = paste0(extdatadir, "/aa/mouse/SBB+06_aa.csv.xz"))
    up2 <- dat$AR > 1
    dat <- cleanup(dat, "PAN", up2)
    pcomp <- protcomp(dat$PAN, aa_file = paste0(extdatadir, "/aa/mouse/SBB+06_aa.csv.xz"))
  } else if(study=="BCMS20") {
    # 20200118 MCF-7 breast cancer cells, Bush et al., 2020
    dat <- read.csv(paste0(datadir, "BCMS20.csv.xz"), as.is = TRUE)
    description <- "MCF-7 breast cancer cells"
    up2 <- dat$log.ratio.M.L.. > 0
    dat <- cleanup(dat, "Entry", up2)
    pcomp <- protcomp(dat$Entry)
  } else if(study=="RVN+20") {
    # 20200403 prostate cancer cell hypoxia with different drug treatments and radiation, Ross et al., 2020
    # RVN+20_DMSO, RVN+20_NO.sul, RVN+20_sul, RVN+20_DMSO.4Gy, RVN+20_NO.sul.4Gy, RVN+20_sul.4Gy
    dat <- read.csv(paste0(datadir, "RVN+20.csv.xz"), as.is = TRUE)
    description <- paste("PC-3 prostate cancer cells in", stage)
    # use data for this treatment
    icol <- grep(paste0("^", stage, "_"), colnames(dat))
    dat <- dat[, c(1, 2, icol)]
    # calculate fold change and keep differential proteins
    ihypoxic <- grep("hypoxic", colnames(dat))
    inormoxic <- grep("normoxic", colnames(dat))
    FC <- dat[, ihypoxic] / dat[, inormoxic]
    dat <- cbind(dat, FC = FC)
    dat <- dat[!is.na(dat$FC), ]
    dat <- dat[dat$FC > 3/2 | dat$FC < 2/3, ]
    up2 <- dat$FC > 3/2
    pcomp <- protcomp(dat$UniProt)
  } else if(study=="SPJ+20") {
    # 20200405 canine OS cells, Song et al., 2020
    # SPJ+20_POS, SPJ+20_HMPOS
    dat <- read.csv(paste0(datadir, "SPJ+20.csv.xz"), as.is = TRUE)
    description <- paste("canine", stage, "cells")
    icol <- match(stage, colnames(dat))
    dat <- dat[!is.na(dat[, icol]), ]
    dat <- check_IDs(dat, "Accession", aa_file = paste0(extdatadir, "/aa/dog/SPJ+20_aa.csv.xz"))
    up2 <- dat[, icol] > 0
    pcomp <- protcomp(dat$Accession, aa_file = paste0(extdatadir, "/aa/dog/SPJ+20_aa.csv.xz"))
  } else if(study=="LLL+19") {
    # 20200406 human periodontal ligament cells, Li et al., 2019
    dat <- read.csv(paste0(datadir, "LLL+19.csv.xz"), as.is = TRUE)
    description <- "human periodontal ligament cells"
    dat <- check_IDs(dat, "Accession")
    up2 <- up2 <- dat$FC > 1
    pcomp <- protcomp(dat$Accession)
  } else if(study=="GPT+19") {
    # 20200406 MIAPaCa-2 pancreatic cancer cells pulse/trace, Gupta et al., 2019
    # GPT+19_Light.S, GPT+19_Heavy.S, GPT+19_Light.SF, GPT+19_Heavy.SF
    dat <- read.csv(paste0(datadir, "GPT+19.csv.xz"), as.is = TRUE)
    dtxt <- gsub(".SF", " serum-free", stage)
    dtxt <- gsub(".S", " serum replete", dtxt)
    description <- paste("MIAPaCa-2 pancreatic cancer cells pulse/trace", dtxt)
    dat <- check_IDs(dat, "Accession")
    # use indicated experiment
    icol <- grep(paste0(stage, "$"), colnames(dat))
    dat <- dat[!is.na(dat[, icol]), ]
    dat <- dat[dat[, icol] > 2 | dat[, icol] < 0.5, ]
    up2 <- dat[, icol] > 2
    dat <- cleanup(dat, "Accession", up2)
    pcomp <- protcomp(dat$Accession)
  } else stop(paste("hypoxia dataset", dataset, "not available"))
  print(paste0("pdat_hypoxia: ", description, " [", dataset, "]"))
  # use the up2 from the cleaned-up data, if it exists 20191120
  if("up2" %in% colnames(dat)) up2 <- dat$up2
  return(list(dataset=dataset, pcomp=pcomp, up2=up2, description=description))
}
