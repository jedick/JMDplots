# JMDplots/pdat_.R
# get protein IDs for differentially expressed proteins or genes
# combined pdat_TCGA.R and pdat_HPA.R 20200301
# pdat_aneuploidy added 20200505

# retrieve protein IDs for differentially expressed genes in TCGA / GTEx [via GEPIA2] 20111123 jmd
pdat_TCGA <- function(dataset = 2020) {
  if(identical(dataset, 2020)) {
    return(c(
             "GEPIA2_ACC", "GEPIA2_BLCA", "GEPIA2_BRCA", "GEPIA2_CESC", "GEPIA2_COAD", 
             "GEPIA2_DLBC", "GEPIA2_ESCA", "GEPIA2_GBM", "GEPIA2_HNSC", "GEPIA2_KICH", 
             "GEPIA2_KIRC", "GEPIA2_KIRP", "GEPIA2_LAML", "GEPIA2_LGG", "GEPIA2_LIHC", 
             "GEPIA2_LUAD", "GEPIA2_LUSC", "GEPIA2_OV", "GEPIA2_PAAD", "GEPIA2_PRAD", 
             "GEPIA2_READ", "GEPIA2_SKCM", "GEPIA2_STAD", "GEPIA2_TGCT", "GEPIA2_THCA", 
             "GEPIA2_THYM", "GEPIA2_UCEC", "GEPIA2_UCS"
             ))
  }
  # remove tags
  dataset <- strsplit(dataset, "=")[[1]][1]
  # get study and stage/condition
  study <- strsplit(dataset, "_")[[1]][1]
  stage <- paste(strsplit(dataset, "_")[[1]][-1], collapse="_")
  extdatadir <- system.file("extdata", package="JMDplots")
  datadir <- paste0(extdatadir, "/diffexpr/pancan/")
  if(study=="GEPIA2") {
    # 20191125 expression levels derived from GEPIA2
    # GEPIA2_ACC, GEPIA2_BLCA, ...
    dat <- read.csv(paste0(datadir, "GEPIA2.csv.xz"), as.is=TRUE)
#    # 20191201 remove HBB, HBA1, HBA2 because they are very high in normal tissue
#    # for THYM and DLBC (producing anomalously low values for phylostrata difference)
#    dat <- dat[!dat$Gene.Symbol %in% c("HBB", "HBA1", "HBA2"), ]
    description <- stage
    icol <- grep(paste0("^", stage), colnames(dat))
    # keep at least 2-fold differentially expressed genes
    ilog2 <- abs(dat[, icol]) >= 1
    ilog2[is.na(ilog2)] <- FALSE
    dat <- dat[ilog2, ]
    up2 <- dat[, icol] > 0
    dat <- cleanup(dat, "Entry", up2)
    pcomp <- protcomp(dat$Entry)
  } else stop(paste("TCGA dataset", dataset, "not available"))
  print(paste0("pdat_TCGA: ", description, " [", dataset, "]"))
  # use the up2 from the cleaned-up data, if it exists 20190429
  if("up2" %in% colnames(dat)) up2 <- dat$up2
  return(list(dataset=dataset, pcomp=pcomp, up2=up2, description=description))
}

# retrieve protein IDs for Human Protein Atlas 20111121 jmd
pdat_HPA <- function(dataset = 2020) {
  if(identical(dataset, 2020)) {
    return(c(
             "HPA19_1", "HPA19_2", "HPA19_3", "HPA19_4", "HPA19_5", "HPA19_6",
             "HPA19_7", "HPA19_8", "HPA19_9", "HPA19_10", "HPA19_11", "HPA19_12",
             "HPA19_13", "HPA19_14", "HPA19_15", "HPA19_16", "HPA19_17", "HPA19_18"
             ))
  }
  # remove tags
  dataset <- strsplit(dataset, "=")[[1]][1]
  # get study and stage/condition
  study <- strsplit(dataset, "_")[[1]][1]
  stage <- paste(strsplit(dataset, "_")[[1]][-1], collapse="_")
  extdatadir <- system.file("extdata", package="JMDplots")
  datadir <- paste0(extdatadir, "/diffexpr/pancan/")
  if(study=="HPA19") {
    # 20191121 expression levels derived from Human Protein Atlas, 2019
    # HPA19_1 .. HPA19_18
    # normal tissues for different cancers -- from cancer chapters on https://www.proteinatlas.org/learn/dictionary/pathology)
    cancer_tissues <- list(
      c("breast", "breast cancer"),
      c("cervix, uterine", "cervical cancer"),
      c("colon", "colorectal cancer"),
      c("endometrium 1", "endometrial cancer"),
      c("cerebral cortex", "glioma"),
      c("salivary gland", "head and neck cancer"),
      c("liver", "liver cancer"),
      c("lung", "lung cancer"),
      c("lymph node", "lymphoma"),
      c("skin 1", "melanoma"),
      c("ovary", "ovarian cancer"),
      c("pancreas", "pancreatic cancer"),
      c("prostate", "prostate cancer"),
      c("kidney", "renal cancer"),
      c("stomach 1", "stomach cancer"),
      c("testis", "testis cancer"),
      c("thyroid gland", "thyroid cancer"),
      c("urinary bladder", "urothelial cancer")
    )
    dat <- read.csv(paste0(datadir, "HPA19.csv.xz"), as.is=TRUE)
    # get name of cancer and normal tissue
    stage <- as.numeric(stage)
    normal <- cancer_tissues[[stage]][1]
    cancer <- cancer_tissues[[stage]][2]
    description <- paste(cancer, "/", normal)
    if(normal == "cervix, uterine") normal <- "cervix"
    inormal <- match(make.names(normal), colnames(dat))
    icancer <- match(make.names(cancer), colnames(dat))
    # keep highly differential proteins
    absdiff <- abs(dat[, icancer] - dat[, inormal])
    isdiff <- absdiff >= 2.5
    isdiff[is.na(isdiff)] <- FALSE
    dat <- dat[isdiff, ]
    # drop missing proteins
    up2 <- dat[, icancer] - dat[, inormal] >= 2.5
    dat <- cleanup(dat, "Entry", up2)
    pcomp <- protcomp(dat$Entry)
  } else stop(paste("HPA dataset", dataset, "not available"))
  print(paste0("pdat_HPA: ", description, " [", dataset, "]"))
  # use the up2 from the cleaned-up data, if it exists 20190429
  if("up2" %in% colnames(dat)) up2 <- dat$up2
  return(list(dataset=dataset, pcomp=pcomp, up2=up2, description=description))
}

pdat_aneuploidy <- function(dataset = 2020) {
  if(identical(dataset, 2020)) {
    return("TNC+19")
  }
  # remove tags
  dataset <- strsplit(dataset, "=")[[1]][1]
  # get study and stage/condition
  study <- strsplit(dataset, "_")[[1]][1]
  stage <- paste(strsplit(dataset, "_")[[1]][-1], collapse = "_")
  extdatadir <- system.file("extdata", package="JMDplots")
  datadir <- paste0(extdatadir, "/diffexpr/aneuploidy/")
  if(study=="TNC+19") {
    # 20200425 yeast aneuploidy, Tsai et al., 2019
    dat <- read.csv(paste0(datadir, "TNC+19.csv.xz"), as.is=TRUE)
    description <- "yeast aneuploidy"
    up2 <- dat$log2FoldChange > 0
    pcomp <- protcomp(dat$Entry, aa_file = paste0(extdatadir, "/aa/yeast/TNC+19_aa.csv.xz"))
  } else stop(paste("aneuploidy dataset", dataset, "not available"))
  print(paste0("pdat_aneuploidy: ", description, " [", dataset, "]"))
  # use the up2 from the cleaned-up data, if it exists 20190407
  if("up2" %in% colnames(dat)) up2 <- dat$up2
  return(list(dataset = dataset, pcomp = pcomp, up2 = up2, description = description))
}

pdat_yeast_stress <- function(dataset = 2020) {
  if(identical(dataset, 2020)) {
    return(c(
      "GSK+00_X1M.sorbitol...5.min", "GSK+00_X1M.sorbitol...15.min", "GSK+00_X1M.sorbitol...30.min", "GSK+00_X1M.sorbitol...45.min",
      "GSK+00_X1M.sorbitol...60.min", "GSK+00_X1M.sorbitol...90.min", "GSK+00_X1M.sorbitol...120.min",
      "GSK+00_Hypo.osmotic.shock...5.min", "GSK+00_Hypo.osmotic.shock...15.min", "GSK+00_Hypo.osmotic.shock...30.min",
      "GSK+00_Hypo.osmotic.shock...45.min", "GSK+00_Hypo.osmotic.shock...60.min"
    ))
  }
  # remove tags
  dataset <- strsplit(dataset, "=")[[1]][1]
  # get study and stage/condition
  study <- strsplit(dataset, "_")[[1]][1]
  stage <- paste(strsplit(dataset, "_")[[1]][-1], collapse = "_")
  extdatadir <- system.file("extdata", package="JMDplots")
  datadir <- paste0(extdatadir, "/diffexpr/yeast_stress/")
  if(study=="GSK+00") {
    # 20200508 yeast gene expression, Gasch et al., 2000
    # GSK+00_X1M.sorbitol...5.min and others
    dat <- read.csv(paste0(datadir, "GSK+00.csv.xz"), as.is=TRUE)
    dtxt <- gsub(".min", " min", stage)
    dtxt <- gsub("X1M.sorbitol...", "1M sorbitol at ", dtxt)
    dtxt <- gsub("Hypo.osmotic.shock...", "hypoosmotic shock at ", dtxt)
    description <- paste("yeast transcriptome in", dtxt)
    icol <- grep(stage, colnames(dat))
    if(length(icol) > 1) stop("multiple columns selected")
    # get genes with at least 2-fold change
    dat <- dat[!is.na(dat[, icol]), ]
    dat <- dat[abs(dat[, icol]) > log2(1.5), ]
    up2 <- dat[, icol] > log2(1.5)
    dat <- cleanup(dat, "Entry", up2)
    pcomp <- protcomp(dat$Entry, aa_file = paste0(extdatadir, "/aa/yeast/GSK+00_aa.csv.xz"))
  } else stop(paste("yeast_stress dataset", dataset, "not available"))
  print(paste0("pdat_yeast_stress: ", description, " [", dataset, "]"))
  # use the up2 from the cleaned-up data, if it exists 20190407
  if("up2" %in% colnames(dat)) up2 <- dat$up2
  return(list(dataset = dataset, pcomp = pcomp, up2 = up2, description = description))
}

# Get differential expression data from Fabre et al., 2019  20210401
pdat_fly <- function(dataset = NULL) {
  if(is.null(dataset)) {
    return(c(
             "FKL+19_mRNA", "FKL+19_protein"
             ))
  }
  # remove tags
  dataset <- strsplit(dataset, "=")[[1]][1]
  # get study and stage/condition
  study <- strsplit(dataset, "_")[[1]][1]
  stage <- paste(strsplit(dataset, "_")[[1]][-1], collapse = "_")
  extdatadir <- system.file("extdata", package = "JMDplots")
  if(study=="FKL+19") {
    # 20200102 Drosophila adult vs embryo, Fabre et al., 2019
    # FKL+19_mRNA, FKL+19_protein
    dat <- read.csv(file.path(extdatadir, "diffexpr/development/FKL+19.csv"), as.is = TRUE)
    description <- paste("Drosophila adult / embryo", stage)
    if(stage == "protein") {
      # get differentially expressed proteins
      dat <- dat[abs(dat$Log2.ratio.Adult.Embryo) >= 1 & dat$X.log10.BH.corrected.p.value.Adult.vs.Embryo >= 2, ]
      up2 <- dat$Log2.ratio.Adult.Embryo > 0
    }
    if(stage == "mRNA") {
      # calculate ratios
      dat <- cbind(dat, ratio = dat$log10.mRNA.adult - dat$log10.mRNA.embryo)
      dat <- dat[abs(dat$ratio) > 1, ]
      up2 <- dat$ratio > 0
    }
    aa_file <- file.path(extdatadir, "aa/fly/FKL+19_aa.csv")
    updates_file <- file.path(extdatadir, "aa/uniprot_updates.csv")
    dat <- check_IDs(dat, "Entry", updates_file = updates_file, aa_file = aa_file)
    pcomp <- protcomp(dat$Entry, aa_file = aa_file)
  } else stop(paste("fly dataset", dataset, "not available"))
  print(paste0("pdat_fly: ", description, " [", dataset, "]"))
  # use the up2 from the cleaned-up data, if it exists 20190407
  if("up2" %in% colnames(dat)) up2 <- dat$up2
  return(list(dataset = dataset, pcomp = pcomp, up2 = up2, description = description))
}
