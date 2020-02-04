# JMDplots/pdat_HPA.R
# retrieve protein IDs for Human Protein Atlas
# 20111121 jmd

pdat_HPA <- function(dataset = 2020, basis = "rQEC") {
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
  datadir <- paste0(extdatadir, "/expression/pancan/")
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
    pcomp <- protcomp(dat$Entry, basis=basis)
  } else stop(paste("HPA dataset", dataset, "not available"))
  print(paste0("pdat_HPA: ", description, " [", dataset, "]"))
  # use the up2 from the cleaned-up data, if it exists 20190429
  if("up2" %in% colnames(dat)) up2 <- dat$up2
  return(list(dataset=dataset, pcomp=pcomp, up2=up2, description=description))
}

