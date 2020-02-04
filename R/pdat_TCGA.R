# JMDplots/pdat_TCGA.R
# retrieve protein IDs for differentially expressed genes in TCGA / GTEx [via GEPIA2]
# 20111123 jmd

pdat_TCGA <- function(dataset = 2020, basis = "rQEC") {
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
  datadir <- paste0(extdatadir, "/expression/pancan/")
  if(study=="GEPIA2") {
    # 20191125 expression levels derived from GEPIA2
    # GEPIA2_BLCA ...
    dat <- read.csv(paste0(datadir, "GEPIA2.csv.xz"), as.is=TRUE)
    # 20191201 remove HBB, HBA1, HBA2 because they are very high in normal tissue
    # for THYM and DLBC (producing anomalously low values for phylostrata difference)
#    dat <- dat[!dat$Gene.Symbol %in% c("HBB", "HBA1", "HBA2"), ]
    description <- stage
    icol <- grep(paste0("^", stage), colnames(dat))
    # keep at least 2-fold differentially expressed genes
    ilog2 <- abs(dat[, icol]) >= 1
    ilog2[is.na(ilog2)] <- FALSE
    dat <- dat[ilog2, ]
    up2 <- dat[, icol] > 0
    dat <- cleanup(dat, "Entry", up2)
    pcomp <- protcomp(dat$Entry, basis=basis)
  } else stop(paste("TCGA dataset", dataset, "not available"))
  print(paste0("pdat_TCGA: ", description, " [", dataset, "]"))
  # use the up2 from the cleaned-up data, if it exists 20190429
  if("up2" %in% colnames(dat)) up2 <- dat$up2
  return(list(dataset=dataset, pcomp=pcomp, up2=up2, description=description))
}

