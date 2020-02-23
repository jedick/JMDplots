# JMDplots/pdat_CH16.R
# common gene expression changes in cancer
# (Chen and He, 2016; doi:10.1093/molbev/msv212)
# 20200223

pdat_CH16 <- function(dataset = 2020, basis = "rQEC") {
  if(identical(dataset, 2020)) {
    return("CH16")
  }
  # remove tags
  dataset <- strsplit(dataset, "=")[[1]][1]
  # get study and stage/condition
  study <- strsplit(dataset, "_")[[1]][1]
  stage <- paste(strsplit(dataset, "_")[[1]][-1], collapse = "_")
  extdatadir <- system.file("extdata", package = "JMDplots")
  datadir <- paste0(extdatadir, "/expression/pancan/")
  if(study=="CH16") {
    # 20200223 common gene expression across cancer types, Chen and He, 2016
    dat <- read.csv(paste0(datadir, "CH16.csv"), as.is=TRUE)
    description <- "common gene expression in cancer"
    up2 <- dat$Up.down.regulated.in.cancer == "Up-regulated"
    dat <- cleanup(dat, "Entry", up2)
    pcomp <- protcomp(dat$Entry, basis)
  } else stop(paste("CH16 dataset", dataset, "not available"))
  print(paste0("pdat_CH16: ", description, " [", dataset, "]"))
  # use the up2 from the cleaned-up data, if it exists 20190407
  if("up2" %in% colnames(dat)) up2 <- dat$up2
  return(list(dataset = dataset, basis = basis, pcomp = pcomp, up2 = up2, description = description))
}
