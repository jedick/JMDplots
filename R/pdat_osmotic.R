# JMDplots/.pdat_osmotic.R
# retrieve protein IDs for hyperosmotic experiments
# 20160926 jmd
# 20200418 renamed to .pdat_osmotic() (for 2017 compilation only)
# 20240426 moved from canprot package

.pdat_osmotic <- function(dataset = 2017) {
  if(identical(dataset, 2017)) {
    return(c(
             "PW08_2h", "PW08_10h", "PW08_12h",
             "WCM+09",
             "OBBH11=ASC",
             "CCC+12_25mM", "CCC+12_100mM",
             "KKG+12_25C_aw0.985", "KKG+12_14C_aw0.985", "KKG+12_25C_aw0.967", "KKG+12_14C_aw0.967",
             "CCCC13_25mM", "CCCC13_100mM", "TSZ+13",
             "GSC14_t30a", "GSC14_t30b", "GSC14_t30c",
             "CLG+15", "KLB+15_trans-suc=transcriptome", "KLB+15_trans-NaCl=transcriptome", 
             "KLB+15_prot-suc", "KLB+15_prot-NaCl",
             "LDB+15_all", "LDB+15_high", "YDZ+15",
             "RBP+16"
             ))
  }
  # remove tags
  dataset <- strsplit(dataset, "=")[[1]][1]
  # get study and stage/condition
  study <- strsplit(dataset, "_")[[1]][1]
  stage <- paste(strsplit(dataset, "_")[[1]][-1], collapse="_")
  extdatadir <- system.file("extdata", package="JMDplots")
  datadir <- paste0(extdatadir, "/diffexpr/osmotic/")
  if(study=="TSZ+13") {
    # 20161113 eel gill (Anguilla japonica), Tse et al., 2013
    dat <- read.csv(paste0(datadir, "TSZ+13.csv.xz"), as.is=TRUE)
    description <- "eel gill"
    up2 <- dat$Fold.Change..FW.SW. > 1
    dat <- cleanup(dat, "Entry", up2)
    pcomp <- protcomp(dat$Entry)
  } else if(study %in% c("PW08", "WCM+09", "CCC+12", "CCCC13", "GSC14", "LDB+15")) {
    # 20200411 datasets from the 2017 compilation that have been moved to pdat_glucose.R
    return(pdat_glucose(dataset))
  } else if(study %in% c("OBBH11", "CLG+15", "YDZ+15", "RBP+16")) {
    # 20200418 datasets from the 2017 compilation that have been moved to pdat_osmotic_euk.R
    return(pdat_osmotic_euk(dataset))
  } else if(study %in% c("KKG+12", "KLB+15")) {
    # 20200418 datasets from the 2017 compilation that have been moved to pdat_osmotic_bact.R
    return(pdat_osmotic_bact(dataset))
  } else stop(paste("osmotic dataset", dataset, "not available"))
  print(paste0(".pdat_osmotic: ", description, " [", dataset, "]"))
  # use the up2 from the cleaned-up data, if it exists 20191120
  if("up2" %in% colnames(dat)) up2 <- dat$up2
  return(list(dataset=dataset, pcomp=pcomp, up2=up2, description=description))
}
