# JMDplots/pdat_multi.R
# retrieve protein IDs from studies that have multiple conditions
# (i.e. both proteome and secretome in hypoxia)
# 20191204 extracted from pdat_secreted.R

.pdat_multi <- function(dataset = 2020) {
  if(identical(dataset, 2020)) {
    return(c("CGH+17_exosomes", "CGH+17_secretome", "CGH+17_whole",
             "CLY+18_proteome", "CLY+18_secretome",
             "KAN+19_proteome", "KAN+19_secretome"
             ))
  }
  # remove tags
  dataset <- strsplit(dataset, "=")[[1]][1]
  # get study and stage/condition
  study <- strsplit(dataset, "_")[[1]][1]
  stage <- paste(strsplit(dataset, "_")[[1]][-1], collapse="_")
  extdatadir <- system.file("extdata", package="JMDplots")
  datadir <- paste0(extdatadir, "/canH2O/multi/")
  if(study=="CGH+17") {
    # 20190324 mouse cardiac fibroblast exosomes, secretome, whole-cell lysate, Cosme et al., 2017
    # CGH+17_exosomes, CGH+17_secretome, CGH+17_whole
    dat <- read.csv(paste0(datadir, "CGH+17.csv.xz"), as.is=TRUE)
    description <- paste("mouse cardiac fibroblasts", stage)
    # use selected dataset
    icol <- grep(stage, colnames(dat))
    dat <- dat[!is.na(dat[, icol[1]]), ]
    pcomp <- protcomp(dat$Entry, aa_file=paste0(extdatadir, "/aa/mouse/CGH+17_aa.csv.xz"))
    up2 <- dat[, icol[1]] > 0
  } else if(study=="CLY+18") {
    # 20190324 HCT116 cells, Chen et al., 2018
    # CLY+18_proteome, CLY+18_secretome
    dat <- read.csv(paste0(datadir, "CLY+18.csv.xz"), as.is=TRUE)
    description <- paste("HCT116 colon cancer cells", stage)
    # use selected dataset
    icol <- grep(stage, colnames(dat))
    dat <- dat[!is.na(dat[, icol]), ]
    pcomp <- protcomp(dat$Accession)
    up2 <- dat[, icol] > 0
  } else if(study=="KAN+19") {
    # 20191226 cancer-associated fibroblasts, Kugeratski et al., 2019
    # KAN+19_proteome, KAN+19_secretome
    dat <- read.csv(paste0(datadir, "KAN+19.csv.xz"), as.is=TRUE)
    description <- paste("cancer-associated fibroblasts", stage)
    if(stage=="proteome") icol <- "Proteome"
    if(stage=="secretome") {
      # for secretome, combine data for Soluble Secretome and EVs
      dat <- cbind(dat, secretome = NA)
      dat <- dat[!is.na(dat$SolubleSecretome) | !is.na(dat$EVs), ]
      # remove proteins with opposite change in Soluble Secretome and EVs
      iambi <- dat$SolubleSecretome != dat$EVs
      iambi[is.na(iambi)] <- FALSE
      dat <- dat[!iambi, ]
      dat$secretome <- dat$SolubleSecretome
      dat$secretome[is.na(dat$secretome)] <- dat$EVs[is.na(dat$secretome)]
      icol <- "secretome"
    }
    dat <- dat[!is.na(dat[, icol]), ]
    dat <- check_IDs(dat, "Protein.IDs..UniProt.")
    up2 <- dat[, icol] == "Up"
    dat <- cleanup(dat, "Protein.IDs..UniProt.", up2)
    pcomp <- protcomp(dat$Protein.IDs..UniProt.)
  } else stop(paste("multi dataset", dataset, "not available"))
  print(paste0(".pdat_multi: ", description, " [", dataset, "]"))
  # use the up2 from the cleaned-up data, if it exists 20191120
  if("up2" %in% colnames(dat)) up2 <- dat$up2
  return(list(dataset=dataset, pcomp=pcomp, up2=up2, description=description))
}
