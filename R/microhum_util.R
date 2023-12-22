# JMDplots/microhum_util.R
# Supporting functions for microhum paper

# Get metadata for a study, appending columns for pch and col 20200914
getmdat_microhum <- function(study, metrics = NULL, dropNA = TRUE, quiet = TRUE) {
  # Read metadata file
  # Remove suffix after underscore 20200929
  studyfile <- gsub("_.*", "", study)
  file <- file.path(getdatadir(), "16S/metadata", paste0(studyfile, ".csv"))
  metadata <- read.csv(file, as.is = TRUE, check.names = FALSE)

  if(dropNA) {
    # Exclude samples with NA name (e.g. outliers?) 20200916
    iname <- match("name", tolower(colnames(metadata)))
    noname <- is.na(metadata[, iname])
    if(any(noname)) {
      if(!quiet) print(paste0("getmdat_microhum [", study, "]: dropping ", sum(noname), " samples with NA name"))
      metadata <- metadata[!is.na(metadata[, iname]), ]
    }
  }
  # Use NULL pch as flag for unavailable dataset 20210820
  pch <- NULL

  # For microhum paper 20220202

  # 20210801 COVID-19 Gut
  if(study == "MMP+21") {
    pch <- sapply(metadata$Severity, switch, Mild = 24, Moderate = 20, Severe = 25)
    col <- sapply(metadata$Severity, switch, Mild = 4, Moderate = 1, Severe = 2)
  }
  # 20220803 COVID-19 Gut
  if(study == "CGC+22") {
    pch <- sapply(metadata$Type, switch, Control = 24, Convalescence = 20, Acute = 25, NA)
    col <- sapply(metadata$Type, switch, Control = 4, Convalescence = 1, Acute = 2, NA)
  }
  # 20220803 COVID-19 Nasopharyngeal
  if(study == "ENJ+21") {
    pch <- sapply(metadata$Description, switch, "COVID-19-NEGATIVE" = 24, "Blank-Control" = 20, "COVID-19-POSITIVE" = 25, NA)
    col <- sapply(metadata$Description, switch, "COVID-19-NEGATIVE" = 4, "Blank-Control" = 1, "COVID-19-POSITIVE" = 2, NA)
  }
  # 20220803 COVID-19 Oral
  if(study == "GBS+22") {
#    pch <- sapply(metadata$Sex, switch, "female" = 24, "male" = 25)
#    col <- sapply(metadata$Sex, switch, "female" = 4, "male" = 2)
    pch <- sapply(metadata$Status, switch, "non-infected" = 24, "infected" = 25)
    col <- sapply(metadata$Status, switch, "non-infected" = 4, "infected" = 2)
  }
  # 20220804 COVID-19 Gut
  if(study == "GCW+20") {
    pch <- sapply(metadata$Status, switch, Control = 24, H1N1 = 20, "COVID-19" = 25)
    col <- sapply(metadata$Status, switch, Control = 4, H1N1 = 1, "COVID-19" = 2)
  }
  # 20220804 COVID-19 Nasopharyngeal
  if(study == "GKJ+22") {
    pch <- sapply(metadata$Status, switch, Control = 24, Covid = 25)
    col <- sapply(metadata$Status, switch, Control = 4, Covid = 2)
  }
  # 20220805 COVID-19 Respiratory samples (oropharyngeal swabs, nasopharyngeal swabs, and tracheal aspirates)
  if(study == "HMH+21") {
    # Include only severe and fatal cases 20220822
    metadata <- metadata[metadata$Status != "mild", ]
    pch <- sapply(metadata$Status, switch, healthy = 24, mild = 25, severe = 25, fatal = 25, NA)
    col <- sapply(metadata$Status, switch, healthy = 4, mild = 2, severe = 2, fatal = 2, NA)
  }
  # 20220805 COVID-19 Oral
  if(study == "MAC+21") {
#    # Remove samples with discord between lab and qPCR results
#    metadata <- metadata[!sapply(metadata$Concord_Discord == "DISCORD", isTRUE), ]
    # covid_status_per_lab is in [0, 1], so add 1 to get [1, 2],
    # which correspond to the unnamed arguments for switch
    pch <- sapply(metadata$covid_status_per_lab + 1, switch, 24, 25)
    col <- sapply(metadata$covid_status_per_lab + 1, switch, 4, 2)
  }
  # 20220805 COVID-19 Respiratory samples (oropharyngeal swabs, nasopharyngeal swabs, and endotracheal aspirates)
  # MLW+21_Oropharyngeal, MLW+21_Nasopharyngeal, MLW+21_Endotracheal
  if(grepl("MLW\\+21", study)) {
    SampleType <- sapply(strsplit(study, "_"), "[", 2)
    if(!is.na(SampleType)) metadata <- metadata[grep(SampleType, metadata$SampleType), ]
    pch <- sapply(metadata$Study_group, switch, "non COVID" = 24, "COVID" = 25, NA)
    col <- sapply(metadata$Study_group, switch, "non COVID" = 4, "COVID" = 2, NA)
  }
  # 20220807 COVID-19 Gut
  if(study == "NGH+21") {
    pch <- sapply(metadata$Status, switch, Control = 24, COVIDrecovered = 20, COVID = 25)
    col <- sapply(metadata$Status, switch, Control = 4, COVIDrecovered = 1, COVID = 2)
  }
  # 20220807 COVID-19 Nasopharyngeal
  if(study == "PMM+22") {
    pch <- sapply(metadata$Status, switch, Control = 24, Asymptomatic = 20, Symptomatic = 25)
    col <- sapply(metadata$Status, switch, Control = 4, Asymptomatic = 1, Symptomatic = 2)
  }
  # 20220807 COVID-19 Gut
  if(study == "RDM+22") {
    pch <- sapply(metadata$Status, switch, other = 24, "MIS-C" = 20, "COVID-19" = 25, NA)
    col <- sapply(metadata$Status, switch, other = 4, "MIS-C" = 1, "COVID-19" = 2, NA)
  }
  # 20220807 COVID-19 Oral and Gut
  # RFH+22_Oral, RFH+22_Gut
  if(grepl("RFH\\+22", study)) {
    Type <- sapply(strsplit(study, "_"), "[", 2)
    if(!is.na(Type)) metadata <- metadata[metadata$Type == Type, ]
    pch <- sapply(metadata$Status, switch, Control = 24, COVID = 25)
    col <- sapply(metadata$Status, switch, Control = 4, COVID = 2)
  }
  # 20220807 COVID-19 Oral and Gut
  # RWC+21_Oral, RWC+21_Gut
  if(grepl("RWC\\+21", study)) {
    Type <- sapply(strsplit(study, "_"), "[", 2)
    if(!is.na(Type)) {
      if(Type == "Oral") metadata <- metadata[metadata$Source == "tongue coating", ]
      if(Type == "Gut") metadata <- metadata[metadata$Source == "feces", ]
    }
    pch <- sapply(metadata$Status, switch, Control = 24, Patient = 25)
    col <- sapply(metadata$Status, switch, Control = 4, Patient = 2)
  }
  # 20220808 COVID-19 Gut
  if(study == "SRK+22") {
    pch <- sapply(metadata$Status, switch, Control = 24, Covid19 = 25)
    col <- sapply(metadata$Status, switch, Control = 4, Covid19 = 2)
  }
  # 20220808 COVID-19 Respiratory
  if(study == "SRS+22") {
    pch <- sapply(metadata$Status, switch, Control = 24, "Negative Control" = 20, COVID = 25)
    col <- sapply(metadata$Status, switch, Control = 4, "Negative Control" = 1, COVID = 2)
  }
  # 20220810 COVID-19 Nasopharyngeal
  # VCV+21_1, VCV+21_2, VCV+21_3
  if(grepl("VCV\\+21", study)) {
    Severity <- sapply(strsplit(study, "_"), "[", 2)
    if(!is.na(Severity)) metadata <- metadata[metadata$Severity %in% c(Severity, 0), ]
    pch <- sapply(as.character(metadata$Severity), switch, "0" = 24, "1" = 25, "2" = 25, "3" = 25)
    col <- sapply(as.character(metadata$Severity), switch, "0" = 4, "1" = 2, "2" = 2, "3" = 2)
  }
  # 20220810 COVID-19 Oral and Gut
  # WCJ+21_Oral, WCJ+21_Gut
  if(grepl("WCJ\\+21", study)) {
    Type <- sapply(strsplit(study, "_"), "[", 2)
    if(!is.na(Type)) {
      if(Type == "Oral") metadata <- metadata[metadata$Source == "throat swab", ]
      if(Type == "Gut") metadata <- metadata[metadata$Source == "fecal sample", ]
    }
    pch <- sapply(metadata$Severity2, switch, Healthy = 24, NonSevere = 20, Severe = 25)
    col <- sapply(metadata$Severity2, switch, Healthy = 4, NonSevere = 1, Severe = 2)
  }
  # 20220810 COVID-19 Respiratory
  if(study == "XLZ+21") {
    # Use airway samples only (gut has no healthy controls)
    metadata <- metadata[metadata$Location == "Airway", ]
    pch <- sapply(metadata$Status, switch, Control = 24, Patient = 25)
    col <- sapply(metadata$Status, switch, Control = 4, Patient = 2)
  }
  # 20220810 COVID-19 Gut
  if(study == "ZZZ+21") {
    pch <- sapply(metadata$Status, switch, Control = 24, Patient = 25)
    col <- sapply(metadata$Status, switch, Control = 4, Patient = 2)
  }
  # 20221012 COVID-19 Nasopharyngeal
  if(study == "CSC+22") {
    pch <- sapply(metadata$SarsCov2, switch, neg = 24, pos = 25)
    col <- sapply(metadata$SarsCov2, switch, neg = 4, pos = 2)
  }
  # 20221012 COVID-19 Gut
  if(study == "FBD+22") {
    pch <- sapply(metadata$Disease, switch, Control = 24, "COVID-19" = 25)
    col <- sapply(metadata$Disease, switch, Control = 4, "COVID-19" = 2)
  }
  # 20221012 COVID-19 Oropharyngeal
  if(study == "GWL+21") {
    pch <- sapply(metadata$Status, switch, Control = 24, "COVID-19" = 25)
    col <- sapply(metadata$Status, switch, Control = 4, "COVID-19" = 2)
  }
  # 20221012 COVID-19 Oropharyngeal
  if(study == "IZC+21") {
    pch <- sapply(metadata$Disease, switch, Control = 24, "COVID-19" = 25)
    col <- sapply(metadata$Disease, switch, Control = 4, "COVID-19" = 2)
  }
  # 20221012 COVID-19 Gut
  if(study == "MIK+22") {
    pch <- sapply(metadata$Disease, switch, "HIV-1 negative" = 24, "COVID-19" = 25, NA)
    col <- sapply(metadata$Disease, switch, "HIV-1 negative" = 4, "COVID-19" = 2, NA)
  }
  # 20221103 COVID-19 Gut
  if(study == "KMG+21") {
    pch <- sapply(metadata$Status, switch, Control = 24, "COVID-19" = 25, NA)
    col <- sapply(metadata$Status, switch, Control = 4, "COVID-19" = 2, NA)
  }
  # 20221105 COVID-19 Nasopharyngeal
  if(study == "SGC+21") {
    pch <- sapply(metadata$Status, switch, Control = 24, "COVID-19" = 25, NA)
    col <- sapply(metadata$Status, switch, Control = 4, "COVID-19" = 2, NA)
  }

  # 20230321 Ulcerative colitis
  if(study == "TWC+22") {
    status <- ifelse(grepl("UC", metadata$Sample), "UC", "Healthy")
    pch <- sapply(status, switch, Healthy = 24, UC = 25, NA)
    col <- sapply(status, switch, Healthy = 4, UC = 2, NA)
  }
  # 20230322 IBD
  if(study == "LAA+19") {
    pch <- sapply(metadata$Disease, switch, nonIBD = 24, UC = 25, CD = 25, NA)
    col <- sapply(metadata$Disease, switch, nonIBD = 4, UC = 2, CD = 2, NA)
  }
  # 20230322 IBD
  if(study == "MDV+22") {
    pch <- sapply(metadata$Subject_ID, switch, Healthy_control = 24, 25)
    col <- sapply(metadata$Subject_ID, switch, Healthy_control = 4, 2)
  }
  # 20230323 Ulcerative colitis
  if(study == "ZTG+21") {
    status <- ifelse(grepl("UC", metadata$SampleName), "UC", "Healthy")
    pch <- sapply(status, switch, Healthy = 24, UC = 25, NA)
    col <- sapply(status, switch, Healthy = 4, UC = 2, NA)
  }
  # 20230323 IBD
  if(study == "ASM+23") {
    pch <- sapply(metadata$Status, switch, Control = 24, Patient = 25)
    col <- sapply(metadata$Status, switch, Control = 4, Patient = 2)
  }
  # 20230323 IBD
  if(study == "AAM+20") {
    pch <- sapply(metadata$Disease, switch, Control = 24, "Ulcerative Colitis" = 25, "Crohn's Disease" = 25)
    col <- sapply(metadata$Disease, switch, Control = 4, "Ulcerative Colitis" = 2, "Crohn's Disease" = 2)
  }
  # 20230323 IBD
  if(study == "RAF+20") {
    status <- ifelse(grepl("HT", metadata$SampleName), "Healthy", "IBD")
    pch <- sapply(status, switch, Healthy = 24, IBD = 25, NA)
    col <- sapply(status, switch, Healthy = 4, IBD = 2, NA)
  }
  # 20230324 IBD
  if(study == "LZD+19") {
    # Stool samples from patients and controls
    idat <- metadata$Tissue == "Stool"
    idat[is.na(idat)] <- FALSE
    metadata <- metadata[idat, ]
    pch <- sapply(metadata$Disease, switch, "Healthy subject" = 24, "Inflammatory Bowel Disease" = 25, NA)
    col <- sapply(metadata$Disease, switch, "Healthy subject" = 4, "Inflammatory Bowel Disease" = 2, NA)
  }
  if(study == "LZD+19_mucosa") {
    # Inflamed and not inflamed tissue from patients and controls
    icontrol <- metadata$Disease == "Healthy subject" & metadata$Tissue == "Not Inflamed"
    ipatient <- metadata$Disease == "Inflammatory Bowel Disease" & metadata$Tissue == "Inflamed"
    idat <- icontrol | ipatient
    idat[is.na(idat)] <- FALSE
    metadata <- metadata[idat, ]
    pch <- sapply(metadata$Disease, switch, "Healthy subject" = 24, "Inflammatory Bowel Disease" = 25, NA)
    col <- sapply(metadata$Disease, switch, "Healthy subject" = 4, "Inflammatory Bowel Disease" = 2, NA)
  }
  # 20230325 IBD
  if(study == "WGL+19") {
    pch <- sapply(metadata$Disease, switch, HC = 24, UC = 25, CD = 25)
    col <- sapply(metadata$Disease, switch, HC = 4, UC = 2, CD = 2)
  }
  # 20230327 IBD
  if(study == "GKD+14") {
    metadata <- metadata[metadata$Location == "stool", ]
    pch <- sapply(metadata$Diagnosis, switch, "Not IBD" = 24, CD = 25)
    col <- sapply(metadata$Diagnosis, switch, "Not IBD" = 4, CD = 2)
  }
  if(study == "GKD+14_Rectum") {
    metadata <- metadata[metadata$Location == "Rectum", ]
    pch <- sapply(metadata$Diagnosis, switch, "Not IBD" = 24, CD = 25)
    col <- sapply(metadata$Diagnosis, switch, "Not IBD" = 4, CD = 2)
  }
  if(study == "GKD+14_Ileum") {
    metadata <- metadata[metadata$Location == "Terminal Ileum", ]
    pch <- sapply(metadata$Diagnosis, switch, "Not IBD" = 24, CD = 25)
    col <- sapply(metadata$Diagnosis, switch, "Not IBD" = 4, CD = 2)
  }

  # 20220814 Body sites
  # BPB+21_NoTreatment, BPB+21_AnyTreatment
  if(grepl("BPB\\+21", study)) {
    Treatment <- sapply(strsplit(study, "_"), "[", 2)
    if(!is.na(Treatment)) {
      if(Treatment == "NoTreatment") metadata <- metadata[metadata$Treatment == "No_treatment", ]
      if(Treatment == "AnyTreatment") metadata <- metadata[metadata$Treatment != "No_treatment", ]
    }
    pch <- sapply(metadata$Site, switch, "Oral cavity" = pch_Oral, "Nasal cavity" = pch_Nasal, "Skin of forearm" = pch_Skin, "feces" = pch_Gut, NA)
    col <- sapply(metadata$Site, switch, "Oral cavity" = adjustcolor(col_Oral, 0.8), "Nasal cavity" = adjustcolor(col_Nasal, 0.9),
      "Skin of forearm" = adjustcolor(col_Skin, 0.8), "feces" = adjustcolor(col_Gut, 0.8), NA)
  }

  # 20231219 Human Microbiome Project
  if(study == "HMP12") {
    pch <- sapply(metadata$"Body site", switch, "Skin" = pch_Skin, "Nasal cavity" = pch_Nasal, "Oral cavity" = pch_Oral, "GI tract" = pch_Gut, "UG tract" = pch_UG)
    col <- sapply(metadata$"Body site", switch, "Skin" = adjustcolor(col_Skin, 0.8), "Nasal cavity" = adjustcolor(col_Nasal, 0.8),
      "Oral cavity" = adjustcolor(col_Oral, 0.8), "GI tract" = adjustcolor(col_Gut, 0.8), "UG tract" = adjustcolor(col_UG, 0.8))
  }

  if(is.null(pch)) stop(paste(study, "metadata file exists, but not set up for processing"))

  metadata <- cbind(metadata, pch, col)
  # Return both metadata and metrics, if provided 20220506
  if(is.null(metrics)) metadata else {
    # Keep metadata only for samples with metrics 20201006
    metadata <- metadata[metadata$Run %in% metrics$Run, ]
    # Put metrics in same order as metadata 20220505
    imet <- match(metadata$Run, metrics$Run)
    metrics <- metrics[imet, ]
    # Insert sample column in metrics
    # Use first column name starting with "sample" or "Sample" 20210818
    sampcol <- grep("^sample", colnames(metadata), ignore.case = TRUE)[1]
    metrics <- cbind(data.frame(Run = metrics$Run, sample = metadata[, sampcol]), metrics[, -1, drop = FALSE])
    list(metadata = metadata, metrics = metrics)
  }
}

# Function to calculate metrics for a given study 20220506
getmetrics_microhum <- function(study, lineage = NULL, mincount = 100, return_AA = FALSE, zero_AA = NULL, quiet = TRUE, ...) {
  # Remove suffix after underscore 20200929
  studyfile <- gsub("_.*", "", study)
  RDPfile <- file.path(getdatadir(), "16S/RDP-GTDB", paste0(studyfile, ".tab.xz"))
  # If there is no .xz file, look for a .tab file 20210607
  if(!file.exists(RDPfile)) RDPfile <- file.path(getdatadir(), "16S/RDP-GTDB", paste0(studyfile, ".tab"))
  RDP <- read_RDP(RDPfile, lineage = lineage, mincount = mincount, quiet = quiet, ...)
  map <- map_taxa(RDP, quiet = quiet)
  get_metrics(RDP, map = map, taxon_AA = taxon_AA[["GTDB"]], return_AA = return_AA, zero_AA = zero_AA)
}

# Function to calculate and plot metrics for a given study 20220506
plotmet_microhum <- function(study, lineage = NULL, quiet = TRUE, ...) {
  metrics <- getmetrics_microhum(study, lineage = lineage, quiet = quiet)
  mdat <- getmdat_microhum(study, metrics)
  pm <- plot_metrics(mdat, xvar = "nO2", yvar = "nH2O", ...)
  # Prepend study column
  cbind(study = study, pm)
}
