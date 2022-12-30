# Process Saanich Inlet metaproteome
# 20220501 jmd version 1
# 20221028 revised for sars16S paper
# 20221227 Use ScanCount for protein abundance

# Get amino acid composition for each protein sequence
# https://ftp.ebi.ac.uk/pride-archive/2017/10/PXD004433/SaanichInlet_LP_ORFs_2015-01-13_Filtered.fasta
refaa <- CHNOSZ::read.fasta("SaanichInlet_LP_ORFs_2015-01-13_Filtered.fasta.xz")
# Remove trailing \r of some IDs
refaa$protein <- gsub("\r", "", refaa$protein)

# Read protein table
# https://ftp.ebi.ac.uk/pride-archive/2017/10/PXD004433/SBI_Metagenome2015_AllProteinsAllExperiments.txt
dat <- read.csv("SBI_Metagenome2015_AllProteinsAllExperiments.txt.xz", sep = "\t")

# Read sample information from Hawley et al. (2017)
# Selected values from Table S2 of https://doi.org/10.1038/sdata.2017.160
samp <- read.csv("HTZ+17_TableS2_selected.csv")

# Sum amino acid compositions for each sample
out <- lapply(1:nrow(samp), function(i) {
  Experiment <- samp$MetaP.Pride.File.Prefix[i]
  print(Experiment)
  idat <- dat$Experiment == Experiment
  thisdat <- dat[idat, ]
  irefaa <- match(thisdat$Reference, refaa$protein)
  aa <- refaa[irefaa, ]
  # Multiply amino acid composition by ScanCount
  aa[, 5:25] <- aa[, 5:25] * thisdat$ScanCount
  # Sum amino acid composition
  CHNOSZ::aasum(aa)
})

out <- do.call(rbind, out)
out$protein <- "metaproteome"
out$organism <- samp$MetaP.Pride.File.Prefix
out$ref <- samp$Cruise.ID
out$abbrv <- samp$Depth..m.
write.csv(out, "HTZ+17_aa.csv", row.names = FALSE, quote = FALSE)
