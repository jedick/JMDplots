# evdevH2O/Main/mkmetrics.R
# Calculate ZC and nH2O and get gene age for each protein in the consensus tables of Liebeskind et al. (2016)
# 20211103 jmd

# Input file (provided here):
# reference_proteomes.csv: UniProt reference proteomes for organisms in Liebeskind et al. (2016)
#   This file was extracted from
#   https://ftp.uniprot.org/pub/databases/uniprot/current%5Frelease/knowledgebase/reference%5Fproteomes/README

# Output files:
# modeAges.csv: names of mode ages for each organism
# In metrics/ directory:
# <OSCODE>.csv: ID, modeAge, ZC, nH2O for each protein (one file for each organism, 31 total)

# NOTE: change to make output <OSCODE>.csv files smaller 20211220:
# ZC and nH2O values are rounded at 3 instead of 6 decimal places
# csv files are xz compressed after running this script

# Other input files (must be downloaded separately):

# In Main/ directory:
# Consensus gene ages tables by Liebskind et al. (2016), available at https://doi.org/10.5281/zenodo.51708
# main_ANOGA.csv  main_CANAL.csv  main_CRYNJ.csv  main_IXOSC.csv  main_MOUSE.csv  main_PANTR.csv  main_SCHPO.csv  main_XENTR.csv
# main_BOVIN.csv  main_CANFA.csv  main_DANRE.csv  main_MACMU.csv  main_NEMVE.csv  main_PHANO.csv  main_SCLS1.csv  main_YARLI.csv
# main_BRAFL.csv  main_CHICK.csv  main_DROME.csv  main_MONBE.csv  main_NEUCR.csv  main_RAT.csv    main_TAKRU.csv  main_YEAST.csv
# main_CAEEL.csv  main_CIOIN.csv  main_HUMAN.csv  main_MONDO.csv  main_ORNAN.csv  main_SCHMA.csv  main_USTMA.csv

# NOTE: Rename one file to match current UniProt organism code (OSCODE):
# mv main_CANFA.csv main_CANLF.csv

# In UniProt/ directory: Reference proteome sequence files for each organism
# UP000000437_7955.fasta.gz    UP000001300_284591.fasta.gz  UP000001940_6239.fasta.gz    UP000002485_284812.fasta.gz  UP000008144_7719.fasta.gz
# UP000000539_9031.fasta.gz    UP000001312_665079.fasta.gz  UP000002149_214684.fasta.gz  UP000002494_10116.fasta.gz   UP000008854_6183.fasta.gz
# UP000000559_237561.fasta.gz  UP000001357_81824.fasta.gz   UP000002254_9615.fasta.gz    UP000005226_31033.fasta.gz   UP000009136_9913.fasta.gz
# UP000000561_237631.fasta.gz  UP000001554_7739.fasta.gz    UP000002277_9598.fasta.gz    UP000005640_9606.fasta.gz
# UP000000589_10090.fasta.gz   UP000001555_6945.fasta.gz    UP000002279_9258.fasta.gz    UP000006718_9544.fasta.gz
# UP000000803_7227.fasta.gz    UP000001593_45351.fasta.gz   UP000002280_13616.fasta.gz   UP000007062_7165.fasta.gz
# UP000001055_321614.fasta.gz  UP000001805_367110.fasta.gz  UP000002311_559292.fasta.gz  UP000008143_8364.fasta.gz

# This script was used to download each proteome:
#dat <- read.csv("reference_proteomes.csv")
#for(i in 1:nrow(dat)) {
#  ID <- dat$Proteome_ID[i]
#  Tax <- dat$Tax_ID[i]
#  URLbase <- "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/"
#  URL <- paste0(URLbase, ID, "/", ID, "_", Tax, ".fasta.gz")
#  cmd <- paste("wget", URL)
#  print(cmd)
#  system(cmd)
#}

# Read data frame with list of UniProt reference proteomes
refprot <- read.csv("reference_proteomes.csv")
# Initialize data frame to hold names of modeAges for different organisms
modeAges <- as.data.frame(t(rep(NA, 9)))[1:nrow(refprot), ]
# Loop over reference proteomes
for(i in 1:nrow(refprot)) {
  # Read amino acid composition from FASTA file
  file <- file.path("UniProt", paste0(refprot$Proteome_ID[i], "_", refprot$Tax_ID[i], ".fasta.gz"))
  aa <- CHNOSZ::read.fasta(file)
  # Get Protein IDs
  aa$protein <- sapply(strsplit(aa$protein, "\\|"), "[", 2)
  # Calculate ZC and nH2O
  ZC <- ZC(protein.formula(aa))
  nH2O <- as.numeric(H2OAA(aa))
  # Read gene age table
  agefile <- file.path("Main", paste0("main_", refprot$OSCODE[i], ".csv"))
  ages <- read.csv(agefile, check.names = FALSE, row.names = 1)
  # Match modeAge to column names to get numeric value
  modeAge <- match(ages$modeAge, colnames(ages))
  # Match protein ID to reference sequence data
  iaa <- match(rownames(ages), aa$protein)
  # Build output data frame
  out <- data.frame(ID = rownames(ages), modeAge = modeAge, ZC = round(ZC[iaa], 3), nH2O = round(nH2O[iaa], 3))
  # Exclude NA values of ZC (i.e. Liebeskind protein ID not in reference proteome)
  out <- out[!is.na(out$ZC), ]
  outfile <- file.path("metrics", paste0(refprot$OSCODE[i], ".csv"))
  write.csv(out, outfile, row.names = FALSE, quote = FALSE)
  # Add row to modeAges data frame
  nAges <- which(colnames(ages) == "modeAge") - 1
  modeAges[i, 1:nAges] <- colnames(ages)[1:nAges]
}
# Write modeAges data frame
colnames(modeAges) <- 1:9
modeAges <- data.frame(OSCODE = refprot$OSCODE, modeAges, check.names = FALSE)
write.csv(modeAges, "modeAges.csv", row.names = FALSE, quote = FALSE)
