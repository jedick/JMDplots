# orp16S/pipeline.R
# Merge and filter 16S sequences from SRA, subsample, find chimeras, and run RDP Classifier

# 20200909 Initial version by Jeffrey M. Dick
# 20210922 Add settings for orp16S datasets
#          Add --skip-technical --clip options to fastq-dump
#          Don't use fastq_filter --stripleft for 454 data
# 20220512 Use GNU parallel in classify()
# 20220521 Replace > in "awk '/>" in subsample() and findchimeras()

## SYSTEM REQUIREMENTS
# (version numbers for information only)

# R 4.0.0              https://www.r-project.org/
# fastq-dump 2.9.0     https://ncbi.github.io/sra-tools/
# vsearch 2.15.0       https://github.com/torognes/vsearch
# seqtk 1.3-r115       https://github.com/lh3/seqtk
# RDP Classifier 2.13  https://sourceforge.net/projects/rdp-classifier/
# Java openjdk version 1.8.0_252 (for RDP); GNU utils (cat, awk, grep)

## USAGE

# Source THIS file. NOTE: first change 'study' (defined below),
# which affects pipeline settings and location of processed files
# > source("pipeline.R") 

# Make sure the working directory is clean (especially no *.fa files)
# To filter one file (generates .fa files in workdir for subsampling)
# > filter("SRR2059382")

# To filter multiple files:
# > RUNID <- c("SRR2059381", "SRR2059380", "SRR2059379")
# > lapply(RUNID, filter)

# To subsample and remove chimeras from the filtered files:
#  - First make sure directories FASTAdir and RDPdir (defined below) exist
# > subsample()
# > findchimeras()

# To run RDP Classifier after chimera removal:
# > classify(RUNID)
# To combine output of RDP Classifier into one file
#   - Result is saved in RDP/<study>.tab below the current directory
# > mkRDP(study, RUNID)

## STUDY SETTINGS

# Change the following line to setup the pipeline for one study
study <- "MTC21"
# Settings for all studies are stored here
file <- tempfile()
# Write spaces here (but don't save them) to make this easier to read
writeLines(con = file, text = gsub(" ", "", c(
  # For 454 experiments, set forwardonly to NA
  "study, forwardonly, trunclen",

  ## Datasets processed for geo16S paper and also used in orp16S paper
  "BCA+21, FALSE, 450",  # PRJNA395513  10.1111/1462-2920.14909
  "HXZ+20, FALSE, 440",  # PRJNA503500  10.1038/s41598-020-62411-2

  ## For orp16S paper 20210922
  "MLL+19, TRUE, 250",
  "RMB+17, FALSE, 250",
  "NTB+21, TRUE, 150",
  "SBP+20, FALSE, 250",
  "MWY+21, FALSE, 250",
  "SAR+13, NA, NA",
  "CTS+17, FALSE, 400",
  "HSF+19, FALSE, 440",
  "HDZ+19, FALSE, 420",
  "ZHZ+19, FALSE, 300",
  "YHK+20, FALSE, 400",
  "SRM+19, FALSE, 250",
  "HLZ+18, FALSE, 420",
  "PSG+20, FALSE, 250",
  "KSR+21, FALSE, 440",
  "ZCZ+21, FALSE, 450",
  "ZZL+21, FALSE, 450",
  "PBU+20, FALSE, 400",
  "MLL+18, TRUE, 250",
  "BWD+19, FALSE, 400",
  "WHL+21, FALSE, 400",
  "ZML+17, TRUE, 400",
  "RBW+14, FALSE, 250", # Winogradsky column
  "DTJ+20, FALSE, 420",
  "WFB+21, FALSE, 440",
  "KLM+16, NA, NA",
  "LMBA21, TRUE, 280",
  "BSPD17, FALSE, 400",
  "CWC+20, FALSE, 440",
  "BMOB18, TRUE, 350",
  "LJC+20, FALSE, 420",
  "DLS21, TRUE, 300",
  "ZZLL21, FALSE, 450",
  "GWS+20, FALSE, 300",
  "CLS+19, TRUE, 150",
  "GZL21, FALSE, 420",
  "APV+20, FALSE, 250",
  "NLE+21, FALSE, 440",
  "LLC+19, TRUE, 400",
  "WHLH21a, TRUE, 440",
  "PCL+18, TRUE, 250",

  # More orp16S 20220504
  "SPA+21, FALSE, 440",
  "GSBT20, FALSE, 400",
  "IBK+22, TRUE, 440",
  "HSF+22, FALSE, 440",
  "HCW+22, TRUE, 240",
  "WKP+22, FALSE, 450",
  "CKB+22, FALSE, 290",
  "WHLH21, FALSE, 450",
  "PSB+21, FALSE, 250",
  "RKSK22, FALSE, 290",
  "DJK+18, TRUE, 250",
  "WLJ+16, TRUE, 200",

  # 20220528-20220610
  "OHL+18, TRUE, 250",
  "MGW+22, FALSE, 290",
  "RKN+17, TRUE, 300",
  "ZLH+22, FALSE, 290",
  "LWJ+21, FALSE, 250",
  "ZDW+19, TRUE, 400",

  # 20220622-20220705
  "WKG+22, FALSE, 440",
  "CLZ+22, TRUE, 290",
  "RARG22, FALSE, 250",
  # 20220917
  "MCR+22, FALSE, 290",
  # 20221020
  "MTC21, FALSE, 250"

  # For RSS+18 (Lake Hazen):
  #  - FASTQ files were downloaded from NCBI cloud
  #  - FASTA files for each site were created with split_libraries.py (from QIIME) and seqtk subseq 20220511
  #  - Start processing with subsample()

)))

# This reads and applies the settings
settings <- read.csv(file, as.is = TRUE)
istudy <- match(study, settings$study)
# Flag for 454 studies
is454 <- is.na(settings$forwardonly[istudy])
# Use forward reads only (for low-quality or missing reverse reads)
forwardonly <- settings$forwardonly[istudy]
# Value of fastq_trunclen
trunclen <- settings$trunclen[istudy]

## FILE LOCATIONS

# Working directory
workdir <- "~/tmp/chem16S"
# Output directory for FASTA files
FASTAdir <- file.path("/home/chem16S", study, "fasta")
# Reference database for chimera detection
refdb <- "/home/chem16S/SILVA/SILVA_138.1_SSURef_NR99_tax_silva.fasta.gz"
#refdb <- "/home/chem16S/chimera/gold.fa"
# RDP jar file
RDPjar <- "/opt/rdp_classifier_2.13/dist/classifier.jar"
# Output directory for RDP Classifier
RDPdir <- file.path("/home/chem16S", study, "RDP")
# RDP Classifier version (java -jar /opt/rdp_classifier_2.13/dist/classifier.jar version)
#Gene:16srrna    Trainset No:18  Taxonomy Version:RDP 16S rRNA training setNo 18 07/2020
#RDP Classifier Version:RDP Naive Bayesian rRNA Classifier Version 2.11, September 2015

filter <- function(RUNID) {
  message("============================================")
  message(paste0("Merging and filtering ", RUNID, " [", study, "]"))
  message("============================================")

  # Change to working directory
  olddir <- setwd(workdir)
  on.exit(setwd(olddir))
  # The output file from this function is a FASTA file with .fa suffix
  outfile <- paste0(RUNID, ".fa")

  if(study %in% c("DJK+18", "WLJ+16")) {
    # For DJK+18, FASTQ files were downloaded from SRA cloud and split with seqtk subseq 20220521
    # For WLJ+16, FASTQ files were downloaded from SRA cloud 20220521
    fqdump <- FALSE
  } else {
    # Generate input FASTQ files with fastq-dump
    fqdump <- TRUE
    cmd <- paste("fastq-dump --split-files --skip-technical --clip", RUNID)
    print(cmd)
    system(cmd)
  }

  if(study %in% c("XXXXXX")) {
    # These are files from MG-RAST (.fasta suffix)
    # Prefix MG-RAST ID to header so we can extract the reads in the subsample() and findchimeras() steps 20211004
    fastafile <- paste0(RUNID, ".fasta")
    lines <- readLines(fastafile)
    ihead <- grep("^>", lines)
    lines[ihead] <- paste0(">", RUNID, ";", gsub("^>", "", lines[ihead]))
    writeLines(lines, outfile)
  } else if(is454) {
    # For 454 studies, fastq-dump --skip-technical puts biological sequences into file with highest-numbered suffix 
    # Use _3.fastq, _2.fastq or _1.fastq if they are present 20200921
    infile <- paste0(RUNID, "_3.fastq")
    if(!file.exists(infile)) infile <- paste0(RUNID, "_2.fastq")
    if(!file.exists(infile)) infile <- paste0(RUNID, "_1.fastq")
    # This is needed for files not generated by fastq-dump
    if(!file.exists(infile)) infile <- paste0(RUNID, ".fastq")
    truncqual <- 15
    if(study %in% c("SAR+13")) truncqual <- 11
    cmd <- paste("vsearch -fastq_filter", infile, "-fastq_minlen 300 -fastq_maxlen 500 -fastq_truncqual", truncqual, "-fastaout", outfile)
    print(cmd)
    system(cmd)
  } else if(forwardonly) {
    # fastq-dump adds _3 suffix for some runs in RKN+17  20220607
    infile <- paste0(RUNID, "_3.fastq")
    if(!file.exists(infile)) infile <- paste0(RUNID, "_2.fastq")
    if(!file.exists(infile)) infile <- paste0(RUNID, "_1.fastq")
    # Use forward reads only
    file.copy(infile, "merged.fastq", overwrite = TRUE)
    # For Ohio Aquifers, file names are different 20220521
    if(study == "DJK+18") file.copy(paste0(RUNID, ".fastq"), "merged.fastq", overwrite = TRUE)
    # For Hetao Plain, file names are different 20220521
    if(study == "WLJ+16") file.copy(paste0(RUNID, "_R1.fastq"), "merged.fastq", overwrite = TRUE)
    nseq <- length(readLines("merged.fastq")) / 4
    print(paste0("Using forward reads only (", nseq, " sequences)"))
  } else {
    # Get file names for forward and reverse reads
    FWDfile <- paste0(RUNID, "_1.fastq")
    REVfile <- paste0(RUNID, "_2.fastq")
    # Check that we have same number of forward and reverse reads 20200917
    # (issue with some datasets: fastq_dump rejects some reads because READLEN < 1)
    ID1 <- sapply(strsplit(gsub("@", "", system(paste("awk '(NR - 1) % 4 == 0'", FWDfile), intern = TRUE)), " "), "[", 1)
    ID2 <- sapply(strsplit(gsub("@", "", system(paste("awk '(NR - 1) % 4 == 0'", REVfile), intern = TRUE)), " "), "[", 1)
    # Forward reads that have matching reverse reads
    forward_matches <- ID1[ID1 %in% ID2]
    if(length(forward_matches) < length(ID1)) {
      writeLines(forward_matches, "forward_matches.txt")
      print(cmd <- paste("seqtk subseq", FWDfile, "forward_matches.txt > new_1.fastq"))
      system(cmd)
      print(cmd <- paste("mv new_1.fastq", FWDfile))
      system(cmd)
    }
    # Reverse reads that have matching forward reads
    reverse_matches <- ID2[ID2 %in% ID1]
    if(length(reverse_matches) < length(ID2)) {
      writeLines(reverse_matches, "reverse_matches.txt")
      print(cmd <- paste("seqtk subseq", REVfile, "reverse_matches.txt > new_2.fastq"))
      system(cmd)
      print(cmd <- paste("mv new_2.fastq", REVfile))
      system(cmd)
    }
    # Make sure forward and reverse reads have same IDs
    ID1 <- sapply(strsplit(gsub("@", "", system(paste("awk '(NR - 1) % 4 == 0'", FWDfile), intern = TRUE)), " "), "[", 1)
    ID2 <- sapply(strsplit(gsub("@", "", system(paste("awk '(NR - 1) % 4 == 0'", REVfile), intern = TRUE)), " "), "[", 1)
    if(!all(ID1 == ID2)) stop(paste("sequence IDs in", FWDfile, "and", REVfile, "don't match!"))
    # Merge paired reads
    print(cmd <- paste("vsearch -fastq_mergepairs", FWDfile, "-reverse", REVfile, "-fastqout merged.fastq"))
    system(cmd)
  }
  if(!is454) {
    # Remove short and high-error reads
    # Use fastq_maxee_rate instead of fastq_maxee 20200920
    maxee_rate <- 0.005
    # Use default -fastq_qmax except for some studies 20210922
    qmax <- 41
    # For Ion Torrent 20210928
    if(study %in% c("PCL+18", "RKN+17")) qmax <- 45
    cmd <- paste("vsearch -fastq_filter merged.fastq -fastq_trunclen", trunclen, "-fastq_qmax", qmax, "-fastq_maxee_rate", maxee_rate, "-fastaout", outfile)
    print(cmd)
    system(cmd)
  }

  # Clean up
  if(fqdump) file.remove(Sys.glob("*.fastq"))
  return()
}

subsample <- function(do.singletons = TRUE) {
  message("=======================")
  message(paste0("Subsampling [", study, "]"))
  message("=======================")

  # Change to working directory
  olddir <- setwd(workdir)
  on.exit(setwd(olddir))

  if(do.singletons) {
    # Make combined FASTA file of filtered sequences 20200914
    print(cmd <- "cat *.fa > filtered.fasta")
    system(cmd)
    # Remove singletons 20200912
    # 1. Dereplicate sequences, keeping those that appear exactly once
    print(cmd <- "vsearch -derep_fulllength filtered.fasta -maxuniquesize 1 -output singletons.fasta")
    system(cmd)
    # 2. Get all sequence IDs
    print(cmd <- 'grep ">" filtered.fasta | sed -e s/\\ .*//g | sed -e s/\\>//g > allIDs.txt')
    system(cmd)
    allIDs <- readLines("allIDs.txt")
    # 3. Get singleton sequence IDs
    print(cmd <- 'grep ">" singletons.fasta | sed -e s/\\ .*//g | sed -e s/\\>//g > singletonIDs.txt')
    system(cmd)
    singletonIDs <- readLines("singletonIDs.txt")
    # 4. Write IDs of sequences that are not singletons
    writeLines(setdiff(allIDs, singletonIDs), "not_singletons.txt")
    # 5. Extract sequences that are not singletons
    print(cmd <- "seqtk subseq filtered.fasta not_singletons.txt > not_singletons.fasta")
    system(cmd)
    singletxt <- paste0("Removed ", length(singletonIDs), " singletons from ", length(allIDs), " sequences (", round(length(singletonIDs) / length(allIDs) * 100, 1), "%)")
    print(singletxt)
  } else {
    # Don't remove singletons, but combine filtered sequences for subsampling 20210821
    print(cmd <- "cat *.fa > not_singletons.fasta")
    system(cmd)
  }

  # Loop over samples
  allRUNID <- gsub("\\.fa$", "", dir(pattern = "\\.fa$"))
  for(thisRUNID in allRUNID) {
    # Extract sequences for this sample
    thisfile <- paste0(thisRUNID, ".fasta")
    # https://stackoverflow.com/questions/26144692/printing-a-sequence-from-a-fasta-file
    # "awk '/>" is used to match sample name at beginning of header line
    # (this is needed for numeric or short sample names that can match other parts of the header)
    print(cmd <- paste0("awk '/>", thisRUNID, "/{p++;print;next} /^>/{p=0} p' not_singletons.fasta > ", thisfile))
    system(cmd)

    # Subsample 10000 sequences - put output in .fa file in FASTAdir
    print(cmd <- paste('grep "^>"', thisfile))
    thisfile.headers <- system(cmd, intern = TRUE)
    nout <- nseq <- length(thisfile.headers)
    outfile <- file.path(FASTAdir, paste0(thisRUNID, ".fa"))
    if(nseq > 10000) {
      print(cmd <- paste("vsearch -fastx_subsample", thisfile, "-sample_size 10000 -randseed 1234 -fastaout", outfile))
      system(cmd)
      nout <- 10000
    } else {
      file.copy(thisfile, outfile, overwrite = TRUE)
    }
    subsamptxt <- paste("Subsampled", nout, "of", nseq, "sequences")
    print(subsamptxt)
  }

  # Clean up
  file.remove(Sys.glob("*.fa"))
  file.remove(Sys.glob("*.fasta"))
  return()
}

# Find chimeras in combined FASTA file, then split result into samples 20200913
# (Saves time by masking and creating k-mer index only once)
findchimeras <- function(threads = 8) {
  message("=========================")
  message(paste0("Finding chimeras [", study, "]"))
  message("=========================")

  # Change to output directory
  olddir <- setwd(FASTAdir)
  on.exit(setwd(olddir))
  # Make combined FASTA file
  print(cmd <- "cat *.fa > combined.fasta")
  system(cmd)

  # Identify chimeras
  print(cmd <- paste("vsearch -threads", threads, "-uchime_ref combined.fasta -nonchimeras nonchimeras.fasta -db", refdb))
  system(cmd)

  # Split result into one file for each sample
  # https://stackoverflow.com/questions/26144692/printing-a-sequence-from-a-fasta-file
  allRUNID <- gsub("\\.fa$", "", dir(pattern = "\\.fa$"))
  for(thisRUNID in allRUNID) {
    print(cmd <- paste0("awk '/>", thisRUNID, "/{p++;print;next} /^>/{p=0} p' nonchimeras.fasta > ", thisRUNID, ".fasta"))
    system(cmd)
  }

  # Clean up
  file.remove("combined.fasta")
  file.remove("nonchimeras.fasta")
}

# Run RDP Classifier on one or multiple files 20200910
# Use GNU parallel 20220509
classify <- function(RUNID, conf = 0.8, jobs = 4) {

  # Make sure RDP output directory exists 20200915
  if(!dir.exists(RDPdir)) stop(paste("directory", RDPdir, "does not exist"))
  # Force evaluation of RUNID before changing the directory
  allRUNID <- RUNID

  # Change to RDP directory
  olddir <- setwd(RDPdir)
  on.exit(setwd(olddir))

  # Loop over RUNIDs
  n <- length(allRUNID)
  for(i in 1:n) {
    thisRUNID <- allRUNID[i]
    # Path to the FASTA file
    FASTAfile <- file.path("../fasta", paste0(thisRUNID, ".fasta"))
    # Copy FASTA file to here without .fasta suffix (so it doesn't get into sample name)
    file.copy(FASTAfile, thisRUNID, overwrite = TRUE)
  }

  # Write input filenames to file
  writeLines(allRUNID, "fastafiles")
  # Write output (.txt and .tab) filenames to files
  txtfiles <- paste0(RUNID, ".txt")
  writeLines(txtfiles, "txtfiles")
  tabfiles <- paste0(RUNID, ".tab")
  writeLines(tabfiles, "tabfiles")

  cmd <- paste("parallel --link --jobs", jobs, "java -jar /opt/rdp_classifier_2.13/dist/classifier.jar classify -c 0.8 -h {3} -o {2} {1} :::: fastafiles tabfiles txtfiles")
  message("=======================================")
  message(paste0("Classifying ", n, " runs with ", jobs, " jobs [", study, "]"))
  message("=======================================")

  system(cmd)

  # Clean up
  file.remove(c("fastafiles", "tabfiles", "txtfiles"))
  file.remove(Sys.glob(paste0("cnadjusted_*")))
  file.remove(Sys.glob(paste0("*.tab")))
  for(i in 1:n) {
    thisRUNID <- allRUNID[i]
    file.remove(thisRUNID)
  }

}

# Combine output files of RDP Classifier into a single CSV file for one study 20200903
# Use RDP Classifier merge-count command 20200910
# Get file list from RDPdir 20210817
mkRDP <- function() {
  # Change to temporary directory
  olddir <- setwd(tempdir())
  on.exit(setwd(olddir))
  # Remove any existing .txt files
  file.remove(Sys.glob("*.txt"))

  # Loop over runs
  runs <- gsub(".txt", "", dir(RDPdir))
  for(run in runs) {
    # The RDP result file
    RDPfile <- file.path(RDPdir, paste0(run, ".txt"))
    print(RDPfile)
    # Copy the RDP result file to here (temporary directory)
    file.copy(RDPfile, ".")
  }

  # Now list all the RDP result files
  files <- dir(pattern = "txt")
  # Name of output file
  outfile <- file.path(olddir, "RDP", paste0(study, ".tab"))
  # Remove output file in case it exists
  if(file.exists(outfile)) file.remove(outfile)
  # RDP jar file
  RDPjar <- "/opt/rdp_classifier_2.13/dist/classifier.jar"
  # Run merge-count command
  print(cmd <- paste("java -jar", RDPjar, "merge-count", outfile, paste(files, collapse = " ")))
  system(cmd)
}
