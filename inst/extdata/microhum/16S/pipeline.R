# microhum/pipeline.R
# Dump 16S rRNA gene sequences from SRA, merge and filter, subsample, find chimeras, and run RDP Classifier

# 20200909 Initial version by Jeffrey M. Dick (for geo16S paper)
# 20210922 Add --skip-technical --clip options to fastq-dump
#          Don't use fastq_filter --stripleft for 454 data
# 20220512 Use GNU parallel in classify()
# 20221013 Add GTDB option (for microhum paper)

## SYSTEM DEPENDENCIES

# R 4.0.0              https://www.r-project.org/
# fastq-dump 2.9.0     https://ncbi.github.io/sra-tools/
# vsearch 2.15.0       https://github.com/torognes/vsearch
# seqtk 1.3-r115       https://github.com/lh3/seqtk
# RDP Classifier 2.13  https://sourceforge.net/projects/rdp-classifier/
# RDP Classifier training files based on GTDB release 207
#                      https://doi.org/10.5281/zenodo.7633100
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
#   - Result is saved in <study>.tab in the RDP/ or RDP-GTDB/ directory
# > mkRDP()

## STUDY SETTINGS

# Change the following line to setup the pipeline for one study
study <- "WZL+23"
# Settings for all studies are stored here
file <- tempfile()
# Write spaces here (but don't save them) to make this easier to read
writeLines(con = file, text = gsub(" ", "", c(
  # For 454 experiments, set forwardonly to NA
  "study, forwardonly, trunclen",
  # COVID 20210731
  "MMP+21, TRUE, 200",
  # COVID 20210717
  "XLZ+21, TRUE, 250",
  "MAC+21, FALSE, 350",
  "NGH+21, TRUE, 300",
  # COVID 20210718
  "WCJ+21, FALSE, 290",
  "MLW+21, FALSE, 300",
  "ZZZ+21, TRUE, 400",
  "ENJ+21, FALSE, 290",
  "GBS+22, FALSE, 290",
  "GKJ+22, FALSE, 290",
  "HMH+21, FALSE, 460",
  # COVID 20210719
  "SRK+22, FALSE, 440",
  "RWC+21, FALSE, 450",
  "PMM+22, TRUE, 600",
  "VCV+21, FALSE, 450",
  "SRS+22, FALSE, 250",
  # COVID 20210720
  "CGC+22, TRUE, 400",
  "GCW+20, FALSE, 450",
  "RDM+22, FALSE, 440",
  "RFH+22, FALSE, 440",
  # Tissue experiments 20220721
  "BPB+21, FALSE, 250",
  # IBD 20220827
  "TWC+22, FALSE, 310",
  # COVID 20220921-20220922
  "FBD+22, TRUE, 305",
  "MIK+22, FALSE, 440",
  "GWL+21, FALSE, 460",
  "CSC+22, FALSE, 440",
  "IZC+21, TRUE, 160",
  # COVID 20221103-20221104
  "KMG+21, FALSE, 440",
  "SGC+21, FALSE, 420",
  # IBD 20230321-20230323
  "LAA+19, FALSE, 250",
  "MDV+22, TRUE, 150",
  "ZTG+21, FALSE, 440",
  "RAF+20, FALSE, 440",
  "ASM+23, TRUE, 290",
  "AAM+20, NA, NA",
  "LZD+19, NA, NA",
  "WGL+19, FALSE, 250",
  "GKD+14, FALSE, 250",
  # Human Microbiome Project 20231214
  "HMP12, TRUE, 400",
  # COVID and IBD gut 20231231
  "AHM+21, FALSE, 250",
  "BKK+17, TRUE, 440",
  "DKK+23, TRUE, 400",
  "HBL+17, TRUE, 99",
  "MLL+16, FALSE, 250",
  "MZW+23, FALSE, 410",
  "PYL+23, FALSE, 440",
  "REP+23, FALSE, 310",
  "WZL+23, FALSE, 410"
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
# Data directory for study
studydir <- file.path("/home/chem16S", study)
# Output directory for FASTA files
FASTAdir <- file.path(studydir, "fasta")
# Reference database for chimera detection
refdb <- "/home/chem16S/SILVA/SILVA_138.1_SSURef_NR99_tax_silva.fasta.gz"
# RDP jar file
RDPjar <- "/opt/rdp_classifier_2.13/dist/classifier.jar"
# Use GTDB <- TRUE and set GTDBpath to use the GTDB training set
#GTDB <- FALSE
GTDB <- TRUE
# Path to rRNAClassifier.properties
GTDBpath <- "/home/chem16S/RDP/GTDB/genus/training_files/rRNAClassifier.properties"

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

  if(study %in% c("TWC+22")) {
    # For TWC+22, FASTQ files were downloaded from Zenodo 20220829
    fqdump <- FALSE
  } else {
    # Generate input FASTQ files with fastq-dump
    fqdump <- TRUE
    cmd <- paste("fastq-dump --split-files --skip-technical --clip", RUNID)
    print(cmd)
    system(cmd)
  }

  if(is454) {
    # For 454 studies, fastq-dump --skip-technical puts biological sequences into file with highest-numbered suffix 
    # Use _3.fastq, _2.fastq or _1.fastq if they are present 20200921
    infile <- paste0(RUNID, "_3.fastq")
    if(!file.exists(infile)) infile <- paste0(RUNID, "_2.fastq")
    if(!file.exists(infile)) infile <- paste0(RUNID, "_1.fastq")
    # This is needed for files not generated by fastq-dump
    if(!file.exists(infile)) infile <- paste0(RUNID, ".fastq")
    minlen <- 300
    truncqual <- 15
    cmd <- paste("vsearch -fastq_filter", infile, "-fastq_minlen", minlen, "-fastq_maxlen 500 -fastq_truncqual", truncqual, "-fastaout", outfile)
    print(cmd)
    system(cmd)
  } else if(forwardonly) {
    # fastq-dump adds _3 suffix for some runs 20220607
    infile <- paste0(RUNID, "_3.fastq")
    if(!file.exists(infile)) infile <- paste0(RUNID, "_2.fastq")
    if(!file.exists(infile)) infile <- paste0(RUNID, "_1.fastq")
    # Use forward reads only
    file.copy(infile, "merged.fastq", overwrite = TRUE)
    # For Human Microbiome Project, reads are in different file 20211216
    if(study == "HMP12") file.copy(paste0(RUNID, "_4.fastq"), "merged.fastq", overwrite = TRUE)
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
    # For Illumina MiSeq 20220717
    if(study %in% c("XLZ+21")) qmax <- 42
    # For Nanopore 20220719
    if(study %in% c("PMM+22")) {
      qmax <- 90
      # No filtering on expected error rate
      maxee_rate <- 1
    }
    # Increase maxee_rate for selected datasets 20220720
    if(study %in% c("RDM+22")) maxee_rate <- 0.01
    cmd <- paste("vsearch -fastq_filter merged.fastq -fastq_trunclen", trunclen, "-fastq_qmax", qmax, "-fastq_maxee_rate", maxee_rate, "-fastaout", outfile)
    print(cmd)
    system(cmd)
  }

  if(study %in% c("TWC+22")) {
    # Prefix Run ID to header so we can extract the reads in the subsample() and findchimeras() steps 20220829
    lines <- readLines(outfile)
    ihead <- grep("^>", lines)
    lines[ihead] <- paste0(">", RUNID, ";", gsub("^>", "", lines[ihead]))
    writeLines(lines, outfile)
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
    # Output .fa file in FASTAdir
    outfile <- file.path(FASTAdir, paste0(thisRUNID, ".fa"))
    # Skip if output .fa file already exists
    if(file.exists(outfile)) {
      print(paste("File exists - skipping subsampling:", outfile))
      next
    }
    # Extract sequences for this sample
    thisfile <- paste0(thisRUNID, ".fasta")
    # https://stackoverflow.com/questions/26144692/printing-a-sequence-from-a-fasta-file
    # "awk '/>" is used to match sample name at beginning of header line
    # (this is needed for numeric or short sample names that can match other parts of the header)
    print(cmd <- paste0("awk '/>", thisRUNID, "/{p++;print;next} /^>/{p=0} p' not_singletons.fasta > ", thisfile))
    system(cmd)
    # Subsample 10000 sequences
    print(cmd <- paste('grep "^>"', thisfile))
    thisfile.headers <- system(cmd, intern = TRUE)
    nout <- nseq <- length(thisfile.headers)
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
classify <- function(RUNID, conf = 0.8, jobs = 4, dryrun = FALSE) {

  # Output directory for RDP Classifier
  if(GTDB) RDPdir <- file.path(studydir, "RDP-GTDB") else RDPdir <- file.path(studydir, "RDP")
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

  if(GTDB) {
    # Command for GTDB training set 20221013
    cmd <- paste("parallel --link --jobs", jobs, "java -jar", RDPjar, "classify -c", conf,
                 "-t", GTDBpath, "-h {3} -o {2} {1} :::: fastafiles tabfiles txtfiles")
  } else {
    cmd <- paste("parallel --link --jobs", jobs, "java -jar", RDPjar, "classify -c", conf,
                 "-h {3} -o {2} {1} :::: fastafiles tabfiles txtfiles")
  }
  message("=======================================")
  message(paste0("Classifying ", n, " runs with ", jobs, " jobs [", study, "]"))
  message("=======================================")

  print(cmd)
  if(dryrun) print("DRY RUN") else system(cmd)

  # Clean up
  file.remove(c("fastafiles", "tabfiles", "txtfiles"))
  file.remove(Sys.glob(paste0("*.tab")))
  if(!GTDB) file.remove(Sys.glob(paste0("cnadjusted_*")))
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

  if(GTDB) RDPdir <- file.path(studydir, "RDP-GTDB") else RDPdir <- file.path(studydir, "RDP")

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
  if(GTDB) outdir <- "RDP-GTDB" else outdir <- "RDP"
  outfile <- file.path(olddir, outdir, paste0(study, ".tab"))
  # Remove output file in case it exists
  if(file.exists(outfile)) file.remove(outfile)
  # Run merge-count command
  print(cmd <- paste("java -jar", RDPjar, "merge-count", outfile, paste(files, collapse = " ")))
  system(cmd)
}
