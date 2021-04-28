# comp16S/pipeline.R
# Merge and filter 16S SRA sequences, subsample, find chimeras, and run RDP Classifier
# 20200909 jmd v1

## SYSTEM REQUIREMENTS
# (version numbers for information only)

# R 4.0.0              https://www.r-project.org/
# fastq-dump 2.9.0     https://ncbi.github.io/sra-tools/
# vsearch 2.15.0       https://github.com/torognes/vsearch
# seqtk 1.3-r115       https://github.com/lh3/seqtk
# RDP Classifier 2.13  https://sourceforge.net/projects/rdp-classifier/
# Java openjdk version 1.8.0_252 (for RDP); GNU utils (cat, awk, grep)

## STUDY SETTINGS

# Change the following line to setup the pipeline for one study
study <- "ZLF+19"
# Settings for all studies are stored here
file <- tempfile()
# Write spaces here (but don't save them) to make this easier to read
writeLines(con = file, text = gsub(" ", "", c(
  # For 454 experiments, set forwardonly to NA
  "system, study, forwardonly, trunclen",

  # Natural environment datasets
  "environmental, BGPF13, NA, NA",
  "environmental, HLA+16, NA, NA",
  "environmental, JHM+16, FALSE, 230",
  "environmental, MPB+17, TRUE, 450",
  "environmental, XDZ+17, TRUE, 250",
  "environmental, VAH+18, FALSE, 250",
  "environmental, SVH+19, NA, NA",
  "environmental, ZLM+16, TRUE, 250",

  # Unconventional oil and gas datasets
  "environmental, UKD+18.sediment, TRUE, 100",
  "environmental, UKD+18.water, TRUE, 100",
  "environmental, CUN+18, TRUE, 250",
  "environmental, MMA+20, FALSE, 250",
  "environmental, CHM+14, NA, NA",
  "environmental, HRR+18, TRUE, 300",
  "environmental, ZLF+19, FALSE, 290"

)))
# This reads and applies the settings
settings <- read.csv(file, as.is = TRUE)
istudy <- match(study, settings$study)
# Different filtering for 454 studies
is454 <- is.na(settings$forwardonly[istudy])
# Use forward reads only (for low-quality or missing reverse reads)
forwardonly <- settings$forwardonly[istudy]
# Value of fastq_trunclen
trunclen <- settings$trunclen[istudy]

## FILE LOCATIONS

# Working directory
workdir <- "~/tmp/comp16S"
# Output directory for FASTA files
FASTAdir <- file.path("/home/comp16S", study, "fasta")
# Reference database for chimera detection
refdb <- "/home/comp16S/SILVA/SILVA_138.1_SSURef_NR99_tax_silva.fasta.gz"
#refdb <- "/home/comp16S/chimera/gold.fa"
# RDP jar file
RDPjar <- "/opt/rdp_classifier_2.13/dist/classifier.jar"
# Output directory for RDP Classifier
RDPdir <- file.path("/home/comp16S", study, "RDP")
# RDP Classifier version (java -jar /opt/rdp_classifier_2.13/dist/classifier.jar version)
#Gene:16srrna    Trainset No:18  Taxonomy Version:RDP 16S rRNA training setNo 18 07/2020
#RDP Classifier Version:RDP Naive Bayesian rRNA Classifier Version 2.11, September 2015

## USAGE

# First, make sure the working directory is clean (especially no *.fa files)
# To filter one file (generates .fa files in workdir for subsampling)
# filter("SRR2059382")

# To filter multiple files:
# SRR <- c("SRR2059381", "SRR2059380", "SRR2059379")
# lapply(SRR, filter)

# To subsample and remove chimeras from the filtered files:
# subsample()
# findchimeras()

# To run RDP Classifier after chimera removal:
# classify(SRR)

filter <- function(SRR) {
  message("============================================")
  message(paste0("Merging and filtering ", SRR, " [", study, "]"))
  message("============================================")

  # Change to working directory
  olddir <- setwd(workdir)
  on.exit(setwd(olddir))
  # Generate input FASTQ files
  cmd <- paste("fastq-dump --split-files", SRR)
  if(study == "CCN+16") cmd <- paste("fastq-dump --split-files -X 500000", SRR)
  print(cmd)
  system(cmd)
  outfile <- paste0(SRR, ".fa")
  if(is454) {
    # For some 454 studies, fastq-dump puts barcodes into _2.fastq, primer and experimental sequence into _3.fastq
    # Use -stripleft to remove primer sequence; minimum and maximum length as described by KGP+12
    # Use _2.fastq or _1.fastq if needed 20200921
    infile <- paste0(SRR, "_3.fastq")
    if(!file.exists(infile)) infile <- paste0(SRR, "_2.fastq")
    if(!file.exists(infile)) infile <- paste0(SRR, "_1.fastq")
    cmd <- paste("vsearch -fastq_filter", infile, "-fastq_stripleft 18 -fastq_minlen 200 -fastq_maxlen 600 -fastq_truncqual 15 -fastaout", outfile)
    print(cmd)
    system(cmd)
  } else if(forwardonly) {
    # Use forward reads only
    file.copy(paste0(SRR, "_1.fastq"), "merged.fastq", overwrite = TRUE)
    nseq <- length(readLines("merged.fastq")) / 4
    print(paste0("Using forward reads only (", nseq, " sequences)"))
  } else {
    # Check that we have same number of forward and reverse reads 20200917
    # (issue with some datasets: fastq_dump rejects some reads because READLEN < 1)
    ID1 <- sapply(strsplit(gsub("@", "", system(paste0("awk '(NR - 1) % 4 == 0' ", SRR, "_1.fastq"), intern = TRUE)), " "), "[", 1)
    ID2 <- sapply(strsplit(gsub("@", "", system(paste0("awk '(NR - 1) % 4 == 0' ", SRR, "_2.fastq"), intern = TRUE)), " "), "[", 1)
    # Forward reads that have matching reverse reads
    forward_matches <- ID1[ID1 %in% ID2]
    if(length(forward_matches) < length(ID1)) {
      writeLines(forward_matches, "forward_matches.txt")
      print(cmd <- paste0("seqtk subseq ", SRR, "_1.fastq forward_matches.txt > new_1.fastq"))
      system(cmd)
      print(cmd <- paste0("mv new_1.fastq ", SRR, "_1.fastq"))
      system(cmd)
    }
    # Reverse reads that have matching forward reads
    reverse_matches <- ID2[ID2 %in% ID1]
    if(length(reverse_matches) < length(ID2)) {
      writeLines(reverse_matches, "reverse_matches.txt")
      print(cmd <- paste0("seqtk subseq ", SRR, "_2.fastq reverse_matches.txt > new_2.fastq"))
      system(cmd)
      print(cmd <- paste0("mv new_2.fastq ", SRR, "_2.fastq"))
      system(cmd)
    }
    # Make sure forward and reverse reads have same IDs
    ID1 <- sapply(strsplit(gsub("@", "", system(paste0("awk '(NR - 1) % 4 == 0' ", SRR, "_1.fastq"), intern = TRUE)), " "), "[", 1)
    ID2 <- sapply(strsplit(gsub("@", "", system(paste0("awk '(NR - 1) % 4 == 0' ", SRR, "_2.fastq"), intern = TRUE)), " "), "[", 1)
    if(!all(ID1 == ID2)) stop(paste0("sequence IDs in ", SRR, "_1.fastq and ", SRR, "_2.fastq don't match!"))

    # Merge paired reads
    print(cmd <- paste0("vsearch -fastq_mergepairs ", SRR, "_1.fastq -reverse ", SRR, "_2.fastq -fastqout merged.fastq"))
    system(cmd)
  }
  if(!is454) {
    # Remove short and high-error reads
    # Use fastq_maxee_rate instead of fastq_maxee 20200920
    cmd <- paste("vsearch -fastq_filter merged.fastq -fastq_trunclen", trunclen, "-fastq_maxee_rate 0.005 -fastaout", outfile)
    # JHM+16: Very high expected error rate 20200923
    if(study == "JHM+16") cmd <- paste("vsearch -fastq_filter merged.fastq -fastq_trunclen", trunclen, "-fastaout", outfile)
    # WYF+20: Ion Torrent needs -fastq_qmax 44 20200928
    if(study == "WYF+20") cmd <- paste("vsearch -fastq_filter merged.fastq -fastq_trunclen", trunclen, "-fastq_maxee_rate 0.005 -fastq_qmax 44 -fastaout", outfile)
    # FPP+18: Ion Torrent needs -fastq_qmax 45 20201104
    if(study == "FPP+18") cmd <- paste("vsearch -fastq_filter merged.fastq -fastq_trunclen", trunclen, "-fastq_maxee_rate 0.005 -fastq_qmax 45 -fastaout", outfile)
    print(cmd)
    system(cmd)
  }

  # Clean up
  file.remove(Sys.glob("*.fastq"))
  return()
}

subsample <- function() {
  message("=======================")
  message(paste0("Subsampling [", study, "]"))
  message("=======================")

  # Change to working directory
  olddir <- setwd(workdir)
  on.exit(setwd(olddir))
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

  # Loop over samples
  allSRR <- gsub("\\.fa$", "", dir(pattern = "\\.fa$"))
  for(thisSRR in allSRR) {
    # Extract sequences for this sample
    thisfile <- paste0(thisSRR, ".fasta")
    # https://stackoverflow.com/questions/26144692/printing-a-sequence-from-a-fasta-file
    print(cmd <- paste0("awk '/>", thisSRR, "/{p++;print;next} /^>/{p=0} p' not_singletons.fasta > ", thisfile))
    system(cmd)

    # Subsample 10000 sequences - put output in .fa file in FASTAdir
    print(cmd <- paste('grep "^>"', thisfile))
    thisfile.headers <- system(cmd, intern = TRUE)
    nout <- nseq <- length(thisfile.headers)
    outfile <- file.path(FASTAdir, paste0(thisSRR, ".fa"))
    if(nseq > 10000) {
      print(cmd <- paste("vsearch -fastx_subsample", thisfile, "-sample_size 10000 -randseed 1234 -fastaout", outfile))
      system(cmd)
      nout <- 10000
    } else {
      file.copy(thisfile, outfile, overwrite = TRUE)
    }
    subsamptxt <- paste("Subsampled", nout, "of", nseq, "sequences")
    print(subsamptxt)
    # Save messages about number of subsampled sequences
    txtfile <- file.path(FASTAdir, paste0(thisSRR, ".txt"))
    writeLines(subsamptxt, txtfile)
  }

  # Clean up
  file.remove(Sys.glob("*.fa"))
  file.remove(Sys.glob("*.fasta"))
  return()
}

# Find chimeras in combined FASTA file, then split result into samples 20200913
# (Saves time by masking and creating k-mer index only once)
findchimeras <- function() {
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
  print(cmd <- paste("vsearch -uchime_ref combined.fasta -nonchimeras nonchimeras.fasta -db", refdb))
  system(cmd)

  # Split result into one file for each sample
  # https://stackoverflow.com/questions/26144692/printing-a-sequence-from-a-fasta-file
  allSRR <- gsub("\\.fa$", "", dir(pattern = "\\.fa$"))
  for(thisSRR in allSRR) {
    print(cmd <- paste0("awk '/>", thisSRR, "/{p++;print;next} /^>/{p=0} p' nonchimeras.fasta > ", thisSRR, ".fasta"))
    system(cmd)
  }

  # Clean up
  file.remove("combined.fasta")
  file.remove("nonchimeras.fasta")
}

# Run RDP Classifier on one or multiple files 20200910
classify <- function(SRR, conf = 0.8) {
  # Make sure RDP output directory exists 20200915
  if(!dir.exists(RDPdir)) stop(paste("directory", RDPdir, "does not exist"))
  # Force evaluation of SRR before changing the directory
  allSRR <- SRR
  # Change to working directory
  olddir <- setwd(workdir)
  on.exit(setwd(olddir))

  # Loop over SRR IDs
  for(thisSRR in allSRR) {
    message("===============================")
    message(paste0("Classifying ", thisSRR, " [", study, "]"))
    message("===============================")

    # Input FASTA file
    FASTAfile <- file.path(FASTAdir, paste0(thisSRR, ".fasta"))

    # Copy FASTA file to here without .fasta suffix (so it doesn't get into sample name)
    file.copy(FASTAfile, thisSRR, overwrite = TRUE)
    # Run classifier
#    print(cmd <- paste0("java -jar ", RDPjar, " classify -o ", thisSRR, ".tab -h ", thisSRR, ".txt ", thisSRR))
    # Add -c (confidence cutoff) 20201007
    print(cmd <- paste0("java -jar ", RDPjar, " classify -c ", conf, " -o ", thisSRR, ".tab -h ", thisSRR, ".txt ", thisSRR))
    system(cmd)

    # Copy results (assignment count for each taxon) to output directories
    results <- paste0(thisSRR, ".txt")
    resfile <- file.path(RDPdir, paste0(thisSRR, ".txt"))
    file.copy(results, resfile, overwrite = TRUE)

    # Clean up
    file.remove(Sys.glob(paste0("*", thisSRR, "*")))
  }
}
