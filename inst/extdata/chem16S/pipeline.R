# chem16S/pipeline.R
# Merge and filter 16S sequences from SRA, subsample, find chimeras, and run RDP Classifier

# 20200909 Initial version by Jeffrey M. Dick
# 20210922 Add settings for orp16S datasets
#          Add --skip-technical --clip options to fastq-dump
#          Don't use fastq_filter --stripleft for 454 data

## SYSTEM REQUIREMENTS
# (version numbers for information only)

# R 4.0.0              https://www.r-project.org/
# fastq-dump 2.9.0     https://ncbi.github.io/sra-tools/
# vsearch 2.15.0       https://github.com/torognes/vsearch
# seqtk 1.3-r115       https://github.com/lh3/seqtk
# RDP Classifier 2.13  https://sourceforge.net/projects/rdp-classifier/
# Java openjdk version 1.8.0_252 (for RDP); GNU utils (cat, awk, grep)

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
# To combine output of RDP Classifier into one file
# -- result is saved in RDP/<study>.tab in the current directory
# mkRDP(study, SRR)

## STUDY SETTINGS

# Change the following line to setup the pipeline for one study
study <- "PCL+18"
# Settings for all studies are stored here
file <- tempfile()
# Write spaces here (but don't save them) to make this easier to read
writeLines(con = file, text = gsub(" ", "", c(
  # For 454 experiments, set forwardonly to NA
  "study, forwardonly, trunclen",

  ## For geo16S paper 20210502

  # Natural environment datasets
  # For HCW+13 (Guerrero Negro), don't use filter(); start with findchimeras()
  #"HCW+13, NA, NA",      # GenBank JN427016-JN539989  10.1038/ismej.2012.79
  "BGPF13, NA, NA",      # BioProject PRJNA207095  10.3389/fmicb.2013.00330
  "HLA+16, NA, NA",      # PRJEB1245    10.3389/fmicb.2016.01883
  "JHM+16, FALSE, 230",  # PRJNA291280  10.1128/AEM.02699-15
  "ZLM+16, TRUE, 250",   # PRJNA294836  10.1128/AEM.03332-15
  "MPB+17, TRUE, 450",   # PRJEB15554   10.1038/ismej.2017.37
  "XDZ+17, TRUE, 250",   # PRJNA388250  10.1038/s41598-017-13608-5
  "SVH+19, NA, NA",      # PRJNA423140  10.1111/gbi.12316

  # Shale gas datasets
  "UKD+18.sediment, TRUE, 100",  # PRJNA394724  10.1038/s41598-018-23679-7
  "UKD+18.water, TRUE, 100",
  "MMA+20, FALSE, 250",  # PRJNA544240  10.1073/pnas.1911458117
  "CHM+14, NA, NA",      # PRJNA229085  10.1021/es501173p
  "HRR+18, TRUE, 300",   # PRJNA438710  10.1016/j.scitotenv.2018.06.067
  "ZLF+19, FALSE, 290",  # PRJNA407226  10.1021/acsearthspacechem.9b00037

  # Stratified water body datasets
  "GBL+15, FALSE, 250",  # PRJNA263621  10.1038/ismej.2015.44
  "BCA+21, FALSE, 450",  # PRJNA395513  10.1111/1462-2920.14909
  "MZG+20, FALSE, 450",  # PRJEB27579   10.1038/s41396-019-0515-8
  "HXZ+20, FALSE, 440",  # PRJNA503500  10.1038/s41598-020-62411-2

  ## For orp16S paper 20210922

  "MLL+19, TRUE, 250",
  "RSJ+21, FALSE, 290",
  "RMB+17, FALSE, 250",
  "NTB+21, TRUE, 150",
  "SBP+20, FALSE, 250",
  "MWY+21, FALSE, 250",
  "SAR+13, NA, NA",
  "CTS+17, FALSE, 400",
  "HSF+19, FALSE, 440",
  "SCM+18, TRUE, 250",
  "HDZ+19, FALSE, 420",
  "BOEM21, FALSE, 400",
  "ZHZ+19, FALSE, 300",
  "YHK+20, FALSE, 400",
  "CNA+20, TRUE, 450",
  "BMJ+19, FALSE, 300",
  "SRM+19, FALSE, 250",
  "HLZ+18, FALSE, 420",
  "XLD+20, TRUE, 250",
  "JHL+12, NA, NA",
  "PSG+20, FALSE, 250",
  "KSR+21, FALSE, 440",
  "ZCZ+21, FALSE, 450",
  "LLZ+20, TRUE, 150",
  "SKP+21, FALSE, 200",
  "ZZL+21, FALSE, 450",
  "PBU+20, FALSE, 400",
  "GWS+19, FALSE, 300",
  "KLY+20, TRUE, 250",
  "SRM+21, FALSE, 250",
  "MLL+18, TRUE, 250",
  "JDP+20, TRUE, 275",
  "BWD+19, FALSE, 400",
  "LXH+20, TRUE, 450",
  "LMG+20, TRUE, 280",
  "WHL+21, FALSE, 400",
  "LLL+21, FALSE, 450",
  "SDH+19, FALSE, 300",
  "GWSS21, FALSE, 300",
  "ZML+17, TRUE, 400",
  "RBW+14, FALSE, 250", # Winogradsky column
  "DTJ+20, FALSE, 420",
  "WFB+21, FALSE, 440",
  "SBW+17, TRUE, 280",
  "KLM+16, NA, NA",
  "LMBA21, TRUE, 280",
  "ZDA+20, TRUE, 450",
  "ZZZ+18, TRUE, 240",
  "BSPD17, FALSE, 400",
  "CWC+20, FALSE, 440",
  "BMOB18, TRUE, 350",
  "JVW+20, FALSE, 250",
  "LJC+20, FALSE, 420",
  "GFE+16, NA, NA",
  "ECS+18, TRUE, 400",
  "FAV+21, FALSE, 400",
  "VMB+19, TRUE, 380",
  "DLS21, TRUE, 300",
  "ZZLL21, FALSE, 450",
  "GWS+20, FALSE, 300",
  "CLS+19, TRUE, 150",
  # For SMS+12 (Bison Pool), start with findchimeras()
  #"SMS+12, NA, NA",
  "BYB+17, NA, NA",
  "MCS+21, FALSE, 400",
  "PMM+20, FALSE, 400",
  "GZL21, FALSE, 420",
  "SCH+16, TRUE, 250",
  "APV+20, FALSE, 250",
  "NLE+21, FALSE, 440",
  "GSY+20, TRUE, 120",
  "WHC+19, FALSE, 420",
  "LZR+17, TRUE, 250",
  "LLC+19, TRUE, 400",
  "YHK+19, TRUE, 250",
  "WHLH21, TRUE, 440",
  "GRG+20, NA, NA",
  "PCL+18, TRUE, 250"

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

filter <- function(SRR) {
  message("============================================")
  message(paste0("Merging and filtering ", SRR, " [", study, "]"))
  message("============================================")

  # Change to working directory
  olddir <- setwd(workdir)
  on.exit(setwd(olddir))
  # The output file from this function is a FASTA file with .fa suffix
  outfile <- paste0(SRR, ".fa")

  if(study %in% c("BYB+17", "GRG+20")) {
    # For BYB+17, FASTQ files are downloaded from SRA cloud 20210917
    fqdump <- FALSE
  } else {
    # Generate input FASTQ files with fastq-dump
    fqdump <- TRUE
    cmd <- paste("fastq-dump --split-files --skip-technical --clip", SRR)
    if(study == "SCH+16") cmd <- paste("fastq-dump --split-files --skip-technical --clip -X 100000", SRR)
    print(cmd)
    system(cmd)
  }

  if(study %in% c("GRG+20")) {
    # These are files from MG-RAST (.fasta suffix)
    # Prefix MG-RAST ID to header so we can extract the reads in the subsample() and findchimeras() steps 20211004
    fastafile <- paste0(SRR, ".fasta")
    lines <- readLines(fastafile)
    ihead <- grep("^>", lines)
    lines[ihead] <- paste0(">", SRR, ";", gsub("^>", "", lines[ihead]))
    writeLines(lines, outfile)
  } else if(is454) {
    # For 454 studies, fastq-dump --skip-technical puts biological sequences into file with highest-numbered suffix 
    # Use _3.fastq, _2.fastq or _1.fastq if they are present 20200921
    infile <- paste0(SRR, "_3.fastq")
    if(!file.exists(infile)) infile <- paste0(SRR, "_2.fastq")
    if(!file.exists(infile)) infile <- paste0(SRR, "_1.fastq")
    # This is needed for files not generated by fastq-dump (BYB+17)
    if(!file.exists(infile)) infile <- paste0(SRR, ".fastq")
    truncqual <- 15
    fastq_minlen <- 200
    if(study %in% c("SAR+13", "GFE+16")) truncqual <- 11
    if(study %in% c("GFE+16")) fastq_minlen <- 100
    cmd <- paste("vsearch -fastq_filter", infile, "-fastq_minlen", fastq_minlen, "-fastq_maxlen 600 -fastq_truncqual", truncqual, "-fastaout", outfile)
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
    maxee_rate <- 0.005
    if(study == "YHK+19") maxee_rate <- 0.05
    # Use default -fastq_qmax except for some studies 20210922
    qmax <- 41
    # For Illumina MiSeq (strange quality encoding?) 20210909
    if(study == "ECS+18") qmax <- 93
    # For Ion Torrent 20210928
    if(study == "PCL+18") qmax <- 45
    cmd <- paste("vsearch -fastq_filter merged.fastq -fastq_trunclen", trunclen, "-fastq_qmax", qmax, "-fastq_maxee_rate", maxee_rate, "-fastaout", outfile)
    print(cmd)
    system(cmd)
  }
  if(study %in% c("BYB+17")) {
    # Append run ID to header so we can extract the reads in the subsample() and findchimeras() steps 20210913
    lines <- readLines(outfile)
    ihead <- grep("^>", lines)
    lines[ihead] <- paste0(lines[ihead], ";", SRR)
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
  allSRR <- gsub("\\.fa$", "", dir(pattern = "\\.fa$"))
  for(thisSRR in allSRR) {
    # Extract sequences for this sample
    thisfile <- paste0(thisSRR, ".fasta")
    # https://stackoverflow.com/questions/26144692/printing-a-sequence-from-a-fasta-file
    # Don't use "awk '/>" because sample names for SBE+17 are at end of header of original FASTQ files, not beginning 20210501
    print(cmd <- paste0("awk '/", thisSRR, "/{p++;print;next} /^>/{p=0} p' not_singletons.fasta > ", thisfile))
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
    # Don't use "awk '/>" because sample names for SBE+17 are at end of header of original FASTQ files, not beginning 20210501
    print(cmd <- paste0("awk '/", thisSRR, "/{p++;print;next} /^>/{p=0} p' nonchimeras.fasta > ", thisSRR, ".fasta"))
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
  n <- length(allSRR)
  for(i in 1:n) {
    thisSRR <- allSRR[i]
    message("=======================================")
    message(paste0("Classifying ", thisSRR, " [", study, "] (", i, "/", n, ")"))
    message("=======================================")

    # Input FASTA file
    FASTAfile <- file.path(FASTAdir, paste0(thisSRR, ".fasta"))

    # Copy FASTA file to here without .fasta suffix (so it doesn't get into sample name)
    file.copy(FASTAfile, thisSRR, overwrite = TRUE)
    # Run classifier with -c (confidence score threshold) 20201007
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
