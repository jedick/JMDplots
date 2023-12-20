# ARAST.R
# Pipeline in R similar to MG-RAST
# (only up to RNA and protein genecalling steps)

# 20180310 v0 by Jeffrey M. Dick
# 20181117 v1 - Used in gradox paper (https://doi.org/10.3389/fmicb.2019.00120)
#   Script is on Zenodo: https://doi.org/10.5281/zenodo.2314933
# 20211216 v2 - Used in geo16S paper (https://doi.org/10.1007/s00248-022-01988-9)
#   Add techtype argument to ARAST() (script no longer has hardcoded list of IDs and techtypes)
# 20221118 v3
#   Remove human DNA using bowtie2 (step 5 added to script)
# 20231216 v4 - Used in microhum paper
#   Output processing statistics

# USAGE EXAMPLES
# To process a batch of files:
#   files <- dir("/home/ARAST/work/", ".fasta", full.names = TRUE)
#   lapply(files, ARAST)
# To stop after the dereplication step and don't gzip the output:
#   lapply(files, ARAST, run = "DNA", run.gzip = FALSE)

# SYSTEM SOFTWARE: modified from MG-RAST's "Processing Information and Downloads" seen for ERX242607_sra_data.fasta
# skewer_0.2.2 (https://github.com/relipmoc/skewer)
# bowtie_v2.5.0 (https://github.com/BenLangmead/bowtie2) [used by autoskewer and for decontamination step]
# autoskewer (https://github.com/MG-RAST/autoskewer)
# ea-utils (https://github.com/ExpressionAnalysis/ea-utils) [fastq-mcf, used by filter_sequences]
# biopython (pip install biopython) [used by seq_length_stats.py]
# pipeline (https://github.com/MG-RAST/pipeline) [seq_length_stats.py, filter_sequences, dereplication.py, parallel_FragGeneScan.py]
# sortmerna (https://github.com/biocore/sortmerna)
# FragGeneScan_FGS1.18 (https://github.com/MG-RAST/FGS)

# DATABASE:
# download md5rna.clust (pipeline-master/CWL/Inputs/DBs/getpredata.sh) [used by sortmerna]
#   curl "http://shock.metagenomics.anl.gov/node/c4c76c22-297b-4404-af5c-8cd98e580f2a?download" > md5rna.clust
# set up index for sortmerna:
#   indexdb --ref md5rna.clust,md5rna.clust.index
# Human genome index for bowtie2
#   wget https://genome-idx.s3.amazonaws.com/bt/GRCh38_noalt_as.zip

ARAST <- function(inputfile, proc = parallel::detectCores(), run = "protein", run.gzip = TRUE, techtype = "illumina", mem_MB = 6144, mem_GB = 16) {

  # This function processes a single FASTA or FASTQ file
  #   'proc' number of threads for bowtie and sortmerna
  #     NOTE: Edit skewoptions in the autoskewer.py script to adjust number of threads for skewer
  #   'techtype' can be sanger, 454, or illumina (sequencing technology needed for FragGeneScan)
  #   'mem_MB' is used by arast_sortme_rna.pl
  #     NOTE: value is limited by free memory resources; tested values include 3941 for 8 GB RAM and 6144 for 16 GB RAM
  #   'mem_GB' is used by dereplication.py
  #   'run' specifies which steps to run:
  #     DNA: trimming, denoising, artifact removal, decontamination
  #     RNA: DNA steps plus RNA genecalling with sortmerna
  #     protein: DNA and RNA steps plus protein genecalling with FragGeneScan
  #     FGS: FragGeneScan only

  ## What steps to run
  steps <- numeric()
  if(run == "DNA") steps <- c(0:5, 23:24)
  if(run == "RNA") steps <- c(0:6, 23:24)
  if(run == "protein") steps <- c(0:9, 23:24)
  if(run == "protein_no_decontamination") steps <- c(0:4, 6:9, 23:24)
  if(run == "FGS") steps <- c(0:1, 9, 23:24)

  if(length(steps) == 0) stop("Invalid value for 'run' argument")

  if(0 %in% steps) {
    print("######## 0. Set up environment ########")
    # Get working directory and run ID from file name
    workdir <- dirname(inputfile)
    ID <- strsplit(basename(inputfile), "_")[[1]][1]
    # File type: fasta or fastq
    if(grepl("\\.fasta$", inputfile)) filetype <- "fasta" else
    if(grepl("\\.fna$", inputfile)) filetype <- "fasta" else
    if(grepl("\\.fastq$", inputfile)) filetype <- "fastq" else {
      print(paste("Not processing", inputfile, "which is not fasta or fastq"))
      return()
    }
    print(paste0("Processing ", filetype, " file for ", ID, " (", techtype, " sequencing)"))
    tmpdir <- file.path(workdir, "tmp")
    if(!dir.exists(tmpdir)) dir.create(tmpdir)
    # Keep list of temporary and output files
    tmpfiles <- outfiles <- character()
  } else stop("Step 0 is required")

  if(1 %in% steps) {
    print("######## 1. Initial sequence statistics ########")
    # Initialize list of statistics
    stats <- list(ID = ID)
    # Read output from fastq-dump, if available
    fastq_dump.logfile <- file.path(workdir, paste0(ID, "_fastq-dump.log"))
    if(file.exists(fastq_dump.logfile)) {
      fastq_dump.log <- readLines(fastq_dump.logfile)
      # Find lines starting with "Read" and "Written"
      iRead <- grep("^Read", fastq_dump.log)
      iWritten <- grep("^Written", fastq_dump.log)
      available_spots <- as.numeric(strsplit(fastq_dump.log[iRead], " ")[[1]][2])
      passed_read_filter <- as.numeric(strsplit(fastq_dump.log[iWritten], " ")[[1]][2])
      # List numbers of available spots and those written by fastq-dump --read-filter pass
      stats <- c(stats, available_spots = available_spots, passed_read_filter = passed_read_filter)
    }
    # Count number of sequences in input file
    input_sequences <- count_sequences(inputfile, filetype)
    stats <- c(stats, inputfile = inputfile, input_sequences = input_sequences)
  } else stop("Step 1 is required")

  if(2 %in% steps) {
    print("######## 2. Adapter Trimming ########")
    scrubbedfile <- file.path(workdir, paste0(ID, "_scrubbed.", filetype))
    autoskewer(inputfile, scrubbedfile, tmpdir)
    tmpfiles <- c(tmpfiles, scrubbedfile)
    scrubbed_sequences <- count_sequences(scrubbedfile, filetype)
    stats <- c(stats, scrubbed_sequencs = scrubbed_sequences)
  } else scrubbedfile <- inputfile

  if(3 %in% steps) {
    print("######## 3. Denoising and normalization ########")
    passedfile <- file.path(workdir, paste0(ID, "_passed.fasta"))
    removedfile <- file.path(workdir, paste0(ID, "_removed.fasta"))
    preprocess(scrubbedfile, passedfile, removedfile, filetype)
    tmpfiles <- c(tmpfiles, passedfile, removedfile)
    passed_sequences <- count_sequences(passedfile, "fasta")
    stats <- c(stats, passed_sequences = passed_sequences)
  } else passedfile <- scrubbedfile

  if(4 %in% steps) {
    print("######## 4. Removal of sequencing artifacts ########")
    derepfile <- file.path(workdir, paste0(ID, "_dereplicated.fasta"))
    dereplication(passedfile, derepfile, tmpdir, mem_GB)
    # Mark this as a temporary file unless this is the final step
    if(5 %in% steps) tmpfiles <- c(tmpfiles, derepfile) else outfiles <- c(outfiles, derepfile)
    dereplicated_sequences <- count_sequences(derepfile, "fasta")
    stats <- c(stats, dereplicated_sequences = dereplicated_sequences)
  } else derepfile <- passedfile

  if(5 %in% steps) {
    print("######## 5. Host DNA contamination removal ########")
    decontamfile <- file.path(workdir, paste0(ID, "_decontaminated.fasta"))
    bowtie2(derepfile, decontamfile, proc, tmpdir)
    if(6 %in% steps) tmpfiles <- c(tmpfiles, decontamfile) else outfiles <- c(outfiles, decontamfile)
    decontaminated_sequences <- count_sequences(decontamfile, "fasta")
    stats <- c(stats, decontaminated_sequences = decontaminated_sequences)
  } else decontamfile <- derepfile

  if(6 %in% steps) {
    print("######## 6. RNA feature identification (rRNA genecalling) ########")
    rna_nr = "/home/ARAST/REFDB/md5rna.clust"
    index <- paste0(rna_nr, ".index")
    rRNAfile <- file.path(workdir, paste0(ID, "_rRNA.fasta"))
    not_rRNAfile <- file.path(workdir, paste0(ID, "_not-rRNA.fasta"))
    sortmerna(decontamfile, rRNAfile, not_rRNAfile, proc, mem_MB, rna_nr, index)
    outfiles <- c(outfiles, rRNAfile)
    if(9 %in% steps) tmpfiles <- c(tmpfiles, not_rRNAfile) else outfiles <- c(outfiles, not_rRNAfile)
    rRNA_sequences <- count_sequences(rRNAfile, "fasta")
    not_rRNA_sequences <- count_sequences(not_rRNAfile, "fasta")
    stats <- c(stats, rRNA_sequences = rRNA_sequences, not_rRNA_sequences = not_rRNA_sequences)
  } else not_rRNAfile <- decontamfile

  # NOTE: Steps labelled 00000000 are in the MG-RAST pipeline but are not implemented here
  if(7 %in% steps) print("00000000 7. RNA clustering (not implemented) 00000000")

  if(8 %in% steps) print("00000000 8. RNA similarity search (not implemented) 00000000")

  if(9 %in% steps) {
    print("######## 9. Identify putative protein coding features (genecalling) ########")
    codingfile <- file.path(workdir, paste0(ID, "_coding"))
    genecalling(not_rRNAfile, codingfile, techtype, tmpdir, proc)
    if(file.exists(paste0(codingfile, ".faa"))) {
      fnafile <- paste0(codingfile, ".fna")
      faafile <- paste0(codingfile, ".faa")
      tmpfiles <- c(tmpfiles, fnafile)
      outfiles <- c(outfiles, faafile)
      coding_sequences <- count_sequences(faafile, "fasta")
      stats <- c(stats, coding_sequences = coding_sequences)
      # Get average protein length
      cmd <- paste('grep -v "^>"', faafile, '| wc -c')
      nbytes <- as.numeric(system(cmd, intern = TRUE))
      # Subtract the number of newlines (one for each sequence) from the number of bytes and divide by number of sequences to get average sequence length
      protein_length <- ( nbytes - coding_sequences ) / coding_sequences
      stats <- c(stats, protein_length = protein_length)
    } else stats <- c(stats, coding_sequences = 0, protein_length = NA)
  }

  if(23 %in% steps) {
    print("######## 23. Summary statistics ########")
    # Remove temporary files and compress output files
    system(paste("rm -f", paste(tmpfiles, collapse = " ")))
    if(run.gzip & length(outfiles) > 0) system(paste("gzip -f", paste(outfiles, collapse = " ")))
    # Save stats file
    statsfile <- file.path(workdir, paste0(ID, "_stats.csv"))
    write.csv(as.data.frame(stats), statsfile, row.names = FALSE, quote = FALSE)
  }

  if(24 %in% steps) print("######## 24. Completed ########")

}

# Run autoskewer.py 20180310
autoskewer <- function(file, scrubbedfile, tmpdir) {
  cmd <- paste("autoskewer.py -v -i", file, "-o", scrubbedfile, "-l ~/tmp/scrubbed.log -t", tmpdir)
  system(cmd)
  # Clean up tmpdir
  tmpfiles <- file.path(tmpdir, "*")
  cmd <- paste("rm", tmpfiles)
  system(cmd)
}

# Run filter_sequences 20180309
# Based on MG-RAST/pipeline/mgcmd/mgrast_preprocess.pl
preprocess <- function(file, passedfile, removedfile, filetype) {
  # Print stats
  cmd <- paste("seq_length_stats.py -i", file, "-t", filetype, "-f")
  system(cmd)
  if(filetype=="fasta") {
    # Get stats
    cmd <- paste("seq_length_stats.py -i", file, "-t", filetype, "-f | cut -f2")
    out <- system(cmd, intern=TRUE)
    # This gives bp_count, sequence_count, average_length, standard_deviation_length, length_min, length_max
    out <- as.numeric(out)
    # Translated from mgrast_preprocess.pl
    mymean <- out[3]
    mystdv <- out[4]
    min_ln <- round(mymean - 2 * mystdv)
    max_ln <- round(mymean + 2 * mystdv)
    if(min_ln < 0) min_ln <- 0
    cmd <- paste("filter_sequences -i", file, "-format fasta -o", passedfile, "-r", removedfile,
                  "-filter_ln -min_ln", min_ln, "-max_ln", max_ln, "-filter_ambig -max_ambig 5")
  }
  if(filetype=="fastq") {
    cmd <- paste("filter_sequences -i", file, "-format fastq -o", passedfile, "-r", removedfile,
                  "-dynamic_trim -min_qual 15 -max_lqb 5")
  }
  system(cmd)
}

# Run dereplication.py 20180310
# Based on MG-RAST/pipeline/mgcmd/mgrast_dereplicate.pl
dereplication <- function(file, derepfile, tmpdir, mem_GB) {
  prefix_size <- 50
  run_dir <- tmpdir
  input <- file
  out_prefix <- tmpdir
  cmd <- paste("dereplication.py -l", prefix_size, "-m", mem_GB, "-d", run_dir, input, out_prefix)
  system(cmd)
  # Rename and clean up files
  system(paste("mv", paste0(tmpdir, ".passed.fna"), derepfile))
  system(paste("rm", paste0(tmpdir, ".*")))
}

# Run sortmerna 20180310
# Based on MG-RAST/pipeline/mgcmd/mgrast_sortme_rna.pl
sortmerna <- function(file, rRNAfile, not_rRNAfile, proc, mem_MB, rna_nr, index) {
  eval <- 0.1
  ident <- 70 # cutoff value listed in the MG-RAST Processing Information page
  fasta <- file
  ## NOTE: MG-RAST picks sequences using the BLAST ouput with an identity threshold;
  ##   here we use sortmerna with command modified for FASTA output (--fastx),
  ##   fast filtering (--num_alignments), and saving both aligned and rejected sequences;
  ##   sortmerna doesn't have an identity threshold option
  #ref <- paste0(rna_nr, ",", index)
  #cmd <- paste("sortmerna -a", proc, "-m", mem, "-e", eval,
  #             "--fastx --num_alignments 1 --ref", ref, "--reads", fasta,
  #             "--aligned", rRNAfile, "--other", not_rRNAfile, "-v")
  # Use a modified version of mgrast_sortme_rna.pl to save both RNA and non-RNA sequences 20180311
  cmd <- paste("./arast_sortme_rna.pl -input", file, "-output", rRNAfile, "-notrna", not_rRNAfile,
               "-rna_nr", rna_nr, "-index", index, "-proc", proc, "-mem", mem_MB, "-eval", eval, "-ident", ident)
  #print(cmd)  
  system(cmd)
  # Clean up files
  WD <- setwd(dirname(file))
  system("rm *.blast *.tab *.id *.sort")
  setwd(WD)
}

# Run FragGeneScan 20180309
# Based on MG-RAST/pipeline/mgcmd/mgrast_genecalling.pl
genecalling <- function(file, codingfile, techtype = c("sanger", "454", "illumina"), tmpdir, proc) {
  out_prefix <- codingfile
  # mgrast_genecalling.pl chooses training files with highest error rate for given techtype
  if(techtype == "sanger") type <- "sanger_10"
  if(techtype == "454") type <- "454_30"
  if(techtype == "illumina") type <- "illumina_10"
  cmd <- paste("parallel_FragGeneScan.py -v -p", proc, "-s 100 -t", type, "-d", tmpdir, file, out_prefix) 
  # If parallel_FragGeneScan.py gives python error (OSError: [Errno 2] No such file or directory),
  # try this instead:
  ##cmd <- paste("run_FragGeneScan.pl -genome", file, "-out", codingfile, "-complete 0 -train", type) 
  system(cmd)
  # Clean up files
  system(paste("rm", paste0(codingfile, ".out")))
  ##system(paste("rm", file.path(tmpdir, "*")))
}

# Run bowtie2 (host sequence removal) 20221117
# Based on MG-RAST/pipeline/mgcmd/mgrast_bowtie_screen.pl
bowtie2 <- function(file, decontamfile, proc, tmpdir) {
  refdb <- "/home/ARAST/genome-idx/GRCh38_noalt_as/GRCh38_noalt_as"
  # 'reorder' option outputs sequences in same order as input file
  cmd <- paste("bowtie2 -f --reorder -p", proc, "--un", tmpdir, "-x", refdb, "-U", file, "> /dev/null")
  print(cmd)
  system(cmd)
  # Rename output file
  file.rename(file.path(tmpdir, "un-seqs"), decontamfile)
}

# Function to count number of sequences in FASTA or FASTQ file 20231217
count_sequences <- function(file, filetype) {
  cmd <- paste("wc -l", file)
  wc_output <- system(cmd, intern = TRUE)
  nlines <- as.numeric(strsplit(wc_output, " ")[[1]][1])
  if(filetype == "fasta") nseq <- nlines / 2
  if(filetype == "fastq") nseq <- nlines / 4
  nseq
}
