# ARAST.R
# pipeline in R similar to MG-RAST
# (up to RNA and protein genecalling steps)

# 20180310 v0
# 20181117 v1:
#   Used in gradox paper: https://doi.org/10.3389/fmicb.2019.00120
#   Script on Zenodo: https://doi.org/10.5281/zenodo.2314933
# 20211216 v2:
#   Remove list of run IDs and techtypes, and add techtype argument to process()
#   To be used in geo16S paper (preprint at https://doi.org/10.1101/2021.05.31.446500)

# usage: to process a batch of files:
# files <- dir("/home/ARAST/work/", ".fasta", full.names=TRUE)
# system.time(lapply(files, process))
# stop after the dereplication step and don't gzip the output:
# system.time(lapply(files, process, run="DNA", run.gzip=FALSE))

# 'run' argument specified which steps to run:
# DNA: trimming, denoising, artifact removal
# RNA: DNA steps plus RNA genecalling (sortmerna)
# all: DNA and RNA steps plus protein genecalling (FragGeneScan)
# FGS: FragGeneScan only

# SOFTWARE: taken from MG-RAST's "Processing Information and Downloads" for ERX242607_sra_data.fasta
# download skewer_0.2.2 from https://github.com/relipmoc/skewer
# download bowtie_v2.3.4.1 from https://github.com/BenLangmead/bowtie2 [used by autoskewer]
# download autoskewer (ZIP) from https://github.com/MG-RAST/autoskewer
# download ea-utils (ZIP) from https://github.com/ExpressionAnalysis/ea-utils [fastq-mcf, used by filter_sequences]
# install biopython [used by seq_length_stats.py]
#   sbopkg: upgrade pip, then run pip install biopython
# download pipeline (ZIP) from https://github.com/MG-RAST/pipeline
#   [seq_length_stats.py, filter_sequences, dereplication.py, parallel_FragGeneScan.py]
# download sortmerna (ZIP) from https://github.com/biocore/sortmerna
# download FragGeneScan_FGS1.18 from https://github.com/MG-RAST/FGS

# DATABASE:
# download md5rna.clust (pipeline-master/CWL/Inputs/DBs/getpredata.sh) [used by sortmerna]
#   curl "http://shock.metagenomics.anl.gov/node/c4c76c22-297b-4404-af5c-8cd98e580f2a?download" > md5rna.clust
# set up index for sortmerna:
#   indexdb --ref md5rna.clust,md5rna.clust.index

# function to process a single file 20180310
# techtype: sequencing technology (needed for FragGeneScan) can be sanger, 454, or illumina
process <- function(file, proc=parallel::detectCores(), rna_nr="/home/ARAST/REFDB/md5rna.clust", run="all", run.gzip=TRUE, techtype = "illumina") {
  ## what steps to run
  if(run=="all") steps <- c(0:9, 23:24)
  if(run=="DNA") steps <- c(0:4, 23:24)
  if(run=="RNA") steps <- c(0:6, 23:24)
  if(run=="FGS") {
    steps <- c(0, 9, 23:24)
    no_rRNAfile <- file
  }
  if(run=="afterDNA") {
    steps <- c(0, 5:9, 23:24)
    derepfile <- file
  }
  if(0 %in% steps) {
    print("######## 0. Set up environment ########")
    # get workdir and ID from file name
    workdir <- dirname(file)
    ID <- strsplit(basename(file), "_")[[1]][1]
    # file type: fasta or fastq
    if(grepl("\\.fasta$", file)) filetype <- "fasta"
    else if(grepl("\\.fna$", file)) filetype <- "fasta"
    else if(grepl("\\.fastq$", file)) filetype <- "fastq"
    else {
      print(paste("not processing", file, "which is not fasta or fastq"))
      return()
    }
    print(paste0("processing ", filetype, " file for ", ID, " (", techtype, " sequencing)"))
    tmpdir <- file.path(workdir, "tmp")
    if(!dir.exists(tmpdir)) dir.create(tmpdir)
    # keep list of temporary (to remove) and output (to gzip) files
    tmpfiles <- outfiles <- character()
  }
  # NOTE: 00000000 steps are omitted from/do not follow MG-RAST pipeline
  if(1 %in% steps) print("00000000 1. Initial sequence statistics 00000000")
  if(2 %in% steps) {
    print("######## 2. Adapter Trimming ########")
    scrubbedfile <- file.path(workdir, paste0(ID, "_scrubbed.", filetype))
    autoskewer(file, scrubbedfile, tmpdir)
    tmpfiles <- c(tmpfiles, scrubbedfile)
  }
  if(3 %in% steps) {
    print("######## 3. Denoising and normalization ########")
    passedfile <- file.path(workdir, paste0(ID, "_passed.fasta"))
    removedfile <- file.path(workdir, paste0(ID, "_removed.fasta"))
    preprocess(scrubbedfile, passedfile, removedfile, filetype)
    tmpfiles <- c(tmpfiles, passedfile, removedfile)
  }
  if(4 %in% steps) {
    print("######## 4. Removal of sequencing artifacts ########")
    derepfile <- file.path(workdir, paste0(ID, "_dereplicated.fasta"))
    dereplication(passedfile, derepfile, tmpdir)
    outfiles <- c(outfiles, derepfile)
  }
  if(5 %in% steps) print("00000000 5. Host DNA contamination removal 00000000")
  if(6 %in% steps) {
    print("######## 6. RNA feature identification (rRNA genecalling) ########")
    index <- paste0(rna_nr, ".index")
    rRNAfile <- file.path(workdir, paste0(ID, "_rRNA.fasta"))
    no_rRNAfile <- file.path(workdir, paste0(ID, "_no-rRNA.fasta"))
    # adjust memory use based on number of processors
    if(proc > 4) mem <- 6144 else mem <- 2048 
    sortmerna(derepfile, rRNAfile, no_rRNAfile, proc, mem, rna_nr, index)
    outfiles <- c(outfiles, rRNAfile, no_rRNAfile)
  }
  if(7 %in% steps) print("00000000 7. RNA clustering 00000000")
  if(8 %in% steps) print("00000000 8. RNA similarity search 00000000")
  if(9 %in% steps) {
    print("######## 9. Identify putative protein coding features (genecalling) ########")
    codingfile <- file.path(workdir, paste0(ID, "_coding"))
    genecalling(no_rRNAfile, codingfile, techtype, tmpdir, proc)
    outfiles <- c(outfiles, paste0(codingfile, ".fna"), paste0(codingfile, ".faa"))
  }
  if(23 %in% steps) {
    print("00000000 23. Summary statistics 00000000")
    # gzip output files and remove intermediate files
    system(paste("rm -f", paste(tmpfiles, collapse=" ")))
    if(run.gzip) system(paste("gzip -f", paste(outfiles, collapse=" ")))
  }
  if(24 %in% steps) print("######## 24. Completed ########")
}

# run autoskewer.py 20180310
autoskewer <- function(file, scrubbedfile, tmpdir) {
  cmd <- paste("autoskewer.py -v -i", file, "-o", scrubbedfile, "-l ~/tmp/scrubbed.log -t", tmpdir)
  system(cmd)
  # clean up tmpdir
  tmpfiles <- file.path(tmpdir, "*")
  cmd <- paste("rm", tmpfiles)
  system(cmd)
}

# run filter_sequences 20180309
# based on MG-RAST/pipeline/mgcmd/mgrast_preprocess.pl
preprocess <- function(file, passedfile, removedfile, filetype) {
  # print stats
  cmd <- paste("seq_length_stats.py -i", file, "-t", filetype, "-f")
  system(cmd)
  if(filetype=="fasta") {
    # get stats
    cmd <- paste("seq_length_stats.py -i", file, "-t", filetype, "-f | cut -f2")
    out <- system(cmd, intern=TRUE)
    # this gives bp_count, sequence_count, average_length, standard_deviation_length, length_min, length_max
    out <- as.numeric(out)
    # translated from mgrast_preprocess.pl
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

# run dereplication.py 20180310
# based on MG-RAST/pipeline/mgcmd/mgrast_dereplicate.pl
dereplication <- function(file, derepfile, tmpdir) {
  prefix_size <- 50
  memory <- 8  # GB
  run_dir <- tmpdir
  input <- file
  out_prefix <- tmpdir
  cmd <- paste("dereplication.py -l", prefix_size, "-m", memory, "-d", run_dir, input, out_prefix)
  system(cmd)
  # rename and clean up files
  system(paste("mv", paste0(tmpdir, ".passed.fna"), derepfile))
  system(paste("rm", paste0(tmpdir, ".*")))
}

# run sortmerna 20180310
# based on MG-RAST/pipeline/mgcmd/mgrast_sortme_rna.pl
sortmerna <- function(file, rRNAfile, no_rRNAfile, proc, mem, rna_nr, index) {
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
  #             "--aligned", rRNAfile, "--other", no_rRNAfile, "-v")
  # 20180311: use a modified version of mgrast_sortme_rna.pl to save both RNA and non-RNA sequences
  cmd <- paste("./arast_sortme_rna.pl -input", file, "-output", rRNAfile, "-notrna", no_rRNAfile,
               "-rna_nr", rna_nr, "-index", index, "-proc", proc, "-mem", mem, "-eval", eval, "-ident", ident)
  system(cmd)
  # clean up files
  WD <- setwd(dirname(file))
  system("rm *.blast *.tab *.id *.sort")
  setwd(WD)
}

# run FragGeneScan 20180309
# based on MG-RAST/pipeline/mgcmd/mgrast_genecalling.pl
genecalling <- function(file, codingfile, techtype=c("sanger", "454", "illumina"), tmpdir, proc) {
  out_prefix <- codingfile
  # sequencing technology: mgrast_genecalling.pl chooses training files with highest error rate
  if(techtype=="sanger") type <- "sanger_10"
  if(techtype=="454") type <- "454_30"
  if(techtype=="illumina") type <- "illumina_10"
  cmd <- paste("parallel_FragGeneScan.py -v -p", proc, "-s 100 -t", type, "-d", tmpdir, file, out_prefix) 
  # sometimes(?) that gives python error:
  # OSError: [Errno 2] No such file or directory
  #cmd <- paste("run_FragGeneScan.pl -genome", file, "-out", codingfile, "-complete 0 -train", type) 
  system(cmd)
  # clean up files
  system(paste("rm", paste0(codingfile, ".out")))
  system(paste("rm", file.path(tmpdir, "*")))
}
