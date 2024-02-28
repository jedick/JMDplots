# microhum/runARAST.R
# Script to run ARAST for particular metagenomes 20220503

# Usage:
# runARAST("LLZ+21") # Process nasopharyngeal metagenomes
# runARAST("ZZL+20") # Process gut metagenomes
# runARAST("CZH+22") # Process oropharyngeal metagenomes
# runARAST("HMP12")  # Process Human Microbiome Project metagenomes
# runARAST("HMP12", screening = FALSE) # Process HMP metagenomes without screening step

# This loads the ARAST() function for the processing pipeline
source("ARAST.R")
# Set working directory
workdir <- "/home/ARAST/work"

runARAST <- function(dataset, screening = TRUE, cleanup = TRUE) {

  # Set techtype for all datasets
  techtype <- "illumina"
  # Start with empty list of SRA Run IDs and sample groups
  Run <- Control <- COVID <- IBD <- UC <- CD <- NULL

  if(dataset == "LLZ+21") {
    # Nasopharyngeal metagenomes 
    Run <- c("SRR12442171", "SRR12442172", "SRR12442173",
            "SRR12442174", "SRR12442175", "SRR12442176",
            "SRR12442177", "SRR12442178", "SRR12442179"
    )
    # List of control samples (from Figure S1 of Liu et al., 2021 and SRA run info PRJNA656660)
    Control <- c("SRR12442173", "SRR12442175", "SRR12442176")
    COVID <- setdiff(Run, Control)
  }

  if(dataset == "ZZL+20") {
    # Gut metagenomes
    Run <- c(
      "SRR12328886", "SRR12328887", "SRR12328888", "SRR12328889", "SRR12328890", "SRR12328891", 
      "SRR12328892", "SRR12328893", "SRR12328894", "SRR12328895", "SRR12328896", "SRR12328897", 
      "SRR12328898", "SRR12328899", "SRR12328900", "SRR12328901", "SRR12328902", "SRR12328903", 
      "SRR12328904", "SRR12328905", "SRR12328906", "SRR12328907", "SRR12328908", "SRR12328909", 
      "SRR12328910", "SRR12328911", "SRR12328912", "SRR12328913", "SRR12328914", "SRR12328915", 
      "SRR12328916", "SRR12328917", "SRR12328918", "SRR12328919", "SRR12328920", "SRR12328921", 
      "SRR12328922", "SRR12328923", "SRR12328924", "SRR12328925", "SRR12328926", "SRR12328927", 
      "SRR12328928", "SRR12328929", "SRR12328930", "SRR12328931", "SRR12328932", "SRR12328933", 
      "SRR12328934", "SRR12328935", "SRR12328936", "SRR12328937", "SRR12328938", "SRR12328939", 
      "SRR12328940", "SRR12328941", "SRR12328942", "SRR12328943", "SRR12328944", "SRR12328945", 
      "SRR12328946", "SRR12328947", "SRR12328948", "SRR12328949", "SRR12328950", "SRR12328951", 
      "SRR12328952", "SRR12328953", "SRR12328954", "SRR12328955", "SRR12328956", "SRR12328957", 
      "SRR12328958"
    )
    # List of control samples (from BioSample metadata PRJNA624223)
    Control <- c(
       "SRR12328893", "SRR12328894", "SRR12328895", "SRR12328897", "SRR12328898",
       "SRR12328899", "SRR12328900", "SRR12328901", "SRR12328902", "SRR12328903",
       "SRR12328904", "SRR12328905", "SRR12328906", "SRR12328908", "SRR12328909"
    )
    COVID <- setdiff(Run, Control)
  }

  if(dataset == "CZH+22") {
    # Oropharyngeal metagenomes
    Run <- c(
      "ERR7672491", "ERR7672496", "ERR7672497", "ERR7672501", "ERR7672502", "ERR7672506", "ERR7672511", "ERR7672513", "ERR7672515", "ERR7672519", 
      "ERR7672520", "ERR7672526", "ERR7672532", "ERR7672536", "ERR7672539", "ERR7672540", "ERR7672545", "ERR7672546", "ERR7672548", "ERR7672551", 
      "ERR7672553", "ERR7672559", "ERR7672564", "ERR7672583", "ERR7672586", "ERR7672587", "ERR7672590", "ERR7672593", "ERR7672594", "ERR7672597", 
      "ERR7672601", "ERR7672602", "ERR7672604", "ERR7672605", "ERR7672608", "ERR7672610", "ERR7672611", "ERR7672613", "ERR7672614", "ERR7672619", 
      "ERR7672627", "ERR7672628", "ERR7672629", "ERR7672630", "ERR7672633", "ERR7672634", "ERR7672635", "ERR7672637", "ERR7672638", "ERR7672639", 
      "ERR7672641", "ERR7672642", "ERR7672645", "ERR7672649", "ERR7672653", "ERR7672655", "ERR7672660", "ERR7672662", "ERR7672663", "ERR7672664", 
      "ERR7672666", "ERR7672667", "ERR7672668", "ERR7672670", "ERR7672675", "ERR7672679", "ERR7672681", "ERR7672686", "ERR7672687", "ERR7672689", 
      "ERR7672695", "ERR7672697", "ERR7672698", "ERR7672700", "ERR7672703", "ERR7672704", "ERR7672705", "ERR7672706", "ERR7672708", "ERR7672710", 
      "ERR7672712", "ERR7672715", "ERR7672717", "ERR7672719", "ERR7672722", "ERR7672730", "ERR7672735", "ERR7672736", "ERR7672738", "ERR7672739", 
      "ERR7672740", "ERR7672741", "ERR7672742", "ERR7672746", "ERR7672749", "ERR7672751", "ERR7672762", "ERR7672763", "ERR7672764", "ERR7672765", 
      "ERR7672766", "ERR7672768", "ERR7672769", "ERR7672770", "ERR7672774", "ERR7672776", "ERR7672778", "ERR7672781", "ERR7672782", "ERR7672784", 
      "ERR7672788", "ERR7672790", "ERR7672792", "ERR7672794", "ERR7672796", "ERR7672797", "ERR7672798", "ERR7672800", "ERR7672802", "ERR7672805", 
      "ERR7672807", "ERR7672808"
    )
    # Couldn't find Control and COVID samples in metadata
  }

  if(dataset == "HMP12") {
    # Various body sites (HMP Consortium, 2012)
    Run <- c(
      # Samples used in Tax4Fun paper
      "SRR059889", "SRR060375", "SRR060433", "SRR060445", "SRR061143", "SRR061144", "SRR061153", "SRR061157", "SRR061164", "SRR061168", 
      "SRR061169", "SRR061171", "SRR061182", "SRR061183", "SRR061191", "SRR061196", "SRR061198", "SRR061202", "SRR061205", "SRR061222", 
      "SRR061226", "SRR061233", "SRR061241", "SRR061255", "SRR061263", "SRR061274", "SRR061326", "SRR061494", "SRR061496", "SRR061557", 
      "SRR061910", "SRR061940", "SRR061960", "SRR061964", "SRR061968", "SRR062004", "SRR062010", "SRR062040", "SRR062063", "SRR062091", 
      "SRR062356", "SRR062364", "SRR062386", "SRR062434", "SRR062506", "SRR062519", "SRR062538", "SRR063917", "SRR346668",
      # Additional samples analyzed by Dick and Tan (2023)
      "SRR059348", "SRR059518", "SRR059872", "SRR060382", "SRR060414", "SRR060436", "SRR061141", "SRR061189", "SRR061210", "SRR061217", 
      "SRR061225", "SRR061239", "SRR061248", "SRR061261", "SRR061264", "SRR061473", "SRR061913", "SRR061921", "SRR061923", "SRR061925", 
      "SRR061927", "SRR061943", "SRR061947", "SRR061965", "SRR062005", "SRR062007", "SRR062015", "SRR062017", "SRR062019", "SRR062025", 
      "SRR062060", "SRR062074", "SRR062076", "SRR062094", "SRR062106", "SRR062335", "SRR062399", "SRR062419", "SRR062429", "SRR062440", 
      "SRR062442", "SRR062458", "SRR062461", "SRR062522", "SRR063912", "SRR063925", "SRR346655", "SRR346690", "SRR346709", "SRR346717", 
      "SRR353633", "SRR514852"
    )
  }

  if(dataset == "LAA+19") {
    # Gut metagenomes (Lloyd-Price et al., 2019)
    Run <- c(
      "SRR5962888", "SRR5962886", "SRR5950683", "SRR5962905", "SRR5962899", "SRR5962887", "SRR5950679", "SRR5950681", "SRR5962892", "SRR5962883", 
      "SRR5962893", "SRR5962891", "SRR5950675", "SRR5962882", "SRR5950687", "SRR5950677", "SRR5962907", "SRR5962906", "SRR5950686", "SRR5962898", 
      "SRR5962901", "SRR5962900", "SRR5962908", "SRR5950690", "SRR5962902", "SRR5950693", "SRR5962897", "SRR5962890", "SRR5950685", "SRR5962904", 
      "SRR5950682", "SRR5950722", "SRR5950714", "SRR5950750", "SRR5950737", "SRR5950710", "SRR5950719", "SRR5950715", "SRR5950729", "SRR5950723", 
      "SRR5950731", "SRR5950751", "SRR5950738", "SRR5950743", "SRR5950708", "SRR5950733", "SRR5950746", "SRR5950705", "SRR5950706", "SRR5950700", 
      "SRR5950713", "SRR5950712", "SRR5950709", "SRR5950762", "SRR5950771", "SRR5950766", "SRR5950790", "SRR5950763", "SRR5950775", "SRR5950760", 
      "SRR5950756", "SRR5950493", "SRR5950786", "SRR5950765", "SRR5950753", "SRR5950781", "SRR5950764", "SRR5950793", "SRR5950497", "SRR5950498", 
      "SRR5950759", "SRR5950783", "SRR5950777", "SRR5950495", "SRR5950543", "SRR5950539", "SRR5950551", "SRR5950536", "SRR5950523", "SRR5950553", 
      "SRR5950533", "SRR5950510", "SRR5950535", "SRR5950537", "SRR5950513", "SRR5950552", "SRR5950520", "SRR5950531", "SRR5950538", "SRR5950541", 
      "SRR5950506", "SRR5950519", "SRR5950518", "SRR5950524", "SRR5950548", "SRR5950560", "SRR5950511", "SRR5950550", "SRR5950525", "SRR5950547", 
      "SRR5950517", "SRR5950614", "SRR5950567", "SRR5950622", "SRR5950589", "SRR5950573", "SRR5950597", "SRR5950571", "SRR5950593", "SRR5950594", 
      "SRR5950615", "SRR5950586", "SRR5950579", "SRR5950570", "SRR5950607", "SRR5950642", "SRR5950635", "SRR5950637", "SRR5950650", "SRR5950666", 
      "SRR5950625", "SRR5950646", "SRR5950638", "SRR5950671", "SRR5950629", "SRR5950657", "SRR5950647", "SRR5950652", "SRR5950651", "SRR5950632", 
      "SRR5950664", "SRR5950667", "SRR5950639", "SRR5950641", "SRR5950624", "SRR5950621", "SRR5950631", "SRR5950661", "SRR5950663"
    )
    UC <- c(
      "SRR5962888", "SRR5962886", "SRR5950683", "SRR5950679", "SRR5962907", "SRR5962898", "SRR5962900", "SRR5962908", "SRR5950685", "SRR5950710", 
      "SRR5950723", "SRR5950751", "SRR5950738", "SRR5950743", "SRR5950746", "SRR5950705", "SRR5950700", "SRR5950712", "SRR5950709", "SRR5950762", 
      "SRR5950771", "SRR5950756", "SRR5950495", "SRR5950551", "SRR5950536", "SRR5950523", "SRR5950510", "SRR5950513", "SRR5950538", "SRR5950560", 
      "SRR5950597", "SRR5950615", "SRR5950570", "SRR5950637", "SRR5950650", "SRR5950671", "SRR5950629", "SRR5950657", "SRR5950647", "SRR5950651", 
      "SRR5950641", "SRR5950624", "SRR5950663"
    )
    CD <- c(
      "SRR5962905", "SRR5962899", "SRR5962887", "SRR5950681", "SRR5962892", "SRR5962891", "SRR5950675", "SRR5950687", "SRR5950677", "SRR5962901", 
      "SRR5962897", "SRR5962890", "SRR5950682", "SRR5950722", "SRR5950714", "SRR5950750", "SRR5950737", "SRR5950731", "SRR5950713", "SRR5950766", 
      "SRR5950790", "SRR5950763", "SRR5950775", "SRR5950760", "SRR5950493", "SRR5950786", "SRR5950765", "SRR5950781", "SRR5950764", "SRR5950793", 
      "SRR5950497", "SRR5950539", "SRR5950553", "SRR5950533", "SRR5950535", "SRR5950537", "SRR5950552", "SRR5950541", "SRR5950548", "SRR5950525", 
      "SRR5950517", "SRR5950614", "SRR5950567", "SRR5950622", "SRR5950589", "SRR5950573", "SRR5950586", "SRR5950579", "SRR5950607", "SRR5950635", 
      "SRR5950666", "SRR5950625", "SRR5950632", "SRR5950667", "SRR5950639", "SRR5950631", "SRR5950661"
    )
    Control <- setdiff(Run, c(UC, CD))
  }

  if(is.null(Run)) stop("invalid dataset")
  # Setting for ARAST
  if(screening) {
    mode <- "protein" 
    aafile <- paste0(dataset, "_aa.csv")
    statsfile <- paste0(dataset, "_stats.csv")
  } else {
    mode <- "protein_no_screening"
    aafile <- paste0(dataset, "_no_screening_aa.csv")
    statsfile <- paste0(dataset, "_no_screening_stats.csv")
  }
  # Names of files to keep processing statistics for each run
  statfiles <- file.path(workdir, paste0(Run, "_stats.csv"))

  # Loop over Runs
  for(i in 1:length(Run)) {
    # If statistics file already exists then assume that this run has been processed
    if(file.exists(statfiles[i])) {
      print(paste("*** Skipping", Run[i], "because statistics file exists ***"))
      next
    }
    # Dump FASTQ files
    # NOTE: If all SRA files have been prefetch'ed, run vdb-config -i and turn off
    # "Enable Remote Access" to prevent fastq-dump from dying when there is no internet
    cmd <- paste("fastq-dump -O", workdir, "--skip-technical --clip --split-3 --read-filter pass", Run[i])
    print(cmd)
    fastq_dump.log <- system(cmd, intern = TRUE)
    # Save output of fastq_dump for putting together processing statistics in ARAST()
    writeLines(fastq_dump.log, file.path(workdir, paste0(Run[i], "_fastq-dump.log")))
    # If fastq_dump produced no output, then it probably didn't work, so skip this run
    if(length(fastq_dump.log) == 0) {
      print(paste("*** Skipping", Run[i], "because fastq_dump produced no output ***"))
      next
    }
    # Use file with forward reads ...
    inputfile <- file.path(workdir, paste0(Run[i], "_pass_1.fastq"))
    # ... or unsplit file, if only that is available
    if(!file.exists(inputfile)) inputfile <- file.path(workdir, paste0(Run[i], "_pass.fastq"))
    if(!file.exists(inputfile)) {
      print(paste("*** Skipping", Run[i], "because no input file could be found ***"))
      next
    }
    # Run ARAST pipeline on 128 GB workstation
    ARAST(inputfile, mode = mode, techtype = techtype, mem_MB = 30000, mem_GB = 100, cleanup = cleanup)
    ## Or use this command for 16 GB laptop
    #ARAST(inputfile, mode = mode, techtype = techtype, mem_MB = 6144, mem_GB = 16, cleanup = cleanup)
    # Remove FASTQ files
    if(cleanup) file.remove(dir(workdir, "*.fastq", full.names = TRUE))
  }

  # Create blank data frame to hold amino acid compositions of proteins in each sample
  aa <- as.data.frame(matrix(ncol = 25, nrow = length(Run)))
  colnames(aa) <- c("protein", "organism", "ref", "abbrv", "chains",
    "Ala", "Cys", "Asp", "Glu", "Phe", "Gly", "His", "Ile", "Lys", "Leu",
    "Met", "Asn", "Pro", "Gln", "Arg", "Ser", "Thr", "Val", "Trp", "Tyr"
  )
  # Put in Run IDs and dataset name
  aa$protein <- Run
  aa$organism <- dataset
  # Loop over files with the amino acid sequences of inferred protein-coding genes
  faafiles <- file.path(workdir, paste0(Run, "_coding.faa.gz"))
  # Only process files that exist
  iexist <- file.exists(faafiles)
  for(i in seq_along(faafiles)[iexist]) {
    ## Slow way: Calculate amino acid composition for each sequence
    #myaa <- canprot::read.fasta(faafiles[i])
    #aa[i, 5:25] <- colSums(aa[, 5:25])
    # Faster way: count total number of each amino acid in file 20231217
    # https://unix.stackexchange.com/questions/5010/how-can-i-count-the-number-of-different-characters-in-a-file
    cmd <- paste('zcat', faafiles[i], '| grep -v "^>" | fold -b -w1 | sort --buffer-size=8G | uniq -c')
    print(cmd)
    uniq_output <- system(cmd, intern = TRUE)
    # Find amino acids in output of uniq command
    iAA <- sapply(c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"), grep, uniq_output)
    uniq_output_iAA <- uniq_output[iAA]
    # Convert output of uniq to numeric counts
    counts <- as.numeric(gsub("[A-Z]", "", gsub(" ", "", uniq_output_iAA)))
    aa[i, 6:25] <- counts
    # Get number of protein sequences
    cmd <- paste("zcat", faafiles[i], "| wc -l")
    wc_output <- system(cmd, intern = TRUE)
    nlines <- as.numeric(strsplit(wc_output, " ")[[1]][1])
    nseq <- nlines / 2
    aa$chains[i] <- nseq
  }

  # Label sample groups in 'abbrv' column 20230731
  aa$abbrv[Run %in% Control] <- "Control"
  aa$abbrv[Run %in% COVID] <- "COVID-19"
  aa$abbrv[Run %in% IBD] <- "IBD"
  aa$abbrv[Run %in% UC] <- "UC"
  aa$abbrv[Run %in% CD] <- "CD"

  # Gather statistics for all runs 20231217
  stats <- do.call(rbind, lapply(statfiles, function(file) read.csv(file)))
  # Put number of input sequences into 'ref' column 20231219
  aa$ref <- stats$input_sequences

  # Write output files
  write.csv(aa, aafile, row.names = FALSE, quote = FALSE)
  write.csv(stats, statsfile, row.names = FALSE, quote = FALSE)

}
