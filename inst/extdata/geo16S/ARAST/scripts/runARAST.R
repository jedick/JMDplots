# geo16S/runARAST.R
# Scripts used for running ARAST for analysis of metagenomes 20211216

# This loads the process() function for the processing pipeline
source("ARAST.R")

# Function to calculate amino acid composition from inferred protein sequences
mkAA <- function(faafiles, environment) {
  # Get Run IDs from file name
  ID <- sapply(strsplit(basename(faafiles), ".fast*"), "[", 1)
  # Create blank data frame for amino acid composition for each sample
  out <- as.data.frame(matrix(ncol = 25, nrow = length(ID)))
  colnames(out) <- c("protein", "organism", "ref", "abbrv", "chains",
    "Ala", "Cys", "Asp", "Glu", "Phe", "Gly", "His", "Ile", "Lys", "Leu",
    "Met", "Asn", "Pro", "Gln", "Arg", "Ser", "Thr", "Val", "Trp", "Tyr"
  )
  # Put in run ID and environment (HMP, Gut, or Soils)
  out$protein <- ID
  out$organism <- environment
  # Read each file and calculate amino acid composition
  for(i in seq_along(faafiles)) {
    aa <- CHNOSZ::read.fasta(faafiles[i])
    out[i, 5:25] <- colSums(aa[, 5:25])
  }
  outfile <- paste0(environment, "_AA.csv")
  write.csv(out, outfile, row.names = FALSE, quote = FALSE)
}

########
## Process HMP metagenomes (HMP Consortium, 2012)
########

# These fasta files initially have first 500000 sequences (1000000 lines) of each SRA run
# (Later iterations use more sequences for specific runs - see Table S3)
ID <- c(
  # Samples used in Tax4Fun paper
  "SRR059889", "SRR060375", "SRR060433", "SRR060445", "SRR061143", 
  "SRR061144", "SRR061153", "SRR061157", "SRR061164", "SRR061168", 
  "SRR061169", "SRR061171", "SRR061182", "SRR061183", "SRR061191", 
  "SRR061196", "SRR061198", "SRR061202", "SRR061205", "SRR061222", 
  "SRR061226", "SRR061233", "SRR061241", "SRR061255", "SRR061263", 
  "SRR061274", "SRR061326", "SRR061494", "SRR061496", "SRR061557", 
  "SRR061910", "SRR061940", "SRR061960", "SRR061964", "SRR061968", 
  "SRR062004", "SRR062010", "SRR062040", "SRR062063", "SRR062091", 
  "SRR062356", "SRR062364", "SRR062386", "SRR062434", "SRR062506", 
  "SRR062519", "SRR062538", "SRR063917", "SRR346668",
  # Additional samples
  "SRR059348", "SRR059518", "SRR059872", "SRR060382", "SRR060414", 
  "SRR060436", "SRR061141", "SRR061189", "SRR061210", "SRR061217", 
  "SRR061225", "SRR061239", "SRR061248", "SRR061261", "SRR061264", 
  "SRR061473", "SRR061913", "SRR061921", "SRR061923", "SRR061925", 
  "SRR061927", "SRR061943", "SRR061947", "SRR061965", "SRR062005", 
  "SRR062007", "SRR062015", "SRR062017", "SRR062019", "SRR062025", 
  "SRR062060", "SRR062074", "SRR062076", "SRR062094", "SRR062106", 
  "SRR062335", "SRR062399", "SRR062419", "SRR062429", "SRR062440", 
  "SRR062442", "SRR062458", "SRR062461", "SRR062522", "SRR063912", 
  "SRR063925", "SRR346655", "SRR346690", "SRR346709", "SRR346717", 
  "SRR353633", "SRR514852"
)
files <- paste0("/home/ARAST/work/", ID, ".fasta")
lapply(files, process, techtype = "illumina")

# The inferred amino acid sequences of protein-coding genes
faafiles <- paste0("/home/ARAST/work/", ID, ".fasta_coding.faa.gz")
# Generate table of summed amino acid frequencies for each sample
mkAA(faafiles, "HMP")

########
## Process Gut metagenomes (Muegge et al., 2011)
########

# These are the complete fasta files downloaded from MG-RAST
ID <- c("mgm4461284.3", "mgm4461285.3", "mgm4461286.3", "mgm4461287.3",
"mgm4461288.3", "mgm4461289.3", "mgm4461290.3", "mgm4461291.3",
"mgm4461292.3", "mgm4461293.3", "mgm4461294.3", "mgm4461295.3",
"mgm4461296.3", "mgm4461297.3", "mgm4461298.3", "mgm4461299.3",
"mgm4461300.3", "mgm4461301.3", "mgm4461341.3", "mgm4461342.3",
"mgm4461343.3", "mgm4461344.3", "mgm4461345.3", "mgm4461346.3",
"mgm4461347.3", "mgm4461348.3", "mgm4461349.3", "mgm4461350.3",
"mgm4461351.3", "mgm4461352.3", "mgm4461353.3", "mgm4461354.3",
"mgm4461355.3", "mgm4461357.3", "mgm4461358.3", "mgm4461360.3",
"mgm4461361.3", "mgm4461362.3", "mgm4461363.3", "mgm4461364.3",
"mgm4461365.3", "mgm4461366.3", "mgm4461367.3", "mgm4461368.3",
"mgm4461369.3", "mgm4461370.3", "mgm4461371.3", "mgm4461372.3",
"mgm4461374.3", "mgm4461375.3", "mgm4461376.3", "mgm4461377.3",
"mgm4461378.3", "mgm4461379.3", "mgm4461380.3", "mgm4461383.3"
)
files <- paste0("/home/ARAST/work/", ID, ".fasta")
lapply(files, process, techtype = "454")

faafiles <- paste0("/home/ARAST/work/", ID, ".fasta_coding.faa.gz")
mkAA(faafiles, "Guts")

########
## Process Soil metagenomes (Fierer et al., 2012)
########

# These fastq files have first 100 MB (1024*1024*100 bytes) of each run
ID <- c("mgm4477803.3", "mgm4477804.3", "mgm4477805.3", "mgm4477807.3", 
"mgm4477872.3", "mgm4477873.3", "mgm4477874.3", "mgm4477875.3", 
"mgm4477876.3", "mgm4477877.3", "mgm4477899.3", "mgm4477900.3", "mgm4477901.3",
"mgm4477902.3", "mgm4477903.3", "mgm4477904.3")
files <- paste0("/home/ARAST/work/", ID, ".fastq")
lapply(files, process, techtype = "illumina")

faafiles <- paste0("/home/ARAST/work/", ID, ".fastq_coding.faa.gz")
mkAA(faafiles, "Soils")

########
## Process Marcellus metagenomes (Daly et al., 2016)
########

# These fastq files have first 500000 sequences (2000000 lines) of each SRA run
# Time points: input, T7, T13, T82, T328
ID <- c("SRR3111417", "SRR3111625", "SRR3111724", "SRR3111729", "SRR3111737")

files <- paste0("/home/ARAST/work/", ID, ".fastq")
lapply(files, process, techtype = "illumina")

faafiles <- paste0("/home/ARAST/work/", ID, ".fastq_coding.faa.gz")
mkAA(faafiles, "Marcellus_Shale")

########
## Process Manus Basin metagenomes (Meier et al., 2017)
########

# These fastq files have first 500000 sequences (2000000 lines) of each SRA run
# Samples: NSu-F2b, NSu-F5, Fw-F1b, Fw-F3, RR-F1b
ID <- c("ERR1679394", "ERR1679395", "ERR1679397", "ERR1679396", "ERR1679398")

files <- paste0("/home/ARAST/work/", ID, ".fastq")
lapply(files, process, techtype = "illumina")

faafiles <- paste0("/home/ARAST/work/", ID, ".fastq_coding.faa.gz")
mkAA(faafiles, "Manus_Basin")

########
## Process Black Sea metagenomes (Villanueva et al., 2021)
########

# These fastq files have first 500000 sequences (2000000 lines) of each SRA run
# Depths: c(50, 70, 80, 85, 90, 95, 100, 105, 110, 130, 170, 250, 500, 1000, 2000)
ID <- c("SRR12347146", "SRR12347145", "SRR12347139", "SRR12347138", "SRR12347137",
        "SRR12347136", "SRR12347135", "SRR12347134", "SRR12347133", "SRR12347132",
        "SRR12347144", "SRR12347143", "SRR12347142", "SRR12347141", "SRR12347140")

files <- paste0("/home/ARAST/work/", ID, ".fastq")
lapply(files, process, techtype = "illumina")

faafiles <- paste0("/home/ARAST/work/", ID, ".fastq_coding.faa.gz")
mkAA(faafiles, "Black_Sea")
