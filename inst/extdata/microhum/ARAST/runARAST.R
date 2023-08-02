# chem16S/runARAST.R
# Scripts used for running ARAST for analysis of metagenomes 20220503

# This loads the process() function for the processing pipeline
source("ARAST.R")

# Function to calculate amino acid composition from inferred protein sequences
mkAA <- function(faafiles, environment, Control = NULL) {
  # Get Run IDs from file name
  ID <- sapply(strsplit(basename(faafiles), ".fast*"), "[", 1)
  # Create blank data frame for amino acid composition for each sample
  out <- as.data.frame(matrix(ncol = 25, nrow = length(ID)))
  colnames(out) <- c("protein", "organism", "ref", "abbrv", "chains",
    "Ala", "Cys", "Asp", "Glu", "Phe", "Gly", "His", "Ile", "Lys", "Leu",
    "Met", "Asn", "Pro", "Gln", "Arg", "Ser", "Thr", "Val", "Trp", "Tyr"
  )
  # Put in run ID and environment
  out$protein <- ID
  out$organism <- environment
  # Read each file and calculate amino acid composition
  for(i in seq_along(faafiles)) {
    aa <- CHNOSZ::read.fasta(faafiles[i])
    out[i, 5:25] <- colSums(aa[, 5:25])
  }
  outfile <- paste0(environment, "_AA.csv")
  # Identify Control and COVID-19 samples in 'ref' column 20230731
  if(!is.null(Control)) out$ref <- ifelse(ID %in% Control, "Control", "COVID-19")
  write.csv(out, outfile, row.names = FALSE, quote = FALSE)
}

########
## Process nasopharyngeal metagenomes [LLZ+21]
########

# These fasta files have first 1000000 sequences (2000000 lines) of each SRA run
ID <- c("SRR12442171", "SRR12442172", "SRR12442173",
        "SRR12442174", "SRR12442175", "SRR12442176",
        "SRR12442177", "SRR12442178", "SRR12442179"
)
files <- paste0("/home/ARAST/work/", ID, ".fasta")
lapply(files, process, techtype = "illumina")

# The inferred amino acid sequences of protein-coding genes
faafiles <- paste0("/home/ARAST/work/", ID, ".fasta_coding.faa.gz")

# List of control samples (from Figure S1 of Liu et al., 2021 and SRA run info PRJNA656660)
Control <- c("SRR12442173", "SRR12442175", "SRR12442176")

# Generate table of summed amino acid frequencies for each sample
mkAA(faafiles, "LLZ+21", Control)

########
## Process gut metagenomes [ZZL+20]
########

ID <- c("SRR12328886", "SRR12328887", "SRR12328888", "SRR12328889", "SRR12328890", "SRR12328891", 
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
files <- paste0("/home/ARAST/work/", ID, ".fasta")
lapply(files, process, techtype = "illumina")

# The inferred amino acid sequences of protein-coding genes
faafiles <- paste0("/home/ARAST/work/", ID, ".fasta_coding.faa.gz")

# List of control samples (from BioSample metadata PRJNA624223)
Control <- c("SRR12328893", "SRR12328894", "SRR12328895", "SRR12328897", "SRR12328898",
             "SRR12328899", "SRR12328900", "SRR12328901", "SRR12328902", "SRR12328903",
             "SRR12328904", "SRR12328905", "SRR12328906", "SRR12328908", "SRR12328909"
)

# Generate table of summed amino acid frequencies for each sample
mkAA(faafiles, "ZZL+20", Control)

########
## Process oropharyngeal metagenomes [CZH+22]
########

ID <- c("ERR7672491", "ERR7672496", "ERR7672497", "ERR7672501", "ERR7672502", 
        "ERR7672506", "ERR7672511", "ERR7672513", "ERR7672515", "ERR7672519", 
        "ERR7672520", "ERR7672526", "ERR7672532", "ERR7672536", "ERR7672539", 
        "ERR7672540", "ERR7672545", "ERR7672546", "ERR7672548", "ERR7672551", 
        "ERR7672553", "ERR7672559", "ERR7672564", "ERR7672583", "ERR7672586", 
        "ERR7672587", "ERR7672590", "ERR7672593", "ERR7672594", "ERR7672597", 
        "ERR7672601", "ERR7672602", "ERR7672604", "ERR7672605", "ERR7672608", 
        "ERR7672610", "ERR7672611", "ERR7672613", "ERR7672614", "ERR7672619", 
        "ERR7672627", "ERR7672628", "ERR7672629", "ERR7672630", "ERR7672633", 
        "ERR7672634", "ERR7672635", "ERR7672637", "ERR7672638", "ERR7672639", 
        "ERR7672641", "ERR7672642", "ERR7672645", "ERR7672649", "ERR7672653", 
        "ERR7672655", "ERR7672660", "ERR7672662", "ERR7672663", "ERR7672664", 
        "ERR7672666", "ERR7672667", "ERR7672668", "ERR7672670", "ERR7672675", 
        "ERR7672679", "ERR7672681", "ERR7672686", "ERR7672687", "ERR7672689", 
        "ERR7672695", "ERR7672697", "ERR7672698", "ERR7672700", "ERR7672703", 
        "ERR7672704", "ERR7672705", "ERR7672706", "ERR7672708", "ERR7672710", 
        "ERR7672712", "ERR7672715", "ERR7672717", "ERR7672719", "ERR7672722", 
        "ERR7672730", "ERR7672735", "ERR7672736", "ERR7672738", "ERR7672739", 
        "ERR7672740", "ERR7672741", "ERR7672742", "ERR7672746", "ERR7672749", 
        "ERR7672751", "ERR7672762", "ERR7672763", "ERR7672764", "ERR7672765", 
        "ERR7672766", "ERR7672768", "ERR7672769", "ERR7672770", "ERR7672774", 
        "ERR7672776", "ERR7672778", "ERR7672781", "ERR7672782", "ERR7672784", 
        "ERR7672788", "ERR7672790", "ERR7672792", "ERR7672794", "ERR7672796", 
        "ERR7672797", "ERR7672798", "ERR7672800", "ERR7672802", "ERR7672805", 
        "ERR7672807", "ERR7672808"
)
files <- paste0("/home/ARAST/work/", ID, ".fasta")
lapply(files, process, techtype = "illumina")

# The inferred amino acid sequences of protein-coding genes
faafiles <- paste0("/home/ARAST/work/", ID, ".fasta_coding.faa.gz")
# Generate table of summed amino acid frequencies for each sample
mkAA(faafiles, "CZH+22")
