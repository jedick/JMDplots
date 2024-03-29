# prepare data files for the paper
# Average oxidation state of carbon in proteins
# by Jeffrey M. Dick, 2013-11-24

# load CHNOSZ
library(CHNOSZ)

# to make ZC_HUMAN.csv or ZC_membrane.csv
mkZC <- function(which="HUMAN") {
  # read fasta files
  # HUMAN.fasta.gz available at
  # ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/proteomes/HUMAN.fasta.gz
  # membrane.fa.gz generated by combining files in Additional File 2 of Almen et al., 2009
  # http://www.biomedcentral.com/content/supplementary/1741-7007-7-50-S2.zip
  if(which=="HUMAN") lines <- readLines("fasta/human/HUMAN.fasta.gz")
  else if(which=="membrane") lines <- readLines("fasta/human/membrane.fa.gz")
  aa <- canprot::read_fasta("", lines=lines)
  pl <- protein.length(aa)
  pf <- protein.formula(aa)
  ZC <- ZC(pf)
  # get UniProt/IPI accessions
  if(which=="HUMAN") acc <- sapply(strsplit(aa$protein, "|", fixed=TRUE), "[", 2)
  else if(which=="membrane") acc <- aa$protein
  # make data frame and csv file
  out <- data.frame(accession=acc, length=pl, ZC=round(ZC, 4))
  write.csv(out, paste("ZC_", which, ".csv", sep=""), row.names=FALSE, quote=FALSE)
}


# to make Sce.csv, amino acid compositions of yeast proteins
# this becomes inst/extdata/protein/Sce.csv.xz in the CHNOSZ package
mkprot <- function() {
  # files under SGD directory were downloaded from the Saccharomyces Genome Database
  # (http://www.yeastgenome.org)
  # read protein_properties.tab
  pp <- read.table("SGD/protein_properties.tab", sep="\t") 
  # give it column names (adapted from protein_properties.README)
  colnames(pp) <- c("ORF", "SGDID", "MW", "PI", "CAI", "LENGTH", "NTERM", "CTERM", "CODON_BIAS",
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    "FOP", "GRAVY", "AROMATICITY", "TYPE"
  )
  # add a column with gene names
  features <- read.table("SGD/SGD_features.tab", sep="\t", quote="")
  iSGD <- match(pp$SGDID, features[, 1])
  pp <- cbind(pp, data.frame(GENE=features[iSGD, 5]))
  # build the output data frame: ORF, SGDID, GENE, amino acid composition
  pp <- pp[, c(1, 2, 34, 10:29)]
  write.csv(pp, "Sce.csv", row.names=FALSE, quote=FALSE)
}

# to summarize top associations; default aspect is cellular component
assoc <- function(n=30, aspect=c("C")) {
  # read the GO information
  ga <- read.table("SGD/gene_association.sgd.gz", sep="\t", skip=8, quote="", comment.char="")
  # go_terms.tab from SGD has been gzipped locally
  gt <- read.table("SGD/go_terms.tab.gz", sep="\t", quote="")
  # count occurrences of each GO term
  ga.counts <- sort(table(ga$V5), decreasing=TRUE)
  # get numeric GO identifiers
  go.numbers <- as.numeric(substr(names(ga.counts), 4, 10))
  igo <- match(go.numbers, gt$V1)
  # prepare output data frame: sort by ga.counts
  out <- gt[igo, ]
  # replace description with counts
  out$V4 <- ga.counts
  # return only top n for the given aspect
  iaspect <- which(out$V3 %in% aspect)
  return(out[iaspect[1:n], ])
}

# to make SGD_associations.csv
mkassoc <- function() {
  # read the GO information
  ga <- read.table("SGD/gene_association.sgd.gz", sep="\t", skip=8, quote="", comment.char="", as.is=TRUE)
  # remove NOT qualifier
  ga <- ga[!ga$V4=="NOT", ]
  # match SGDID to protein properties
  pp <- read.table("SGD/protein_properties.tab", sep="\t") 
  iSGD <- match(ga$V2, pp$V2)
  # remove NA ones (RNA-coding genes etc)
  ga <- ga[!is.na(iSGD), ]
  # the GO terms of interest
  myterms <- c(
    C = "GO:0005737",  # cytoplasm
    N = "GO:0005634",  # nucleus
    M = "GO:0005739",  # mitochondrion
    I = "GO:0016021",  # integral to membrane
    P = "GO:0005886",  # plasma membrane
    E = "GO:0005783",  # endoplasmic reticulum
    U = "GO:0005730",  # nucleolus
    R = "GO:0005743",  # mitochondrial inner membrane
    D = "GO:0005789",  # endoplasmic reticulum membrane
    G = "GO:0005794",  # Golgi apparatus
    K = "GO:0005935",  # cellular bud neck
    T = "GO:0005741",  # mitochondrial outer membrane
    V = "GO:0005773",  # vacuole
    X = "GO:0005576",  # extracellular region
    O = "GO:0005774",  # vacuolar membrane
    B = "GO:0005934",  # cellular bud tip
    A = "GO:0000139",  # Golgi membrane
    L = "GO:0031965"   # nuclear membrane
  )
  # find all SGDID that contain GO terms of interest
  for(i in 1:length(myterms)) {
    assign(names(myterms)[i], unique(ga$V2[grep(myterms[i], ga$V5)]))
  }
  # membranes: keep only proteins that are integral to membrane 
  for(loc in c("P", "R", "D", "O", "T", "A", "L")) {
    length.orig <- length(get(loc))
    assign(loc, intersect(get(loc), I))
    length.new <- length(get(loc))
    print(paste("location", loc, length.new, "/", length.orig, "(", round(length.new/length.orig, 3), ")"))
  }
  # extracellular region: no change
  # assemble the associations for each SGDID
  SGDID <- sort(unique(c(C, N, M, E, U, G, K, V, B, P, R, D, O, T, A, L, X)))
  associations <- character(length(SGDID))
  for(location in c("C", "N", "M", "E", "U", "G", "K", "V", "B", "P", "R", "D", "O", "T", "A", "L", "X")) {
    iassoc <- match(get(location), SGDID)
    associations[iassoc] <- paste(associations[iassoc], location, sep=",")
  }
  # remove a starting ","
  associations <- substr(associations, 2, 99)
  # write output file
  accession <- pp$V1[match(SGDID, pp$V2)]
  out <- data.frame(accession=accession, SGDID=SGDID, associations=associations)
  write.csv(out, "SGD_associations.csv", row.names=FALSE)
}
