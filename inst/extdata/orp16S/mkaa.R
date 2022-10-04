# chem16S/metaproteome/HWLH22/mkaa.R
# Get peptide sequences for each sample
# 20220830 version 1
# 20221003 Process data for two months

for(month in c("Jun", "Sep")) {

  ## Read source data
  # https://download.iprox.cn/IPX0003299000/IPX0003299001/Bh09LQ_all.pep.xml
  if(month == "Jun") dat <- readLines("Bh09LQ_all.pep.xml")
  # https://download.iprox.cn/IPX0003315000/IPX0003315001/Bk12LQ_all.pep.xml
  if(month == "Sep") dat <- readLines("Bk12LQ_all.pep.xml")

  ## Extract data from XML files
  # Get line numbers of spectrum query results
  iquery <- grep("<spectrum_query", dat, fixed = TRUE)
  # Get sample names
  samples <- substr(gsub('.*spectrum\\=\\"', "", dat[iquery]), 8, 9)
  # Get line numbers of first search hits for each query
  isearch <- iquery + 2
  # Get peptide sequences
  peptides <- gsub('\\".*', "", gsub('.*peptide\\=\\"', "", dat[isearch]))
  # Make sure we have one peptide for each spectrum
  stopifnot(length(samples) == length(peptides))

  ## Sum amino acid composition of peptides for each sample
  usamp <- unique(samples)
  aa <- lapply(usamp, function(sample){
    isamp <- samples == sample
    aa <- lapply(peptides[isamp], seq2aa, protein = "metaproteome")
    aa <- do.call(rbind, aa)
    CHNOSZ::aasum(aa)
  })

  ## Prepare output
  aa <- do.call(rbind, aa)
  aa$organism <- usamp
  aa$ref <- "HWLH22"
  aa$abbrv <- paste0(month, "_18")
  write.csv(aa, paste0("HWLH22_", month, "_2018_aa.csv"), row.names = FALSE, quote = FALSE)

}
