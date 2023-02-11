# KWL22/mkaa.R
# Get amino acid compositions from protein sequence files
# 20221028 jmd

files <- dir("prodigal", full.names = TRUE)
out <- lapply(1:length(files), function(i) {
  if(i %% 100 == 0) print(i)
  file <- files[i]
  aa <- suppressMessages(CHNOSZ::read.fasta(file))
  aa <- CHNOSZ::aasum(aa)
  # Get run and bin ID from file name
  run_bin <- gsub(".fa.faa.xz", "", gsub("prodigal/", "", file, fixed = TRUE))
  aa$protein <- strsplit(run_bin, "_")[[1]][1]
  aa$organism <- strsplit(run_bin, "_")[[1]][2]
  aa
})

# Save results
out <- do.call(rbind, out)
write.csv(out, "KWL22_MAGs_prodigal_aa.csv", row.names = FALSE, quote = FALSE)
