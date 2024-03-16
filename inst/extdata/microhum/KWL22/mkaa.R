# microhum/KWL22/mkaa.R
# Get amino acid compositions from protein sequence files
# 20221028 jmd

# REQUIRED FILES:
# prodigal/*.faa
#   - Download MAG.zip from https://figshare.com/s/a426a12b463758ed6a54
#     - This zip contains 5403 .fa files
#   - Obtain prodigal from https://github.com/hyattpd/Prodigal
#   - Run prodigal on each .fa file:
#     - for i in `ls *.fa`; do prodigal -i $i -o /dev/null -a $i.faa; done
#   - Move the resulting fa.faa files to the prodigal/ directory referenced below  

files <- dir("prodigal", full.names = TRUE)
out <- lapply(1:length(files), function(i) {
  if(i %% 100 == 0) print(i)
  file <- files[i]
  aa <- suppressMessages(canprot::read_fasta(file))
  aa <- canprot::sum_aa(aa)
  # Get run and bin ID from file name
  run_bin <- gsub(".fa.faa.xz", "", gsub("prodigal/", "", file, fixed = TRUE))
  aa$protein <- strsplit(run_bin, "_")[[1]][1]
  aa$organism <- strsplit(run_bin, "_")[[1]][2]
  aa
})

# Save results
out <- do.call(rbind, out)
write.csv(out, "KWL22_MAGs_prodigal_aa.csv", row.names = FALSE, quote = FALSE)
