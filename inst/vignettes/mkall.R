# JMDplots/inst/vignettes/mkall.R
# 20200509 This compiles each of the vignettes for compositional analysis of differential expression,
#   which updates the CSV files that are used in the gradH2O and canH2O projects.
# 20200802 Move html files to ../doc
# 20200817 Render vignettes for various basis species (save CSV files for all; HTML for last one only)

files <- dir(pattern = "Rmd")
print(system.time(
  for(basis in c("QCa")) {
    options(basis = basis)
    for(f in files) {
      sep <- paste(rep("=", nchar(f) + 4), collapse = "")
      message()
      print(sep, quote = FALSE)
      print(paste("|", f, "|"), quote = FALSE)
      print(sep, quote = FALSE)
      render(f)
      from <- gsub("Rmd", "html", f)
      to <- file.path("../doc", from)
      file.rename(from, to)
    }
  }
))
