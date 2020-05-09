# JMDplots/inst/vignettes/mkall.R
# This compiles each of the vignettes for compositional analysis of differential expression,
# which also updates the CSV files that are used in the gradH2O and canH2O projects.

files <- dir(pattern = "Rmd")
print(system.time(
  for(f in files) {
    sep <- paste(rep("=", nchar(f) + 4), collapse = "")
    message()
    print(sep, quote = FALSE)
    print(paste("|", f, "|"), quote = FALSE)
    print(sep, quote = FALSE)
    render(f)
  }
))
