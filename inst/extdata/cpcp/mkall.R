# JMDplots/inst/extdata/cpcp/mkall.R
# process all vignettes 20191207
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
