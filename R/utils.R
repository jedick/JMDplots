# JMDplots/R/utils.R

# Add exif metadata to PDF file 20191027
addexif <- function(name, title, paperref) {
  titlearg <- paste0('-Title="', title, '"')
  authorarg <- '-Author="Jeffrey M. Dick"'
  creatorarg <- paste0('-Creator="R (JMDplots::', name, ')"')
  subjectarg <- paste0('-Subject="Plot based on ', paperref, '"')
  allargs <- paste(titlearg, authorarg, creatorarg, subjectarg)
  file <- paste0(name, ".pdf")
  cmd <- paste("exiftool -overwrite_original", allargs, file)
  # this will do nothing on systems that don't have exiftool (instead of producing an error)
  tryCatch(system(cmd), error = function(e) {})
}

# To get hyphen instead of minus sign 20220630
# https://stackoverflow.com/questions/10438398/any-way-to-disable-the-minus-hack-in-pdf-poscript-output
hyphen.in.pdf <- function(x) {
  # We only want to make the substitution in a pdf device (won't work in png, e.g. for knitr vignettes)
  if(identical(names(dev.cur()), "pdf")) gsub("-", "\uad", x, fixed = TRUE) else x
}
