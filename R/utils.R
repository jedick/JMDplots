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

# Add transparency to given color (moved from geo16S.R and exported on 20220223)
addalpha <- function(col, alpha) {
  x <- col2rgb(col)
  newcol <- rgb(x[1], x[2], x[3], maxColorValue = 255)
  newcol <- paste0(newcol, alpha)
  newcol
}
