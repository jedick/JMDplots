# JMDplots/R/utils.R

# add exif metadata to PDF file 20191027
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

# run obigt or OBIGT depending on CHNOSZ version 20200704
OBIGT <- function(...) {
  if(packageVersion("CHNOSZ") > "1.3.6") fun <- get("OBIGT", asNamespace("CHNOSZ"))
  else fun <- get("obigt", asNamespace("CHNOSZ"))
  fun(...)
}

# run add.obigt or add.OBIGT depending on CHNOSZ version 20200704
add.OBIGT <- function(...) {
  if(packageVersion("CHNOSZ") > "1.3.6") fun <- get("add.OBIGT", asNamespace("CHNOSZ"))
  else fun <- get("add.obigt", asNamespace("CHNOSZ"))
  fun(...)
}

# run mod.obigt or mod.OBIGT depending on CHNOSZ version 20200704
mod.OBIGT <- function(...) {
  if(packageVersion("CHNOSZ") > "1.3.6") fun <- get("mod.OBIGT", asNamespace("CHNOSZ"))
  else fun <- get("mod.obigt", asNamespace("CHNOSZ"))
  fun(...)
}
