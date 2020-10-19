# JMDplots/makevig.R
# compile and view vignettes from command line
# 20200421 jmd

makevig <- function(vig = NULL) {
  vig.allowed <- gsub("pdat_", "", grep("pdat_", ls("package:JMDplots"), value = TRUE))
  vig.allowed <- vig.allowed[!vig.allowed == "aneuploidy"]
  isnull <- is.null(vig)
  toomany <- length(vig) > 1
  notallowed <- !any(vig %in% vig.allowed)
  if(isnull | toomany | notallowed) stop("'vig' should be one of: ", paste(vig.allowed, collapse = ", "))
  # location of the vignette directory
  vigdir <- system.file("vignettes", package = "JMDplots")
  # names of the vignette source and html output files
  vigfile <- file.path(vigdir, paste0(vig, ".Rmd"))
  htmlfile <- tempfile(pattern = paste0(vig, "_"), fileext = ".html")
  # compile the vignette and open it in the browser
  rmarkdown::render(vigfile, output_file = htmlfile, knit_root_dir = vigdir)
  # set 'browser' option so vignettes open under Rstudio 20201019
  # https://stackoverflow.com/questions/62536479/the-command-exams2html-does-not-generate-html-page-when-it-is-run-from-rstudio
  if(.Platform$OS.type == "windows") oldopt <- options(browser = NULL)
  browseURL(htmlfile)
  if(.Platform$OS.type == "windows") options(oldopt)
  # return the path of html file
  htmlfile
}
