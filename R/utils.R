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

# Function to add significant difference letters to boxplot 20220531
# Moved from utogig.R and exported 20220609
cldfun <- function(Zclist, bp, dy) {
  # Add names if missing 20240202
  if(is.null(names(Zclist))) names(Zclist) <- 1:length(Zclist)
  # Remove empty elements 20240202
  Zclist_1 <- Zclist[sapply(Zclist, length) != 0]
  # Turn the list into a data frame with "group" column taken from names of the list elements
  Zcdat <- do.call(rbind, sapply(1:length(Zclist_1), function(i) data.frame(group = names(Zclist_1)[i], Zc = Zclist_1[[i]]), simplify = FALSE))
  # One-way ANOVA and Tukey's Honest Significant Differences
  # Adapted from https://statdoe.com/one-way-anova-and-box-plot-in-r/
  anova <- aov(Zc ~ group, data = Zcdat)
  tukey <- TukeyHSD(anova)
  # Compact letter display
  cld <- multcompLetters4(anova, tukey, reversed = TRUE)$group$Letters
  # Get into same order as data
  cld <- cld[match(names(Zclist), names(cld))]
  # Add to plot
  n <- length(Zclist)
  text((1:n) + 0.35, bp$stats[4, ] + dy, cld)
}

# To get hyphen instead of minus sign 20220630
# https://stackoverflow.com/questions/10438398/any-way-to-disable-the-minus-hack-in-pdf-poscript-output
hyphen.in.pdf <- function(x) {
  # We only want to make the substitution in a pdf device (won't work in png, e.g. for knitr vignettes)
  if(identical(names(dev.cur()), "pdf")) gsub("-", "\uad", x, fixed = TRUE) else x
}
