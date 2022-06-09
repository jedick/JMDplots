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

# Function to add significant difference letters 20220531
# Moved from utegig.R and exported 20220609
cldfun <- function(ZClist, bp, dy) {
  # Ugly one-liner to turn a list into a data frame with "group" column taken from names of the list elements
  ZCdat <- do.call(rbind, sapply(1:length(ZClist), function(i) data.frame(group = names(ZClist)[i], ZC = ZClist[[i]]), simplify = FALSE))
  # One-way ANOVA and Tukey's Honest Significant Differences
  # Adapted from https://statdoe.com/one-way-anova-and-box-plot-in-r/
  anova <- aov(ZC ~ group, data = ZCdat)
  tukey <- TukeyHSD(anova)
  # Compact letter display
  cld <- multcompLetters4(anova, tukey, reversed = TRUE)$group$Letters
  # Get into same order as data
  cld <- cld[match(names(ZClist), names(cld))]
  # Add to plot
  n <- length(ZClist)
  text((1:n) + 0.35, bp$stats[4, ] + dy, cld)
}

