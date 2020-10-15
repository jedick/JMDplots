# create csv files in extdata/summary
# (used by groupplots())
for(what in c("colorectal", "pancreatic", "hypoxia", "osmotic")) {
  pdat_fun <- paste0("pdat_", what)
  if(what == "osmotic") pdat_fun <- ".pdat_osmotic"
  datasets <- get(pdat_fun)(2017)
  comptab <- lapply(datasets, function(dataset) {
    pdat <- get(pdat_fun)(dataset, basis = "QEC")
    get_comptab(pdat, mfun = "mean", oldstyle = TRUE)
  })
  # write summary table
  comptab <- do.call(rbind, comptab)
  comptab <- cbind(set = c(letters, LETTERS)[1:nrow(comptab)], comptab)
  comptab[, 6:15] <- signif(comptab[, 6:15], 4)
  filename <- paste0("summary_", what, ".csv")
  write.csv(comptab, filename, row.names = FALSE, quote = 3)
}
