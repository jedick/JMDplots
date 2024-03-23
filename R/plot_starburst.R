# Plot chemical metrics for specific taxa and their children 20200911
# Adapted from taxacomp() previously in geo16S_util.R 20231229
plot_starburst <- function(
  taxa = c("Bacteria", "Archaea"), metrics = c("Zc", "nH2O"), refdb = "RefSeq",
  remove_species_20000 = TRUE, terminal_H2O = 0,
  xlim = NULL, ylim = NULL, pch = NULL, col = seq_along(taxa), lcol = NULL,
  lwd = rep(1, length(taxa)), hline = NULL, legend.x = NA, identify = FALSE) {

  # Compute chemical metrics of all taxa in reference database
  datadir <- system.file(file.path("extdata", refdb), package = "chem16S")
  if(refdb == "UHGG") datadir <- system.file(file.path("extdata/RefDB", refdb), package = "JMDplots")
  aa_refdb <- read.csv(file.path(datadir, "taxon_AA.csv.xz"))
  refdb_metrics <- data.frame(
    rank = aa_refdb$protein,
    taxon = aa_refdb$organism,
    ntaxa = aa_refdb$ref,
    parent = aa_refdb$abbrv,
    # Use get() to call canprot::Zc, canprot::nO2, or canprot::nH2O
    Xvals = get(metrics[1])(aa_refdb, terminal_H2O = terminal_H2O),
    Yvals = get(metrics[2])(aa_refdb, terminal_H2O = terminal_H2O)
  )
  # Default point symbols
  if(is.null(pch)) pch <- 21
  pch <- rep(pch, length.out = length(taxa))

  # Use semi-transparent colors for lines 20210518
  if(is.null(lcol)) lcol <- adjustcolor(col, alpha.f = 0.5)

  # Store all values to compute plot limits and for identify()
  Xvals <- Yvals <- numeric()
  names <- character()

  # Loop 1: Get values to plot
  vals <- list()
  # Where to keep amino acid compositions of species
  aa_species <- NULL

  for(i in seq_along(taxa)) {

    taxon <- taxa[i]
    # Get the chemical metrics for this taxon
    itaxon <- refdb_metrics$taxon == taxon
    if(sum(itaxon) > 1) warning(paste0("found more than one ", taxon, " (", paste(refdb_metrics$rank[itaxon], collapse = ", "), "); using the first"))
    taxon_metrics <- refdb_metrics[which(itaxon)[1], ]
    if(identical(taxon_metrics$rank, "genus")) {

      # For a genus, look for children (species) in reference database 20210603
      if(is.null(aa_species)) {

        # Read amino acid compositions and taxonomy for the specified reference database
        aa_refdb_all <- read.csv(system.file(file.path("extdata/RefDB", refdb, "genome_AA.csv.xz"), package = "JMDplots"))
        taxonomy <- read.csv(system.file(file.path("extdata/RefDB", refdb, "taxonomy.csv.xz"), package = "JMDplots"))

        if(refdb == "RefSeq" & remove_species_20000) {
          # Take out species with > 20000 sequences (biased to high Zc/low nH2O in RefSeq) 20210604
          is_species <- !is.na(taxonomy$species)
          is_highseq <- aa_refdb_all$chains > 20000
          is_species_20000 <- is_species & is_highseq
          # NOTE: the following genera are completely removed: Buchnera, Sorangium, Dolosigranulum, Enhygromyxa, Ruthenibacterium
          message(paste("plot_taxa_children: removing", sum(is_species_20000), "species with > 20000 sequences"))
          aa_refdb_all <- aa_refdb_all[!is_species_20000, ]
          taxonomy <- taxonomy[!is_species_20000, ]
        }

        # Keep species-level taxa that have a genus name
        is_not_orphan_species <- !is.na(taxonomy$species) & !is.na(taxonomy$genus)
        aa_species <- aa_refdb_all[is_not_orphan_species, ]
        # Put genus name in "abbrv" column
        aa_species$abbrv <- taxonomy$genus[is_not_orphan_species]
        # Keep species with at least 500 sequences
        aa_species <- aa_species[aa_species$chains >= 500, ]

      }

      # Find children (species) of this genus and calculate chemical metrics
      is_child <- aa_species$abbrv == taxon
      aa_children <- aa_species[is_child, ]
      children_metrics <- data.frame(
        taxon = aa_children$ref,
        chains = aa_children$chains,
        Xvals = get(metrics[1])(aa_children, terminal_H2O = terminal_H2O),
        Yvals = get(metrics[2])(aa_children, terminal_H2O = terminal_H2O)
      )

    } else {

      # Get the chemical metrics for all children
      ichildren <- refdb_metrics$parent == taxon
      children_metrics <- refdb_metrics[ichildren, ]

    }
    # Store the values to make the plot
    vals[[i]] <- list(taxon = taxon_metrics, children = children_metrics)

    # Keep values for identification
    Xvals <- c(Xvals, taxon_metrics$Xvals, children_metrics$Xvals)
    Yvals <- c(Yvals, taxon_metrics$Yvals, children_metrics$Yvals)
    names <- c(names, taxon_metrics$taxon, children_metrics$taxon)

  }
  # Initialize plot
  if(is.null(xlim)) xlim <- extendrange(Xvals)
  if(is.null(ylim)) ylim <- extendrange(Yvals)
  plot(xlim, ylim, xlab = chemlab(metrics[1]), ylab = chemlab(metrics[2]), type = "n", xaxs = "i", yaxs = "i")
  # Add horizontal lines to show range of following plot 20200925
  if(!is.null(hline)) abline(h = hline, lty = 2, col = "gray40")

  # Loop 2: Plot lines from parents to all children 20200925
  for(i in seq_along(taxa)) {
    taxon_metrics <- vals[[i]]$taxon
    children_metrics <- vals[[i]]$children
    for(j in seq_along(children_metrics$Xvals)) lines(c(taxon_metrics$Xvals, children_metrics$Xvals[j]),
      c(taxon_metrics$Yvals, children_metrics$Yvals[j]), col = lcol[i], lwd = lwd[i])
  }

  # Loop 3: Plot points for parents
  for(i in seq_along(taxa)) {
    taxon_metrics <- vals[[i]]$taxon
    children <- vals[[i]]$children
    points(taxon_metrics$Xvals, taxon_metrics$Yvals, pch = pch[i], cex = 1.5, col = col[i], bg = col[i])
  }

  # Loop 4: Plot points for children
  for(i in seq_along(taxa)) {
    taxon_metrics <- vals[[i]]$taxon
    children_metrics <- vals[[i]]$children
    # Use white outline for black points 20210518
    pt.col <- 1
    if(is.numeric(col[i])) {
      if((col[i] - 1) %% 8 == 0) pt.col <- "white"
    }
    points(children_metrics$Xvals, children_metrics$Yvals, pch = pch[i], cex = 0.7, col = pt.col, bg = col[i], lwd = 0.5)
  }

  if(!is.null(legend.x) & !identical(legend.x, NA)) legend(legend.x, taxa, pch = pch, col = seq_along(taxa), pt.bg = seq_along(taxa), cex = 0.9, bg = "white")
  if(identify) identify(Xvals, Yvals, names)
  # Return values invisibly 20210603
  names(vals) <- taxa
  invisible(vals)

}
