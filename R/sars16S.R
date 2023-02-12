# JMDplots/sars16S.R
# Plots for chem16S paper 20220202
# Renamed to sars16S 20220907
# Moved to JMDplots 20230211

# Plot symbols and colors for body sites 20221125
pch_Oral <- 21
pch_Nasal <- 22
pch_Skin <- 23
pch_Gut <- 24
col_Oral <- "#D62728"
col_Nasal <- "#56B4E9"
col_Skin <- "#9467BD"
col_Gut <- "#E69F00"

##########################
### Plotting Functions ###
##########################

# Chemical metrics of reference proteomes and comparison with metaproteomes 20221125
sars16S_1 <- function(pdf = FALSE) {

  # Setup plot
  if(pdf) pdf("Figure_1.pdf", width = 7, height = 5)
  mat <- matrix(c(1,1, 2,2, 5,5, 3,3,3, 4,4,4), nrow = 2, byrow = TRUE)
  layout(mat)
  par(mar = c(3, 4, 3, 1))
  par(mgp = c(2.5, 1, 0))

  ## Panel A: Differences of ZC and nH2O between obligate anaerobic and aerotolerant genera 20221017

  # Use GTDB-based reference proteomes
  refdb <- "GTDB"

  # Read table from Million and Raoult, 2018
  dat <- read.csv(system.file("extdata/orp16S/MR18_Table_S1.csv", package = "JMDplots"))
  # Clean up names
  dat$Genus.name <- gsub(" ", "", dat$Genus.name)

  # Loop over ZC and nH2O
  for(metric in c("ZC", "nH2O")) {

      # Get amino acid compositions for genera
      aa <- taxon_AA[[refdb]]
      # Calculate ZC or nH2O
      if(metric == "ZC") {
        values <- ZCAA(aa)
        ylim <- c(-0.24, -0.1)
      }
      if(metric == "nH2O") {
        values <- H2OAA(aa)
        ylim <- c(-0.82, -0.65)
      }
      # Match genus names to RefSeq
      iref <- match(dat$Genus.name, aa$organism)
      values <- values[iref]
      # Get values for obligate anaerobes and aerotolerant genera
      values <- list(
        Anaerobe = values[dat$Obligate.anerobic.prokaryote == 1],
        Aerotolerant = values[dat$Obligate.anerobic.prokaryote == 0]
      )
      boxplot(values, ylab = cplab[[metric]], ylim = ylim)
      # Calculate p-value
      pval <- wilcox.test(values$Anaerobe, values$Aerotolerant)$p.value
      ptext <- bquote(italic(p) == .(format(signif(pval, 1), scientific = 2)))
      legend("bottomright", legend = ptext, bty = "n", inset = c(0, -0.03))
      # Show mean difference
      md <- round(diff(sapply(values, mean, na.rm = TRUE)), 4)
      mdtext <- paste("MD =", md)
      legend("topleft", legend = mdtext, bty = "n", inset = c(-0.05, -0.02))
      # Print coverage
      coverage <- round(sum(!is.na(iref)) / length(iref) * 100, 1)
      ctext <- paste0("Coverage = ", coverage, "%")
      if(metric == "ZC") print(ctext)
      if(metric == "ZC") label.figure("A. Genus reference proteomes", cex = 1.5, font = 2, adj = 0, xfrac = 0.02)

  }

  ## Panel B: ZC and nH2O from metaproteomes and reference proteomes 20220827
  par(mar = c(4, 4, 2.5, 1))
  par(mgp = c(2.5, 1, 0))

  # Function to calculate RMSD 20221018
  RMSD <- function(x, y) sqrt( mean( (y - x) ^ 2 ) )

  # Calculate chemical metrics
  dat <- getMP_sars16S(refdb)

  # Make ZC plot
  ZClim <- c(-0.20, -0.05)
  plot(ZClim, ZClim, xlab = quote(italic(Z)[C]~"of metaproteome"), ylab = quote(italic(Z)[C]~"of reference proteome"), type = "n")
  lines(ZClim, ZClim, lty = 2, col = 8)
  points(dat$ZC_MP, dat$ZC_16S, pch = dat$pch, bg = dat$bg, col = dat$col)
  # Show RMSD in legend 20221018
  RMSDtxt <- paste("RMSD =", round(RMSD(dat$ZC_MP, dat$ZC_16S), 3))
  legend("topleft", RMSDtxt, bty = "n")

  # Label figure
  label.figure("B. Community reference proteomes", font = 2, cex = 1.5, adj = 0, xfrac = 0.01)

  # Make nH2O plot
  nH2Olim <- c(-0.80, -0.59)
  plot(nH2Olim, nH2Olim, xlab = quote(italic(n)[H[2]*O]~"of metaproteome"), ylab = quote(italic(n)[H[2]*O]~"of reference proteome"), type = "n")
  lines(nH2Olim, nH2Olim, lty = 2, col = 8)
  points(dat$nH2O_MP, dat$nH2O_16S, pch = dat$pch, bg = dat$bg, col = dat$col)
  RMSDtxt <- paste("RMSD =", round(RMSD(dat$nH2O_MP, dat$nH2O_16S), 3))
  legend("topleft", RMSDtxt, bty = "n")

  # Make legend for datasets
  idup <- duplicated(dat$Name)
  # Get dataset names and number of samples
  names <- dat$Name[!idup]
  nsamp <- sapply(lapply(names, `==`, dat$Name), sum)
  names <- paste0(names, " (", nsamp, ")")
  dat <- dat[!idup, ]
  opar <- par(mar = c(0, 0, 4, 0))
  plot.new()
  legend("bottom", sapply(names, hyphen.in.pdf), pch = dat$pch, col = dat$col, pt.bg = dat$bg, bty = "n", cex = 1.2, inset = -0.18, xpd = NA)
  par(opar)

  if(pdf) dev.off()

}

# Chemical metrics of communities in body sites; viral inactivation and cross-omic comparison 20221125
sars16S_2 <- function(pdf = FALSE, refdb = "GTDB") {

  if(refdb == "GTDB") PDFfile <- "Figure_2.pdf" else PDFfile <- "Figure_2_RefSeq.pdf"
  if(pdf) pdf(PDFfile, width = 8, height = 6)
  mat <- matrix(c(
    1,1,1, 2,2, 3,3, 4,4,
    1,1,1, 7,7, 5,5, 6,6,
    8,8,8, 9,9,9, 10,10,10
  ), nrow = 3, byrow = TRUE)
  layout(mat, heights = c(1, 1, 1.5))
  par(mgp = c(2.5, 1, 0))

  ## Panels A-B: Community reference proteomes for body sites and viral inactivation
  ## based on 16S rRNA gene sequences from Boix-Amoros et al. (2021) 20220814

  # Panel A: nH2O-ZC plot for all samples
  par(mar = c(4, 4, 3, 1))
  # Plot small symbols: Any treatment
  Any <- plotmet_sars16S("BPB+21_AnyTreatment", extracolumn = c("Subject", "Site", "Treatment"),
    refdb = refdb, title = FALSE, pt.open.col = NA, plot.bg = FALSE, xlim = c(-0.20, -0.12),
    xlab = canprot::cplab$ZC, ylab = canprot::cplab$nH2O)
  # Plot large symbols: No treatment
  No <- plotmet_sars16S("BPB+21_NoTreatment", add = TRUE, cex = 2, extracolumn = c("Subject", "Site"),
    refdb = refdb, title = FALSE, pt.open.col = NA, plot.bg = FALSE)
  # Draw convex hull around all samples 20221125
  AnyNo <- list(Any, No)
  AnyNo[[2]] <- cbind(AnyNo[[2]], Treatment = NA)
  AnyNo <- do.call(rbind, AnyNo)
  AnyNo <- AnyNo[AnyNo$Site != "Control", ]
  ANhull <- chull(AnyNo$ZC, AnyNo$nH2O)
  polygon(AnyNo$ZC[ANhull], AnyNo$nH2O[ANhull], border = 8, lty = 2)

  # Plot p-values 20230204
  # Use metrics for untreated samples
  met <- getmetrics_sars16S("BPB+21_NoTreatment", refdb = refdb)
  met <- met[match(No$Run, met$Run), ]
  # Use metrics for gut and oral samples
  igut <- No$Site == "feces"
  ioral <- No$Site == "Oral cavity"
  plot.p.values(met$ZC[igut], met$ZC[ioral], met$nH2O[igut], met$nH2O[ioral], ypos = "top")

  legend("bottomleft", legend = c("", "  No treatment", "", "  Viral inactivation"), inset = c(-0.05, 0), bg = "white", bty = "n", cex = 0.9)
  legend("bottomleft", legend = c("Large symbols", "", "Small symbols", ""), text.font = 2, inset = c(-0.05, 0), bty = "n", cex = 0.9)
  title(hyphen.in.pdf("Community reference proteomes\n(Boix-Amor\u00f3s et al., 2021)"), font.main = 1)
  label.figure("A", font = 2, cex = 1.8, yfrac = 0.97)

  # Panel B: Viral inactivation
  par(mar = c(4, 4, 2, 1))
  # Keep samples for subjects (not controls)
  Any <- Any[grep("^S", Any$sample), ]
  No <- No[grep("^S", No$sample), ]

  plottreated <- function(Treatment) {
    # Get chemical metrics for treated samples
    Treated <- Any[Any$Treatment == Treatment, ]
    # Put in same order as untreated samples
    No_subject_site <- paste(No$Subject, No$Site)
    Treated_subject_site <- paste(Treated$Subject, Treated$Site)
    iTreated <- match(No_subject_site, Treated_subject_site)
    Treated <- Treated[iTreated, ]
    # Get colors and symbols
    pch <- sapply(Treated$Site, switch, "Oral cavity" = pch_Oral, "Nasal cavity" = pch_Nasal, "Skin of forearm" = pch_Skin, "feces" = pch_Gut, NA)
    col <- sapply(Treated$Site, switch, "Oral cavity" = col_Oral, "Nasal cavity" = col_Nasal, "Skin of forearm" = col_Skin, "feces" = col_Gut, NA)
    col <- sapply(col, add.alpha, "d0")
    # Calculate D_nH2O and D_ZC
    D_nH2O <- Treated$nH2O - No$nH2O
    D_ZC <- Treated$ZC - No$ZC
    # Plot D_nH2O and D_ZC
    plot(D_ZC, D_nH2O, xlab = canprot::cplab$DZC, ylab = canprot::cplab$DnH2O, xlim = c(-0.03, 0.03), ylim = c(-0.03, 0.015), type = "n")
    abline(h = 0, v = 0, lty = 2, col = 8)
    points(D_ZC, D_nH2O, pch = pch, bg = col, col = NA)
    # Plot p-values 20230204
    plot.p.values(No$ZC, Treated$ZC, No$nH2O, Treated$nH2O, paired = TRUE)
    title(Treatment, font.main = 1)
  }

  plottreated("Ethanol")
  label.figure("B", font = 2, cex = 1.8, yfrac = 0.94)
  plottreated("Formaldehyde")
  plottreated("Heat")
  plottreated("Psoralen")
  plottreated("Trizol")

  # Make legend
  plot.new()
  legend("center", legend = c("Nasal", "Skin", "Oral", "Gut"), pch = c(pch_Nasal, pch_Skin, pch_Oral, pch_Gut),
    pt.bg = c(col_Nasal, col_Skin, col_Oral, col_Gut), col = NA, bty = "n", cex = 1.5)

  ### Panels C-#: nH2O-ZC for reference proteomes, metagenomes and metaproteomes from various body sites 20221118
  par(mar = c(4, 4, 3, 1))
  xlim <- c(-0.2, -0.10)
  ylim <- c(-0.84, -0.60)

  ## Panel C: Plot 16S-based reference proteomes for controls in COVID-19 datasets 20220822
  # Setup plot
  plot(xlim, ylim, xlab = canprot::cplab$ZC, ylab = canprot::cplab$nH2O, type = "n")
  # Colors and point symbols for sample types
  col <- list(oro = col_Oral, naso = col_Nasal, gut = col_Gut)
  col <- lapply(col, add.alpha, "d0")
  pch <- list(oro = pch_Oral, naso = pch_Nasal, gut = pch_Gut)
  # Read precomputed mean values 20220823
  datadir <- system.file("extdata/sars16S", package = "JMDplots")
  means <- read.csv(file.path(datadir, "COVID_means_GTDB.csv"))
  # Loop over body sites
  for(type in c("gut", "oro", "naso")) {
    itype <- means$type == type
    # Plot points for control samples
    points(means$ZC_dn[itype], means$nH2O_dn[itype], pch = pch[[type]], bg = col[[type]], col = NA)
  }
  # Plot convex hull from Panel A
  polygon(AnyNo$ZC[ANhull], AnyNo$nH2O[ANhull], border = 8, lty = 2)
  # Plot p-values 20230204
  gut <- subset(means, type == "gut")
  oro <- subset(means, type == "oro")
  plot.p.values(gut$ZC_dn, oro$ZC_dn, gut$nH2O_dn, oro$nH2O_dn)
  legend("topleft", legend = c("Nasopharyngeal", "Oropharyngeal", "Gut"), pch = c(pch_Nasal, pch_Oral, pch_Gut),
    pt.bg = c(col_Nasal, col_Oral, col_Gut), col = NA, bty = "n")
  legend("topright", c("Various", "studies -", "see Table 1"), bty = "n")
  title(hyphen.in.pdf("Community reference proteomes\n(controls in COVID-19 studies)"), font.main = 1)
  label.figure("C", font = 2, cex = 1.8, yfrac = 0.96)

  ## Panel D: Metagenomes from different body sites 20221124
  # Setup plot
  plot(xlim, ylim, xlab = canprot::cplab$ZC, ylab = canprot::cplab$nH2O, type = "n")
  # Studies are for gut, oral, nasal
  studies <- c("ZZL+20", "CZH+22", "LLZ+21")
  pchs <- c(pch_Gut, pch_Oral, pch_Nasal)
  cols <- c(col_Gut, col_Oral, col_Nasal)
  cols <- sapply(cols, add.alpha, "b0")
  for(i in 1:length(studies)) {
    datadir <- system.file("extdata/sars16S", package = "JMDplots")
    file <- paste0(datadir, "/ARAST/", studies[i], "_AA.csv")
    dat <- read.csv(file)
    ZC <- ZCAA(dat)
    nH2O <- H2OAA(dat)
    points(ZC, nH2O, pch = pchs[i], bg = cols[i], col = NA)
    if(studies[i] == "ZZL+20") gut <- list(ZC = ZC, nH2O = nH2O)
    if(studies[i] == "CZH+22") oro <- list(ZC = ZC, nH2O = nH2O)
  }
  # Plot convex hull from Panel A
  polygon(AnyNo$ZC[ANhull], AnyNo$nH2O[ANhull], border = 8, lty = 2)
  # Plot p-values 20230204
  plot.p.values(gut$ZC, oro$ZC, gut$nH2O, oro$nH2O)
  # Add legend
  legend <- c("Nasopharyngeal", "Oropharyngeal", "Gut")
  legend("topleft", legend, pch = rev(pchs), pt.bg = c(col_Nasal, col_Oral, col_Gut), col = NA, bty = "n")
  legend <- c("Liu'21", "de Castilhos'22", "Zuo'20")
  legend("topright", legend, bty = "n")
  title(hyphen.in.pdf("Proteins from metagenomes\n(controls and COVID-19 patients)"), font.main = 1)
  label.figure("D", font = 2, cex = 1.8, yfrac = 0.96)

  ## Panel E: Metaproteomes from various body sites 20221114
  # Setup plot
  plot(xlim, ylim, xlab = canprot::cplab$ZC, ylab = canprot::cplab$nH2O, type = "n")
  # Define studies and point symbols
  studies_MP <- c(
    "TWC+22", "MLL+17", # Gut
    "GNT+21_cells", "JZW+22" # Oral
  )
  pchs <- c(
    pch_Gut, 25,
    pch_Oral, pch_Oral
  )
  cols <- c(
    col_Gut, col_Gut,
    col_Oral, col_Oral
  )
  cols <- sapply(cols, add.alpha, "d0")
  cexs <- ifelse(studies_MP %in% c("JZW+22"), 0.7, 1)
  # Store all data to calculate p-values 20230203
  gut.ZC <- gut.nH2O <- oral.ZC <- oral.nH2O <- numeric()
  # Loop over studies
  for(i in 1:length(studies_MP)) {
    # Get amino acid composition from metaproteome
    studydir <- strsplit(studies_MP[i], "_")[[1]][1]
    datadir <- system.file("extdata/sars16S", package = "JMDplots")
    aa <- read.csv(paste0(datadir, "/metaproteome/", studydir, "/", studies_MP[i], "_aa.csv"))
    # Calculate ZC and nH2O
    ZC <- ZCAA(aa)
    nH2O <- H2OAA(aa)
    # Add points
    points(ZC, nH2O, pch = pchs[i], bg = cols[i], cex = cexs[i], col = NA)
    if(studies_MP[i] %in% c("TWC+22", "MLL+17")) gut.ZC <- c(gut.ZC, ZC)
    if(studies_MP[i] %in% c("TWC+22", "MLL+17")) gut.nH2O <- c(gut.nH2O, nH2O)
    if(studies_MP[i] %in% c("GNT+21_cells", "JZW+22")) oral.ZC <- c(oral.ZC, ZC)
    if(studies_MP[i] %in% c("GNT+21_cells", "JZW+22")) oral.nH2O <- c(oral.nH2O, nH2O)
  }
  # Plot convex hull from Panel A
  polygon(AnyNo$ZC[ANhull], AnyNo$nH2O[ANhull], border = 8, lty = 2)
  # Plot p-values 20230204
  plot.p.values(gut.ZC, oral.ZC, gut.nH2O, oral.nH2O)
  # Add legend
  legend("topright", c("Jiang'22", "Granato'21"), title = "Oral",
    pch = c(pch_Oral, pch_Oral), pt.bg = col_Oral, col = NA,
    pt.cex = c(0.7, 1), cex = 0.9, bg = "transparent")
  legend("topleft", c(hyphen.in.pdf("Thuy-Boun'22"), "Maier'17"), title = "Gut", title.adj = 0.4,
    pch = c(pch_Gut, 25), pt.bg = col_Gut, col = NA,
    cex = 0.9, bg = "transparent", inset = c(0, 0.2))
  title(hyphen.in.pdf("Metaproteomes\n(gut and oral, not COVID-19)"), font.main = 1)
  label.figure("E", font = 2, cex = 1.8, yfrac = 0.96)

  if(pdf) dev.off()

}

# Differences of chemical metrics of oropharyngeal, nasopharyngeal, and
# gut microbiomes between control subjects and COVID-19 positive patients 20220806
sars16S_3 <- function(pdf = FALSE, refdb = "GTDB") {

  # Start plot
  if(refdb == "GTDB") PDFfile <- "Figure_3.pdf" else PDFfile <- "Figure_3_RefSeq.pdf"
  if(pdf) pdf(PDFfile, width = 15, height = 12)
  mat <- matrix(c(1,1,1,1, 2,2,2,2, 3,3,3,3, 4,4,4, 6,6,6, 5,5,5, 7,7,7, 8,8,8, 9,9,9, 10,10,10, 11,11,11), nrow = 3, byrow = TRUE)
  layout(mat)

  # Define plot settings
  par(cex = 1.2)
  par(mar = c(4, 4, 3, 1))
  startplot <- function() {
    plot(c(-0.008, 0.010), c(-0.012, 0.020), type  = "n", pch = ".", xlab = cplab$DZC, ylab = cplab$DnH2O)
    abline(h = 0, v = 0, lty = 2, col = 8)
  }
  # Colors and point symbols for sample types
  col <- list(naso = col_Nasal, oro = col_Oral, gut = col_Gut)
  pch <- list(naso = pch_Nasal, oro = pch_Oral, gut = pch_Gut)
  # Read precomputed mean values 20220823
  datadir <- system.file("extdata/sars16S", package = "JMDplots")
  means <- read.csv(file.path(datadir, "COVID_means_GTDB.csv"))

  ## Panel A: nH2O-ZC plots for 16S data for nasopharyngeal, oral/oropharyngeal, and gut datasets
  for(type in c("naso", "oro", "gut")) {
    startplot()
    itype <- means$type == type
    label <- 1:sum(itype)

    # Add points
    points(means$D_ZC[itype], means$D_nH2O[itype], pch = pch[[type]], col = col[[type]], bg = col[[type]])
    dx <- rep(0, sum(itype))
    dy <- rep(0.0018, sum(itype))
    if(type == "gut") {
      dy[c(1, 3, 7)] <- -0.0018
      dx[1] <- -0.0002
      dx[3] <- 0.0002
      dx[10] <- 0.0003
    }
    # Label points
    text(means$D_ZC[itype] + dx, means$D_nH2O[itype] + dy, label, cex = 0.8)
    # For gut, plot dropline at mean difference for bacterial metaproteome ZC 20220902
    at <- c(-0.0027, -0.0019)
    if(type == "gut") axis(1, at = at, labels = FALSE, tcl = -3, col = 8, lty = 2, lwd = 1.5)
    # Plot p-values 20230204
    plot.p.values(means$ZC_dn[itype], means$ZC_up[itype], means$nH2O_dn[itype], means$nH2O_up[itype], paired = TRUE)
    # Add plot title
    titles <- c(naso = "Nasopharyngeal", oro = "Oropharyngeal", gut = "Gut")
    title(titles[type], font.main = 1, line = 0.5)
    # Add panel title
    start <- hyphen.in.pdf("A. Community reference proteomes (")
    covid <- hyphen.in.pdf("COVID-19")
    label <- bquote(bold(.(start) * Delta == .(covid) ~ "minus control)"))
    if(type == "naso") label.figure(label, font = 2, xfrac = 0.695, yfrac = 0.95, cex = 1.1)
  }

  # Common nH2O and ZC limits for boxplots in Panels B and C
  nH2Olim <- c(-0.83, -0.68)
  ZClim <- c(-0.23, -0.09)

  ## Panel B: boxplots of ZC and nH2O of bacterial metaproteome in control and COVID-19 patients 20220830
  par(mar = c(3, 4, 4, 1))
  par(mgp = c(2.5, 1, 0))
  # Loop over taxonomy
  taxonomies <- c("Bacteria", "All")
  for(i in 1:2) {
    taxonomy <- taxonomies[i]
    datadir <- system.file("extdata/sars16S", package = "JMDplots")
    if(taxonomy == "All") aa <- read.csv(file.path(datadir, "metaproteome/HZX+21/HZX+21_aa.csv"))
    if(taxonomy == "Bacteria") aa <- read.csv(file.path(datadir, "metaproteome/HZX+21/HZX+21_bacteria_aa.csv"))
    # Identify control and COVID-19 patients
    icontrol <- grep("Ctrl", aa$organism)
    icovid <- grep("P", aa$organism)
    # Calculate ZC and nH2O
    ZC <- ZCAA(aa)
    nH2O <- H2OAA(aa)
    ZC_list <- list(Control = ZC[icontrol], "COVID-19" = ZC[icovid])
    nH2O_list <- list(Control = nH2O[icontrol], "COVID-19" = nH2O[icovid])
    names(ZC_list)[2] <- hyphen.in.pdf(names(ZC_list)[2])
    names(nH2O_list)[2] <- hyphen.in.pdf(names(nH2O_list)[2])

    # Add number of samples to group names
    len <- sapply(nH2O_list, length)
    labels <- paste0(names(nH2O_list), " (", len, ")")
    names(nH2O_list) <- ""
    # Make nH2O plot
    boxplot(nH2O_list, ylab = cplab$nH2O,  ylim = nH2Olim)
    # Make rotated labels (modified from https://www.r-bloggers.com/rotated-axis-labels-in-r-plots/)
    text(x = (1:2)+0.5, y = par()$usr[3] - 1.5 * strheight("A"), labels = labels, srt = 15, adj = 1, xpd = TRUE)
    # Add p-value
    nH2O_pvalue <- wilcox.test(nH2O_list[[1]], nH2O_list[[2]])$p.value
    legend <- bquote(italic(p) == .(format(signif(nH2O_pvalue, 1), scientific = 2)))
    legend("bottomleft", legend = legend, bty = "n", inset = c(-0.1, -0.05))
    # Add title
    main <- ifelse(taxonomy == "All", "All UniProt sequences", "Bacterial sequences")
    title(main, font.main = 1, line = 0.5)
    # Show mean difference
    nH2O_diff <- mean(nH2O_list[[2]]) - mean(nH2O_list[[1]])
    diffval <- signif(nH2O_diff, 2)
    difftxt <- bquote(Delta*italic(n)[H[2]*O] == .(diffval))
    title(difftxt)
    if(taxonomy == "Bacteria") label.figure("B. Gut metaproteomes (He et al., 2021)", font = 2, xfrac = 0.55, cex = 1.1)

    # Add number of samples to group names
    len <- sapply(ZC_list, length)
    labels <- paste0(names(ZC_list), " (", len, ")")
    names(ZC_list) <- ""
    # Make ZC plot
    boxplot(ZC_list, ylab = cplab$ZC, ylim = c(-0.16, -0.10))
    # Make rotated labels (modified from https://www.r-bloggers.com/rotated-axis-labels-in-r-plots/)
    text(x = (1:2)+0.5, y = par()$usr[3] - 1.5 * strheight("A"), labels = labels, srt = 15, adj = 1, xpd = TRUE)
    # Add p-value
    ZC_pvalue <- wilcox.test(ZC_list[[1]], ZC_list[[2]])$p.value
    legend <- bquote(italic(p) == .(format(signif(ZC_pvalue, 1), scientific = 2)))
    legend("topleft", legend = legend, bty = "n", inset = c(-0.1, -0.05))
    # Add title
    title(main, font.main = 1, line = 0.5)
    # Show mean difference
    ZC_diff <- mean(ZC_list[[2]]) - mean(ZC_list[[1]])
    diffval <- signif(ZC_diff, 2)
    # Make sure we plotted the ticks in Panel A at the correct locations
    stopifnot(diffval == at[i])
    difftxt <- bquote(Delta*italic(Z)[C] == .(diffval))
    title(difftxt)
  }
  # Draw arrows to ticks in Panel A
  arrows(-0.6, -0.083, 0.3, -0.065, length = 0.2, col = 8, lwd = 1.5, xpd = NA)
  arrows(0.7, -0.083, 0.5, -0.065, length = 0.2, col = 8, lwd = 1.5, xpd = NA)

  ## Panel C: Boxplots for nH2O and ZC in fecal MAGs 20221029
  # Read amino acid compositions of proteins predicted from MAGs
  datadir <- system.file("extdata/sars16S", package = "JMDplots")
  aa <- read.csv(file.path(datadir, "KWL22/KWL22_MAGs_prodigal_aa.csv.xz"))
  # https://github.com/Owenke247/COVID-19/blob/main/Pre-processed_Files/COVID19_metadata.txt
  dat <- read.csv(file.path(datadir, "KWL22/COVID19_metadata.txt"), sep = "\t")

  # Loop over nH2O and ZC
  for(metric in c("nH2O", "ZC")) {
    # Loop over SRA run prefixes:
    # Zuo et al. (PRJNA624223) and Yeoh et al. (PRJNA650244)
    for(SRAprefix in c("SRR1232", "SRR1307")) {
      # Get amino acid compositions for this BioProject
      iaa <- grep(SRAprefix, aa$protein)
      thisaa <- aa[iaa, ]
      # Calculate ZC or nH2O
      if(metric == "ZC") x <- ZCAA(thisaa) else x <- H2OAA(thisaa)
      if(metric == "ZC") ylim <- ZClim else ylim <- nH2Olim
      if(metric == "ZC") legend.x <- "topleft" else legend.x <- "bottomleft"
      # Get names of groups
      idat <- match(thisaa$protein, dat$Dat)
      group <- dat$Group[idat]
      # Make list of values in each group
      if(SRAprefix == "SRR1232") x_list <- list(Control = x[group == "Healthy_controls"], "COVID-19" = x[group == "COVID19"])
      if(SRAprefix == "SRR1307") x_list <- list(Control = x[group == "Healthy_control"], "COVID-19" = x[group == "COVID19"])
      # Add number of samples to group names
      len <- sapply(x_list, length)
      labels <- paste0(names(x_list), " (", len, ")")
      names(x_list) <- ""
      # Make boxplot
      boxplot(x_list, ylim = ylim, ylab = cplab[[metric]], xlab = "")
      # Make rotated labels (modified from https://www.r-bloggers.com/rotated-axis-labels-in-r-plots/)
      text(x = (1:2)+0.5, y = par()$usr[3] - 1.5 * strheight("A"), labels = labels, srt = 15, adj = 1, xpd = TRUE)
      # Add p-value
      x_pvalue <- wilcox.test(x_list[[1]], x_list[[2]])$p.value
      legend <- bquote(italic(p) == .(format(signif(x_pvalue, 1), scientific = 2)))
      legend(legend.x, legend = legend, bty = "n", inset = c(-0.1, -0.05))
      # Add title
      if(SRAprefix == "SRR1232") main <- "Zuo et al." else main <- "Yeoh et al."
      title(main, font.main = 1, line = 0.5)
      # Show mean difference
      x_diff <- mean(x_list[[2]]) - mean(x_list[[1]])
      diffval <- signif(x_diff, 2)
      if(metric == "ZC") difftxt <- bquote(Delta*italic(Z)[C] == .(diffval))
      if(metric == "nH2O") difftxt <- bquote(Delta*italic(n)[H[2]*O] == .(diffval))
      title(difftxt)
      figlab <- hyphen.in.pdf("C. Gut metagenome-assembled genomes (Ke et al., 2022)")
      if(metric == "nH2O" & SRAprefix == "SRR1232") label.figure(figlab, font = 2, xfrac = 0.814, cex = 1.1)
    }
  }

  if(pdf) dev.off()
}

# Summary of evidence for chemical variation of the human microbiome 20230112
sars16S_4 <- function(pdf = FALSE) {

  if(pdf) pdf("Figure_4.pdf", width = 8, height = 6)

  ylim <- c(0, 12)
  xlim <- c(0, 16)

  par(mar = c(0.1, 0.1, 0.1, 0.1))
  par(xaxs = "i", yaxs = "i")
  plot(xlim, ylim, type = "n", axes = FALSE, xlab = "", ylab = "")

  # -dx for adjusting math text with phantom(.)
  # +dx for adjusting columns after "Reference proteome"
  dx <- 0.2

  ## Header text
  y <- 12
  text(8, y, "Evidence for chemical\nvariation of microbiome", adj = c(0, 1))
  y <- 10.75
  text(4, y, "Body\nSite", adj = c(0, 1))
  text(5.5, y, "\nComparison", adj = c(0, 1))
  text(8, y, "Reference\nproteome", adj = c(0, 1), cex = 0.8, font = 2)
  text(8 + 4/3 + dx, y, hyphen.in.pdf("Meta-\ngenome/\nMAG"), adj = c(0, 1), cex = 0.8, font = 2)
  text(8 + 8/3 + dx, y, hyphen.in.pdf("Meta-\nproteome"), adj = c(0, 1), cex = 0.8, font = 2)
  y <- 10.1
  text(1.9, y, "Chemical\nprocesses", adj = c(0, 1), cex = 0.9, font = 3)

  # Dividing line
  y <- 9.5
  lines(c(4, 12 + dx), c(y, y))

  ## Normal body sites
  # Nasal
  y <- 9
  text(4, y, "Nasal", adj = c(0, 1))
  text(5.5, y, "vs other sites", adj = c(0, 1))
  text(8-dx, y, quote(phantom(.) %up% italic(Z)[C]), adj = c(0, 1))
  # Oral
  y <- 8
  text(4, y, "Oral", adj = c(0, 1))
  text(5.5, y, "vs other sites", adj = c(0, 1))
  text(8-dx, y, quote(phantom(.) %down% italic(Z)[C]), adj = c(0, 1))
  text(8 + 8/3, y, quote(phantom(.) %up% italic(Z)[C]), adj = c(0, 1))
  # Gut
  y <- 7
  text(4, y, "Gut", adj = c(0, 1))
  text(5.5, y, "vs other sites", adj = c(0, 1))
  text(8-dx, y, quote(phantom(.) %down% italic(n)[H[2]*O]), adj = c(0, 1))
  text(8 + 4/3, y, quote(phantom(.) %down% italic(n)[H[2]*O]), adj = c(0, 1))
  text(8 + 8/3, y, quote(phantom(.) %down% italic(Z)[C]), adj = c(0, 1))
  # Figure numbers
  y <- 6
  text(6.5, y, "Figures", adj = c(0, 1), cex = 0.8, font = 2)
  text(8, y, "2a,c", adj = c(0, 1), cex = 0.8, font = 2)
  text(8 + 4/3 + dx, y, "2d", adj = c(0, 1), cex = 0.8, font = 2)
  text(8 + 8/3 + dx, y, "2e", adj = c(0, 1), cex = 0.8, font = 2)

  # Dividing line
  y <- 5
  lines(c(4, 12 + dx), c(y, y))

  ## Ethanol treatment
  y <- 4.5
  text(4, y, "Multiple", adj = c(0, 1))
  text(5.5, y, "Ethanol\ntreatment\nvs untreated", adj = c(0, 1))
  text(8-dx, y, quote(phantom(.) %down% italic(n)[H[2]*O]), adj = c(0, 1))
  # Figure number
  y <- 3.5
  text(8, y, "2b", adj = c(0, 1), cex = 0.8, font = 2)

  # Dividing line
  y <- 2.5
  lines(c(4, 12 + dx), c(y, y))

  ## COVID
  y <- 1.5
  text(4, y, "Gut", adj = c(0, 1))
  text(5.5, y, hyphen.in.pdf("COVID-19\nvs control"), adj = c(0, 1))
  text(8-dx, y, quote(phantom(.) %down% italic(Z)[C]), adj = c(0, 1))
  text(8 + 4/3, y, quote(phantom(.) %down% italic(Z)[C]), adj = c(0, 1))
  text(8 + 8/3, y, quote(phantom(.) %down% italic(Z)[C]), adj = c(0, 1))
  # Figure numbers
  y <- 0.5
  text(8, y, "3a", adj = c(0, 1), cex = 0.8, font = 2)
  text(8 + 4/3 + dx, y, "3c", adj = c(0, 1), cex = 0.8, font = 2)
  text(8 + 8/3 + dx, y, "3b", adj = c(0, 1), cex = 0.8, font = 2)

  ## Draw images
  # Aspect ratio of image is ca. 4:9
  datadir <- system.file("extdata/sars16S", package = "JMDplots")
  body <- readPNG(file.path(datadir, "images/body.png"))
  rasterImage(body, 0, 5.6, 1.6, 9.2)
  rasterImage(body, 16, 0, 14.4, 3.6)
  # Aspect ratio of image is 1:1
  testtube <- readPNG(file.path(datadir, "images/testtube.png"))
  rasterImage(testtube, 11.5, 3.8, 12.5, 4.8)
  text(12.2, 3.95, "EtOH", cex = 0.8)

  ## Add lines and text
  # Nasal
  lines(c(1.15, 3.8), c(8.75, 8.9), col = 8)
  text(2.55, 8.825, "Oxidation", col = col_Nasal, adj = c(0.5, -0.1), cex = 0.9, font = 3, srt = 3.3)
  # Oral
  lines(c(1.1, 3.8), c(8.6, 7.9), col = 8)
  text(2.5, 8.24, "Reduction", col = col_Oral, adj = c(0.5, -0.1), cex = 0.9, font = 3, srt = -14.5)
  text(2.5, 8.24, "Oxidation", col = col_Nasal, adj = c(0.5, 1.1), cex = 0.9, font = 3, srt = -14.5)
  # Gut
  lines(c(0.8, 3.8), c(7.2, 6.9), col = 8)
  text(2.6, 7.02, "Dehydration", col = col_Gut, adj = c(0.5, -0.1), cex = 0.9, font = 3, srt = -5.6)
  text(2.6, 7.02, "Reduction", col = col_Oral, adj = c(0.5, 1.1), cex = 0.9, font = 3, srt = -5.6)
  # Ethanol treatment
  lines(c(9.2, 11.4), c(4.3, 4.3), col = 8)
  text(10.3, 4.3, "Dehydration", col = col_Gut, adj = c(0.5, -0.1), cex = 0.9, font = 3)
  # COVID-19
  lines(c(11.7, 15.2), c(1.3, 1.6), col = 8)
  text(13.4, 1.44, "Reduction", col = col_Oral, adj = c(0.5, -0.1), cex = 0.9, font = 3, srt = 4.8)

  if(pdf) dev.off()

}

##################################
### Data Processing Functions  ###
##################################

# Gather values of ZC and nH2O for metaproteomes and 16S-based estimates 20220828
# Add refdb argument 20221017
# Add zero_AA argument 20221018
getMP_sars16S <- function(refdb = "RefSeq", zero_AA = NULL) {

  studies_MP <- c("MLL+17", "TWC+22", "MPB+19", "PMM+18", 
    "KTS+17", "KTS+17.mock", "HTZ+17")
  studies_16S <- c("MLL+17", "TWC+22", "MPB+19", "RYP+14", 
    "KTS+17", "KTS+17.mock", "HTZ+17")
  n <- length(studies_16S)
  pchs <- c(21, 22, 24, 25, 12, 18, 20)
  bgs <- sapply(c(4, 8, 5, 6, NA, NA, NA), add.alpha, alpha = "d0")
  cols <- c(NA, NA, 1, 1, 2, "#ff000080", "#00000080")

  out <- lapply(seq_along(studies_16S), function(i) {
    # Get 16S estimates
    metrics <- getmetrics_sars16S(studies_16S[i], refdb = refdb, zero_AA = zero_AA)
    mdat <- getmdat_sars16S(studies_16S[i], metrics)
    # Get metaproteome values - look in orp16S then sars16S directory
    datadir <- system.file("extdata/orp16S", package = "JMDplots")
    file <- paste0(datadir, "/metaproteome/", studies_MP[i], "/", studies_MP[i], "_aa.csv")
    if(!file.exists(file)) {
      datadir <- system.file("extdata/sars16S", package = "JMDplots")
      file <- paste0(datadir, "/metaproteome/", studies_MP[i], "/", studies_MP[i], "_aa.csv")
    }
    aa <- read.csv(file)
    # Set abundances of selected amino acids to zero 20221018
    if(!is.null(zero_AA)) aa[, zero_AA] <- 0
    # Match 16S samples to metaproteome
    iaa <- match(mdat$metadata$Metaproteome, aa$organism)
    # Drop non-matching samples
    ina <- is.na(iaa)
    aa <- aa[iaa[!ina], ]
    metadata <- mdat$metadata[!ina, , drop = FALSE]
    metrics <- mdat$metrics[!ina, , drop = FALSE]
    stopifnot(all(metadata$Metaproteome == aa$organism))
    # Make data frame
    data.frame(Name = metadata$Name,
      Study_16S = studies_16S[i], Run_16S = metadata$Run, Sample_16S = metadata$Sample,
      Study_MP = studies_MP[i], Sample_MP = aa$organism,
      ZC_16S = round(metrics$ZC, 6), ZC_MP = round(ZCAA(aa), 6),
      nH2O_16S = round(metrics$nH2O, 6), nH2O_MP = round(H2OAA(aa), 6),
      #pch = metadata$pch, col = metadata$col
      pch = pchs[i], bg = bgs[i], col = cols[i]
    )
  })

  # Return value instead of saving to file 20221018
  do.call(rbind, out)

}

# Calculate mean values of chemical metrics for patients and controls in COVID-19 datasets 20220823
COVID_means <- function(refdb = "GTDB") {

  # Function to calculate mean ZC and nH2O for patients and controls
  getmeans <- function(study) {
    print(study)
    metrics <- getmetrics_sars16S(study, refdb = refdb)
    mdat <- getmdat_sars16S(study, metrics)
    # "up" for disease/positive, "down" for control/negative
    iup <- sapply(mdat$metadata$pch == 25, isTRUE)
    idn <- sapply(mdat$metadata$pch == 24, isTRUE)
    # Calculate means of chemical metrics
    ZC_dn <- mean(mdat$metrics$ZC[idn])
    nH2O_dn <- mean(mdat$metrics$nH2O[idn])
    ZC_up <- mean(mdat$metrics$ZC[iup])
    nH2O_up <- mean(mdat$metrics$nH2O[iup])
    # Calculate p-values 20220905
    ZC_pvalue <- wilcox.test(mdat$metrics$ZC[idn], mdat$metrics$ZC[iup])$p.value
    nH2O_pvalue <- wilcox.test(mdat$metrics$nH2O[idn], mdat$metrics$nH2O[iup])$p.value
    # Return values
    # Include number of samples and p-values 20220905
    list(n_dn = sum(idn), n_up = sum(iup), ZC_dn = ZC_dn, nH2O_dn = nH2O_dn, ZC_up = ZC_up, nH2O_up = nH2O_up, ZC_pvalue = ZC_pvalue, nH2O_pvalue = nH2O_pvalue)
  }

  # List datasets
  studies <- list(
    # Nasopharyngeal datasets
    naso = c("CSC+22", "ENJ+21", "GKJ+22", "HMH+21", "MLW+21_Nasopharyngeal", "PMM+22", "SGC+21", "SRS+22", "VCV+21"),
    # Oral/oropharyngeal datasets
    oro = c("GBS+22", "GWL+21", "IZC+21", "MAC+21", "MLW+21_Oropharyngeal", "RFH+22_Oral", "RWC+21_Oral", "WCJ+21_Oral", "XLZ+21"),
    # Gut datasets
    gut = c("CGC+22", "FBD+22", "GCW+20", "KMG+21", "MIK+22", "MMP+21", "NGH+21", "RFH+22_Gut", "RWC+21_Gut", "RDM+22", "SRK+22", "WCJ+21_Gut", "ZZZ+21")
  )

  out <- list()
  # Calculate means for each sample type
  for(i in 1:3) {
    means <- lapply(studies[[i]], getmeans)
    n_dn <- sapply(means, "[[", "n_dn")
    n_up <- sapply(means, "[[", "n_up")
    ZC_dn <- sapply(means, "[[", "ZC_dn")
    nH2O_dn <- sapply(means, "[[", "nH2O_dn")
    ZC_up <- sapply(means, "[[", "ZC_up")
    nH2O_up <- sapply(means, "[[", "nH2O_up")
    ZC_pvalue <- sapply(means, "[[", "ZC_pvalue")
    nH2O_pvalue <- sapply(means, "[[", "nH2O_pvalue")
    thisout <- data.frame(type = names(studies)[i], study = studies[[i]], n_dn = n_dn, n_up = n_up,
      ZC_dn = ZC_dn, nH2O_dn = nH2O_dn, ZC_up = ZC_up, nH2O_up = nH2O_up, ZC_pvalue = ZC_pvalue, nH2O_pvalue = nH2O_pvalue)
    out[[i]] <- thisout
  }

  # Put together data frames
  out <- do.call(rbind,out)
  # Calculate mean differences of chemical metrics
  D_ZC <- out$ZC_up - out$ZC_dn
  D_nH2O <- out$nH2O_up - out$nH2O_dn
  out <- cbind(out, D_ZC, D_nH2O)
  file <- paste0("COVID_means_", refdb, ".csv")
  write.csv(out, file, row.names = FALSE, quote = FALSE)

}

# Function to add p-values to x and y axes 20230204
plot.p.values <- function(ZC.1, ZC.2, nH2O.1, nH2O.2, paired = FALSE, ypos = "bottom") {
  p.ZC <- format(signif(wilcox.test(ZC.1, ZC.2, paired = paired)$p.value, 1), scientific = 2)
  p.nH2O <- format(signif(wilcox.test(nH2O.1, nH2O.2, paired = paired)$p.value, 1), scientific = 2)
  p.ZC.txt <- bquote(italic(p) == .(p.ZC))
  p.nH2O.txt <- bquote(italic(p) == .(p.nH2O))
  pu <- par("usr")
  dx <- (pu[2] - pu[1]) / 30
  dy <- (pu[4] - pu[3]) / 30
  text(pu[2] - dx, pu[3] + dy/2, p.ZC.txt, adj = c(1, 0), cex = 0.9)
  if(ypos == "bottom") text(pu[1] + dx/2, pu[3] + dy, p.nH2O.txt, srt = 90, adj = c(0, 1), cex = 0.9)
  if(ypos == "top") text(pu[1] + dx/2, pu[4] - dy, p.nH2O.txt, srt = 90, adj = c(1, 1), cex = 0.9)
}
