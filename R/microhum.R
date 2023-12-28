# JMDplots/microhum.R
# Plots for chem16S paper 20220202
# Renamed to sars16S 20220907
# Moved to JMDplots 20230211
# Renamed to hum16S 20230530
# Renamed to microhum 20230723

# Plot symbols and colors for body sites 20221125
pch_Oral <- 21
pch_Nasal <- 22
pch_Skin <- 23
pch_Gut <- 24
pch_IBD <- 25
pch_UG <- 25
col_Oral <- "#D62728"
col_Nasal <- "#56B4E9"
col_Skin <- "#9467BD"
col_Gut <- "#E69F00"
col_IBD <- "#009E73"
col_UG <- "#9E9E9E"

# Location of data files
getdatadir <- function() {
  system.file("extdata/microhum", package = "JMDplots")
}

##########################
### Plotting Functions ###
##########################

# Compatibility of inferences from shotgun metagenomes and community reference proteomes
# 20211218 First version: comparison of Zc used in geo16S paper
#          Based on samples used by Aßhauer et al. (2015) (Tax4Fun paper) with additions by Dick and Tan (2023)
# 20231218 Added human DNA screening and nH2O-nO2 plots for microhum paper
microhum_1 <- function(pdf = FALSE) {

  if(pdf) pdf("Figure_1.pdf", width = 10, height = 6)
  mat <- matrix(c(1,1,2,2,3,3,4,4, 0,5,5,5,6,6,6,0), nrow = 2, byrow = TRUE)
  layout(mat)
  par(mar = c(3.5, 3.5, 2, 1), mgp = c(2.5, 1, 0))
  par(cex = 0.8)

  # Get chemical metrics for community reference proteomes
  metrics <- getmetrics_microhum("HMP12")
  # Get sample metadata
  metadata <- getmdat_microhum("HMP12")
  # Define colors and symbols
  bg <- sapply(metadata$"Body site", switch, "Skin" = col_Skin, "Nasal cavity" = col_Nasal, "Oral cavity" = col_Oral, "GI tract" = col_Gut, "UG tract" = col_UG)
  pch <- sapply(metadata$"Body site", switch, "Skin" = pch_Skin, "Nasal cavity" = pch_Nasal, "Oral cavity" = pch_Oral, "GI tract" = pch_Gut, "UG tract" = pch_UG)
  # Use semi-transparent colors for symbol outline 20220122
  c1 <- addalpha(1, "80")
  c1_light <- addalpha(1, "40")

  # Define cutoff for number of proteins (% of MG reads)
  cutoff <- 40

  ## Top row: nO2 and nH2O of community reference proteome vs shotgun metagenome
  # Loop over metrics
  for(metric in c("nO2", "nH2O")) {
    
    if(metric == "Zc") {
      xlim <- c(-0.2, -0.08)
      ylim <- c(-0.2, -0.12)
      main <- "Carbon oxidation state (Zc)"
    }

    if(metric == "nO2") {
      xlim <- c(-0.8, -0.55)
      ylim <- c(-0.8, -0.6)
      main <- quote("Stoichiometric oxidation state (" * italic(n)[O[2]] * ")")
    }

    if(metric == "nH2O") {
      xlim <- c(-1.05, -0.72)
      ylim <- c(-0.8, -0.72)
      main <- quote("Stoichiometric hydration state (" * italic(n)[H[2]*O] * ")")
    }

    # Compute min/max limits for 1:1 line
    xylim <- c(min(xlim, ylim), max(xlim, ylim))

    # Loop over MG pipelines
    for(pipeline in c("no_screening", "screening")) {

      # Read MG datasets processed without or with contamination removal
      if(pipeline=="no_screening") {
        MG_aafile <- file.path(getdatadir(), "ARAST/HMP12_no_screening_aa.csv")
        xlab <- "Metagenomes - NO screening"
      }
      if(pipeline=="screening") {
        MG_aafile <- file.path(getdatadir(), "ARAST/HMP12_aa.csv")
        xlab <- "Metagenomes - WITH screening"
      }
      MG_aa <- read.csv(MG_aafile)

      # Put MG and 16S samples in same order
      i16S <- match(metadata$Run, metrics$Run)
      metrics <- metrics[i16S, ]
      iMG <- match(metadata$Metagenome, MG_aa$protein)
      MG_aa <- MG_aa[iMG, ]
      # Make sure the 16S and metagenomes are paired correctly
      stopifnot(all(na.omit(metrics$Run == metadata$Run)))
      stopifnot(all(MG_aa$protein == metadata$Metagenome))

      # Get chemical metric for community reference proteomes
      metric_16S <- metrics[, metric]
      # Get chemical metric for metagenomes
      # This executes canprot::Zc(), canprot::nO2(), or canprot::nH2O()
      metric_MG <- get(metric)(MG_aa)

      # Start plot
      if(metric == "nO2" & pipeline == "no_screening") ylab <- "Community reference proteomes" else ylab <- ""
      plot(xlim, ylim, type = "n", xlab = xlab, ylab = ylab)
      lines(xylim, xylim, lty = 2, col = "gray40")
      # Use open symbols (no color) for MG with low numbers of protein fragments 20231218
      ilow <- MG_aa$chains / MG_aa$ref * 100 < cutoff
      ilow[is.na(ilow)] <- TRUE
      mybg <- bg
      mybg[ilow] <- "transparent"
      col <- ifelse(ilow, c1_light, c1)
      points(metric_MG, metric_16S, pch = pch, bg = mybg, col = col)

      # Show R-squared values excluding MG with low numbers of protein fragments 20231219
      mylm <- lm(metric_16S[!ilow] ~ metric_MG[!ilow])
      R2 <- round(summary(mylm)$r.squared, 2)
      R2txt <- bquote(italic(R)^2 == .(R2))
      legend("bottomleft", legend = R2txt, bty = "n", cex = 0.85)

      if(pipeline == "no_screening") {
        if(metric == "nO2") {
          legend("topleft", c("Skin", "Nasal cavity", "Oral cavity", "GI tract", "UG tract"),
            pch = c(pch_Skin, pch_Nasal, pch_Oral, pch_Gut, pch_UG),
            pt.bg = c(col_Skin, col_Nasal, col_Oral, col_Gut, col_UG), col = c1, bty = "n", cex = 0.7)
          # Use mtext to put title over two plots
          mtext(main, at = -0.5, line = 0.5)
          label.figure("A", font = 2, cex = 1.8, yfrac = 0.97)
        }
        if(metric == "nH2O") {
          mtext(main, at = -0.65, line = 0.5)
          label.figure("B", font = 2, cex = 1.8, yfrac = 0.97)
        }
      }

      if(pipeline == "screening" & metric == "nO2") {
        cutoff_txt <- c(paste("<", cutoff, "% of reads"), paste(">", cutoff, "% of reads"))
        legend("bottomright", cutoff_txt, pch = pch_UG, pt.bg = c("transparent", col_UG), col = c(c1_light, c1),
          bty = "n", cex = 0.7, title = "Protein sequences in MG")
        legend("bottomright", c("", ""), pch = pch_Nasal, pt.bg = c("transparent", col_Nasal), col = c(c1_light, c1),
          bty = "n", cex = 0.7, inset = c(0.5, 0))
      }

    }
  }

  ## Bottom row: nH2O vs nO2 for community reference proteomes and metagenomes
  xlim <- c(-0.77, -0.63)
  ylim <- c(-0.82, -0.72)
  par(mar = c(4, 4, 3, 1))
  par(cex.lab = 1.2)
  plotmet_microhum("HMP12", title = FALSE, pt.open.col = c1, xlim = xlim, ylim = ylim)
  legend("bottomright", c("Skin", "Nasal cavity", "Oral cavity", "GI tract", "UG tract"),
    pch = c(pch_Skin, pch_Nasal, pch_Oral, pch_Gut, pch_UG),
    pt.bg = c(col_Skin, col_Nasal, col_Oral, col_Gut, col_UG), col = c1, bty = "n", cex = 0.9)
  title("Community reference proteomes", font.main = 1)
  label.figure("C", font = 2, cex = 1.8, xfrac = 0.12, yfrac = 0.92)

  # Read amino acid composition of proteins inferred from metagenome
  MG_aafile <- file.path(getdatadir(), "ARAST/HMP12_aa.csv")
  MG_aa <- read.csv(MG_aafile)
  # Exclude samples with low numbers of proteins
  ilow <- MG_aa$chains / MG_aa$ref * 100 < cutoff
  ilow[is.na(ilow)] <- TRUE
  MG_aa[ilow, ] <- NA
  # Put into order of metadata table
  iMG <- match(metadata$Metagenome, MG_aa$protein)
  MG_aa <- MG_aa[iMG, ]
  plot(nO2(MG_aa), nH2O(MG_aa), xlim = xlim, ylim = ylim, xlab = chemlab("nO2"), ylab = chemlab("nH2O"), pch = pch, bg = bg, col = col)
  title("Metagenomes - WITH screening", font.main = 1)

  if(pdf) dev.off()

}

# Chemical variation of microbial proteins across body sites and multi-omics comparison 20221125
microhum_2 <- function(pdf = FALSE) {

  if(pdf) pdf("Figure_2.pdf", width = 7, height = 6)
  par(mfrow = c(2, 2))
  par(mgp = c(2.5, 1, 0))

  ## Panel A: Community reference proteomes for body sites
  ## based on 16S rRNA gene sequences from Boix-Amoros et al. (2021) 20220814
  xlim <- c(-0.77, -0.59)
  ylim <- c(-0.8, -0.72)

  # nH2O-nO2 plot for all samples
  par(mar = c(4, 4, 3, 1))
  # Plot data for "no treatment" samples
  BodySites <- plotmet_microhum("BPB+21_NoTreatment", extracolumn = c("Subject", "Site"),
    title = FALSE, pt.open.col = NA, xlim = xlim, ylim = ylim)
  # Keep samples for subjects (not controls)
  BodySites <- BodySites[BodySites$Site != "Control", ]
  # Draw convex hull around all samples 20221125
  BShull <- chull(BodySites$nO2, BodySites$nH2O)
  polygon(BodySites$nO2[BShull], BodySites$nH2O[BShull], border = 8, lty = 2)

  # Plot p-values 20230204
  # Use metrics for untreated samples
  met <- getmetrics_microhum("BPB+21_NoTreatment")
  met <- met[match(BodySites$Run, met$Run), ]
  # Use metrics for gut and oral samples
  igut <- BodySites$Site == "feces"
  ioral <- BodySites$Site == "Oral cavity"
  plot.p.values(met$nO2[igut], met$nO2[ioral], met$nH2O[igut], met$nH2O[ioral], ypos = "bottom")

  # Make legend
  legend("topright", legend = c("Skin", "Nasal", "Oral", "Gut"), pch = c(pch_Skin, pch_Nasal, pch_Oral, pch_Gut),
    pt.bg = c(col_Skin, col_Nasal, col_Oral, col_Gut), col = NA, bty = "n", cex = 0.8)
  title(hyphen.in.pdf("Community reference proteomes\n(data from Boix-Amor\u00f3s et al., 2021)"), font.main = 1)
  label.figure("A", font = 2, cex = 1.8, yfrac = 0.97)

  ## Panel B: Community reference proteomes for controls in COVID-19 datasets 20220822
  # Setup plot
  plot(xlim, ylim, xlab = canprot::cplab$nO2, ylab = canprot::cplab$nH2O, type = "n")
  # Colors and point symbols for sample types
  col <- list(oro = col_Oral, naso = col_Nasal, gut = col_Gut)
  col <- lapply(col, add.alpha, "d0")
  pch <- list(oro = pch_Oral, naso = pch_Nasal, gut = pch_Gut)
  # Read precomputed mean values 20220823
  means <- read.csv(file.path(getdatadir(), "16S/dataset_metrics.csv"))
  # Loop over body sites
  for(type in c("gut", "oro", "naso")) {
    itype <- means$type == type
    # Plot points for control samples
    points(means$nO2_dn[itype], means$nH2O_dn[itype], pch = pch[[type]], bg = col[[type]], col = NA)
  }
  # Plot convex hull from Panel A
  polygon(BodySites$nO2[BShull], BodySites$nH2O[BShull], border = 8, lty = 2)
  # Plot p-values 20230204
  gut <- subset(means, type == "gut")
  oro <- subset(means, type == "oro")
  plot.p.values(gut$nO2_dn, oro$nO2_dn, gut$nH2O_dn, oro$nH2O_dn, ypos = "bottom")
  legend("topright", legend = c("Nasopharyngeal", "Oropharyngeal", "Gut"), pch = c(pch_Nasal, pch_Oral, pch_Gut),
    pt.bg = c(col_Nasal, col_Oral, col_Gut), col = NA, bty = "n", cex = 0.8)
  title(hyphen.in.pdf("Community reference proteomes\n(controls in COVID-19 studies)"), font.main = 1)
  label.figure("B", font = 2, cex = 1.8, yfrac = 0.96)

  ## Panel C: Metagenomes from different body sites 20221124
  ylim <- c(-0.84, -0.60)

  # Setup plot
  plot(xlim, ylim, xlab = canprot::cplab$nO2, ylab = canprot::cplab$nH2O, type = "n")
  # Studies are for gut, oral, nasal
  studies <- c("ZZL+20", "CZH+22", "LLZ+21")
  pchs <- c(pch_Gut, pch_Oral, pch_Nasal)
  cols <- c(col_Gut, col_Oral, col_Nasal)
  cols <- sapply(cols, add.alpha, "b0")
  for(i in 1:length(studies)) {
    file <- file.path(getdatadir(), "ARAST", paste0(studies[i], "_aa.csv"))
    dat <- read.csv(file)
    nO2 <- nO2(dat)
    nH2O <- nH2O(dat)
    points(nO2, nH2O, pch = pchs[i], bg = cols[i], col = NA)
    if(studies[i] == "ZZL+20") gut <- list(nO2 = nO2, nH2O = nH2O)
    if(studies[i] == "CZH+22") oro <- list(nO2 = nO2, nH2O = nH2O)
  }
  # Plot convex hull from Panel A
  polygon(BodySites$nO2[BShull], BodySites$nH2O[BShull], border = 8, lty = 2)
  # Plot p-values 20230204
  plot.p.values(gut$nO2, oro$nO2, gut$nH2O, oro$nH2O, ypos = "bottom")
  # Add legend
  legend <- c("Nasopharyngeal", "Oropharyngeal", "Gut")
  legend("topleft", legend, pch = rev(pchs), pt.bg = c(col_Nasal, col_Oral, col_Gut), col = NA, bty = "n", cex = 0.8)
  legend <- c("Liu'21", "de Castilhos'22", "Zuo'20")
  legend("topright", legend, bty = "n", cex = 0.8)
  title(hyphen.in.pdf("Proteins from metagenomes\n(controls and COVID-19 patients)"), font.main = 1)
  label.figure("C", font = 2, cex = 1.8, yfrac = 0.96)

  ## Panel D: Metaproteomes from various body sites 20221114
  # Setup plot
  plot(xlim, ylim, xlab = canprot::cplab$nO2, ylab = canprot::cplab$nH2O, type = "n")
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
  gut.nO2 <- gut.nH2O <- oral.nO2 <- oral.nH2O <- numeric()
  # Loop over studies
  for(i in 1:length(studies_MP)) {
    # Get amino acid composition from metaproteome
    studydir <- strsplit(studies_MP[i], "_")[[1]][1]
    aa <- read.csv(file.path(getdatadir(), "metaproteome", studydir, paste0(studies_MP[i], "_aa.csv")))
    # Calculate nO2 and nH2O
    nO2 <- nO2(aa)
    nH2O <- nH2O(aa)
    # Add points
    points(nO2, nH2O, pch = pchs[i], bg = cols[i], cex = cexs[i], col = NA)
    if(studies_MP[i] %in% c("TWC+22", "MLL+17")) gut.nO2 <- c(gut.nO2, nO2)
    if(studies_MP[i] %in% c("TWC+22", "MLL+17")) gut.nH2O <- c(gut.nH2O, nH2O)
    if(studies_MP[i] %in% c("GNT+21_cells", "JZW+22")) oral.nO2 <- c(oral.nO2, nO2)
    if(studies_MP[i] %in% c("GNT+21_cells", "JZW+22")) oral.nH2O <- c(oral.nH2O, nH2O)
  }
  # Plot convex hull from Panel A
  polygon(BodySites$nO2[BShull], BodySites$nH2O[BShull], border = 8, lty = 2)
  # Plot p-values 20230204
  plot.p.values(gut.nO2, oral.nO2, gut.nH2O, oral.nH2O, ypos = "bottom")
  # Add legend
  legend("topright", c("Jiang'22", "Granato'21"), title = "Oral",
    pch = c(pch_Oral, pch_Oral), pt.bg = col_Oral, col = NA,
    pt.cex = c(0.7, 1), cex = 0.8, bg = "transparent")
  legend("topleft", c(hyphen.in.pdf("Thuy-Boun'22"), "Maier'17"), title = "Gut", title.adj = 0.4,
    pch = c(pch_Gut, 25), pt.bg = col_Gut, col = NA,
    cex = 0.8, bg = "transparent", inset = c(0, 0.2))
  title(hyphen.in.pdf("Metaproteomes\n(controls and patients, not COVID-19)"), font.main = 1)
  label.figure("D", font = 2, cex = 1.8, yfrac = 0.96)

  if(pdf) dev.off()

}

# Differences of chemical metrics between controls and COVID-19/IBD patients 20220806
microhum_3 <- function(pdf = FALSE) {

  # Start plot
  if(pdf) pdf("Figure_3.pdf", width = 13, height = 12)
  mat <- matrix(c(1,1,1,1, 2,2,2,2, 3,3,3,3, 4,4,4, 5,5,5, 6,6,6, 7,7,7, 0,0,0,0, 8,8,8,8, 0,0,0,0), nrow = 3, byrow = TRUE)
  layout(mat)

  # Define plot settings
  par(cex = 1.2)
  par(mgp = c(2.5, 1, 0))
  startplot <- function() {
    plot(c(-0.010, 0.010), c(-0.012, 0.020), type  = "n", pch = ".", xlab = cplab$DnO2, ylab = cplab$DnH2O)
    abline(h = 0, v = 0, lty = 2, col = 8)
  }
  # Colors and point symbols for sample types
  col <- list(naso = col_Nasal, oro = col_Oral, gut = col_Gut, IBD = col_IBD)
  pch <- list(naso = pch_Nasal, oro = pch_Oral, gut = pch_Gut, IBD = pch_IBD)
  # Read precomputed mean values 20220823
  means <- read.csv(file.path(getdatadir(), "16S/dataset_metrics.csv"))

  ## Panel A: nH2O-nO2 plots for nasopharyngeal, oral/oropharyngeal, and gut communities
  par(mar = c(4, 4, 3, 1))
  for(type in c("naso", "oro", "gut")) {

    startplot()
    itype <- means$type == type
    label <- 1:sum(itype)

    # Add points
    points(means$D_nO2[itype], means$D_nH2O[itype], pch = pch[[type]], col = col[[type]], bg = col[[type]])
    dx <- rep(0, sum(itype))
    dy <- rep(0.0018, sum(itype))
    if(type == "gut") {
      dy[c(5, 6, 7)] <- -0.0018
      dx[5] <- -0.0002
      dx[6] <- 0.0003
      dx[9] <- 0.0003
    }
    # Label points
    text(means$D_nO2[itype] + dx, means$D_nH2O[itype] + dy, label, cex = 0.8)
    # For gut, plot dropline at mean difference for bacterial metaproteome nO2 20220902
    at <- -0.0052
    if(type == "gut") axis(1, at = at, labels = FALSE, tcl = -3, col = 8, lty = 2, lwd = 2)
    # Plot p-values 20230204
    plot.p.values(means$nO2_dn[itype], means$nO2_up[itype], means$nH2O_dn[itype], means$nH2O_up[itype], paired = TRUE)
    # Add plot title
    titles <- c(naso = "Nasopharyngeal", oro = "Oropharyngeal", gut = hyphen.in.pdf("Gut (COVID-19)"))
    title(titles[type], font.main = 1, line = 0.5)
    # Add panel title
    start <- hyphen.in.pdf("A. Community reference proteomes (")
    covid <- hyphen.in.pdf("COVID-19")
    label <- bquote(bold(.(start) * Delta == .(covid) ~ "minus control)"))
    if(type == "naso") label.figure(label, font = 2, xfrac = 0.8, yfrac = 0.95, cex = 1.1)
  }

  # Common nH2O and nO2 limits for boxplots in Panels B and C
  nH2Olim <- c(-0.83, -0.68)
  nO2lim <- c(-0.85, -0.55)

  ## Panel B: Boxplots for nH2O and nO2 in fecal MAGs 20221029
  par(mar = c(4, 4, 3, 1))
  # Read amino acid compositions of proteins predicted from MAGs
  aa <- read.csv(file.path(getdatadir(), "KWL22/KWL22_MAGs_prodigal_aa.csv.xz"))
  # https://github.com/Owenke247/COVID-19/blob/main/Pre-processed_Files/COVID19_metadata.txt
  dat <- read.csv(file.path(getdatadir(), "KWL22/COVID19_metadata.txt"), sep = "\t")

  # Set SRA run prefix to choose study:
  # SRR1232: Zuo et al. (PRJNA624223)
  # SRR1307: Yeoh et al. (PRJNA650244)
  SRAprefix <- "SRR1307"
  # Loop over nH2O and nO2
  for(metric in c("nH2O", "nO2")) {
    # Get amino acid compositions for this BioProject
    iaa <- grep(SRAprefix, aa$protein)
    thisaa <- aa[iaa, ]
    # Calculate nO2 or nH2O
    if(metric == "nO2") x <- nO2(thisaa) else x <- nH2O(thisaa)
    if(metric == "nO2") ylim <- nO2lim else ylim <- nH2Olim
    if(metric == "nO2") legend.x <- "topleft" else legend.x <- "bottomleft"
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
    if(x_pvalue < 0.05) legend <- bquote(bolditalic(p) == bold(.(format(signif(x_pvalue, 1), scientific = 2))))
    legend(legend.x, legend = legend, bty = "n", inset = c(-0.1, -0.05))
    # Show mean difference
    x_diff <- mean(x_list[[2]]) - mean(x_list[[1]])
    diffval <- signif(x_diff, 2)
    if(metric == "nO2") difftxt <- bquote(Delta*italic(n)[O[2]] == .(diffval))
    if(metric == "nH2O") difftxt <- bquote(Delta*italic(n)[H[2]*O] == .(diffval))
    title(difftxt, line = 0.7)
    figlab <- hyphen.in.pdf("B. Gut MAGs (Ke et al., 2022; Yeoh et al., 2021)")
    if(metric == "nH2O" & SRAprefix == "SRR1307") label.figure(figlab, font = 2, xfrac = 0.77, cex = 1.1)
  }

  ## Panel C: boxplots of nO2 and nH2O of bacterial metaproteome in control and COVID-19 patients 20220830
  # Limit to bacterial taxonomy
  taxonomy <- "Bacteria"
  if(taxonomy == "All") aa <- read.csv(file.path(getdatadir(), "metaproteome/HZX+21/HZX+21_aa.csv"))
  if(taxonomy == "Bacteria") aa <- read.csv(file.path(getdatadir(), "metaproteome/HZX+21/HZX+21_bacteria_aa.csv"))
  # Identify control and COVID-19 patients
  icontrol <- grep("Ctrl", aa$organism)
  icovid <- grep("P", aa$organism)
  # Calculate nO2 and nH2O
  nO2 <- nO2(aa)
  nH2O <- nH2O(aa)
  nO2_list <- list(Control = nO2[icontrol], "COVID-19" = nO2[icovid])
  nH2O_list <- list(Control = nH2O[icontrol], "COVID-19" = nH2O[icovid])
  names(nO2_list)[2] <- hyphen.in.pdf(names(nO2_list)[2])
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
  if(nH2O_pvalue < 0.05) legend <- bquote(bolditalic(p) == bold(.(format(signif(nH2O_pvalue, 1), scientific = 2))))
  legend("bottomleft", legend = legend, bty = "n", inset = c(-0.1, -0.05))
  # Show mean difference
  nH2O_diff <- mean(nH2O_list[[2]]) - mean(nH2O_list[[1]])
  diffval <- signif(nH2O_diff, 2)
  difftxt <- bquote(Delta*italic(n)[H[2]*O] == .(diffval))
  title(difftxt, line = 0.7)
  label.figure("C. Gut metaproteomes (He et al., 2021)", font = 2, xfrac = 0.6, cex = 1.1)

  # Add number of samples to group names
  len <- sapply(nO2_list, length)
  labels <- paste0(names(nO2_list), " (", len, ")")
  names(nO2_list) <- ""
  # Make nO2 plot
  boxplot(nO2_list, ylab = cplab$nO2, ylim = nO2lim)
  # Make rotated labels (modified from https://www.r-bloggers.com/rotated-axis-labels-in-r-plots/)
  text(x = (1:2)+0.5, y = par()$usr[3] - 1.5 * strheight("A"), labels = labels, srt = 15, adj = 1, xpd = TRUE)
  # Add p-value
  nO2_pvalue <- wilcox.test(nO2_list[[1]], nO2_list[[2]])$p.value
  legend <- bquote(italic(p) == .(format(signif(nO2_pvalue, 1), scientific = 2)))
  if(nO2_pvalue < 0.05) legend <- bquote(bolditalic(p) == bold(.(format(signif(nO2_pvalue, 1), scientific = 2))))
  legend("topleft", legend = legend, bty = "n", inset = c(-0.1, -0.05))
  # Show mean difference
  nO2_diff <- mean(nO2_list[[2]]) - mean(nO2_list[[1]])
  diffval <- signif(nO2_diff, 2)
  # Make sure we plotted the ticks in Panel A at the correct locations
  stopifnot(diffval == at)
  difftxt <- bquote(Delta*italic(n)[O[2]] == .(diffval))
  title(difftxt, line = 0.7)
  # Draw arrows to ticks in Panel A
  arrows(0.65, -0.5, 0.2, -0.42, length = 0.2, col = 8, lwd = 1.5, xpd = NA)

  ## Panel D: nH2O-nO2 plots for 16S data for IBD 20230723
  par(mar = c(4, 4, 3, 1))
  type <- "IBD"
  startplot()
  itype <- means$type == type
  label <- 1:sum(itype)
  # Add points
  points(means$D_nO2[itype], means$D_nH2O[itype], pch = pch[[type]], col = col[[type]], bg = col[[type]])
  dx <- rep(0, sum(itype))
  dy <- rep(0.0018, sum(itype))
  # Label points
  text(means$D_nO2[itype] + dx, means$D_nH2O[itype] + dy, label, cex = 0.8)
  # Plot p-values 20230204
  plot.p.values(means$nO2_dn[itype], means$nO2_up[itype], means$nH2O_dn[itype], means$nH2O_up[itype], paired = TRUE)
  # Add plot title
  title("Gut (IBD)", font.main = 1, line = 0.5)
  # Add panel title
  start <- hyphen.in.pdf("D. Community reference proteomes (")
  ibd <- hyphen.in.pdf("IBD")
  label <- bquote(bold(.(start) * Delta == .(ibd) ~ "minus control)"))
  label.figure(label, font = 2, xfrac = 0.52, yfrac = 0.95, cex = 1.1)

  if(pdf) dev.off()
}

# Differences of relative abundances of genera between controls and patients 20231227
microhum_4 <- function(pdf = FALSE) {

  if(pdf) pdf("Figure_4.pdf", width = 10, height = 8.5)

  # First, get the genera with largest abundance differences in individual COVID-19 gut and IBD datasets
  covid <- get_abundance("gut")
  ibd <- get_abundance("IBD")
  # Get the unique genus names
  genus_names <- unique(c(colnames(covid), colnames(ibd)))
  # Then, get the abundance differences for these genera in COVID-19 gut and IBD datasets
  COVID <- get_abundance("gut", genus_names)
  IBD <- get_abundance("IBD", genus_names)

  # Setup plot
  layout(matrix(c(3, 1, 2)), heights = c(0.3, 1, 1))
  par(mar = c(7, 5, 0, 4))
  par(tcl = -0.3)
  par(las = 1)
  par(cex.lab = 1.2)
  par(mgp = c(3.2, 1, 0))
  # Define colors and breaks for heatmap
  col <- function(n) hcl.colors(n, "RdYlBu")
  breaks <- c(-.25, -.10, -.05, 0, .05, .10, .25)

  for(disease in c("COVID", "IBD")) {

    # Get differential abundance data
    D_abundance <- get(disease)
    # Convert to matrix to use plot.matrix
    D_abundance <- as.matrix(D_abundance)
    # Get amino acid compositions of reference proteomes for genera and calculate nO2
    AAcomp <- taxon_AA[["GTDB"]]
    AAcomp <- AAcomp[match(colnames(D_abundance), AAcomp$organism), ]
    nO2 <- calc_metrics(AAcomp, "nO2")[, 1]
    # Reorder genera from most reduced to most oxidized
    onO2 <- order(nO2)
    D_abundance <- D_abundance[, onO2]
    nO2 <- nO2[onO2]

    # Truncate values to the range for the color scale
    islo <- D_abundance < min(breaks)
    D_abundance[islo] <- min(breaks)
    ishi <- D_abundance > max(breaks)
    D_abundance[ishi] <- max(breaks)

    # Plot heatmap
    ## This isn't needed because import("plot.matrix") has been added to NAMESPACE
    #requireNamespace("plot.matrix")
    # Temporarily suppress the x-axis
    opar <- par(xaxt = "n")
    plot(D_abundance, col = col, breaks = breaks, main = "", xlab = "", ylab = paste(disease, "Dataset"))
    par(opar)

    # Add triangles to indicate values beyond the color scale
    pch <- c(6, 2)
    islh <- list(islo, ishi)
    for(i in 1:2) {
      isi <- islh[[i]]
      if(any(isi)) {
        xy <- which(isi, arr.ind = TRUE)
        x <- xy[, 2]
        y <- nrow(D_abundance) + 1 - xy[, 1]
        points(x, y, pch = pch[i], col = "white")
      }
    }

    # Find oxygen tolerance of genera
    genus <- colnames(D_abundance)
    oxygen.tolerance <- get.oxytol(genus)
    # Add symbols to indicate obligate anaerobes
    ianaerobe <- oxygen.tolerance == "obligate anaerobe"
    genus[ianaerobe] <- paste(genus[ianaerobe], "*")
    genus[genus == "Bifidobacterium"] <- "Bifidobacterium +"
    labels <- genus
    # Make rotated labels (modified from https://www.r-bloggers.com/rotated-axis-labels-in-r-plots/)
    text(x = seq_along(labels), y = par()$usr[3] - strheight("A"), labels = labels, srt = 40, adj = 1, xpd = TRUE)
    # Add tick marks
    axis(1, at = seq_along(labels), labels = FALSE)
    if(disease == "COVID") {
      # Add legend title
      text(29.5, -1.5, "Change in\nrelative abundance\n(patient - control)", xpd = NA, cex = 1.2)
    }

  }

  # Plot nO2 of genera at top
  par(mar = c(1, 5, 1, 4))
  # Calculate x-axis limits to account for width of boxes and legend in heatmap
  xlim <- c(1 - 0.5, ncol(D_abundance) + 2)
  plot(xlim, range(nO2), xaxs = "i", xaxt = "n", xlab = "", yaxs = "i", ylab = quote(italic(n)[O[2]]~"of genus RP"), type = "n", bty = "n")
  for(i in 1:ncol(D_abundance)) lines(c(i-0.5, i+0.5), rep(nO2[i], 2), lwd = 2)


  if(pdf) dev.off()

}


# Oxygen tolerance of genera in body sites, COVID-19, and IBD 20230726
microhum_5 <- function(pdf = FALSE) {

  # Start plot
  if(pdf) pdf("Figure_5.pdf", width = 10, height = 9)
  par(mfrow = c(3, 4))
  par(mar = c(4, 4, 2.8, 1), mgp = c(2.5, 1, 0))

  # Panels A and B: oxygen tolerance of genera in body sites / in selected COVID-19 and IBD studies
  for(site in c("Nasal", "Oral", "Skin", "Feces", "COVID_control", "COVID", "IBD_control", "IBD")) {
    dat <- calc.oxytol(site)
    plot.oxytol(dat)
    title(site, font.main = 1, line = 0.5)
    btitle <- "A. Oxygen tolerance of genera in body sites (data from Boix-Amor\u00f3s et al., 2021)"
    if(site == "Nasal") label.figure(btitle, font = 2, xfrac = 1.29, yfrac = 0.97, cex = 1.5)
    ctitle <- "B. Oxygen tolerance of genera in selected COVID-19 and IBD studies (data from Schult et al., 2022 and Weng et al., 2019)"
    if(site == "COVID_control") label.figure(hyphen.in.pdf(ctitle), font = 2, xfrac = 1.92, yfrac = 0.97, cex = 1.5)
  }

  # Panel C: Percent of aerotolerant genera in COVID-19 or IBD vs controls 20230726
  startplot <- function(ylab) {
    plot(c(0, 100), c(0, 100), type = "n", xlab = "Controls (% aerotolerant)", ylab = ylab)
    lines(c(0, 100), c(0, 100), lty = 2, col = 8)
  }
  # Colors and point symbols for sample types
  col <- list(naso = col_Nasal, oro = col_Oral, gut = col_Gut, IBD = col_IBD)
  pch <- list(naso = pch_Nasal, oro = pch_Oral, gut = pch_Gut, IBD = pch_IBD)
  # Read precomputed values
  metrics <- read.csv(file.path(getdatadir(), "16S/dataset_metrics.csv"))

  # Calculate percentage of aerotolerant genera among those with known oxygen tolerance 20230726
  control <- metrics$control_aerotolerant / (metrics$control_aerotolerant + metrics$control_anaerobe) * 100
  disease <- metrics$disease_aerotolerant / (metrics$disease_aerotolerant + metrics$disease_anaerobe) * 100

  # Loop over sample groups
  for(type in c("naso", "oro", "gut", "IBD")) {
    if(type == "IBD") startplot("IBD (% aerotolerant)") else startplot(hyphen.in.pdf("COVID-19 (% aerotolerant)"))
    itype <- metrics$type == type
    label <- 1:sum(itype)

    # Add points
    points(control[itype], disease[itype], pch = pch[[type]], col = col[[type]], bg = col[[type]])
    dx <- rep(0, sum(itype))
    dy <- rep(4, sum(itype))
    if(type == "naso") {
      dy[9] <- -4
    }
    if(type == "oro") {
      dy[9] <- 0
      dx[9] <- -4
    }
    if(type == "gut") {
      dy[c(3, 9)] <- -4
    }
    if(type == "IBD") {
      dy[7] <- -4
      dy[10] <- 0
      dx[10] <- 4
      dx[4] <- 2
      dx[3] <- -2
    }
    # Label points
    text(control[itype] + dx, disease[itype] + dy, label)
    # Add plot title
    titles <- c(naso = "Nasopharyngeal", oro = "Oropharyngeal", gut = hyphen.in.pdf("Gut (COVID-19)"), IBD = "Gut (IBD)")
    title(titles[type], font.main = 1, line = 0.5)
    # Add panel title
    start <- hyphen.in.pdf("A. Community reference proteomes (")
    covid <- hyphen.in.pdf("COVID-19")
    label <- bquote(bold(.(start) * Delta == .(covid) ~ "minus control)"))

    if(type == "naso") label.figure(
      hyphen.in.pdf("C. Cumulative abundance of aerotolerant genera in COVID-19 and IBD (data sources listed in Table 1)"),
      font = 2, xfrac = 1.6, yfrac = 0.97, cex = 1.5
    )
  }

  if(pdf) dev.off()

}

# Amount of putative human DNA removed from HMP metagenomes in screening step 20231222
microhum_S1 <- function(pdf = FALSE) {
  if(pdf) pdf("Figure_S1.pdf", width = 10, height = 8)
  # Get sequence processing statistics
  statsfile <- file.path(getdatadir(), "ARAST/HMP12_stats.csv")
  stats <- read.csv(statsfile)
  # Get sample metadata and put in same order as processed sequences
  mdat <- getmdat_microhum("HMP12")
  imdat <- match(stats$ID, mdat$Metagenome)
  mdat <- mdat[imdat, ]
  # Calculate number of sequences removed by screening
  n_removed <- stats$scrubbed_sequences - stats$screened_sequences
  # Calculate number of removed sequences as percentage of input sequences
  perc_removed <- n_removed / stats$input_sequences * 100
  # Setup plot
  par(mfrow = c(3, 2))
  par(cex = 0.8)
  # Loop over body sites
  sites <- c("Skin", "Nasal cavity", "Oral cavity", "GI tract", "UG tract")
  col <- c(col_Skin, col_Nasal, col_Oral, col_Gut, col_UG)
  for(i in seq_along(sites)) {
    hist(perc_removed[mdat$"Body site" == sites[i]], breaks = seq(0, 100, 20), col = col[i],
      xlab = "Percentage of sequences removed by human DNA screening", main = sites[i])
  }
  if(pdf) dev.off()
}

# Changes of chemical metrics for microbiomes associated with viral inactivation 20221125
microhum_S2 <- function(pdf = FALSE) {

  if(pdf) pdf("Figure_S2.pdf", width = 6, height = 4)
  par(mfrow = c(2, 3))
  par(mgp = c(2.5, 1, 0))
  par(mar = c(4, 4, 2, 1))

  ## Community reference proteomes for body sites and viral inactivation
  ## based on 16S rRNA gene sequences from Boix-Amoros et al. (2021) 20220814

  # Get data for any treatment
  Any <- plotmet_microhum("BPB+21_AnyTreatment", extracolumn = c("Subject", "Site", "Treatment"), plot.it = FALSE)
  # Get data for no treatment
  No <- plotmet_microhum("BPB+21_NoTreatment", extracolumn = c("Subject", "Site"), plot.it = FALSE)

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
    # Calculate D_nH2O and D_nO2
    D_nH2O <- Treated$nH2O - No$nH2O
    D_nO2 <- Treated$nO2 - No$nO2
    # Plot D_nH2O and D_nO2
    plot(D_nO2, D_nH2O, xlab = canprot::cplab$DnO2, ylab = canprot::cplab$DnH2O, xlim = c(-0.03, 0.03), ylim = c(-0.03, 0.015), type = "n")
    abline(h = 0, v = 0, lty = 2, col = 8)
    points(D_nO2, D_nH2O, pch = pch, bg = col, col = NA)
    # Plot p-values 20230204
    plot.p.values(No$nO2, Treated$nO2, No$nH2O, Treated$nH2O, paired = TRUE)
    title(Treatment, font.main = 1)
  }

  plottreated("Ethanol")
  plottreated("Formaldehyde")
  plottreated("Heat")
  # Make legend
  plot.new()
  legend("center", legend = c("Skin", "Nasal", "Oral", "Gut"), pch = c(pch_Skin, pch_Nasal, pch_Oral, pch_Gut),
    pt.bg = c(col_Skin, col_Nasal, col_Oral, col_Gut), col = NA, bty = "n", cex = 1.5)
  plottreated("Psoralen")
  plottreated("Trizol")

  if(pdf) dev.off()

}

# Differences of Zc and nH2O between obligate anaerobic and aerotolerant genera 20221017
# Changed Zc to nO2 20230729
microhum_S3 <- function(pdf = FALSE) {

  # Setup plot
  if(pdf) pdf("Figure_S3.pdf", width = 7, height = 4)
  par(mfrow = c(1, 2))
  par(mar = c(3, 4, 1, 1))
  par(mgp = c(2.5, 1, 0))

  # Read table from Million and Raoult, 2018
  dat <- read.csv(system.file("extdata/microhum/MR18_Table_S1_modified.csv", package = "JMDplots"))
  # Drop organisms with unknown oxygen tolerance
  dat <- dat[dat$Obligate.anerobic.prokaryote %in% c(0, 1, 2), ]
  # Find obligate anaerobes and aerotolerant genera
  isAnaerobe <- dat$Obligate.anerobic.prokaryote %in% c(1, 2)
  isAerotolerant <- dat$Obligate.anerobic.prokaryote == 0
  # Make sure no genus has conflicting assignments
  stopifnot(length(intersect(dat$Genus.name[isAnaerobe], dat$Genus.name[isAerotolerant])) == 0)

  # Use GTDB-based reference proteomes
  refdb <- "GTDB"
  # Get amino acid compositions for genera
  aa <- taxon_AA[[refdb]]
  aa <- aa[aa$protein == "genus", ]
  ngenera_GTDB_all <- nrow(aa)
  # Remove alphabetical suffixes
  aa$organism <- sapply(strsplit(aa$organism, "_"), "[", 1)
  # Keep genera with known oxygen tolerance
  aa <- aa[aa$organism %in% dat$Genus.name, ]
  ngenera_GTDB_known <- nrow(aa)
  # Make data frames for obligate anaerobes and aerotolerant genera
  aaAnaerobe <- aa[aa$organism %in% dat$Genus.name[isAnaerobe], ]
  aaAerotolerant <- aa[aa$organism %in% dat$Genus.name[isAerotolerant], ]

  # Print numbers of genera
  print(paste("Reference proteomes are available for", sum(dat$Genus.name[isAnaerobe] %in% aaAnaerobe$organism),
    "of", sum(isAnaerobe), "obligately anaerobic genera in the list"))
  print(paste("Reference proteomes are available for", sum(dat$Genus.name[isAerotolerant] %in% aaAerotolerant$organism),
    "of", sum(isAerotolerant), "aerotolerant genera in the list"))
  print(paste("Oxygen tolerance is listed for", ngenera_GTDB_known, "of", ngenera_GTDB_all, "genera in GTDB"))

  # Loop over nO2 and nH2O
  for(metric in c("nO2", "nH2O")) {

    # Calculate nO2 or nH2O
    if(metric == "nO2") {
      # Get values for obligate anaerobes and aerotolerant genera
      values <- list(
        Anaerobe = nO2(aaAnaerobe),
        Aerotolerant = nO2(aaAerotolerant)
      )
    }
    if(metric == "nH2O") {
      values <- list(
        Anaerobe = nH2O(aaAnaerobe),
        Aerotolerant = nH2O(aaAerotolerant)
      )
    }
    boxplot(values, ylab = cplab[[metric]], ylim = c(-0.83, -0.59), col = c("palevioletred", "lightblue"))
    # Calculate p-value
    pval <- wilcox.test(values$Anaerobe, values$Aerotolerant)$p.value
    ptext <- bquote(italic(p) == .(format(signif(pval, 1), scientific = 2)))
    legend("bottomright", legend = ptext, bty = "n", inset = c(0, -0.03))
    # Show mean difference
    md <- format(round(diff(sapply(values, mean, na.rm = TRUE)), 4), scientific = 1)
    if(metric == "nO2") difftxt <- bquote(Delta*italic(n)[O[2]] == .(md))
    if(metric == "nH2O") difftxt <- bquote(Delta*italic(n)[H[2]*O] == .(md))
    legend("topleft", legend = difftxt, bty = "n", inset = c(-0.05, -0.02))

  }

  if(pdf) dev.off()

}

##################################
### Data Processing Functions  ###
##################################

# Calculate mean values of chemical metrics for patients and controls 20220823
# Include abundance of genera grouped according to oxygen tolerance 20230725
dataset_metrics <- function() {

  # Function to calculate mean values of metrics for patients and controls
  getmeans <- function(study) {
    print(study)
    metrics <- getmetrics_microhum(study)
    mdat <- getmdat_microhum(study, metrics)
    # "up" for disease/positive, "down" for control/negative
    iup <- sapply(mdat$metadata$pch == 25, isTRUE)
    idn <- sapply(mdat$metadata$pch == 24, isTRUE)
    # Calculate means of chemical metrics
    Zc_dn <- mean(mdat$metrics$Zc[idn])
    Zc_up <- mean(mdat$metrics$Zc[iup])
    nO2_dn <- mean(mdat$metrics$nO2[idn])
    nO2_up <- mean(mdat$metrics$nO2[iup])
    nH2O_dn <- mean(mdat$metrics$nH2O[idn])
    nH2O_up <- mean(mdat$metrics$nH2O[iup])
    # Calculate p-values 20220905
    Zc_pvalue <- wilcox.test(mdat$metrics$Zc[idn], mdat$metrics$Zc[iup])$p.value
    nO2_pvalue <- wilcox.test(mdat$metrics$nO2[idn], mdat$metrics$nO2[iup])$p.value
    nH2O_pvalue <- wilcox.test(mdat$metrics$nH2O[idn], mdat$metrics$nH2O[iup])$p.value
    # Calculate oxygen tolerance 20230725
    disease <- calc.oxytol(study = study)
    disease_anaerobe <- sum(disease$abundance[disease$oxygen.tolerance == "obligate anaerobe"])
    disease_aerotolerant <- sum(disease$abundance[disease$oxygen.tolerance == "aerotolerant"])
    disease_unknown <- sum(disease$abundance[disease$oxygen.tolerance == "unknown"])
    control <- calc.oxytol("control", study = study)
    control_anaerobe <- sum(control$abundance[control$oxygen.tolerance == "obligate anaerobe"])
    control_aerotolerant <- sum(control$abundance[control$oxygen.tolerance == "aerotolerant"])
    control_unknown <- sum(control$abundance[control$oxygen.tolerance == "unknown"])
    # Include number of samples 20220905
    data.frame(n_dn = sum(idn), n_up = sum(iup),
         Zc_dn = Zc_dn, Zc_up = Zc_up, Zc_pvalue = Zc_pvalue, 
         nO2_dn = nO2_dn, nO2_up = nO2_up, nO2_pvalue = nO2_pvalue, 
         nH2O_dn = nH2O_dn, nH2O_up = nH2O_up, nH2O_pvalue = nH2O_pvalue,
         control_anaerobe = control_anaerobe, disease_anaerobe = disease_anaerobe, 
         control_aerotolerant = control_aerotolerant, disease_aerotolerant = disease_aerotolerant,
         control_unknown = control_unknown, disease_unknown = disease_unknown
    )
  }

  out <- list()

  # List COVID-19 and IBD datasets
  microhum_studies <- list(
    # COVID-19 nasopharyngeal
    naso = c("ENJ+21", "PMM+22", "SGC+21", "HMH+21", "VCV+21", "SRS+22", "CSC+22", "GKJ+22", "MLW+21_Nasopharyngeal"),
    # COVID-19 oral/oropharyngeal
    oro = c("RFH+22_Oral", "IZC+21", "GBS+22", "WCJ+21_Oral", "XLZ+21", "MAC+21", "MLW+21_Oropharyngeal", "GWL+21", "RWC+21_Oral"),
    # COVID-19 gut
    gut = c("ZZZ+21", "RFH+22_Gut", "KMG+21", "WCJ+21_Gut", "CGC+22", "GCW+20", "MMP+21", "NGH+21", "RDM+22", "MIK+22", "FBD+22", "RWC+21_Gut", "SRK+22"),
    # IBD gut
    IBD = c("TWC+22", "ZTG+21", "ASM+23", "AAM+20", "LZD+19", "MDV+22", "LAA+19", "GKD+14", "WGL+19", "RAF+20")
  )
  # Loop over groups of datasets
  for(i in 1:4) {
    # Calculate means for each dataset
    means <- lapply(microhum_studies[[i]], getmeans)
    means <- do.call(rbind, means)
    means <- cbind(type = names(microhum_studies)[i], study = microhum_studies[[i]], means)
    out[[i]] <- means
  }
  # Put together data frames
  out <- do.call(rbind, out)

  # Calculate mean differences of chemical metrics
  D_Zc <- out$Zc_up - out$Zc_dn
  D_nO2 <- out$nO2_up - out$nO2_dn
  D_nH2O <- out$nH2O_up - out$nH2O_dn
  out <- cbind(out, D_Zc, D_nO2, D_nH2O)

  # Round values 20230212
  out[, 5:22] <- signif(out[, 5:22], 6)
  file <- "dataset_metrics.csv"
  write.csv(out, file, row.names = FALSE, quote = FALSE)

}

#############################
### Unexported Functions  ###
#############################

# Function to add p-values to x and y axes 20230204
plot.p.values <- function(X.1, X.2, Y.1, Y.2, paired = FALSE, ypos = "top") {
  p.value.X <- wilcox.test(X.1, X.2, paired = paired)$p.value
  p.value.Y <- wilcox.test(Y.1, Y.2, paired = paired)$p.value
  p.X <- format(signif(p.value.X, 1), scientific = 2)
  p.Y <- format(signif(p.value.Y, 1), scientific = 2)
  p.X.txt <- bquote(italic(p) == .(p.X))
  p.Y.txt <- bquote(italic(p) == .(p.Y))
  # Use bold text for significant p-values 20230723
  if(p.value.X < 0.05) p.X.txt <- bquote(bolditalic(p) == bold(.(p.X)))
  if(p.value.Y < 0.05) p.Y.txt <- bquote(bolditalic(p) == bold(.(p.Y)))
  pu <- par("usr")
  dx <- (pu[2] - pu[1]) / 30
  dy <- (pu[4] - pu[3]) / 30
  text(pu[2] - dx, pu[3] + dy/2, p.X.txt, adj = c(1, 0), cex = 0.9)
  if(ypos == "bottom") text(pu[1] + dx/2, pu[3] + dy, p.Y.txt, srt = 90, adj = c(0, 1), cex = 0.9)
  if(ypos == "top") text(pu[1] + dx/2, pu[4] - dy, p.Y.txt, srt = 90, adj = c(1, 1), cex = 0.9)
}

# Get oxygen tolerance for specified genera 20231227
get.oxytol <- function(genus) {

  # Read the table of oxygen tolerance for genera
  dat <- read.csv(system.file("extdata/microhum/MR18_Table_S1_modified.csv", package = "JMDplots"))
  obligate.anaerobe <- dat$Genus.name[dat$Obligate.anerobic.prokaryote %in% c(1, 2)]
  aerotolerant <- dat$Genus.name[dat$Obligate.anerobic.prokaryote == 0]
  # Remove suffixes used in GTDB (_A, _B, etc.)
  genus <- sapply(strsplit(genus, "_"), "[", 1)
  # Assign oxygen tolerance
  oxygen.tolerance <- rep("unknown", length(genus))
  oxygen.tolerance[genus %in% obligate.anaerobe] <- "obligate anaerobe"
  oxygen.tolerance[genus %in% aerotolerant] <- "aerotolerant"
  oxygen.tolerance

}

# Summarize abundances of genera in segments (body sites, disease, or control) and list their oxygen tolerance 20230725
calc.oxytol <- function(segment = "Feces", study = NULL) {
  
  if(tolower(segment) %in% c("feces", "nasal", "oral", "skin")) {
    # Identify RDP file
    RDPfile <- file.path(getdatadir(), "16S/RDP-GTDB/BPB+21.tab.xz")
    # Read metadata
    mdat <- getmdat_microhum("BPB+21_NoTreatment")
    # Get run IDs for this body segment
    Run <- mdat$Run[grep(tolower(segment), tolower(mdat$Site))]
  }

  if(tolower(segment) %in% c("covid", "covid_control")) {
    RDPfile <- file.path(getdatadir(), "16S/RDP-GTDB/SRK+22.tab.xz")
    mdat <- getmdat_microhum("SRK+22")
    if(tolower(segment) == "covid") Run <- mdat$Run[grep("Covid19", mdat$Status)]
    if(tolower(segment) == "covid_control") Run <- mdat$Run[grep("Control", mdat$Status)]
  }

  if(tolower(segment) %in% c("ibd", "ibd_control")) {
    RDPfile <- file.path(getdatadir(), "16S/RDP-GTDB/WGL+19.tab.xz")
    mdat <- getmdat_microhum("WGL+19")
    if(tolower(segment) == "ibd") Run <- mdat$Run[mdat$Disease != "HC"]
    if(tolower(segment) == "ibd_control") Run <- mdat$Run[mdat$Disease == "HC"]
  }

  # Get classifications and metadata for any COVID-19 or IBD study
  if(!is.null(study)) {
    # Remove suffix after underscore 20200929
    studyfile <- gsub("_.*", "", study)
    RDPfile <- file.path(getdatadir(), "16S/RDP-GTDB", paste0(studyfile, ".tab.xz"))
    mdat <- getmdat_microhum(study)
    if(tolower(segment) == "control") Run <- mdat$Run[mdat$pch == 24] else Run <- mdat$Run[mdat$pch == 25]
  }

  # Read RDP file
  RDP <- read_RDP(RDPfile, quiet = TRUE)
  # Keep only genus-level classifications
  RDP <- RDP[RDP$rank == "genus", ]
  # Use genus names for row names
  rownames(RDP) <- RDP$name
  # Keep only selected samples
  Run <- Run[Run %in% colnames(RDP)]
  RDP <- RDP[, Run]
  # Remove genera with zero counts
  RDP <- RDP[rowSums(RDP) > 0, ]
  # Sum abundances across samples and normalize
  abundance <- rowSums(RDP)
  abundance <- abundance / sum(abundance) * 100

  # Get amino acid compositions of reference proteomes for genera and calculate nO2
  AAcomp <- taxon_AA[["GTDB"]]
  AAcomp <- AAcomp[match(names(abundance), AAcomp$organism), ]
  nO2 <- calc_metrics(AAcomp, "nO2")[, 1]
  # Get oxygen tolerance
  genus <- rownames(RDP)
  oxygen.tolerance <- get.oxytol(genus)

  # Make summary table
  data.frame(oxygen.tolerance, abundance, nO2)

}

# Line plot: genus abundance vs nO2, grouped by oxygen tolerance 20230725
plot.oxytol <- function(dat) {
  # Start with 0 total abundance
  total.abundance <- 0
  plot(c(-0.81, -0.595), c(0, 100), xlab = chemlab("nO2"), ylab = "Abundance (%)", type = "n")
  # List classifications and colors
  oxytols <- c("obligate anaerobe", "aerotolerant", "unknown")
  cols <- c(2, 4, 8)
  # Loop over classifications
  for(i in 1:3) {
    # Get genera with this oxygen tolerance
    thisdat <- dat[dat$oxygen.tolerance == oxytols[i], ]
    # Order by nO2
    thisdat <- thisdat[order(thisdat$nO2), ]
    # Get y values
    y2 <- total.abundance + cumsum(thisdat$abundance)
    y1 <- head(c(total.abundance, y2), -1)
    # Plot a rectangle to encompass entire range
    rect(head(thisdat$nO2, 1), head(y1, 1), tail(thisdat$nO2, 1), tail(y2, 1), col = adjustcolor(cols[i], alpha.f = 0.3), border = NA)
    # Plot a line for the mean
    nO2.mean <- sum(thisdat$nO2 * thisdat$abundance) / sum(thisdat$abundance)
    lines(c(nO2.mean, nO2.mean), c(head(y1, 1), tail(y2, 1)), lwd = 1.5, col = "white")
    # Loop over genera
    for(j in 1:nrow(thisdat)) {
      # Plot a line with length corresponding to abundance of this genus
      lines(c(thisdat$nO2[j], thisdat$nO2[j]), c(y1[j], y2[j]), col = cols[i], lwd = 2)
      # Label lines for most-abundant genera
      genus.abundance <- y2[j] - y1[j]
      if(genus.abundance > 3) {
        ymid <- (y1[j] + y2[j]) / 2
        if(thisdat$nO2[j] > -0.7) adj <- 1 else adj <- 0
        text(thisdat$nO2[j], ymid, paste0(" ", rownames(thisdat)[j], " "), adj = adj)
      }
    }
    # Update the total abundance
    total.abundance <- tail(y2, 1)
 }
}

# Get differences of relative abundance of genera in COVID and IBD
# 20231225 first version jmd
get_abundance <- function(study_type, genus_names = NULL) {

  # List studies with datasets
  metrics <- read.csv(file.path(getdatadir(), "16S/dataset_metrics.csv"))
  studies <- metrics$study[metrics$type == study_type]

  # Loop over studies
  for(istudy in seq_along(studies)) {

    # Remove suffix after underscore
    studyname <- strsplit(studies[istudy], "_")[[1]][1]
    RDPfile <- file.path(getdatadir(), "16S/RDP-GTDB", paste0(studyname, ".tab.xz"))
    mdat <- getmdat_microhum(studies[istudy])

    # Read RDP file
    RDP <- read_RDP(RDPfile, quiet = TRUE)
    # Keep only genus-level classifications
    RDP <- RDP[RDP$rank == "genus", ]
    # Use genus names for row names
    rownames(RDP) <- RDP$name

    # Loop over control and disease
    for(segment in c("control", "disease")) {
      if(segment == "control") in.segment <- sapply(mdat$pch == 24, isTRUE)
      if(segment == "disease") in.segment <- sapply(mdat$pch == 25, isTRUE)
      myRuns <- mdat$Run[in.segment]
      # Remove runs without taxonomic classifications
      # (possibly due to low-count samples discarded by read_RDP)
      myRuns <- myRuns[myRuns %in% colnames(RDP)]
      # Get genus-level classifications for these runs
      myRDP <- RDP[, myRuns]
      # Sum genus abundances for all samples in this segment
      myabundance <- t(rowSums(myRDP))
      if(istudy == 1) abundance <- myabundance else {
        abundance <- get(segment)
        abundance <- merge(abundance, myabundance, all = TRUE, sort = FALSE)
      }
      assign(segment, abundance)
    }

  }

  # Make sure genus names are the same for control and disease
  stopifnot(all(colnames(control) == colnames(disease)))
  # Set NA counts to 0 and normalize counts to sum to 1
  control[is.na(control)] <- 0
  control <- control / rowSums(control)
  disease[is.na(disease)] <- 0
  disease <- disease / rowSums(disease)
  # Calculate differences of abundance of genera between disease and control
  D_abundance_all <- disease - control

  if(!is.null(genus_names)) {
    return(D_abundance_all[, genus_names])
  } else {
    # Find the genera with the largest change
    cutoff <- 0.05
    is_top <- apply(abs(D_abundance_all) > cutoff, 2, any)
    return(D_abundance_all[, is_top])
  }

}
