# JMDplots/MaximAct.R

# 20201218 Calculate optimal logaH2O and logfO2 for target proteins
# 20210215 Use target proteins given in AA_target argument
#   - also add xlab, O2, H2O, nbackground, arguments
# 20210307-08 Add pH and names arguments
# 20210401 Add plot argument
# 20210718 Allow H2O = 0 (unit activity of H2O)
MaximAct <- function(AA_target, seed = 1:100, nbackground = 2000, plot.it = TRUE,
  xlab = "sample", names = NULL, O2 = c(-72, -67), H2O = c(-2, 6), pH = NULL,
  AA_background = NULL) {

  # Load target proteins
  iptarget <- add.protein(AA_target)
  if(is.null(names)) names <- AA_target$protein

  if(is.null(AA_background)) {
    # Get background proteins:
    ## (UniProt reference human proteome) 20210215
    #AA_background <- readRDS(system.file("/extdata/protein/human_base.rds", package = "canprot"))
    # Only use proteins that have phylostrata assignments in both Trigos and Liebeskind datasets 20210402
    TPPG17 <- read.csv(system.file(paste0("extdata/phylostrata/TPPG17.csv.xz"), package = "canprot"), as.is = TRUE)
    LMM16 <- read.csv(system.file(paste0("extdata/phylostrata/LMM16.csv.xz"), package = "canprot"), as.is = TRUE)
    # Take the intersection of lists of proteins from the two sources
    Entry <- na.omit(intersect(TPPG17$Entry, LMM16$UniProt))
    # Get amino acid compositions of the proteins
    AA_background <- protcomp(Entry)$aa
  }

  # Set up system
  if(is.null(pH)) basis("QEC") else basis("QEC+")
  # Initialize output values and plot
  outO2 <- outH2O <- outpH <- list()
  if(plot.it) {
    logH2Olab <- quote(bold(log)*bolditalic(a)[bold(H[2]*O)])
    logO2lab <- quote(bold(log)*bolditalic(f)[bold(O[2])])
    pHlab <- "pH"
    if(is.null(pH)) split.screen(c(2, 1)) else split.screen(c(3, 1))
    screen(1)
    par(mar = c(4, 4, 1, 1), mgp = c(2.5, 1, 0), las = 1)
    plot(range(1:length(names)), O2[1:2], xlab = xlab, ylab = logO2lab, type = "n", xaxt = "n", xaxs = "i", yaxs = "i", font.lab = 2)
    axis(1, 1:length(names), names)
    screen(2)
    par(mar = c(4, 4, 1, 1), mgp = c(2.5, 1, 0), las = 1)
    plot(range(1:length(names)), H2O[1:2], xlab = xlab, ylab = logH2Olab, type = "n", xaxt = "n", xaxs = "i", yaxs = "i", font.lab = 2)
    axis(1, 1:length(names), names)
    abline(h = 0, lty = 3)
    if(!is.null(pH)) {
      screen(3)
      par(mar = c(4, 4, 1, 1), mgp = c(2.5, 1, 0), las = 1)
      plot(range(1:length(names)), pH, xlab = xlab, ylab = pHlab, type = "n", xaxt = "n", xaxs = "i", yaxs = "i", font.lab = 2)
      axis(1, 1:length(names), names)
      abline(h = 7, lty = 3)
    }
  }

  # Prevent multicore usage (runs out of memory)
  thermo("opt$paramin" = 10000)
  # Loop over random seeds
  for(iseed in seq_along(seed)) {
    # Calculate affinities for target proteins and a sample of human proteins (background)
    print(paste("seed is", seed[iseed]))
    set.seed(seed[iseed])
    iback <- 1:nrow(AA_background)
    # Sample background proteins only if number of available proteins is greater than sample size 20210712
    if(nrow(AA_background) > nbackground) iback <- sample(1:nrow(AA_background), nbackground)
    ipback <- add.protein(AA_background[iback, ])
    if(identical(H2O, 0)) {
      # Use unit activity of H2O 20210718
      if(is.null(pH)) a <- suppressMessages(affinity(O2 = O2, iprotein = c(iptarget, ipback)))
      if(!is.null(pH)) a <- suppressMessages(affinity(O2 = c(O2[1:2], 40), pH = c(pH, 40), iprotein = c(iptarget, ipback)))
    } else {
      if(is.null(pH)) a <- suppressMessages(affinity(O2 = O2, H2O = H2O, iprotein = c(iptarget, ipback)))
      if(!is.null(pH)) a <- suppressMessages(affinity(O2 = c(O2[1:2], 40), H2O = c(H2O[1:2], 40), pH = c(pH, 40), iprotein = c(iptarget, ipback)))
    }
    # Equilibrate and find maximum activity for each target protein
    e <- suppressMessages(equilibrate(a, as.residue = TRUE, loga.balance = 0))
    optO2 <- optH2O <- optpH <- numeric()
    for(i in seq_along(names)) {
      imax <- arrayInd(which.max(e$loga.equil[[i]]), dim(e$loga.equil[[i]]))
      optO2 <- c(optO2, e$vals$O2[imax[1]])
      if(identical(H2O, 0)) optH2O <- c(optH2O, 0) else optH2O <- c(optH2O, e$vals$H2O[imax[2]])
      if(!is.null(pH)) optpH <- c(optpH, e$vals$pH[imax[3]])
    }
    if(plot.it) {
      # Plot logfO2 and logaH2O values
      # FIXME: need two screen() calls to make this work 20201218
      screen(1, FALSE); screen(1, FALSE)
      lines(1:length(names), optO2, lwd = 0.5, col = "gray")
      screen(2, FALSE); screen(2, FALSE)
      lines(1:length(names), optH2O, lwd = 0.5, col = "gray")
      if(!is.null(pH)) {
        screen(3, FALSE); screen(3, FALSE)
        lines(1:length(names), optpH, lwd = 0.5, col = "gray")
      }
    }
    # Store results
    outO2[[iseed]] <- optO2
    outH2O[[iseed]] <- optH2O
    if(!is.null(pH)) outpH[[iseed]] <- optpH
  }
  # Gather values
  outO2 <- do.call(rbind, outO2)
  outH2O <- do.call(rbind, outH2O)
  if(plot.it) {
    # Plot mean values
    meanO2 <- colMeans(outO2)
    meanH2O <- colMeans(outH2O)
    screen(1, FALSE); screen(1, FALSE)
    lines(1:length(names), meanO2, col = 2, lwd = 2)
    screen(2, FALSE); screen(2, FALSE)
    lines(1:length(names), meanH2O, col = 2, lwd = 2)
  }

  if(!is.null(pH)) {
    outpH <- do.call(rbind, outpH)
    if(plot.it) {
      meanpH <- colMeans(outpH)
      screen(3, FALSE); screen(3, FALSE)
      lines(1:length(names), meanpH, col = 2, lwd = 2)
    }
  }

  if(plot.it) {
    # Finalize plots
    close.screen(all.screens = TRUE)
  }

  # Create output values
  outO2 <- data.frame(outO2)
  colnames(outO2) <- names
  outO2 <- cbind(seed = seed, outO2)
  outH2O <- data.frame(outH2O)
  colnames(outH2O) <- names
  outH2O <- cbind(seed = seed, outH2O)
  out <- list(O2 = outO2, H2O = outH2O)
  # Handle pH option
  if(!is.null(pH)) {
    outpH <- data.frame(outpH)
    colnames(outpH) <- names
    outpH <- cbind(seed = seed, outpH)
    out <- c(out, list(pH = pH))
  }

  out

}

# Get mean amino acid composition for each phylostratum 20201219
getphyloaa <- function(PS_source) {
  dat <- read.csv(system.file(paste0("extdata/phylostrata/", PS_source, ".csv.xz"), package = "canprot"), as.is = TRUE)
  if(PS_source == "LMM16") {
    colnames(dat)[c(1,3)] <- c("Entry", "Phylostrata")
    # remove entries that have ENSP instead of UniProt IDs
    dat <- dat[!grepl("^ENSP", dat$Entry), ]
  }
  dat <- check_IDs(dat, "Entry")
  dat <- cleanup(dat, "Entry")
  pcomp <- protcomp(dat$Entry)
  # Set up blank amino acid data frame
  PS <- sort(unique(dat$Phylostrata))
  aa <- thermo()$protein[rep(1, length(PS)), ]
  aa$protein <- PS
  aa$organism <- PS_source
  aa$ref <- aa$abbrv <- NA
  aa$chains <- 1
  aa[, 6:25] <- 0
  # Loop over phylostrata
  for(i in seq_along(PS)) {
    iPS <- dat$Phylostrata == PS[i]
    aaPS <- pcomp$aa[iPS, ]
    aamean <- colMeans(aaPS[, 6:25])
    aa[i, 6:25] <- aamean
  }
  # Return both the mean compositions (aa) and all proteins (pcomp)
  list(aa = aa, pcomp = pcomp)
}
