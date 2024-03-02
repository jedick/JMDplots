# JMDplots/canprot.R
# This file has functions previously in the canprot package
# that are used for chemical analysis of differential expression data
# in the gradH2O and canH2O papers.
# 20240426 Functions moved from canprot package

# Compile and view vignettes from command line
# 20200414 jmd
mkvig <- function(vig = NULL) {
  vig.allowed <- gsub("pdat_", "", grep("pdat_", ls("package:JMDplots"), value = TRUE))
  vig.allowed <- vig.allowed[!vig.allowed %in% c("aneuploidy", "fly")]
  isnull <- is.null(vig)
  toomany <- length(vig) > 1
  notallowed <- !any(vig %in% vig.allowed)
  if(isnull | toomany | notallowed) stop("'vig' should be one of: ", paste(vig.allowed, collapse = ", "))
  # Location of the vignette directory
  vigdir <- system.file("vignettes", package = "JMDplots")
  # Names of the vignette source and html output files
  vigfile <- file.path(vigdir, paste0(vig, ".Rmd"))
  htmlfile <- tempfile(pattern = paste0(vig, "_"), fileext = ".html")
  # Compile the vignette and open it in the browser
  rmarkdown::render(vigfile, output_file = htmlfile, knit_root_dir = vigdir)
  # Set 'browser' option so vignettes open under Rstudio 20201017
  # https://stackoverflow.com/questions/62536479/the-command-exams2html-does-not-generate-html-page-when-it-is-run-from-rstudio
  if(.Platform$OS.type == "windows") oldopt <- options(browser = NULL)
  browseURL(htmlfile)
  if(.Platform$OS.type == "windows") options(oldopt)
  # Return the path of html file
  htmlfile
}

# Function to identify known UniProt IDs
# 20160703 jmd
check_IDs <- function(dat, IDcol, aa_file = NULL, updates_file = NULL) {
  # The input IDs that are NA
  input.NA <- is.na(dat[, IDcol]) | dat[, IDcol] == ""
  # The candidate IDs separated into a list
  ID_list <- strsplit(dat[, IDcol], ";")
  # The list of IDs as a vector
  ID <- unlist(ID_list)
  # Get the UniProt ID in case we have e.g. sp|P62308|RUXG_HUMAN
  # for NJVS19 dataset 20191226
  if(any(grepl("\\|", ID))) ID <- sapply(strsplit(ID, "\\|"), "[", 2)
  # Human proteins
  aa <- get("human_aa", canprot)
  # Add amino acid compositions from external file if specified
  if(!is.null(aa_file)) {
    aa_dat <- read.csv(aa_file, as.is=TRUE)
    aa <- rbind(aa_dat, aa)
  }
  # Assemble known IDs
  knownIDs <- sapply(strsplit(aa$protein, "|", fixed = TRUE), "[", 2)
  # If that is NA (i.e. no | separator is present) use the entire string
  ina <- is.na(knownIDs)
  knownIDs[ina] <- aa$protein[ina]
  # Also include obsolete UniProt ID
  updates <- get("uniprot_updates", canprot)
  if(!is.null(updates_file)) {
    updates_dat <- read.csv(updates_file, as.is = TRUE)
    updates <- rbind(updates_dat, updates)
  }
  knownIDs <- c(knownIDs, updates$old)
  # Check if the candidate IDs are known
  known <- match(ID, knownIDs)
  known_IDs <- ID[known > 0]
  # Get the IDs back into list form
  known_IDs <- relist(known_IDs, ID_list)
  # Take the first (non-NA) match 20191119
  ID <- sapply(sapply(known_IDs, na.omit), "[", 1)
  # The output IDs that are NA
  output.NA <- is.na(ID) | ID == ""
  # Print which IDs became NA 20191127
  if(sum(output.NA) > sum(input.NA)) {
    new.NA <- output.NA & !input.NA
    NA.IDs <- dat[new.NA, IDcol]
    if(sum(new.NA)==1) IDtxt <- "ID:" else IDtxt <- "IDs:"
    print(paste("check_IDs:", sum(new.NA), "unavailable UniProt", IDtxt, paste(NA.IDs, collapse = " ")))
  }
  # Now apply the updates 20191119
  iold <- match(ID, updates$old)
  if(any(!is.na(iold))) {
    oldIDs <- updates$old[na.omit(iold)]
    newIDs <- updates$new[na.omit(iold)]
    if(sum(!is.na(iold))==1) IDtxt <- "ID:" else IDtxt <- "IDs:"
    print(paste("check_IDs: updating", sum(!is.na(iold)), "old UniProt", IDtxt, paste(oldIDs, collapse = " ")))
    ID[!is.na(iold)] <- newIDs
  }
  dat[, IDcol] <- ID
  dat
}

# Function to clean up data (remove proteins with unavailable IDs,
# ambiguous expression ratios, and duplicated IDs)
# 20190407 jmd
cleanup <- function(dat, IDcol, up2 = NULL) {

  # Utility function to remove selected entries and print a message 20160713
  remove_entries <- function(dat, irm, description) {
    if(sum(irm) > 0) {
      dat <- dat[!irm, , drop = FALSE]
      if(sum(irm)==1) ptxt <- "protein" else ptxt <- "proteins"
      print(paste("cleanup: removing", sum(irm), description, ptxt))
    }
    return(dat)
  }

  if(is.logical(IDcol)) {
    dat <- remove_entries(dat, IDcol, "selected")
  } else {
    # Add up2 column to dat
    if(!is.null(up2)) dat <- cbind(dat, up2 = up2)
    # Which column has the IDs
    if(is.character(IDcol)) IDcol <- match(IDcol, colnames(dat))
    # Drop NA or "" IDs
    unav <- is.na(dat[, IDcol]) | dat[, IDcol] == ""
    dat <- remove_entries(dat, unav, "unavailable")
    if(!is.null(up2)) {
      # Drop unquantified proteins 20191120
      dat <- remove_entries(dat, is.na(dat$up2), "unquantified")
      # Drop proteins with ambiguous expression ratios
      up <- dat[dat$up2, IDcol]
      down <- dat[!dat$up2, IDcol]
      iambi <- dat[, IDcol] %in% intersect(up, down)
      dat <- remove_entries(dat, iambi, "ambiguous")
    }
    # Drop duplicated proteins
    dat <- remove_entries(dat, duplicated(dat[, IDcol]), "duplicated")
  }
  dat
}

# canprot/R/get_comptab.R
# Merge old Zc_nH2O() and CNS() functions, and add volume 20170718
# CNS: elemental abundance (C, N, S) per residue 20170124
# Zc_nH2O: plot and summarize Zc and nH2O/residue of proteins 20160706
get_comptab <- function(pdat, var1="Zc", var2="nH2O", plot.it=FALSE, mfun="median", oldstyle = FALSE) {
  # Define functions for the possible variables of interest
  # Calculate metrics with canprot functions, not CHNOSZ 20201015
  Zc <- function() canprot::Zc(pdat$pcomp$aa)
  nH2O <- function() canprot::nH2O(pdat$pcomp$aa)
  nO2 <- function() canprot::nO2(pdat$pcomp$aa)
  # Calculate protein length 20201015
  #AA3 <- CHNOSZ::aminoacids(3)
  AA3 <- c("Ala", "Cys", "Asp", "Glu", "Phe", "Gly", "His", "Ile", "Lys", "Leu",
           "Met", "Asn", "Pro", "Gln", "Arg", "Ser", "Thr", "Val", "Trp", "Tyr")
  # The columns for the amino acids
  icol <- match(AA3, colnames(pdat$pcomp$aa))
  nAA <- function() rowSums(pdat$pcomp$aa[, icol])
  # Memoize the volumes so we don't depend on the CHNOSZ database 20200509
  # Changed to canprot::V0 20240301
  V0 <- function() canprot::V0(pdat$pcomp$aa)
  # GRAVY and pI added 20191028
  # canprot:: is used to access the functions in the package namespace, not the ones defined here
  GRAVY <- function() canprot::GRAVY(pdat$pcomp$aa)
  pI <- function() canprot::pI(pdat$pcomp$aa)
  # MW (molecular weight) added 20200501
  MW <- function() canprot::MW(pdat$pcomp$aa)

  # Get the values of the variables using the functions
  val1 <- get(var1)()
  val2 <- get(var2)()
  if(plot.it) {
    # Set symbol shape and color
    # Hollow red square for up, filled blue circle for down in cancer
    col <- ifelse(pdat$up2, 2, 4)
    pch <- ifelse(pdat$up2, 0, 19)
    cex <- ifelse(pdat$up2, 1, 0.8)
    # Shuffle the order of points to mitigate overplotting effects
    i <- sample(1:length(val1))
    plot(val1[i], val2[i], xlab=cplab[[var1]], ylab=cplab[[var2]], col=col[i], pch=pch[i], cex=cex[i])
    title(pdat$description)
    mtext(pdat$dataset, side=4, cex=0.85, las=0, adj=0, line=-0.1)
  }
  # Calculate and print sample size
  val1_dn <- val1[!pdat$up2]
  val1_up <- val1[pdat$up2]
  val2_dn <- val2[!pdat$up2]
  val2_up <- val2[pdat$up2]
  message(paste0(pdat$dataset, " (", pdat$description ,"): n1 ", length(val1_dn), ", n2 ", length(val1_up)))
  # Calculate difference of means/medians, CLES, p-value
  mfun1_dn <- get(mfun)(val1_dn)
  mfun1_up <- get(mfun)(val1_up)
  mfun2_dn <- get(mfun)(val2_dn)
  mfun2_up <- get(mfun)(val2_up)
  val1.diff <- mfun1_up - mfun1_dn
  val2.diff <- mfun2_up - mfun2_dn
  out <- data.frame(dataset=pdat$dataset, description=pdat$description,
    n1=length(val1_dn), n2=length(val1_up),
    val1.median1=mfun1_dn, val1.median2=mfun1_up, val1.diff,
    val2.median1=mfun2_dn, val2.median2=mfun2_up, val2.diff, stringsAsFactors=FALSE)
  if(oldstyle) {
    val1.CLES <- 100*CLES(val1_dn, val1_up, distribution = NA)
    val2.CLES <- 100*CLES(val2_dn, val2_up, distribution = NA)
    val1.p.value <- val2.p.value <- NA
    if(!any(is.na(val1_dn)) & !any(is.na(val1_up))) val1.p.value <- stats::wilcox.test(val1_dn, val1_up)$p.value
    if(!any(is.na(val2_dn)) & !any(is.na(val2_up))) val2.p.value <- stats::wilcox.test(val2_dn, val2_up)$p.value
    # print summary messages
    nchar1 <- nchar(var1)
    nchar2 <- nchar(var2)
    start1 <- paste0(var1, substr("      MD ", nchar1, 10))
    start2 <- paste0(var2, substr("      MD ", nchar2, 10))
    message(paste0(start1, format(round(val1.diff, 3), nsmall=3, width=6),
                 ", CLES ", round(val1.CLES), "%",
                 ", p-value ", format(round(val1.p.value, 3), nsmall=3)))
    message(paste0(start2, format(round(val2.diff, 3), nsmall=3, width=6),
                 ", CLES ", round(val2.CLES), "%",
                 ", p-value ", format(round(val2.p.value, 3), nsmall=3), "\n"))
    out <- data.frame(dataset=pdat$dataset, description=pdat$description,
      n1=length(val1_dn), n2=length(val1_up),
      val1.median1=mfun1_dn, val1.median2=mfun1_up, val1.diff, val1.CLES, val1.p.value,
      val2.median1=mfun2_dn, val2.median2=mfun2_up, val2.diff, val2.CLES, val2.p.value, stringsAsFactors=FALSE)
  }
  # Convert colnames to use names of variables 
  colnames(out) <- gsub("val1", var1, colnames(out))
  colnames(out) <- gsub("val2", var2, colnames(out))
  # Convert colnames
  if(mfun == "mean") colnames(out) <- gsub("median", "mean", colnames(out))
  return(invisible(out))
}

# canprot/R/diffplot.R
# Plot mean or median differences of Zc and nH2O, or other variables
# 20160715 jmd
# 20190329 add oldstyle = FALSE (no drop lines; show kernel density)
diffplot <- function(comptab, vars=c("Zc", "nH2O"), col="black", plot.rect=FALSE, pt.text=c(letters, LETTERS),
                     cex.text = 0.85, oldstyle = FALSE, pch = 1, cex = 2.1, contour = TRUE, col.contour = par("fg"),
                     probs = 0.5, add = FALSE, labtext = NULL, ...) {

  # Convert to data frame if needed
  if(!is.data.frame(comptab)) comptab <- do.call(rbind, comptab)
  # Which columns we're using
  stats <- c("diff", "CLES", "p.value")
  iX <- unlist(sapply(paste(vars[1], stats, sep=".*"), grep, colnames(comptab)))
  iY <- unlist(sapply(paste(vars[2], stats, sep=".*"), grep, colnames(comptab)))
  # Get mean/median difference, common language effect size and p-value
  X_d <- comptab[, iX[1]]
  Y_d <- comptab[, iY[1]]

  # Figure out axis labels
  # Only use part before underscore 20191207
  Dx <- paste0("D", strsplit(vars[1], "_")[[1]][1])
  Dy <- paste0("D", strsplit(vars[2], "_")[[1]][1])
  if(oldstyle) {
    # "oldstyle" labels including overbar
    cplabbar <- cplab
    cplabbar$nH2O <- quote(bar(italic(n))[H[2]*O])
    cplabbar$DnH2O <- quote(Delta*bar(italic(n))[H[2]*O])
    xvar <- cplabbar[[Dx]]
    yvar <- cplabbar[[Dy]]
    # For oldstyle plots, also get common language effect size and p-value
    X_e <- signif(comptab[, iX[2]], 2)
    X_p <- comptab[, iX[3]]
    Y_e <- signif(comptab[, iY[2]], 2)
    Y_p <- comptab[, iY[3]]
  } else {
    xvar <- cplab[[Dx]]
    yvar <- cplab[[Dy]]
  }
  # Use colnames to figure out whether the difference is of the mean or median
  if(is.null(labtext)) {
    # Treat the x- and y-variables separately in case one is median and one is mean 20191127
    xfun <- gsub("1", "", strsplit(grep(vars[1], colnames(comptab), value = TRUE)[1], "\\.")[[1]][2])
    yfun <- gsub("1", "", strsplit(grep(vars[2], colnames(comptab), value = TRUE)[1], "\\.")[[1]][2])
    xparen <- paste0("(", xfun, " difference)")
    yparen <- paste0("(", yfun, " difference)")
  } else {
    labtext <- rep(labtext, length.out = 2)
    if(is.list(labtext)) lt1 <- labtext[[1]] else lt1 <- labtext[1]
    if(is.list(labtext)) lt2 <- labtext[[2]] else lt2 <- labtext[2]
    xparen <- bquote("("*.(lt1)*")")
    yparen <- bquote("("*.(lt2)*")")
  }
  if(identical(labtext[1], NA)) xlab <- xvar else xlab <- bquote(.(xvar) ~ .(xparen))
  if(identical(labtext[2], NA)) ylab <- yvar else ylab <- bquote(.(yvar) ~ .(yparen))

  # Initialize plot: add a 0 to make sure we can see the axis
  # Prevent NA values from influencing the scale of the plot 20200103
  ina <- is.na(X_d) | is.na(Y_d)
  if(!add) plot(type="n", c(X_d[!ina], 0), c(Y_d[!ina], 0), xlab=xlab, ylab=ylab, ...)
  # Add a reference rectangle
  if(plot.rect) rect(-0.01, -0.01, 0.01, 0.01, col="grey80", lwd=0)
  # show axis lines
  abline(h=0, lty=3, col = "gray30")
  abline(v=0, lty=3, col = "gray30")
  if(oldstyle) {
    # Show drop lines: dotted/solid if p-value/effect size meet criteria
    lty.X <- ifelse(abs(X_e - 50) >= 10, 1, ifelse(X_p < 0.05, 2, 0))
    lty.Y <- ifelse(abs(Y_e - 50) >= 10, 1, ifelse(Y_p < 0.05, 2, 0))
    for(i in seq_along(X_d)) {
      lines(rep(X_d[i], 2), c(0, Y_d[i]), lty=lty.X[i])
      lines(c(0, X_d[i]), rep(Y_d[i], 2), lty=lty.Y[i])
    } 
    # Point symbols: open circle, filled circle, filled square (0, 1 or 2 vars with p-value < 0.05)
    p_signif <- rowSums(data.frame(X_p < 0.05, Y_p < 0.05))
    pch <- ifelse(p_signif==2, 15, ifelse(p_signif==1, 19, 21))
  }

  # Plot points with specified color (and point style, only for oldstyle = FALSE)
  col <- rep(col, length.out=nrow(comptab))
  pch <- rep(pch, length.out=nrow(comptab))
  points(X_d, Y_d, pch=pch, col=col, bg="white", cex=cex)
  if(!identical(pt.text, NA) | !identical(pt.text, FALSE)) {
    # Add letters; for dark background, use white letters
    col[pch %in% c(15, 19)] <- "white"
    text(X_d, Y_d, pt.text[seq_along(X_d)], col=col, cex=cex.text)
  }

  # Contour 2-D kernel density estimate 20190329
  # https://stats.stackexchange.com/questions/31726/scatterplot-with-contour-heat-overlay
  if(!oldstyle & any(contour)) {
    # Include points specified by 'contour' 20191102
    Xcont <- X_d[contour]
    Ycont <- Y_d[contour]
    # Remove NA points 20191127
    iNA <- is.na(Xcont) | is.na(Ycont)
    Xcont <- Xcont[!iNA]
    Ycont <- Ycont[!iNA]
    if(length(Xcont) > 0) {
      dens <- kde2d(Xcont, Ycont, n = 200)
      # Add contour around 50% of points (or other fractions specified by 'probs') 20191126
      # https://stackoverflow.com/questions/16225530/contours-of-percentiles-on-level-plot
      # (snippet from emdbook::HPDregionplot from @benbolker)
      dx <- diff(dens$x[1:2])
      dy <- diff(dens$y[1:2])
      sz <- sort(dens$z)
      c1 <- cumsum(sz) * dx * dy
      levels <- sapply(probs, function(x) {
        approx(c1, sz, xout = 1 - x)$y
      })
      # Use lty = 2 and lwd = 1 if points are being plotted, or lty = 1 and lwd = 2 otherwise 20191126
      if(identical(pch[1], NA)) lty <- 1 else lty <- 2
      if(identical(pch[1], NA)) lwd <- 2 else lwd <- 1
      # Don't try to plot contours for NA levels 20191207
      if(!any(is.na(levels))) contour(dens, drawlabels = FALSE, levels = levels, lty = lty, lwd = lwd, add = TRUE, col = col.contour)
    }
  }

}

# Quantile distribution for up- and down-regulated proteins
# First version 20200428
# Use ecdf() to calculate knots 20200506
qdist <- function(pdat, vars = c("Zc", "nH2O"), show.steps = FALSE) {
  # Initialize plot
  if(length(vars)==2) par(mfrow = c(2, 1))
  par(yaxs = "i")
  for(var in vars) {
    if(var=="Zc") {
      X <- Zc(pdat$pcomp$aa)
      xlab <- cplab$Zc
    }
    if(var=="nH2O") {
      X <- nH2O(pdat$pcomp$aa)
      xlab <- cplab$nH2O
    }
    up <- X[pdat$up2]
    dn <- X[!pdat$up2]
    # Start the plot
    plot(range(up, dn), c(0, 1), type = "n", xlab = xlab, ylab = "Quantile point")
    # A function to plot the points and lines
    pfun <- function(x, ...) {
      Fn <- ecdf(x)
      # Plot the values (knots) and verticals
      if(show.steps) plot(Fn, add = TRUE, cex = 0.25, col.01line = NA, verticals = TRUE, col.hor = "gray70", ...)
      # Add lines that bisect the verticals (to intersect the quantile points)
      x <- knots(Fn)
      y <- sort(Fn(x)) - 0.5 / length(x)
      lines(x, y, ...)
    }
    pfun(dn)
    pfun(up, col = 2, lty = 2)
    # Draw a line at the 0.5 quantile
    lines(c(median(dn), median(up)), c(0.5, 0.5), lwd = 2, col = "slategray3")
  }
}

# Format the summary table using xtable
# 20160709 jmd
xsummary <- function(comptab, vars=c("Zc", "nH2O")) {
  # Convert to data frame if needed
  if(!is.data.frame(comptab)) comptab <- do.call(rbind, comptab)
  # Create letter labels
  rownames(comptab) <- c(letters, LETTERS)[1:nrow(comptab)]
  # Return this data frame when the function exits
  comptab.out <- comptab
  # Get the publication key from the dataset name
  publication <- sapply(strsplit(comptab$dataset, "_"), "[", 1)
  # Format the publication key (monospaced font)
  publication <- paste0("<code>", publication, "</code>")
  # Combine the publication and description
  comptab$description <- paste0(publication, " (", comptab$description, ")")
  # Datasets have letter labels
  comptab$dataset <- c(letters, LETTERS)[1:nrow(comptab)]
  # To save column width, change "dataset" to "set"
  colnames(comptab)[1] <- "set"
  # Select the columns to print
  comptab <- comptab[, c(1:4, 7:9, 12:14)]
  # Round values in some columns
  comptab[, c(5, 8)] <- round(comptab[, c(5, 8)], 3)    # mean difference
  comptab[, c(6, 9)] <- signif(comptab[, c(6, 9)], 2)   # CLES
  # Place a marker around high effect size and low p-value
  for(k in vars) {
    # Effect size
    jes <- paste0(k, ".CLES")
    ihigh <- abs(comptab[, jes] - 50) >= 10
    comptab[ihigh, jes] <- paste("**", comptab[ihigh, jes], "**")
    # p-value
    jpv <- paste0(k, ".p.value")
    ilow <- comptab[, jpv] < 0.05
    comptab[, jpv] <- format(comptab[, jpv], digits=1)
    comptab[ilow, jpv] <- paste("**", comptab[ilow, jpv], "**")
    # Bold mean difference if both effect size and p-value are highlighted
    jmd <- paste0(k, ".diff")
    comptab[, jmd] <- format(comptab[, jmd])
    comptab[ihigh & ilow, jmd] <- paste("**", comptab[ihigh & ilow, jmd], "**") 
    # Underline mean difference if only p-value or only effect size is highlighted
    comptab[xor(ilow, ihigh), jmd] <- paste("++", comptab[xor(ilow, ihigh), jmd], "++") 
  }
  # Create xtable
  x <- xtable::xtable(comptab, align=c("c", "l", "l", rep("r", 8)))
  x <- capture.output(xtable::print.xtable(x, type="html",
                                   include.rownames=FALSE,
                                   math.style.exponents=TRUE,
                                   sanitize.text.function=function(x){x}))
  # Bold the indicated effect size, p-values and mean differences
  x <- gsub("> ** ", "> <B>", x, fixed=TRUE)
  x <- gsub(" ** <", "</B> <", x, fixed=TRUE)
  # Underline the indicated mean differences
  x <- gsub("> ++ ", "> <U>", x, fixed=TRUE)
  x <- gsub(" ++ <", "</U> <", x, fixed=TRUE)
  # Add headers that span multiple columns
  span_empty <- "<th colspan=\"4\"></th>"
  span_Zc <- "<th colspan=\"3\"><i>Z</i><sub>C</sub></th>"
  span_nH2O <- "<th colspan=\"3\"><i>n</i><sub>H<sub>2</sub>O</sub></sub></th>"
  span_nO2 <- "<th colspan=\"3\"><i>n</i><sub>O<sub>2</sub></sub></sub></th>"
  span_nAA <- "<th colspan=\"3\"><i>n</i><sub>AA</sub></th>"
  span_var1 <- get(paste0("span_", vars[1]))
  span_var2 <- get(paste0("span_", vars[2]))
  x <- gsub("<table border=1>",
            paste("<table border=1> <tr>", span_empty, span_var1, span_var2, "</tr>"),
            x, fixed=TRUE)
  # More formatting of the headers
  x <- gsub("description", "reference (description)", x, fixed=TRUE)
  x <- gsub("n1", "<i>n</i><sub>1</sub>", x, fixed=TRUE)
  x <- gsub("n2", "<i>n</i><sub>2</sub>", x, fixed=TRUE)
  x <- gsub(paste0(vars[1], ".diff"), "MD", x, fixed=TRUE)
  x <- gsub(paste0(vars[2], ".diff"), "MD", x, fixed=TRUE)
  x <- gsub(paste0(vars[1], ".CLES"), "ES", x, fixed=TRUE)
  x <- gsub(paste0(vars[2], ".CLES"), "ES", x, fixed=TRUE)
  x <- gsub(paste0(vars[1], ".p.value"), "<i>p</i>-value", x, fixed=TRUE)
  x <- gsub(paste0(vars[2], ".p.value"), "<i>p</i>-value", x, fixed=TRUE)
  # Take out extraneous spaces (triggers pre-formatted text - in markdown?)
  x <- gsub("   ", " ", x, fixed=TRUE)
  # Done!
  write.table(x, "", row.names=FALSE, col.names=FALSE, quote=FALSE)
  return(invisible(comptab.out))
}

# Make table for new vignettes
# 20191206
xsummary2 <- function(comptab1, comptab2) {
  ct1 <- do.call(rbind, comptab1)
  ct2 <- do.call(rbind, comptab2)
  # Get all data
  # Include medians or means for up and down groups for making summary .csv files 20200125
  out <- data.frame(
    dataset = ct1$dataset,
    description = ct2$description,
    n1 = ct1$n1,
    n2 = ct1$n2,

    Zc.down = ct1$Zc.median1,
    Zc.up = ct1$Zc.median2,
    Zc.diff = ct1$Zc.diff,

    nH2O.down = ct1$nH2O.median1,
    nH2O.up = ct1$nH2O.median2,
    nH2O.diff = ct1$nH2O.diff,

    nAA.down = ct2$nAA.median1,
    nAA.up = ct2$nAA.median2,
    nAA.diff = ct2$nAA.diff,

    MW.down = ct2$MW.median1,
    MW.up = ct2$MW.median2,
    MW.diff = ct2$MW.diff,

    stringsAsFactors = FALSE
  )

  # Prepare table
  x <- out[, c("dataset", "description", "n1", "n2", "Zc.diff", "nH2O.diff",
               "nAA.diff", "MW.diff")]
  # Get the publication key from the dataset name
  publication <- sapply(strsplit(x$dataset, "_"), "[", 1)
  # Format the publication key (monospaced font)
  publication <- paste0("<code>", publication, "</code>")
  # Italicize species names (first words in description, surrounded by underscore)
  x$description <- gsub("^_", "<i>", x$description)
  x$description <- gsub("_", "</i>", x$description)
  # Combine the publication and description
  x$description <- paste0(publication, " (", x$description, ")")
  # Datasets have letter labels
  x$dataset <- c(letters, LETTERS)[1:nrow(x)]
  # To save column width, change "dataset" to "set"
  colnames(x)[1] <- "set"
  # Multiply values of Zc and nH2O by 1000
  x[, 5:6] <- x[, 5:6] * 1000
  # Multiply values of MW by 100
  x[, 8] <- x[, 8] * 100
  # Round values
  x[, 5:8] <- round(x[, 5:8])

  # Put markers around negative values
  for(icol in 5:8) {
    ineg <- x[, icol] < 0
    ineg[is.na(ineg)] <- FALSE
    x[ineg, icol] <- paste("**", x[ineg, icol], "**")
  }

  # Create xtable
  x <- xtable::xtable(x, align=c("c", "l", "l", rep("r", ncol(x) - 2)))
  x <- capture.output(xtable::print.xtable(x, type="html",
                                   include.rownames=FALSE,
                                   math.style.exponents=TRUE,
                                   sanitize.text.function=function(x){x}))

  # Make the marked negative values bold
  x <- gsub("> ** ", "> <B>", x, fixed=TRUE)
  x <- gsub(" ** <", "</B> <", x, fixed=TRUE)
  # Change "NaN" to "NA"
  x <- gsub("NaN", "NA", x, fixed=TRUE)

  # Add headers that span multiple columns
  span_empty6 <- "<td align=\"center\" colspan=\"6\"></td>"
  span_empty2 <- "<td align=\"center\" colspan=\"2\"></td>"
  border <- paste("<table border=1> <tr>", span_empty6, span_empty2, "</tr>")
  x <- gsub("<table border=1>", border, x, fixed=TRUE)

  # More formatting of the headers
  x <- gsub("description", "reference (description)", x, fixed=TRUE)
  x <- gsub("n1", "<i>n</i><sub>down</sub>", x, fixed=TRUE)
  x <- gsub("n2", "<i>n</i><sub>up</sub>", x, fixed=TRUE)
  x <- gsub("Zc.diff", "&Delta;<i>Z</i><sub>C</sub>", x, fixed=TRUE)
  x <- gsub("nH2O.diff", "&Delta;<i>n</i><sub>H<sub>2</sub>O</sub>", x, fixed=TRUE)
  x <- gsub("nAA.diff", "&Delta;<i>n</i><sub>AA</sub>", x, fixed=TRUE)
  x <- gsub("MW.diff", "&Delta;MW", x, fixed=TRUE)

  # Done!
  write.table(x, "", row.names=FALSE, col.names=FALSE, quote=FALSE)
  return(invisible(out))
}

# Make table for new vignettes (with GRAVY and pI)
# 20200418
xsummary3 <- function(comptab1, comptab2, comptab3) {
  ct1 <- do.call(rbind, comptab1)
  ct2 <- do.call(rbind, comptab2)
  ct3 <- do.call(rbind, comptab3)
  # Get all data
  # Include medians or means for up and down groups for making summary .csv files 20200125
  out <- data.frame(
    dataset = ct1$dataset,
    description = ct2$description,
    n1 = ct1$n1,
    n2 = ct1$n2,

    Zc.down = ct1$Zc.median1,
    Zc.up = ct1$Zc.median2,
    Zc.diff = ct1$Zc.diff,

    nH2O.down = ct1$nH2O.median1,
    nH2O.up = ct1$nH2O.median2,
    nH2O.diff = ct1$nH2O.diff,

    pI.down = ct2$pI.median1,
    pI.up = ct2$pI.median2,
    pI.diff = ct2$pI.diff,

    GRAVY.down = ct2$GRAVY.median1,
    GRAVY.up = ct2$GRAVY.median2,
    GRAVY.diff = ct2$GRAVY.diff,

    nAA.down = ct3$nAA.median1,
    nAA.up = ct3$nAA.median2,
    nAA.diff = ct3$nAA.diff,

    MW.down = ct3$MW.median1,
    MW.up = ct3$MW.median2,
    MW.diff = ct3$MW.diff,

    stringsAsFactors = FALSE
  )

  # Prepare table
  x <- out[, c("dataset", "description", "n1", "n2", "Zc.diff", "nH2O.diff",
               "pI.diff", "GRAVY.diff", "nAA.diff", "MW.diff")]
  # Get the publication key from the dataset name
  publication <- sapply(strsplit(x$dataset, "_"), "[", 1)
  # Format the publication key (monospaced font)
  publication <- paste0("<code>", publication, "</code>")
  # Italicize species names (first words in description, surrounded by underscore)
  x$description <- gsub("^_", "<i>", x$description)
  x$description <- gsub("_", "</i>", x$description)
  # Combine the publication and description
  x$description <- paste0(publication, " (", x$description, ")")
  # Datasets have letter labels
  x$dataset <- c(letters, LETTERS)[1:nrow(x)]
  # To save column width, change "dataset" to "set"
  colnames(x)[1] <- "set"
  # Multiply values of Zc, nH2O and GRAVY by 1000
  x[, c(5:6, 8)] <- x[, c(5:6, 8)] * 1000
  # Multiply values of pI and MW by 100
  x[, c(7, 10)] <- x[, c(7, 10)] * 100
  # Round values
  x[, 5:10] <- round(x[, 5:10])

  # Put markers around negative values
  for(icol in 5:10) {
    ineg <- x[, icol] < 0
    ineg[is.na(ineg)] <- FALSE
    x[ineg, icol] <- paste("**", x[ineg, icol], "**")
  }

  # Create xtable
  x <- xtable::xtable(x, align=c("c", "l", "l", rep("r", ncol(x) - 2)))
  x <- capture.output(xtable::print.xtable(x, type="html",
                                   include.rownames=FALSE,
                                   math.style.exponents=TRUE,
                                   sanitize.text.function=function(x){x}))

  # Make the marked negative values bold
  x <- gsub("> ** ", "> <B>", x, fixed=TRUE)
  x <- gsub(" ** <", "</B> <", x, fixed=TRUE)
  # Change "NaN" to "NA"
  x <- gsub("NaN", "NA", x, fixed=TRUE)

  # More formatting of the headers
  x <- gsub("description", "reference (description)", x, fixed=TRUE)
  x <- gsub("n1", "<i>n</i><sub>down</sub>", x, fixed=TRUE)
  x <- gsub("n2", "<i>n</i><sub>up</sub>", x, fixed=TRUE)
  x <- gsub("Zc.diff", "&Delta;<i>Z</i><sub>C</sub>", x, fixed=TRUE)
  x <- gsub("nH2O.diff", "&Delta;<i>n</i><sub>H<sub>2</sub>O</sub>", x, fixed=TRUE)
  x <- gsub("pI.diff", "&Delta;pI", x, fixed=TRUE)
  x <- gsub("GRAVY.diff", "&Delta;GRAVY", x, fixed=TRUE)
  x <- gsub("nAA.diff", "&Delta;<i>n</i><sub>AA</sub>", x, fixed=TRUE)
  x <- gsub("MW.diff", "&Delta;MW", x, fixed=TRUE)

  # Done!
  write.table(x, "", row.names=FALSE, col.names=FALSE, quote=FALSE)
  return(invisible(out))
}
