# canprot/R/Ehplot.R
# show values of Eh as a function of logaH2O/logfO2
# 20160710 jmd
# 20200829 add 'xlim' and 'ylim' arguments
Ehplot <- function(T=37, pH=7.4, xlim = c(-75, -55), ylim = c(-10, 10), dy = 0.2) {
  logfO2 <- xlim
  # logK for the reaction H2O(liq) = 2H+ + 2e- + 0.5O2(g)
  logK <- subcrt(c("H2O", "H+", "e-", "oxygen"), c(-1, 2, 2, 0.5), T=T)$out$logK
  # to calculate logaH2O at a given logfO2 and Eh
  logaH2O <- function(logfO2, Eh) {
    pe <- convert(Eh, "pe", T=convert(T, "K"))
    return(0.5*logfO2 - 2*pH - 2*pe - logK)
  }
  plot(0, 0, xlim=xlim, ylim=ylim, xlab=quote(log~italic("f")[O[2]]), ylab=quote(log~italic("a")[H[2]*O]), type="n", xaxs="i", yaxs="i")
  for(Eh in seq(-0.8, 0.2, by=dy)) {
    lines(logfO2, logaH2O(logfO2, Eh)) 
    text(-61+19*Eh, logaH2O(-61+19*Eh, Eh) + 1, Eh)
  }
  # add legend 20200829
  legend("topright", describe.property(c("T", "pH"), c(T, pH), digits = 1), bg = "white")
  title(main="Eh (volt)", cex.main=1.1)
}

# canprot/R/rankdiff.R
# calculate weighted difference of sums of ranks
# 20160710 jmd
rankdiff <- function(rank1, rank2, n1=NULL, n2=NULL, as.fraction=TRUE) {
  if(!is.null(n1)) sum1 <- rank1
  else {
    sum1 <- sum(rank1)
    n1 <- length(rank1)
  }
  if(!is.null(n2)) sum2 <- rank2
  else {
    sum2 <- sum(rank2)
    n2 <- length(rank2)
  }
  # scaling to account for different numbers of proteins
  # up-expressed in normal (n1) and cancer (n2) tissue
  sum1 <- 2 * sum1 * n2 / (n1 + n2)
  sum2 <- 2 * sum2 * n1 / (n1 + n2)
  sumdiff <- sum2 - sum1
  if(!as.fraction) return(sumdiff)
  else{
    # what is the maximum possible rank difference? e.g. 20 for H.H.H.H.H.C.C.C.C
    # use numeric values to avoid maxC * n1 giving NAs
    # produced by integer overflow for large data sets (e.g. WDO+15)
    min1 <- as.numeric(sum(1:n1))
    maxC <- as.numeric(sum((1:n2) + n1))
    maxdiff <- 2 * (maxC * n1 - min1 * n2) / (n1 + n2)
    # express the rank difference as a fraction of maximum
    rankdiff <- sumdiff / maxdiff
    return(rankdiff)
  }
}

# canprot/R/rankplot.R
# make logaH2O - logfO2 diagrams of ranking of chemical affinity
# 20160710 jmd
rankplot <- function(pdat, T=37, what="rankdiff", main=NULL, res=300, plot.it=TRUE, xlim=c(-75, -55), ylim=c(-10, 10), basis = "QEC+") {
  # get protein data
  aa <- pdat$pcomp$aa
  ip <- add.protein(aa)
  # assemble limits of logfO2, logaH2O
  # use res+1 here to ensure that x- and y-directions aren't accidentally transposed
  # (would cause an error in image())
  H2O <- c(ylim, res)
  O2 <- c(xlim, res+1)
  ys <- seq(H2O[1], H2O[2], length.out=H2O[3])
  xs <- seq(O2[1], O2[2], length.out=O2[3])
  # what is the log activity of protein corresponding to unit activity of residues?
  plength <- protein.length(pdat$pcomp$aa)
  loga.protein <- log10(1/plength)
  # calculate affinities
  basis(basis)
  a <- affinity(O2=O2, H2O=H2O, T=T, iprotein=ip, loga.protein=loga.protein)
  if(identical(what, "affinity")) {
    col <- ifelse(pdat$up2, 2, 4)
    d <- diagram(a, names=pdat$names, fill=col, as.residue=TRUE, tplot=FALSE, format.names=FALSE)
    # redraw box because it gets obscured by the fill
    box()
  } else if(identical(what, "rankdiff")) {
    # plot affinity difference (healthy - cancer)
    Aarr <- list2array(a$values)
    # grid points on the plot
    grid <- expand.grid(seq_along(xs), seq_along(ys))
    rank_C <- palply("", 1:nrow(grid), function(i) {
      # the ranks at this grid point
      r <- rank(Aarr[grid[i, 1], grid[i, 2], ]/plength)
      # the sum of ranks for cancer proteins at this grid point
      sum_C <- sum(r[pdat$up2])
      return(sum_C)
    })
    # numbers of cancer, healthy, and total proteins
    n_C <- sum(pdat$up2)
    n_H <- sum(!pdat$up2)
    # cancer and healthy rank sums in matrix format
    rank_C <- matrix(unlist(rank_C), nrow=length(xs), ncol=length(ys))
    rank_H <- sum(1:(n_C + n_H)) - rank_C
    # calculate weighted rank difference in percent
    rankdiff <- 100*rankdiff(rank_H, rank_C, n_H, n_C)
    print(paste("weighted rank difference % range", paste(range(rankdiff), collapse=" ")))
    if(!plot.it) return(list(xs=xs, ys=ys, rankdiff=rankdiff, xlab=quote(log~italic("f")[O[2]]), ylab=quote(log~italic("a")[H[2]*O])))
    # display greater healthy and cancer ranks by colored zones on image
    col <- get_colors(rankdiff, max50=TRUE)
    image(xs, ys, rankdiff, col=col, useRaster=TRUE, xlab=quote(log~italic("f")[O[2]]), ylab=quote(log~italic("a")[H[2]*O]))
    # show equal-rank line
    contour(xs, ys, rankdiff, levels=0, add=TRUE, drawlabels=FALSE)
    contour(xs, ys, rankdiff, levels=seq(-100, 100, by=10), lty=3, add=TRUE, drawlabels=FALSE)
    box()
  }
  if(is.null(main)) main <- pdat$description
  title(main=main, cex.main=1.1)
  mtext(pdat$dataset, side=4, cex=0.85, las=0, adj=0, line=-0.1)
}

# canprot/potential.R
# plotting potential diagrams
# 20170611

# calculate rank potential diagrams for each dataset in a group
groupplots <- function(group="hypoxia_ZC_down", each100=FALSE, res=50, plot.it=TRUE) {
  idat <- get_idat(group)
  # make figures and return coordinates of equipotential lines
  if(plot.it) {
    # number of plots - 1 extra for title
    np <- sum(idat$idat) + 1
    if(np <= 28) mfrow <- c(4, 7)
    if(np <= 24) mfrow <- c(4, 6)
    if(np <= 20) mfrow <- c(4, 5)
    if(np <= 18) mfrow <- c(3, 6)
    if(np <= 16) mfrow <- c(4, 4)
    if(np <= 15) mfrow <- c(3, 5)
    if(np <= 12) mfrow <- c(3, 4)
    if(np <= 9) mfrow <- c(3, 3)
    if(np <= 8) mfrow <- c(2, 4)
    if(np <= 6) mfrow <- c(2, 3)
    if(np <= 4) mfrow <- c(2, 2)
    par(mfrow=mfrow)
    par(mar=c(3.5, 4, 2, 1.5))
    par(mgp=c(2.5, 1, 0))
    # first plot the title
    plot.new()
    opar <- par(xpd=TRUE)
    text(0.5, 0.8, idat$group, cex=2)
    text(0.5, 0.6, get_main(group), cex=2)
    par(opar)
  }
  # m=1 here to calculate lines for each dataset individually
  # expand x- and y-axes (48 log units)
  cpeach <- calcpot(what=idat$what, datasets=idat$datasets, res=res, plot.it=plot.it, each100=each100, xlim=c(-88, -40), ylim=c(-24, 24), m=1)
  return(list(group=group, cpeach=cpeach))
}

# plot merged ranked potential diagrams for groups of datasets
mergedplot <- function(gpresult, each100=FALSE, res=50) {
  idat <- get_idat(gpresult$group)
  cpresult <- calcpot(what=idat$what, datasets=idat$datasets, res=res, each100=each100)
  # the merged potential diagram
  plotit(cpresult[[1]][[1]], cpresult[[2]][[1]])
  # median and interquartile range of the individual equipotential lines
  plot_median(gpresult$cpeach, add=TRUE, transpose=grepl("ZC", gpresult$group))
  # plot the merged equipotential line again (for emphasis)
  plotit(cpresult[[1]][[1]], cpresult[[2]][[1]], add=TRUE)
}

### Internal Functions ###

get_idat <- function(group) {
  what <- strsplit(group, "_")[[1]][1]
  metric <- strsplit(group, "_")[[1]][2]
  direction <- strsplit(group, "_")[[1]][3]
  extdatadir <- system.file("extdata", package="JMDplots")
  file <- paste0(extdatadir, "/cpcp/summary/summary_", what, ".csv")
  dat <- read.csv(file, as.is=TRUE)
  cnames <- colnames(dat)
  # the difference for this metric
  diff <- dat[, grepl(metric, cnames) & grepl("diff", cnames)]
  # the p-value and CLES for the other metric
  p.value <- dat[, !grepl(metric, cnames) & grepl("p.value", cnames)]
  CLES <- dat[, !grepl(metric, cnames) & grepl("CLES", cnames)]
  # the rows that have the difference in the specified direction
  # and the magnitude of the difference is at least 0.01
  # and for which the other metric doesn't have a large or significant change
  if(direction=="up") idat <- diff > 0.01 & !(p.value < 0.05 | abs(signif(CLES, 2) - 50) >= 10)
  if(direction=="down") idat <- diff < -0.01 & !(p.value < 0.05 | abs(signif(CLES, 2) - 50) >= 10)
  return(list(what=what, idat=idat, datasets=dat$dataset[idat]))
}

# function to plot values with colors and labels
plotit <- function(values, rankdat, add=FALSE) {
  col <- get_colors(values)
  with(rankdat, {
    if(!add) {
      image(xs, ys, values, col=col, useRaster=TRUE, xlab=xlab, ylab=ylab, axes=FALSE)
      if(diff(range(xs)) > 9) {
        axis(1, at=seq(-88, -40, by=12))
        axis(2, at=seq(-24, 24, by=12), las=1)
      } else {
        axis(1, at=seq(-70, -62, by=2))
        axis(2, at=seq(-6, 2, by=2), las=1)
      }
      box()
    }
    # show equal-rank line
    contour(xs, ys, values, levels=0, add=TRUE, drawlabels=FALSE, lwd=1.5, col="white")
  })
}

get_main <- function(group, n=NA) {
  what <- strsplit(group, "_")[[1]][1]
  metric <- strsplit(group, "_")[[1]][2]
  direction <- strsplit(group, "_")[[1]][3]
  if(direction=="up") than <- "> 0.01" else than <- "< -0.01"
  if(metric=="ZC") main <- substitute(Delta*italic(Z)[C]~than, list(than=than))
  if(metric=="H2O") main <- substitute(Delta*italic(bar(n))[H[2]*O]~than, list(than=than))
  if(!is.na(n)) {
    # make long titles to be placed on two lines
    if(what %in% c("hypoxia", "osmotic")) {
      if(what=="hypoxia") what <- "hypoxia or 3D culture"
      if(what=="osmotic") what <- "hyperosmotic stress"
      main <- c(what, substitute(main~~"("*x*")", list(what=what, main=main, x=n)))
    } else main <- substitute(what~~main~~"("*x*")", list(what=what, main=main, x=n))
  }
  return(main)
}

# plot median and 1st and 3rd quartiles of equipotential lines 20161021
plot_median <- function(eachplotres, each100=FALSE, transpose=FALSE, add=FALSE, type=7) {
  # show scatterplot of points on all lines
  if(!add) plot(eachplotres$x, eachplotres$y, pch=".")
  if(transpose) {
    ys <- seq(-24, 24, by=0.2)
    xs <- sapply(unique(eachplotres$comb), function(comb) {
      idat <- eachplotres$comb == comb
      # use loess to get smooth lines
      indlo <- loess(x~y, data=eachplotres[idat, ], span=0.2)
      spx <- predict(indlo, newdata=ys)
      # show loess lines
      if(!add) lines(spx, ys, col="gray")
      return(spx)
    })
    # calculate and plot median
    xmedian <- apply(xs, 1, median, na.rm=TRUE)
    lines(xmedian, ys, lwd=2)
    # calculate and plot 1st and 3rd quartiles
    xq1 <- apply(xs, 1, quantile, probs=0.25, na.rm=TRUE, type=type)
    xq3 <- apply(xs, 1, quantile, probs=0.75, na.rm=TRUE, type=type)
    lines(xq1, ys, lty=2)
    lines(xq3, ys, lty=2)
  } else {
    xs <- seq(-88, -40, by=0.2)
    ys <- sapply(unique(eachplotres$comb), function(comb) {
      idat <- eachplotres$comb == comb
      # use loess to get smooth lines
      indlo <- loess(y~x, data=eachplotres[idat, ], span=0.2)
      spy <- predict(indlo, newdata=xs)
      # show loess lines
      if(!add) lines(xs, spy, col="gray")
      return(spy)
    })
    # calculate and plot median
    ymedian <- apply(ys, 1, median, na.rm=TRUE)
    lines(xs, ymedian, lwd=1.5)
    # calculate and plot 1st and 3rd quartiles
    yq1 <- apply(ys, 1, quantile, probs=0.25, na.rm=TRUE, type=type)
    yq3 <- apply(ys, 1, quantile, probs=0.75, na.rm=TRUE, type=type)
    lines(xs, yq1, lty=2)
    lines(xs, yq3, lty=2)
  }
}

# merging potential diagrams for pancreatic cancer 20161012
calcpot <- function(what="pancreatic", datasets=c("LHE+04", "PCS+11_PDAC"),
                            res=50, plot.it=FALSE,
                            each100=FALSE, xlim=c(-70.2, -61.8), ylim=c(-6.2, 2.2), m=0) {
  # get the data
  pdat_fun <- paste0("pdat_", what)
  # Special for .pdat_osmotic (2017 osmotic stress datasets) 20201015
  if(what == "osmotic") pdat_fun <- paste0(".pdat_", what)
  rankdat <- lapply(datasets, function(dataset) {
    pdat <- get(pdat_fun)(dataset)
    rankplot(pdat, res=res, plot.it=FALSE, xlim=xlim, ylim=ylim)
  })
  rankdiff <- lapply(rankdat, "[[", "rankdiff")
  extdatadir <- system.file("extdata", package="JMDplots")
  file <- paste0(extdatadir, "/cpcp/summary/summary_", what, ".csv")
  dat <- read.csv(file, as.is=TRUE)
  # function to scale all values by same factor such that maximum is 100
  scale100 <- function(valuelist, each100=FALSE) {
    if(each100) lapply(valuelist, function(x) x * 100 / max(abs(range(x))))
    else {
      totalrange <- range(valuelist)
      lapply(valuelist, function(x) x * 100 / max(abs(totalrange)))
    }
  }
  # scale and display the rank differences
  rankdiff <- scale100(rankdiff, each100)
  # calculate the percentages of total values
  rabs <- lapply(rankdiff, abs)
  rsum <- sapply(rabs, sum)
  rperc <- round(100 * rsum / sum(rsum), 1)
  if(plot.it) for(i in seq_along(rankdat)) {
    plotit(rankdiff[[i]], rankdat[[i]])
    mtext(datasets[i], side=4, cex=0.85, las=0, adj=0, line=0.1)
    # show the dataset description and letter
    idat <- match(datasets[i], dat$dataset)
    label.figure(dat$set[idat], xfrac=0.25, cex=1.5)
    label.figure(paste0("     (", rperc[i], "%)"), adj=0, xfrac=0.25, cex=1.5)
    # put a circle around the letter
    opar <- par(xpd=NA)
    points(grconvertX(0.25, "nfc"), grconvertY(0.95, "nfc"), cex=3.5)
    par(opar)
  }
  if(m) {
    # merge the rank differences for all combinations of 'm' number of datasets
    combs <- combn(length(rankdat), m)
    comb <- x <- y <- numeric()
    print(paste("calcpot: calculating equipotential lines for", ncol(combs), "combinations"))
    for(i in 1:ncol(combs)) {
      merged <- Reduce("+", rankdiff[combs[, i]])
      cl <- contourLines(rankdat[[1]]$xs, rankdat[[1]]$ys, merged, levels=0)
      if(length(cl) > 0) {
        # in case the line has multiple segments
        clx <- do.call(c, lapply(cl, "[[", 2))
        cly <- do.call(c, lapply(cl, "[[", 3))
        x <- c(x, clx)
        y <- c(y, cly)
        comb <- c(comb, rep(i, length(clx)))
      }
    }
    out <- data.frame(comb=comb, x=x, y=y)
    return(out)
  } else {
    # merge the rank differences for all datasets
    merged <- Reduce("+", rankdiff) / length(rankdat)
    merged <- scale100(list(merged))
    return(invisible(list(merged=merged, rankdat=rankdat)))
  }
}

# canprot/R/get_colors.R
# get colors for rank-difference (potential) diagrams
# 20160710 jmd
get_colors <- function(x, max50=FALSE) {
  # diverging (blue - light grey - red) palette
  # max50: values over 50% are all deepest color (red or blue)
  # read precomputed colors:
  # colorspace::diverge_hcl(1000, c = 100, l = c(50, 90), power = 1)
  # colorspace::diverge_hcl(2000, c = 100, l = c(50, 90), power = 1)
  if(max50) file <- system.file("extdata/misc/bluered1000.txt", package = "JMDplots")
  else file <- system.file("extdata/misc/bluered2000.txt", package = "JMDplots")
  dcol <- read.table(file, as.is=TRUE)[[1]]
  # the range of values
  xrange <- range(x)
  # select range of colors corresponding to values
  if(any(xrange < 0)) {
    # white to deep blue
    if(max50) blues <- rev(c(rep(dcol[1], 500), dcol[1:500]))
    else blues <- rev(dcol[1:1000])
    blues <- blues[1:-round(10*xrange[1])]
  } else blues <- character()
  if(any(xrange > 0)) {
    # white to deep red
    if(max50) reds <- c(dcol[501:1000], rep(dcol[1000], 500))
    else reds <- dcol[1001:2000]
    reds <- reds[1:round(10*xrange[2])]
  } else reds <- character()
  col <- c(rev(blues), reds)
  return(col)
}
