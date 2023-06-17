# JMDplots/grad_util.R
# utility functions (not exported) for gradox and gradH2O papers

# get metadata (location names and sequencing IDs) for a study
# extracted from mprep/mplot 20180312
get.mdata <- function(studies, study, seqtype, remove.NA = TRUE) {
  samples <- studies[[study]][[1]]
  if(is.null(samples)) stop("metadata for ", study, " study not available")
  xlabels <- studies[[study]]$xlabels
  if(is.null(xlabels)) xlabels <- samples
  group <- studies[[study]][["group"]]
  seqtype.for.ID <- seqtype
  # change e.g. SRA_MGP to SRA_MG
  seqtype.for.ID <- gsub("_MG.$", "_MG", seqtype.for.ID)
  seqtype.for.ID <- gsub("_MT.$", "_MT", seqtype.for.ID)
  IDs <- studies[[study]][[seqtype.for.ID]]
  if(is.null(IDs)) stop(seqtype.for.ID, " IDs not available for ", study, " study")
  # remove NA IDs and corresponding samples, xlabels, and groups
  if(remove.NA) {
    samples[is.na(IDs)] <- NA
    samples <- na.omit(samples)
    xlabels[is.na(IDs)] <- NA
    xlabels <- na.omit(xlabels)
    if(length(group) > 1) {
      group[is.na(IDs)] <- NA
      group <- na.omit(group)
    }
    IDs <- na.omit(IDs)
  }
  abbrev <- studies[[study]][["abbrev"]]
  techtype <- studies[[study]][["techtype"]]
  dx <- studies[[study]][["dx"]]
  dy <- studies[[study]][["dy"]]
  return(list(samples=samples, xlabels=xlabels, IDs=IDs, group=group, abbrev=abbrev, techtype = techtype, dx = dx, dy = dy))
}

# function to plot ZC of metagenomic DNA 20180215
# or ZC of metagenomic proteins 20180228
# NOTE: dataset should end in "_MG" or "_MT" (plot DNA and RNA compositions)
# or "_MGP" or "_MTP" (plot protein compositions)
## optional: vioplotx package is needed for split violin plots
plotMG <- function(dataset="Guerrero_Negro_IMG_MG", plottype="bars",
  samples=formatC(10:1, width=2, flag="0"), labels=formatC(10:1, width=2, flag="0"),
  group="mat", xlab="layer", ylim=NULL, abbrev=NULL, dsDNA=TRUE, plot.RNA=TRUE,
  taxid=NULL, lwd=1, lty=2, lwd.bars=2, col=NULL, extendrange=FALSE, add.label=TRUE,
  plot_real_x=FALSE, maxdepth=NULL, H2O=FALSE, plot.it = TRUE, add.title = TRUE, yline = 2,
  basis = getOption("basis"), techtype = NULL, dx = NULL, dy = NULL, datadir = NULL,
  add = FALSE, all.labels = NULL, pch = 19, var = NULL, srt = 45, ilabel = NULL) {
  # samples: (used for suffixes on file names)
  # labels: (used for labeling x-axis ticks)
  # xlab: axis label: "layer", "depth", ...
  isprotein <- grepl("_MGP$", dataset) | grepl("_MTP$", dataset)
  # where to keep mean and high/lo (+/- SD) of ZC at each site
  pI <- pI.SD <- GRAVY <- GRAVY.SD <- GC <- Xlo <- Xmean <- Xhi <- numeric()
  # are we using user-supplied data?
  user_data <- TRUE
  if(is.null(datadir)) {
    # which paper is this dataset from?
    if(any(sapply(names(gradox), grepl, dataset))) paper <- "gradox"
    if(any(sapply(names(gradH2O), grepl, dataset))) paper <- "gradH2O"
    # gradox or gradH2O data location in JMDplots package 20190928
    datadir <- system.file(paste0("extdata/", paper), package = "JMDplots")
    user_data <- FALSE
  }
  # set up for proteins or DNA
  if(isprotein) {
    if(user_data) filestart <- paste0(datadir, "/MGP/", dataset)
    else {
      rdsdata <- paste0(paper, "_MGP")
      filestart <- dataset
    }
    if(H2O) Xfun <- nH2O else Xfun <- Zc
    if(is.null(col)) col <- "darkgreen"
  } else {
    # data directory for DNA
    if(user_data) filestart <- paste0(datadir, "/MGD/", dataset, "D")
    else {
      rdsdata <- paste0(paper, "_MGD")
      filestart <- paste0(dataset, "D")
    }
    Xfun <- ZCnuc
    if(is.null(col)) col <- "red"
  }
  # keep sample values for violin plot 20180515
  DNA <- data.frame(ZC=numeric(length(samples)*100), sample=character(length(samples)*100), stringsAsFactors=FALSE)
  meancomp <- NULL
  # read the data from the RDS file 20191022
  if(!user_data) RDS <- get(rdsdata, JMDplots)
  for(i in 1:length(samples)) {
    sample <- samples[i]
    if(!is.null(taxid)) file <- paste0(filestart, "_", sample, "-", taxid, ".csv")
    else file <- paste0(filestart, "_", sample, ".csv")
    # use NA if the file is missing 20180529
    if(user_data) fexists <- file.exists(file)
    else fexists <- file %in% names(RDS)
    if(!fexists) {
      Xmean <- c(Xmean, NA)
      Xlo <- c(Xlo, NA)
      Xhi <- c(Xhi, NA)
      message(paste("missing file", file))
    } else {
      if(user_data) mycomp <- read.csv(file, as.is=TRUE)
      else mycomp <- RDS[[file]]
      if(isprotein) {
        # calculate GRAVY 20191024
        myGRAVY <- GRAVY(mycomp)
        GRAVY <- c(GRAVY, mean(myGRAVY))
        GRAVY.SD <- c(GRAVY.SD, sd(myGRAVY))
        # calculate isoelectric point 20191027
        mypI <- pI(mycomp)
        pI <- c(pI, mean(mypI))
        pI.SD <- c(pI.SD, sd(mypI))
        # add basis argument (ie QEC or QCa) here
        myX <- Xfun(mycomp, basis)
      } else {
        # use base-paired (double-stranded) DNA
        if(dsDNA) mycomp <- make_dsDNA(mycomp)
        # calculate GC ratio 20180309
        if(!isprotein) GC <- c(GC, GCnuc(mycomp))
        myX <- Xfun(mycomp, "deoxyribose")
      }
      Xmean <- c(Xmean, mean(myX))
      Xlo <- c(Xlo, mean(myX) - sd(myX))
      Xhi <- c(Xhi, mean(myX) + sd(myX))
      # initialize data frame for mean compositions
      if(is.null(meancomp)) {
        meancomp <- mycomp[1:length(samples), ]
        meancomp[] <- NA
        rownames(meancomp) <- samples
      }
      # get mean compositions 20180505
      meancomp[i, ] <- colMeans(mycomp)
      # keep sample values for violin plot 20180515
      istart <- (i-1) * 100 + 1
      iend <- istart + 99
      DNA$ZC[istart:iend] <- myX
      DNA$sample[istart:iend] <- paste0("X", letters[i])
    }
  }
  # make plot
  if(plot.it) {
    # replace "X" (ZC or nH2O) with GRAVY or pI 20191028
    plotXmean <- Xmean
    plotXhi <- Xhi
    plotXlo <- Xlo
    if(identical(var, "GRAVY")) {
      plotXmean <- GRAVY
      plotXhi <- GRAVY + GRAVY.SD
      plotXlo <- GRAVY - GRAVY.SD
    }
    if(identical(var, "pI")) {
      plotXmean <- pI
      plotXhi <- pI + pI.SD
      plotXlo <- pI - pI.SD
    }
    isdeep <- logical(length(labels))
    if(is.null(all.labels)) all.labels <- labels
    if(plot_real_x) {
      if(!is.null(maxdepth)) isdeep <- labels > maxdepth
      xlim <- range(labels[!isdeep])
      if(extendrange) xlim <- extendrange(xlim)
      if(!grepl("Bison_Pool", dataset) & !grepl("Menez_Gwen", dataset)) xlim <- rev(xlim)
      atx <- labels
    } else {
      nsamp <- length(all.labels)
      xlim <- c(1, nsamp)
      if(extendrange) xlim <- extendrange(xlim)
      if(plottype=="violin") xlim[1] <- xlim[1] - 0.5
      if(plottype=="violin") xlim[length(xlim)] <- xlim[length(xlim)] + 0.5
      # for polygon plots: find the correct x-values (may not be 1:nsamp because of missing data) 20191005
      atx <- match(labels, all.labels)
    }
    if(!add & (is.null(taxid) | identical(taxid, 0))) {
      plot(0, 0, xlim=xlim, ylim=ylim, xlab=xlab, ylab=NA, xaxt="n")
      plot.labels <- all.labels
      if(!is.null(ilabel)) plot.labels[-ilabel] <- ""
      if(!is.numeric(plot.labels)) {
        # rotate labels 20190113
        if(!is.na(srt)) {
          # https://www.r-bloggers.com/rotated-axis-labels-in-r-plots/
          # modified to use offset calculated with strheight 20191004
          text(x=1:nsamp, y=par()$usr[3]-1.5*strheight("A"), labels=plot.labels, srt=srt, adj=1, xpd=TRUE)
          # add tick marks 
          axis(1, at=1:nsamp, labels=NA)
        } else {
          # loop over labels so that R doesn't omit any (because they're crowded)
          for(i in 1:nsamp) axis(1, at=i, labels=plot.labels[i])
        }
      } else {
        # for numeric labels, lower gap.axis to show more labels 20181215 (requires R 3.6.0)
        if(getRversion() >= "3.6.0") axis(1, at=atx, labels=plot.labels, gap.axis=0.02)
        else axis(1, at=atx, labels=plot.labels)
      }
      # add y-axis: ZC (or nH2O 20181231)
      if(identical(var, "GRAVY")) mtext("GRAVY", side = 2, line = yline, las = 0)
      else if(identical(var, "pI")) mtext("pI", side = 2, line = yline, las = 0)
      else if(H2O) mtext(quote(italic(n)[H[2]*O]), side=2, line=yline, las=0)
      else mtext(quote(italic(Z)[C]), side=2, line=yline, las=0)
    }
    if(plottype=="lines") {
      lines(atx[!isdeep], plotXlo[!isdeep], col=col, lty=3)
      lines(atx[!isdeep], plotXmean[!isdeep], col=col)
      lines(atx[!isdeep], plotXhi[!isdeep], col=col, lty=3)
    }
#    if(plottype=="violin") {
#      vioplotx(ZC~sample, DNA, add=TRUE, col = "palevioletred", plotCentre = "line", side = "left", pchMed = 21, colMed = "palevioletred4", colMed2 = "palevioletred2")
#    }
    if(plottype=="bars") {
      # error bar plots 20180515
      # apply small offset to x-position to separate DNA and RNA
      if(!is.null(taxid)) xx <- 0
      else xx <- abs(diff(par("usr")[1:2])) / 120
      # bars (whiskers) at one SD from mean (arrows trick from https://stackoverflow.com/questions/13032777/scatter-plot-with-error-bars)
      arrows(atx-xx, plotXlo, atx-xx, plotXhi, length = 0.03, angle = 90, code = 3, col=col, lwd=lwd.bars)
      # remove NAs so we can draw lines between all sites 20180529
      iNA <- is.na(plotXmean)
      lines((atx)[!iNA & !isdeep]-xx, plotXmean[!iNA & !isdeep], col=col, lty=lty, lwd=lwd)
    }
    if(substr(plottype, 1, 1)=="#") {
      # polygon plots 20191005
      # assemble top and bottom coordinates of polygon
      polygon(c(atx, rev(atx)), c(plotXhi, rev(plotXlo)), col = plottype, border = NA)
      # remove NAs so we can draw lines between all sites 20180529
      iNA <- is.na(plotXmean)
      lines((atx)[!iNA & !isdeep], plotXmean[!iNA & !isdeep], col=col, lty=lty, lwd=lwd)
      points((atx)[!iNA & !isdeep], plotXmean[!iNA & !isdeep], pch = pch, col=col)
    }
    # add points to show > 1% species abundance 20181118
    if(!is.null(taxid) & !identical(taxid, 0)) {
      # gradox data location in JMDplots package 20190928
      datadir <- system.file("extdata/gradox", package = "JMDplots")
      file <- paste0(datadir, "/one_percent/", dataset, ".csv")
      dat <- read.csv(file)
      # which samples have this taxid with at least 1% abundance?
      dat <- dat[dat$taxid==taxid, ]
      dat <- dat[dat$percentage >= 1, ]
      ioneperc <- samples %in% dat$sample
      points(atx[ioneperc], plotXmean[ioneperc], pch=pch, cex=0.5, col=col)
    }
  }
  # now do it for RNA
  if(!isprotein) {
    RNA_Xlo <- RNA_Xmean <- RNA_Xhi <- numeric()
    # keep sample values for violin plot 20180515
    RNA <- data.frame(ZC=numeric(length(samples)*100), sample=character(length(samples)*100), stringsAsFactors=FALSE)
    # gradox or gradH2O data location in JMDplots package 20190928
    if(!user_data) {
      datadir <- system.file(paste0("extdata/", paper), package = "JMDplots")
      rdsdata <- paste0(paper, "_MGR")
      RDS <- get(rdsdata, JMDplots)
      filestart <- paste0(dataset, "R")
    }
    for(i in 1:length(samples)) {
      sample <- samples[i]
      if(user_data) {
        file <- paste0(datadir, "/MGR/", dataset, "R_", sample, ".csv")
        fexists <- file.exists(file)
      } else {
        file <- paste0(filestart, "_", sample, ".csv")
        fexists <- file %in% names(RDS)
      }
      if(fexists) {
        if(user_data) myRNA <- read.csv(file, as.is=TRUE)
        else myRNA <- RDS[[file]]
        myX <- Xfun(myRNA, "ribose")
        RNA_Xmean <- c(RNA_Xmean, mean(myX))
        RNA_Xlo <- c(RNA_Xlo, mean(myX) - sd(myX))
        RNA_Xhi <- c(RNA_Xhi, mean(myX) + sd(myX))
      }
      # keep sample values for violin plot 20180515
      istart <- (i-1) * 100 + 1
      iend <- istart + 99
      RNA$ZC[istart:iend] <- myX
      RNA$sample[istart:iend] <- paste0("X", letters[i])
    }
    if(plot.RNA & length(RNA_Xmean) > 0 & plot.it) {
      # subtract offset for ZC of RNA
      dZC <- -0.28
      if(plottype=="lines") {
        lines(atx, RNA_Xlo + dZC, col="blue", lty=3)
        lines(atx, RNA_Xmean + dZC, col="blue")
        lines(atx, RNA_Xhi + dZC, col="blue", lty=3)
      }
#      if(plottype=="violin") {
#        RNA$ZC <- RNA$ZC + dZC
#        vioplotx(ZC~sample, RNA, add=TRUE, col = "lightblue", plotCentre = "line", side = "right", pchMed = 21, colMed = "lightblue4", colMed2 = "lightblue2")
#      }
      if(plottype=="bars") {
        # apply small offset to x-position to separate DNA and RNA
        xx <- abs(diff(par("usr")[1:2])) / 200
        arrows(atx+xx, RNA_Xlo + dZC, atx+xx, RNA_Xhi + dZC, length = 0.03, angle = 90, code = 3, col="blue", lwd=lwd.bars)
        lines(atx+xx, RNA_Xmean + dZC, col="blue", lty=2, lwd=lwd)
      }
    }
    # return ZC values
    outval <- list(DNA=Xmean, RNA=RNA_Xmean, GC=GC, group=group, meancomp=meancomp, abbrev=abbrev, techtype = techtype, dx = dx, dy = dy, H2O = H2O)
  } else {
    outval <- list(AA=Xmean, GRAVY=GRAVY, pI=pI, group=group, meancomp=meancomp, abbrev=abbrev, techtype = techtype, dx = dx, dy = dy, H2O = H2O)
  }
  # add title 20180225
  if((is.null(taxid) | identical(taxid, 0)) & plot.it & add.title & !add) {
    main <- dataset2main(dataset, abbrev)
    title(main=main, font.main=1)
    # add in-plot label: MG or MT 20180829
    if(add.label) {
      label <- "MG"
      if(grepl("_MT.*", dataset)) label <- "MT"
      # change position for some datasets
      xfrac <- 0.1
      if(grepl("OMZ", dataset)) xfrac <- 0.9
      if(grepl("Negro", dataset)) xfrac <- 0.9
      if(grepl("ALOHA", dataset)) xfrac <- 0.9
      CHNOSZ::label.plot(label, xfrac=xfrac, yfrac=0.9)
    }
  }
  # done!
  return(outval)
}

# parse dataset name to plot title 20180505
dataset2main <- function(dataset, abbrev=NULL) {
  main <- dataset
  main <- gsub("_SRA", "", main)
  main <- gsub("_IMG", "", main)
  main <- gsub("_MGRAST", "", main)
  main <- gsub("_GenBank", "", main)
  main <- gsub("_iMicrobe", "", main)
  # strip _MG or _MGP
  main <- gsub("_MG.*", "", main)
  main <- gsub("_MT.*", "", main)
  # replace underscore with space, and put space in Baltic Sea
  main <- gsub("_", " ", main)
  main <- gsub("BalticSea", "Baltic Sea", main)
  main <- gsub("Mud Volcano", "SYNH Mud Volcano", main)
  # add dataset abbreviation at end 20180829
  if(identical(abbrev, "GS")) {
    # replace "GS" (Global Survey of Baltic Sea) by "MG" or "MT"
    if(grepl("_MTP", dataset)) abbrev <- "MT"
    if(grepl("_MGP", dataset)) abbrev <- "MG"
  }
  if(!is.null(abbrev)) main <- paste0(main, " (", abbrev, ")")
  main
}

### utility functions for ZC, GC and nH2O calculations ###

# calculate GC ratio for given nucleobase compositions 20180309
GCnuc <- function(nuccomp) {
  # find columns with names for the nucleobases
  isbase <- colnames(nuccomp) %in% c("A", "C", "G", "T", "U")
  # keep only these columns
  nuccomp <- nuccomp[, isbase]
  # find GC columns and not-GC columns
  isGC <- colnames(nuccomp) %in% c("G", "C")
  isnotGC <- colnames(nuccomp) %in% c("A", "T", "U")
  # calculate GC ratio
  sum(nuccomp[, isGC]) / sum(nuccomp)
}

# calculate ZC for given nucleobase compositions 20171223
# changed to nucleosides 20180512
ZCnuc <- function(nuccomp, sugar="deoxyribose") {
  ## ZC and nC of nucleobases
  #ZC_nuc <- c(A=2, C=1.5, G=2.4, T=0.8, U=1.5)
  #nC_nuc <- c(A=5, C=4, G=5, T=5, U=4)
  # the number of carbons of the nucleosides (same in DNA and RNA)
  nC_nuc <- c(A=10, C=9, G=10, T=10, U=9)
  if(sugar=="deoxyribose") {
    # ZC of the nucleosides in DNA
    # DNAsides <- c("deoxyadenosine", "deoxycytidine", "deoxyguanosine", "deoxythymidine", "deoxyuridine")
    # ZC_nuc <- ZC(info(info(DNAsides))$formula)
    ZC_nuc <- c(A=0.8, C=4/9, G=1, T=0.2, U=4/9)
  }
  if(sugar=="ribose") {
    # ZC of the nucleosides in RNA
    # RNAsides <- c("adenosine", "cytidine", "guanosine", "thymidine", "uridine")
    # ZC_nuc <- ZC(info(info(RNAsides))$formula)
    ZC_nuc <- c(A=1, C=2/3, G=1.2, T=0.4, U=2/3)
  }
  # find columns with names for the nucleosides
  isbase <- colnames(nuccomp) %in% c("A", "C", "G", "T", "U")
  ibase <- match(colnames(nuccomp)[isbase], names(ZC_nuc))
  # calculate the nC from the counts
  multC <- t(t(nuccomp[, isbase]) * nC_nuc[ibase])
  # multiply nC by ZC to get total formal charge on carbon (Z)
  multZ <- t(t(multC) * ZC_nuc[ibase])
  # calculate the Z and nC in each composition (row)
  Ztot <- rowSums(multZ)
  nCtot <- rowSums(multC)
  # the ZC in each composition (row)
  Ztot / nCtot
}

# function to convert ssDNA base counts to dsDNA 20180318
make_dsDNA <- function(nuccomp) {
  nuccomp[, c("A", "C", "G", "T")] <- nuccomp[, c("A", "C", "G", "T")] + nuccomp[, c("T", "G", "C", "A")]
  nuccomp
}

