# JMDplots/grad_mplot.R
# plot chemical compositions of metagenomes
# 20180214 first version (MG.R)
# 20180831 updated for paper submission
# 20181215 updated for paper revision
# 20181222 gradH2O branched from gradox/MG.R
# 20190711 add plot.it argument to ppage, mplot, plotMG
# 20190930 merged with JMDplots package

# usage:
# mplot("BalticSea_Sediment", "SRA_MT")   # plot DNA and RNA ZC
# mplot("BalticSea_Sediment", "SRA_MTP")  # plot protein ZC
# mout <- mpage()    # plot DNA and RNA ZC for all datasets
# pout <- ppage()    # plot protein ZC for all datasets
# mcomp(mout)        # compare DNA ZC with RNA ZC or GC content
# pcomp(mout, pout)  # compare DNA ZC with protein ZC

# usage for sequences extracted for species classified by Kraken
# mplot("Menez_Gwen", "SRA_MG", taxid=c(1427364, 1420917))

# naming conventions:
# dataset = study_seqtype
# study = basestudy[|-substudy]
# basestudy = name1_name2
# seqtype = [SRA|IMG|MGRAST|GenBank]_[MG|MT]
# append [D|R|P] to seqtype to indicate processed data (DNA, RNA, protein)

# identify the datasets used for different figures in the papers
usedin <- list(
  # for figures in gradox paper
  gradoxMS = c("BalticSea_Sediment", "Shimokita_Peninsula",  # sediment
               "Diffuse_Vents", "Menez_Gwen",                # hydrothermal
               "ETNP_OMZ", "ETSP_OMZ",                       # ocean
               "Mono_Lake", "Organic_Lake",                  # hypersaline
               "Bison_Pool", "Guerrero_Negro"),              # microbial mat
  # for supporting information of gradox paper
  gradoxSI = names(gradox),
  # for goldschmidt poster 20190713
  gradoxGS = c("Diffuse_Vents", "Bison_Pool", "Guerrero_Negro"),              
  balticsurface = c("Baltic_Sea-0.1s"),
  balticdeep = c("Baltic_Sea-0.1d"),
  # all metagenomic data sets for Baltic Sea 20191004
  balticMG = c("Baltic_Sea-0.1s", "Baltic_Sea-0.8s", "Baltic_Sea-3.0s",
               "Baltic_Sea-0.1d", "Baltic_Sea-0.8d", "Baltic_Sea-3.0d"),
  balticMT = c("Baltic_Sea-0.1s", "Baltic_Sea-0.8s", "Baltic_Sea-3.0s",
               "Baltic_Sea-0.1d", "Baltic_Sea-0.8d", "Baltic_Sea-3.0d")
)

# function to plot sampled compositions for indicated study and sequence type 20180222
# e.g. mplot("Columbia_River", "IMG_MT")
mplot <- function(study, seqtype, plottype = "bars", ylim = NULL, plot.RNA = TRUE, taxid = NULL,
  dsDNA = TRUE, abbrev = NULL, col = NULL, add.label = TRUE, maxdepth = NULL, H2O = FALSE,
  plot.it = TRUE, add.title = TRUE, yline = 2, basis = "QEC", datadir = NULL, mdata = studies) {
  md <- get.mdata(mdata, study, seqtype)
  # get labels, groups, and abbreviation
  samples <- md$samples
  xlabels <- md$xlabels
  group <- md$group
  if(missing(abbrev)) abbrev <- md$abbrev
  techtype <- md$techtype
  dx <- md$dx
  dy <- md$dy
  xlab <- names(mdata[[study]][1])
  # get ZC range
  if(is.null(ylim)) {
    if(grepl("_MGP", seqtype)) ylim <- mdata[[study]][["MGP_range"]]
    else if(grepl("_MTP", seqtype)) ylim <- mdata[[study]][["MTP_range"]]
    else if(grepl("_MG", seqtype) & is.null(taxid)) ylim <- mdata[[study]][["MG_range"]]
    else if(grepl("_MG", seqtype) & !is.null(taxid)) ylim <- mdata[[study]][["MG_srange"]]
    else if(grepl("_MT", seqtype) & is.null(taxid)) ylim <- mdata[[study]][["MT_range"]]
    else if(grepl("_MT", seqtype) & !is.null(taxid)) ylim <- mdata[[study]][["MT_srange"]]
    else stop("invalid seqtype:", seqtype)
    if(is.null(ylim)) stop("ylim range not available for ", study, " ", seqtype)
  }
  # reverse the scale?
  if(diff(ylim) < 0) {
    samples <- rev(md$samples)
    xlabels <- rev(md$xlabels)
    ylim <- rev(ylim)
    group <- rev(group)
  }
  # if taxids are given, plot total ZC first (DNA, not RNA) and add ZC for each species
  if(!is.null(taxid)) {
    for(i in 1:length(taxid)) {
      # color lines for species common to multiple datasets (Marinobacter salarius, Thioglobus singularis, Pelagibacter ubique)
      col <- "black"
      if(taxid[i]==1420917) col <- "blue"
      if(taxid[i]==1427364) col <- "red"
      if(taxid[i]==1977865) col <- "darkorange"
      if(taxid[i]==1898749) col <- "purple"
      if(taxid[i]==1410606) col <- "darkgreen"
      # for taxid 0 (all sequences after dereplication), use thick solid lines 20181117
      lwd <- ifelse(taxid[i]==0, 1.8, 1.2)
      lty <- ifelse(taxid[i]==0, 1, 2)
      lwd.bars <- ifelse(taxid[i]==0, 2, 1)
      add.label <- ifelse(taxid[i]==0, FALSE, add.label)
      plotMG(paste0(study, "_", seqtype), plottype, samples, xlabels, group, xlab, ylim, abbrev, dsDNA,
             plot.RNA = FALSE, taxid = taxid[i], lwd = lwd, lty = lty, lwd.bars = lwd.bars, col = col, extendrange = TRUE,
             add.label = add.label, plot_real_x = mdata[[study]][["plot_real_x"]], maxdepth = maxdepth, H2O = H2O,
             plot.it = plot.it, add.title = add.title, yline = yline, basis = basis, techtype = techtype, dx = dx, dy = dy, datadir = datadir)
    }
  } 
  else plotMG(paste0(study, "_", seqtype), plottype, samples, xlabels, group, xlab, ylim, abbrev, dsDNA,
              plot.RNA, taxid, col = col,
              add.label = add.label, plot_real_x = mdata[[study]][["plot_real_x"]], maxdepth = maxdepth, H2O = H2O,
              plot.it = plot.it, add.title = add.title, yline = yline, basis = basis, techtype = techtype, dx = dx, dy = dy, datadir = datadir)
}

# make page of plots for MG/MT 20180225
# add subset and H2O arguments 20181231
mpage <- function(subset="gradoxSI", H2O=FALSE, plottype="bars", dsDNA=TRUE, set.par=TRUE, add.label = TRUE, mfrow = NULL) {
  if(is.list(subset)) {
    # when subset is a list, it gives the studies
    studies <- subset
    # use user's directory for data
    datadir <- "data/"
  } else {
    # extract only the datasets used for the paper
    studies <- studies[match(usedin[[subset]], names(studies))]
    # if studies is empty, we have a problem
    if(length(studies)==0) stop("no data for ", subset, "; available groups are ", paste(names(usedin), collapse = " "))
    # use NULL datadir for mplot() (gets data from JMDplots/extdata)
    datadir <- NULL
  }
  # plot setup
  if(is.null(mfrow)) {
    mfrow <- c(1, 1)
    if(length(studies) > 1) mfrow <- c(1, 2)
    if(length(studies) > 2) mfrow <- c(2, 2)
    if(length(studies) > 4) mfrow <- c(2, 3)
    if(length(studies) > 6) mfrow <- c(3, 3)
    if(length(studies) > 9) mfrow <- c(3, 4)
    if(length(studies) > 12) mfrow <- c(3, 5)
    if(length(studies) > 15) mfrow <- c(4, 5)
  }
  if(set.par) opar <- par(mfrow=mfrow, mar=c(4, 3.5, 2, 1), mgp=c(2.5, 1, 0))
  # variable for labeling plots
  iletter <- 1
  # initialize output of ZC values
  mout <- list()
  for(i in 1:length(studies)) {
    # name of the study without seqtype (e.g. Guerrero_Negro)
    study <- names(studies[i])
    # skip HOT_ALOHA-2010
    if(grepl("-", study)) next
    for(j in 1:length(studies[[i]])) {
      # the seqtype (e.g. IMG_MG)
      seqtype <- names(studies[[i]][j])
      # the seqtype isn't the samples (first position) or the xlabels, ZC range, abbreviation, group, dx, or dx
      if(j==1 | grepl("xlabels", seqtype) | grepl("range", seqtype) |
         seqtype=="abbrev" | seqtype=="group" | seqtype=="dx" | seqtype=="dy" | seqtype=="plot_real_x" | seqtype == "techtype") next
      # for the figure in the paper, take only one dataset for each study (i.e. metagenome, except for Mono Lake)
      if(identical(subset, "gradoxMS") & (grepl("_MT", seqtype) & study!="Mono_Lake")) next
      ZC <- list(mplot(study, seqtype, plottype, dsDNA=dsDNA, add.label=add.label, H2O=H2O, datadir = datadir, mdata = studies))
      names(ZC) <- paste0(study, "_", seqtype)
      mout <- c(mout, ZC)
    }
    # add figure label 20181210
    if(identical(subset, "gradoxMS")) {
      label.figure(LETTERS[iletter], cex=1.6, font=2, yfrac=0.936)
      iletter <- iletter + 2
    }
  }
  if(set.par) par(opar)
  mout
}

# scatterplot of DNA and RNA ZC 20180307
# mout <- mpage(); mcomp(mout)
# yvar can be RNA or GC
mcomp <- function(mout, yvar="RNA") {
  # set up plot
  if(yvar=="RNA") plot(c(0.57, 0.65), c(0.22, 0.32), xlab="ZC(DNA)", ylab="DZC(RNA - DNA)", type="n")
  if(yvar=="GC") plot(c(0.57, 0.65), c(0.3, 0.65), xlab="ZC(DNA)", ylab="GC", type="n")
  for(i in 1:length(mout)) {
    DNA <- mout[[i]]$DNA
    if(yvar=="RNA") Y <- mout[[i]]$RNA
    if(yvar=="GC") Y <- mout[[i]]$GC
    # order points by increasing DNA value
    ord <- order(DNA)
    DNA <- DNA[ord]
    Y <- Y[ord]
    # for RNA, plot difference from DNA
    if(yvar=="RNA") Y <- Y - DNA
    # color: red for MG, blue for MT
    if(grepl("_MG", names(mout[i]))) col <- "red"
    if(grepl("_MT", names(mout[i]))) col <- "blue"
    # plot lines
    lines(DNA, Y, col=col)
  }
}

# make page of plots for MGP/MTP 20180225
# add subset and H2O arguments 20181231
# add plot.it argument 20190711
ppage <- function(subset = "gradoxSI", H2O = FALSE, set.par = TRUE, plot.it = TRUE, basis = "QEC", add.label = TRUE, mfrow = NULL) {
  if(is.list(subset)) {
    # when subset is a list, it gives the studies
    studies <- subset
    # use user's directory for data
    datadir <- "data/"
  } else {
    # extract only the datasets used for the paper
    studies <- studies[match(usedin[[subset]], names(studies))]
    # if studies is empty, we have a problem
    if(length(studies)==0) stop("no data for ", subset, "; available groups are ", paste(names(usedin), collapse = " "))
    # use NULL datadir for mplot() (gets data from JMDplots/extdata)
    datadir <- NULL
  }
  # plot setup
  if(is.null(mfrow)) {
    mfrow <- c(1, 1)
    if(length(studies) > 1) mfrow <- c(1, 2)
    if(length(studies) > 2) mfrow <- c(2, 2)
    if(length(studies) > 4) mfrow <- c(2, 3)
    if(length(studies) > 6) mfrow <- c(3, 3)
    if(length(studies) > 9) mfrow <- c(3, 4)
    if(length(studies) > 12) mfrow <- c(3, 5)
    if(length(studies) > 15) mfrow <- c(4, 5)
  }
  if(set.par & plot.it) opar <- par(mfrow = mfrow, mar = c(4, 3.5, 2, 1), mgp = c(2.5, 1, 0))
  # variable for labeling plots (starts at B for protein plots in gradox paper SI figure)
  iletter <- 2
  # initialize output of ZC or H2O values
  pout <- list()
  for(i in 1:length(studies)) {
    # name of the study without seqtype (e.g. Columbia_River)
    study <- names(studies[i])
    # we only handle studies with proteins (MGP or MTP)
    if(!any(c("MGP_range", "MTP_range") %in% names(studies[[i]]))) next
    for(j in 1:length(studies[[i]])) {
      # the seqtype
      seqtype <- names(studies[[i]][j])
      # the seqtype isn't the samples (first position) or the xlabels, ZC or H2O range, abbreviation, group, dx, or dy
      if(j==1 | grepl("xlabels", seqtype) | grepl("range", seqtype) |
         seqtype=="abbrev" | seqtype=="group" | seqtype=="dx" | seqtype=="dy" | seqtype=="plot_real_x" | seqtype=="techtype") next
      # for the figure in the paper, take only one dataset for each study (i.e. metagenome, except for Mono Lake)
      if(identical(subset, "gradoxMS") & (grepl("_MT", seqtype) & study!="Mono_Lake")) next
      # for Baltic Sea, use either MG or MT 20190713
      if(identical(subset, "balticMG") & !grepl("_MG", seqtype)) next
      if(identical(subset, "balticMT") & !grepl("_MT", seqtype)) next
      # for comparison with Eiler et al. data, use only metagenomes 20190720
      if(identical(subset, "eiler") & !grepl("_MG", seqtype)) next
      # add "P" for proteins
      seqtype <- paste0(seqtype, "P")
      X <- list(mplot(study, seqtype, add.label = add.label, H2O = H2O, plot.it = plot.it, basis = basis, datadir = datadir, mdata = studies))
      names(X) <- paste0(study, "_", seqtype)
      pout <- c(pout, X)
    }
    # add figure label 20181210
    if(identical(subset, "gradoxMS") & plot.it) {
      label.figure(LETTERS[iletter], cex = 1.6, font = 2, yfrac = 0.936)
      iletter <- iletter + 2
    }
  }
  if(set.par & plot.it) par(opar)
  pout
}

# scatterplot of DNA and protein ZC 20180307
# mout <- mpage(); pout <- ppage(); pcomp(mout, pout, type="ZC")
# mout <- ppage(); pout <- ppage(H2O=TRUE); pcomp(mout, pout, type="both")
# mout <- mpage(H2O=TRUE); pout <- ppage(H2O=TRUE); pcomp(mout, pout, type="H2O")
pcomp <- function(mout, pout, seqtype="MG", type="ZC", parts=c("plot", "legend"), yline = 2,
                  xlim = NULL, ylim = NULL, reorder = TRUE, plot.techtype = FALSE, add = FALSE, pch = NULL) {
  if("plot" %in% parts) {
    # set up plot
    if(type=="ZC") {
      if(seqtype=="MG" & is.null(xlim)) xlim <- c(0.58, 0.65)
      if(seqtype=="MT" & is.null(xlim)) xlim <- c(0.585, 0.63)
      xlab <- quote(italic(Z)[C]~of~DNA)
      if(seqtype=="MG" & is.null(ylim)) ylim <- c(-0.22, -0.098)
      if(seqtype=="MT" & is.null(ylim)) ylim <- c(-0.21, -0.14)
      ylab <- quote(italic(Z)[C]~of~protein)
    }
    if(type=="both") {
      if(is.null(xlim)) xlim <- c(-0.22, -0.098)
      xlab <- quote(italic(Z)[C]~of~protein)
      if(is.null(ylim)) ylim <- c(-0.78, -0.72)
      ylab <- quote(italic(n)[H[2]*O])
    }
#    if(xaxis=="RNA") {xlim <- c(1.85, 1.92); xlab <- "ZC(RNA)"; dx <- -0.004}
#    if(xaxis=="CM") {xlim <- c(0.024, 0.037); xlab <- "Cys+Met fraction"; dx <- -0.0003}
#    if(yaxis=="CM") {ylim <- c(0.024, 0.037); ylab <- "Cys+Met fraction"}
    if(!add) {
      plot(xlim, ylim, xlab=xlab, ylab=NA, type="n")
      mtext(ylab, side=2, line=yline, las=0, cex=par("cex"))
    }
    for(ipout in 1:length(pout)) {
      # get study name and change MGP or MTP to MG or MT
      study <- names(pout[ipout])
      study <- gsub("MGP", "MG", study)
      study <- gsub("MTP", "MT", study)
      # select MG or MT
      if(seqtype=="MG" & !grepl("_MG$", study)) next
      if(seqtype=="MT" & !grepl("_MT$", study)) next
      # assemble ZC or nH2O values and get environment group
      if(type=="ZC") {
        # find this study in mout
        imout <- match(study, names(mout))
        xvals <- mout[[imout]]$DNA
        yvals <- pout[[ipout]]$AA
      }
      if(type=="both") {
        imout <- ipout
        xvals <- mout[[imout]]$AA
        yvals <- pout[[ipout]]$AA
      }
#      if(xaxis=="RNA") xvals <- mout[[imout]]$RNA
#      if(xaxis=="CM") xvals <- pout[[i]]$CM
#      if(yaxis=="CM") yvals <- pout[[i]]$CM
      group <- rep(mout[[imout]]$group, length.out = length(xvals))
      if(reorder) {
        # order points by increasing DNA/RNA ZC value
        ord <- order(xvals)
        xvals <- xvals[ord]
        yvals <- yvals[ord]
        group <- group[ord]
      }
      # color: by group
      col <- rep("black", length(group))
      col[group %in% c("yellowstone", "yellowstone1")] <- "orange"
      col[group %in% c("rock", "rock0")] <- "brown"
      col[group %in% c("vent", "vent0")] <- "red"
      col[group %in% c("mat", "mat1")] <- "green3"
      col[group %in% c("oxic", "oxic0")] <- "blue"
      col[group %in% c("hypersaline", "hypersaline0")] <- "turquoise3"
      col[group %in% c("OMZ", "OMZ0")] <- "black"
      col[group %in% c("plume", "plume0", "plumePA")] <- "purple1"
      col[group %in% c("sediment", "sediment0")] <- "slategrey"
      col[group %in% c("river", "riverPA")] <- "green3"
      col[group %in% c("lake")] <- "yellowgreen"
      col[group %in% c("HSsediment", "HSsediment0")] <- "slategrey"
      # plot lines and points
      # for lines, use the color of most of the points 20180501
      # change this to gray 20181114
      lines(xvals, yvals, col="dimgray", lwd=0.8, lty=2)
      if(grepl("Baltic_Sea", study)) {
        # filled triangle, circle, square for 0.1, 0.8, 3.0 size fractions
        if(grepl("0.1", study)) mypch <- 17
        if(grepl("0.8", study)) mypch <- 16
        if(grepl("3.0", study)) mypch <- 15
      } else {
        # filled circles for marine, filled squares for terrestrial 20181114
        mypch <- rep(19, length(xvals))
        mypch[group %in% c("yellowstone", "yellowstone1", "rock", "rock0", "mat", "mat1", "hypersaline", "hypersaline0")] <- 15
        # filled squares for Amazon River particle associated 20190723
        mypch[group %in% c("riverPA", "plumePA")] <- 15
      }
      studyname <- paste(strsplit(study, "_")[[1]][1:2], collapse="_")
      if(plot.techtype) {
        # use open symbols for 454 or Sanger techtypes 20190723
        techtype <- pout[[ipout]][["techtype"]]
        techtype <- rep(techtype, length(xvals))
        mypch[mypch==19 & techtype %in% c("454", "Sanger")] <- 1
        mypch[mypch==15 & techtype %in% c("454", "Sanger")] <- 0
      }
      # use pch from argument, otherwise mypch as determined here
      if(is.null(pch)) thispch <- mypch else thispch <- pch[ipout]
      points(xvals, yvals, col=col, pch=thispch)
      # outline circle or square for groups with "0" (upper part of OMZ etc)
      pch0 <- rep(1, length(xvals))
      pch0[group %in% c("yellowstone1", "rock0", "mat1", "hypersaline0")] <- 0
      cex0 <- rep(1.8, length(xvals))
      cex0[group %in% c("yellowstone1", "rock0", "mat1", "hypersaline0")] <- 1.6
      i0 <- grepl("0", group)
      if(any(i0)) points(xvals[i0], yvals[i0], pch=pch0[i0], col=col[i0], cex=cex0[i0])
      # outline green circle or square for groups with "1" (phototrophic mat in Yellowstone / upper part of Guerrero Negro mat)
      i1 <- grepl("1", group)
      if(any(i1)) points(xvals[i1], yvals[i1], pch=pch0[i1], col="green3", cex=cex0[i1])
      # get abbreviation for this study
      abbrev <- pout[[ipout]][["abbrev"]]
      dx <- pout[[ipout]][["dx"]][[seqtype]]
      dy <- pout[[ipout]][["dy"]][[seqtype]]
      if(is.null(dx) | type!="ZC") dx <- 0
      if(is.null(dy) | type!="ZC") {
        if(type=="both") dy <- 0.002
        else dy <- 0.003
      }
      # add text label
      #text(xvals[1] + dx, yvals[1], abbrev, cex=0.7)
      imax <- which.max(yvals)
      if(abbrev != "GS") text(xvals[imax] + dx, yvals[imax] + dy, abbrev, cex=0.7)
    }
  }
  # add legend
  if(type=="ZC" & "legend" %in% parts) {
    if("plot" %in% parts) {
      marine.x <- "topleft"
      terrestrial.x <- "bottomright"
      bty <- "o"
    } else {
      marine.x <- "left"
      terrestrial.x <- "right"
      bty <- "n"
    }
    # marine
    legend(-0.063, 1.014, c("", "", "", "", "", "", "", "", ""), lty=2, bty="n", col=c(NA, rep("dimgray", 8)))
    legend(marine.x, lty=2, lwd=NA, bty=bty, pch=19,
           legend=c("MARINE", "OMZ", "ocean surface", "oxic", "vent", "plume", "background SW", "sediment", "sediment surface"),
           col=c(NA, "black", "black", "blue", "red", "purple1", "purple1", "slategrey", "slategrey"))
    # overlay ocean surface, seawater, and sediment surface points (blue circles)
    legend(marine.x, lty=0, lwd=0, bty=bty, pch=1, pt.cex=1.8, pt.lwd=1,
           legend=c("", "", "", "", "", "", "", "", ""),
           col=c(NA, NA, "black", NA, NA, NA, "purple1", NA, "slategrey"))
    # terrestrial
    legend(0.47, 1.014, c("", "", "", "", "", "", "", "", ""), lty=2, bty="n", col=c(NA, rep("dimgray", 8)))
    legend(terrestrial.x, lty=2, lwd=NA, pch=15,
           legend=c("TERRESTRIAL", "hypersaline", "lake surface", "mat", "upper mat", "rock-derived", "surface mixing", "Yellowstone", "phototrophic mat"),
           col=c(NA, "turquoise3", "turquoise3", "green3", "green3", "brown", "brown", "orange", "orange"), bty=bty)
    # overlay lake surface and phototrophic mat points (blue and green)
    legend(0.492, 1.014, lty=0, lwd=0, bty=bty, pch=0, pt.cex=1.6, pt.lwd=1,
           legend=c("", "", "", "", "", "", "", "", ""),
           col=c(NA, NA, "turquoise3", NA, "green3", NA, "brown", NA, "green3"))
  }
}

