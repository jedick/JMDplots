# JMDplots/grad_util.R
# utility functions for gradox and gradH2O papers

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
  basis = "rQEC", techtype = NULL, dx = NULL, dy = NULL, datadir = NULL,
  add = FALSE, all.labels = NULL, pch = 19) {
  # samples: (used for suffixes on file names)
  # labels: (used for labeling x-axis ticks)
  # xlab: axis label: "layer", "depth", ...
  isprotein <- grepl("_MGP$", dataset) | grepl("_MTP$", dataset)
  # where to keep mean and high/lo (+/- SD) of ZC at each site
  CM <- GC <- ZClo <- ZCmean <- ZChi <- numeric()
  # which paper is this dataset from?
  if(any(sapply(names(gradox), grepl, dataset))) paper <- "gradox"
  if(any(sapply(names(gradH2O), grepl, dataset))) paper <- "gradH2O"
  # gradox or gradH2O data location in JMDplots package 20190928
  if(is.null(datadir)) datadir <- system.file(paste0("extdata/", paper), package = "JMDplots")
  # set up for proteins or DNA
  if(isprotein) {
    filestart <- paste0(datadir, "/MGP/", dataset)
    if(H2O) ZCfun <- H2OAA else ZCfun <- ZCAA
    if(is.null(col)) col <- "darkgreen"
  } else {
    # data directory for DNA
    filestart <- paste0(datadir, "/MGD/", dataset, "D")
    ZCfun <- ZCnuc
    if(is.null(col)) col <- "red"
  }
  # keep sample values for violin plot 20180515
  DNA <- data.frame(ZC=numeric(length(samples)*100), sample=character(length(samples)*100), stringsAsFactors=FALSE)
  meancomp <- NULL
  for(i in 1:length(samples)) {
    sample <- samples[i]
    if(!is.null(taxid)) file <- paste0(filestart, "_", sample, "-", taxid, ".csv")
    else file <- paste0(filestart, "_", sample, ".csv")
    # use NA if the file is missing 20180529
    if(!file.exists(file)) {
      ZCmean <- c(ZCmean, NA)
      ZClo <- c(ZClo, NA)
      ZChi <- c(ZChi, NA)
      message(paste("missing file", file))
    } else {
      mycomp <- read.csv(file, as.is=TRUE)
      if(isprotein) {
        # calculate Cys+Met fraction 20180324
        CM <- c(CM, CMAA(mycomp))
        # add basis argument (QEC or rQEC) here
        myZC <- ZCfun(mycomp, basis)
      } else {
        # use base-paired (double-stranded) DNA
        if(dsDNA) mycomp <- make_dsDNA(mycomp)
        # calculate GC ratio 20180309
        if(!isprotein) GC <- c(GC, GCnuc(mycomp))
        myZC <- ZCfun(mycomp, "deoxyribose")
      }
      ZCmean <- c(ZCmean, mean(myZC))
      ZClo <- c(ZClo, mean(myZC) - sd(myZC))
      ZChi <- c(ZChi, mean(myZC) + sd(myZC))
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
      DNA$ZC[istart:iend] <- myZC
      DNA$sample[istart:iend] <- paste0("X", letters[i])
    }
  }
  # make plot
  if(plot.it) {
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
      srt <- 45
      if(!is.numeric(all.labels)) {
        # rotate labels 20190113
        if(!is.na(srt)) {
          # tweak to remove "GS" from Baltic Sea labels 20190715
          if(grepl("^GS", all.labels[1])) all.labels <- gsub("^GS", "", all.labels)
          # https://www.r-bloggers.com/rotated-axis-labels-in-r-plots/
          # modified to use offset calculated with strheight 20191004
          text(x=1:nsamp, y=par()$usr[3]-1.5*strheight("A"), labels=all.labels, srt=srt, adj=1, xpd=TRUE)
          # add tick marks 
          axis(1, at=1:nsamp, labels=NA)
        } else {
          # loop over labels so that R doesn't omit any (because they're crowded)
          for(i in 1:nsamp) axis(1, at=i, labels=all.labels[i])
        }
      } else {
        # for numeric labels, lower gap.axis to show more labels 20181215 (requires R 3.6.0)
        if(getRversion() >= "3.6.0") axis(1, at=atx, labels=all.labels, gap.axis=0.02)
        else axis(1, at=atx, labels=all.labels)
      }
      # add y-axis: ZC (or nH2O 20181231)
      if(H2O) mtext(quote(italic(n)[H[2]*O]), side=2, line=yline, las=0)
      else mtext(quote(italic(Z)[C]), side=2, line=yline, las=0)
    }
    if(plottype=="lines") {
      lines(atx[!isdeep], ZClo[!isdeep], col=col, lty=3)
      lines(atx[!isdeep], ZCmean[!isdeep], col=col)
      lines(atx[!isdeep], ZChi[!isdeep], col=col, lty=3)
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
      arrows(atx-xx, ZClo, atx-xx, ZChi, length = 0.03, angle = 90, code = 3, col=col, lwd=lwd.bars)
      # remove NAs so we can draw lines between all sites 20180529
      iNA <- is.na(ZCmean)
      lines((atx)[!iNA & !isdeep]-xx, ZCmean[!iNA & !isdeep], col=col, lty=lty, lwd=lwd)
    }
    if(substr(plottype, 1, 1)=="#") {
      # polygon plots 20191005
      # assemble top and bottom coordinates of polygon
      polygon(c(atx, rev(atx)), c(ZChi, rev(ZClo)), col = plottype, border = NA)
      # remove NAs so we can draw lines between all sites 20180529
      iNA <- is.na(ZCmean)
      lines((atx)[!iNA & !isdeep], ZCmean[!iNA & !isdeep], col=col, lty=lty, lwd=lwd)
      points((atx)[!iNA & !isdeep], ZCmean[!iNA & !isdeep], pch = pch, col=col)
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
      points(atx[ioneperc], ZCmean[ioneperc], pch=pch, cex=0.5, col=col)
    }
  }
  # now do it for RNA
  if(!isprotein) {
    RNA_ZClo <- RNA_ZCmean <- RNA_ZChi <- numeric()
    # keep sample values for violin plot 20180515
    RNA <- data.frame(ZC=numeric(length(samples)*100), sample=character(length(samples)*100), stringsAsFactors=FALSE)
    for(i in 1:length(samples)) {
      sample <- samples[i]
      # gradox or gradH2O data location in JMDplots package 20190928
      if(is.null(datadir)) datadir <- system.file(paste0("extdata/", paper), package = "JMDplots")
      file <- paste0(datadir, "/MGR/", dataset, "R_", sample, ".csv")
      if(file.exists(file)) {
        myRNA <- read.csv(file, as.is=TRUE)
        myZC <- ZCfun(myRNA, "ribose")
        RNA_ZCmean <- c(RNA_ZCmean, mean(myZC))
        RNA_ZClo <- c(RNA_ZClo, mean(myZC) - sd(myZC))
        RNA_ZChi <- c(RNA_ZChi, mean(myZC) + sd(myZC))
      }
      # keep sample values for violin plot 20180515
      istart <- (i-1) * 100 + 1
      iend <- istart + 99
      RNA$ZC[istart:iend] <- myZC
      RNA$sample[istart:iend] <- paste0("X", letters[i])
    }
    if(plot.RNA & length(RNA_ZCmean) > 0 & plot.it) {
      # subtract offset for ZC of RNA
      dZC <- -0.28
      if(plottype=="lines") {
        lines(atx, RNA_ZClo + dZC, col="blue", lty=3)
        lines(atx, RNA_ZCmean + dZC, col="blue")
        lines(atx, RNA_ZChi + dZC, col="blue", lty=3)
      }
#      if(plottype=="violin") {
#        RNA$ZC <- RNA$ZC + dZC
#        vioplotx(ZC~sample, RNA, add=TRUE, col = "lightblue", plotCentre = "line", side = "right", pchMed = 21, colMed = "lightblue4", colMed2 = "lightblue2")
#      }
      if(plottype=="bars") {
        # apply small offset to x-position to separate DNA and RNA
        xx <- abs(diff(par("usr")[1:2])) / 200
        arrows(atx+xx, RNA_ZClo + dZC, atx+xx, RNA_ZChi + dZC, length = 0.03, angle = 90, code = 3, col="blue", lwd=lwd.bars)
        lines(atx+xx, RNA_ZCmean + dZC, col="blue", lty=2, lwd=lwd)
      }
    }
    # return ZC values
    outval <- list(DNA=ZCmean, RNA=RNA_ZCmean, GC=GC, group=group, meancomp=meancomp, abbrev=abbrev, techtype = techtype, dx = dx, dy = dy)
  } else {
    outval <- list(AA=ZCmean, CM=CM, group=group, meancomp=meancomp, abbrev=abbrev, techtype = techtype, dx = dx, dy = dy)
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

# calculate Cys+Met fraction of amino acid compositions 20180324
CMAA <- function(AAcomp) {
  # find columns with names for the amino acids
  AA <- c("Ala", "Cys", "Asp", "Glu", "Phe", "Gly", "His", "Ile", "Lys",
    "Leu", "Met", "Asn", "Pro", "Gln", "Arg", "Ser", "Thr", "Val", "Trp", "Tyr")
  isAA <- colnames(AAcomp) %in% AA
  # columns for Cys and Met
  isCM <- colnames(AAcomp) %in% c("Cys", "Met")
  # calculate (Cys+Met) / (total AA)
  sum(AAcomp[, isCM]) / sum(AAcomp[, isAA])
}

#########################################
### newly exported functions 20191005 ###
#########################################

# calculate ZC for amino acid compositions 20180228
ZCAA <- function(AAcomp, nothing=NULL) {
  # a dummy second argument is needed because of how this function is used in plotMG()
  # the number of carbons of the amino acids
  nC_AA <- c(Ala = 3, Cys = 3, Asp = 4, Glu = 5, Phe = 9, Gly = 2, His = 6, 
    Ile = 6, Lys = 6, Leu = 6, Met = 5, Asn = 4, Pro = 5, Gln = 5, 
    Arg = 6, Ser = 3, Thr = 4, Val = 5, Trp = 11, Tyr = 9)
  # the Ztot of the amino acids == CHNOSZ::ZC(info(info(aminoacids("")))$formula) * nC_AA
  Ztot_AA <- c(Ala = 0, Cys = 2, Asp = 4, Glu = 2, Phe = -4, Gly = 2, His = 4, 
    Ile = -6, Lys = -4, Leu = -6, Met = -2, Asn = 4, Pro = -2, Gln = 2, 
    Arg = 2, Ser = 2, Thr = 0, Val = -4, Trp = -2, Tyr = -2)
  # the ZC of the amino acids == CHNOSZ::ZC(info(info(aminoacids("")))$formula)
  ZC_AA <- Ztot_AA / nC_AA
  # find columns with names for the amino acids
  isAA <- colnames(AAcomp) %in% c("Ala", "Cys", "Asp", "Glu", "Phe", "Gly", "His", "Ile", "Lys", 
    "Leu", "Met", "Asn", "Pro", "Gln", "Arg", "Ser", "Thr", "Val", "Trp", "Tyr")
  iAA <- match(colnames(AAcomp)[isAA], names(ZC_AA))
  # calculate the nC for all occurrences of each amino acid
  multC <- t(t(AAcomp[, isAA]) * nC_AA[iAA])
  # multiply nC by ZC
  multZC <- t(t(multC) * ZC_AA[iAA])
  # calculate the total ZC and nC, then the overall ZC
  ZCtot <- rowSums(multZC)
  nCtot <- rowSums(multC)
  ZCtot / nCtot
}

# calculate nH2O for amino acid compositions 20181228
H2OAA <- function(AAcomp, basis = "rQEC") {
  # how to use CHNOSZ to get the number of H2O in reactions
  # to form amino acid residues from the "QEC" basis:
  ## basis("QEC")
  ## species(aminoacids(3))
  ## nH2O_AA <- species()[["H2O"]]
  # subtract one H2O to make residues
  ## nH2O_AA <- nH2O_AA - 1
  ## names(nH2O_AA) <- aminoacids(3)
  ## dput(nH2O_AA)
  if(basis == "QEC") {
    nH2O_AA <- c( Ala = -0.4, Cys =   -1, Asp = -1.2, Glu =   -1, Phe = -3.2, Gly = -0.6, His = -2.8,
      Ile =  0.2, Lys =  0.2, Leu =  0.2, Met = -0.6, Asn = -1.2, Pro =   -1, Gln =   -1,
      Arg = -0.8, Ser = -0.4, Thr = -0.2, Val =    0, Trp = -4.8, Tyr = -3.2)
  }
  # residual water content with QEC basis
  ## round(residuals(lm(nH2O_AA ~ ZC(species()$ispecies))), 3)
  if(basis == "rQEC") {
    nH2O_AA <- c(Ala = 0.724, Cys = 0.33, Asp = 0.233, Glu = 0.248, Phe = -2.213,
      Gly = 0.833, His = -1.47, Ile = 1.015, Lys = 1.118, Leu = 1.015,
      Met = 0.401, Asn = 0.233, Pro = 0.001, Gln = 0.248, Arg = 0.427,
      Ser = 0.93, Thr = 0.924, Val = 0.877, Trp = -3.732, Tyr = -2.144)
  }
  # find columns with names for the amino acids
  isAA <- colnames(AAcomp) %in% names(nH2O_AA)
  iAA <- match(colnames(AAcomp)[isAA], names(nH2O_AA))
  # calculate total number of H2O in reactions to form proteins
  nH2O <- rowSums(t(t(AAcomp[, isAA]) * nH2O_AA[iAA]))
  # add one to account for terminal groups
  nH2O <- nH2O + 1
  # divide by number of residues (length of protein)
  nH2O / rowSums(AAcomp[, isAA])
  # to check this function:
  #  basis("QEC")
  #  H2O.ref <- protein.basis(1:6)[, "H2O"] / protein.length(1:6)
  #  AAcomp <- thermo()$protein[1:6, ]
  #  H2O.fun <- H2OAA(AAcomp, "QEC")
  #  stopifnot(H2O.ref == H2O.fun)
}

# calculate ZC and nH2O of proteomes encoding different Nif homologs (Poudel et al., 2018)
# 20191014
NifProteomes <- function() {
  # read file with Nif genome classifications and taxids
  Niffile <- system.file("extdata/gradH2O/Nif_homolog_genomes.csv", package = "JMDplots")
  Nif <- read.csv(Niffile, as.is = TRUE)
  # drop NA taxids
  Nif <- Nif[!is.na(Nif$taxid), ]
  # read refseq data
  RSfile <- system.file("extdata/refseq/protein_refseq.csv.xz", package = "JMDplots")
  refseq <- read.csv(RSfile, as.is = TRUE)
  # the Nif types, arranged from anaerobic to aerobic
  types <- c("Nif-D", "Nif-C", "Nif-B", "Nif-A")
  # assemble the compositional metrics
  ZC <- ZC.SD <- nH2O <- nH2O.SD <- numeric()
  for(type in types) {
    # get the taxids for genomes with this type of Nif
    iNif <- Nif$Type == type
    taxid <- Nif$taxid[iNif]
    # get the row number in the refseq data frame
    irefseq <- match(taxid, refseq$organism)
    # include only organisms with at least 1000 protein sequences
    i1000 <- refseq$chains[irefseq] >= 1000
    irefseq <- irefseq[i1000]
    print(paste(type, "represented by", length(irefseq), "genomes with at least 1000 protein sequences"))
    # get the amino acid composition from refseq
    AAcomp <- refseq[irefseq, ]
    # calculate ZC and nH2O
    ZC <- c(ZC, mean(ZCAA(AAcomp)))
    ZC.SD <- c(ZC.SD, sd(ZCAA(AAcomp)))
    nH2O <- c(nH2O, mean(H2OAA(AAcomp)))
    nH2O.SD <- c(nH2O.SD, sd(H2OAA(AAcomp)))
  }
  # return values
  list(types = types, ZC = ZC, ZC.SD = ZC.SD, nH2O = nH2O, nH2O.SD = nH2O.SD)
}

