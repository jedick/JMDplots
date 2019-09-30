# JMDplots/grad_util.R
# utility functions for gradox and gradH2O papers
# these are exported (in NAMESPACE) to also be available for user scripts
# mplot extracted from MG.R 20190930

# function to plot sampled compositions for indicated study and sequence type 20180222
# e.g. mplot("Columbia_River", "IMG_MT")
mplot <- function(study, seqtype, plottype="bars", ylim=NULL, plot.RNA=TRUE, taxid=NULL, dsDNA=TRUE, abbrev=NULL, col=NULL, add.label=TRUE, maxdepth=NULL) {
  md <- mdata(study, seqtype)
  # get labels, groups, and abbreviation
  samples <- md$samples
  xlabels <- md$xlabels
  group <- md$group
  if(missing(abbrev)) abbrev <- md$abbrev
  xlab <- names(studies[[study]][1])
  # get ZC range
  if(is.null(ylim)) {
    if(grepl("_MGP", seqtype)) ylim <- studies[[study]][["MGP_range"]]
    else if(grepl("_MTP", seqtype)) ylim <- studies[[study]][["MTP_range"]]
    else if(grepl("_MG", seqtype) & is.null(taxid)) ylim <- studies[[study]][["MG_range"]]
    else if(grepl("_MG", seqtype) & !is.null(taxid)) ylim <- studies[[study]][["MG_srange"]]
    else if(grepl("_MT", seqtype) & is.null(taxid)) ylim <- studies[[study]][["MT_range"]]
    else if(grepl("_MT", seqtype) & !is.null(taxid)) ylim <- studies[[study]][["MT_srange"]]
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
             plot.RNA=FALSE, taxid=taxid[i], lwd=lwd, lty=lty, lwd.bars=lwd.bars, col=col, extendrange=TRUE,
             add.label=add.label, plot_real_x=studies[[study]][["plot_real_x"]], maxdepth=maxdepth)
    }
  } 
  else plotMG(paste0(study, "_", seqtype), plottype, samples, xlabels, group, xlab, ylim, abbrev, dsDNA,
              plot.RNA, taxid, col=col, add.label=add.label, plot_real_x=studies[[study]][["plot_real_x"]], maxdepth=maxdepth)
}

