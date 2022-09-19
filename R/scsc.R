# JMDplots/scsc.R
# Two plots from the paper by Dick (2009)
# previously in CHNOSZ/demo/yeastgfp.R; added to JMDplots 20191019

# Localizations and abundances of proteins from YeastGFP are used
# to calculate an abundance-weighted average of amino acid compositions of proteins
# in different subcellular compartments of yeast.

scsc2 <- function() {
  ## Oxygen fugacity - activity of H2O predominance 
  ## diagrams for average protein composition in 23 YeastGFP localizations
  # use superseded properties of [Met], [Gly], and [UPBB] (Dick et al., 2006)
  if(packageVersion("CHNOSZ") > "1.4.3") OldAAfile <- "extdata/OBIGT/OldAA.csv" else "extdata/OBIGT/OldAA_old.csv"
  add.OBIGT(system.file(OldAAfile, package = "JMDplots"))
  # arranged by decreasing metastability:
  # order of this list of locations is based on the 
  # (dis)appearance of species on the current set of diagrams
  names <- c("vacuole", "early.Golgi", "ER", "lipid.particle",
    "cell.periphery", "ambiguous", "Golgi", "mitochondrion",
    "bud", "actin", "cytoplasm", "late.Golgi",
    "endosome", "nucleus", "vacuolar.membrane", "punctate.composite",
    "peroxisome", "ER.to.Golgi", "nucleolus", "spindle.pole",
    "nuclear.periphery", "bud.neck", "microtubule")
  nloc <- c(4, 5, 3, 4, 4, 3)
  # define the system
  basis("CHNOS+")
  # get protein names and abundances in each location
  gfp <- yeastgfp(names)
  # get amino acid compositions of proteins
  aa <- yeast.aa(gfp$protein)
  # calculate average amino acid compositions 
  for(i in 1:length(names)) {
    avgaa <- aasum(aa[[i]], gfp$abundance[[i]], average=TRUE, protein=names[i])
    add.protein(avgaa)
  }
  species(names, "Sce")
  a <- affinity(H2O=c(-5, 0, 256), O2=c(-80, -66, 256))
  # setup the plot
  opar <- par(no.readonly = TRUE)
  layout(matrix(c(1, 1,2:7), byrow=TRUE, nrow=4), heights=c(0.7, 3, 3, 3))
  par(mar=c(0, 0, 0, 0))
  plot.new()
  text(0.5, 0.7, expression("Proteins in subcellular locations of"~italic("S. cerevisiae")~"(Dick, 2009)"), cex=1.5)
  text(0.5, 0.2, describe.basis(ibasis=c(1, 3, 4, 6), oneline=TRUE), cex=1.5)
  par(mar=c(3, 4, 1, 1), xpd=TRUE)
  fill <- heat.colors(length(names))
  inames <- 1:length(names)
  for(i in 1:length(nloc)) {
    diagram(a, normalize=TRUE, names=names[inames], groups=as.list(inames),
      fill=fill[inames], cex.axis=0.75, cex.names=1.2, format.names = FALSE)
    label.plot(letters[i], xfrac=0.95, yfrac=0.9, paren=TRUE, italic=TRUE)
    title(main=paste(length(inames), "locations"))
    # take out the stable species
    inames <- inames[-(1:nloc[i])]
  }
  # return to plot defaults
  layout(matrix(1))
  par(xpd=FALSE)
  par(opar)
  # reset thermodynamic database
  reset()
}

scsc3 <- function() {
  ## This figure is similar to Fig. 3 of Dick (2009). 
  if(packageVersion("CHNOSZ") > "1.4.3") OldAAfile <- "extdata/OBIGT/OldAA.csv" else "extdata/OBIGT/OldAA_old.csv"
  add.OBIGT(system.file(OldAAfile, package = "JMDplots"))
  locations <- yeastgfp()
  gfp <- yeastgfp(locations)
  aa <- yeast.aa(gfp$protein)
  for(i in 1:length(locations)) {
    avgaa <- aasum(aa[[i]], gfp$abundance[[i]], average=TRUE, protein=locations[i])
    add.protein(avgaa)
  }
  basis("CHNOS+")
  species(locations, "Sce")
  a <- affinity(O2=c(-82, -65))
  e <- equilibrate(a, loga.balance=0, normalize=TRUE)
  mycolor <- topo.colors(length(locations))
  diagram(e, names=locations, ylim=c(-5, -3), col=mycolor, lwd=2, format.names = FALSE)
  title(main=expression("Proteins in subcellular locations of"~italic("S. cerevisiae")))

  # reset thermodynamic database
  reset()
}

###########################
### UNEXPORTED FUNCTION ###
###########################

# aasum - combine amino acid counts (sum, average, or weighted sum by abundance)
# Moved here from CHNOSZ 20220620
aasum <- function(aa, abundance=1, average=FALSE, protein=NULL, organism=NULL) {
  # returns the sum of the amino acid counts in aa,
  # multiplied by the abundances of the proteins
  abundance <- rep(abundance, length.out=nrow(aa))
  # drop any NA rows or abundances
  ina.aa <- is.na(aa$chains)
  ina.ab <- is.na(abundance)
  ina <- ina.aa | ina.ab
  if(any(ina)) {
    aa <- aa[!ina, ]
    abundance <- abundance[!ina]
    message("aasum: dropped ", sum(ina), " proteins with NA composition and/or abundance")
  }
  # multiply
  aa[, 6:25] <- aa[, 6:25] * abundance
  # sum
  out <- aa[1, ]
  out[, 5:25] <- colSums(aa[, 5:25])
  # average if told to do so
  if(average) {
    # polypeptide chains by number of proteins, residues by frequence
    out[, 5] <- out[, 5]/nrow(aa)
    out[, 6:25] <- out[, 6:25]/sum(abundance)
  }
  # add protein and organism names if given
  if(!is.null(protein)) out$protein <- protein
  if(!is.null(organism)) out$organism <- organism
  return(out)
}

