# gradox/plot.R
# make plots for paper
# Changes in carbon oxidation state of metagenomes along geochemical redox gradients
# 20180428 initial version
# 20180831 updated for paper submission
# 20181215 updated for paper revision
# 20190928 first addition to JMDplots package

# USAGE (R commands):
# source("plot.R")        - source this file
# Fig1(TRUE)              - make Figure_1.pdf
# Fig2(TRUE)              - make Figure_2.pdf
# mout <- FigS1(TRUE)     - make Figure_S1.pdf
# pout <- FigS2(TRUE)     - make Figure_S2.pdf
# Fig3(mout, pout, TRUE)  - make Figure_3.pdf
# Fig4(mout, pout, TRUE)  - make Figure_4.pdf
# Fig5(TRUE)              - make Figure_5.pdf
# Fig5(TRUE, NULL)        - make Figure_S3.pdf

# NOTE: Fig5() generates some expected messages about
# missing files, indicating samples where data are not available.

## requires CHNOSZ 1.1.3 and R 3.6.0 (for 'gap.axis' argument in axis())
#library(CHNOSZ)
#data(thermo)
#
## code for mpage/ppage/pcomp/mplot (Figures 2, 3, 4, 5)
#source("MG.R")

# general characteristics of ZC of DNA, RNA and proteins
gradox1 <- function(pdf=FALSE) {
  if(pdf) pdf("gradox1.pdf", width=8, height=6)
  par(mfrow=c(2, 2), mgp=c(2.5, 1, 0), las=1)
  ZC_base()
  label.figure("A", xfrac=0.035, yfrac=0.93, cex=1.6, font=2)
  ZC_GC()
  label.figure("B", xfrac=0.035, yfrac=0.93, cex=1.6, font=2)
  AA_codon(146891)
  label.figure("C", xfrac=0.035, yfrac=0.93, cex=1.6, font=2)
  AA_codon(300852)
  label.figure("D", xfrac=0.035, yfrac=0.93, cex=1.6, font=2)
  if(pdf) invisible(dev.off())
}

# selected plots of DNA, RNA, and protein ZC in one figure 20181113
gradox2 <- function(pdf=FALSE) {
  if(pdf) pdf("gradox2.pdf", width=9, height=9)
  mat <- matrix(c(0, 1, 2, 3, 4, 5, 10, 20, 11, 21, 6, 12, 22, 13, 23, 7, 14, 24, 15, 25, 8, 16, 26, 17, 27, 9, 18, 28, 19, 29), byrow=TRUE, nrow=6)
  layout(mat, widths=c(0.6, 3, 3, 3, 3), heights=c(0.6, 3, 3, 3, 3, 3))
  # add column and row titles
  par(mar=c(0, 0, 0, 0))
  plot.new()
  text(0.5, 0.5, "     DNA and RNA", cex=1.5, font=2)
  plot.new()
  text(0.5, 0.5, "     Protein", cex=1.5, font=2)
  plot.new()
  text(0.5, 0.5, "     DNA and RNA", cex=1.5, font=2)
  plot.new()
  text(0.5, 0.5, "     Protein", cex=1.5, font=2)
  plot.new()
  text(0.5, 0.5, "    sediment", srt=90, cex=1.5, font=2)
  plot.new()
  text(0.5, 0.5, "      hydrothermal vent", srt=90, cex=1.5, font=2)
  plot.new()
  text(0.5, 0.5, "     ocean", srt=90, cex=1.5, font=2)
  plot.new()
  text(0.5, 0.5, "     hypersaline", srt=90, cex=1.5, font=2)
  plot.new()
  text(0.5, 0.5, "     microbial mat", srt=90, cex=1.5, font=2)
  # plot the data
  par(mar=c(4, 3.5, 2, 1), mgp=c(2.5, 1, 0))
  mpage("gradoxMS", set.par=FALSE)
  ppage("gradoxMS", set.par=FALSE)
  if(pdf) invisible(dev.off())
}

# all plots of DNA and RNA ZC on a single page
gradoxS1 <- function(pdf=FALSE) {
  if(pdf) pdf("gradoxS1.pdf", width=6.5, height=10)
  mout <- mpage()
  if(pdf) dev.off()
  invisible(mout)
}

# all plots of protein ZC on a single page
gradoxS2 <- function(pdf=FALSE) {
  if(pdf) pdf("gradoxS2.pdf", width=6.5, height=10)
  pout <- ppage()
  if(pdf) dev.off()
  invisible(pout)
}

# ZC of proteins vs DNA (metagenomes and metatranscriptomes)
gradox3 <- function(mout, pout, pdf=FALSE) {
  if(pdf) pdf("gradox3.pdf", width=8, height=5)
  # setup plot
  mat <- matrix(1, nrow=5, ncol=8)
  mat[1:2, 6:8] <- 2
  mat[3:5, 6:8] <- 3
  layout(mat)
  par(mar=c(4, 5, 2, 1), cex=0.8, las=1)
  # metagenomes
  pcomp(mout, pout, parts="plot")
  title(main="Metagenomes")
  label.figure("A", font=2, cex=1.7, yfrac=0.97)
  # legend
  opar <- par(mar=c(1, 0.2, 1, 0.2))
  plot.new()
  pcomp(parts="legend")
  par(opar)
  # metatranscriptomes
  pcomp(mout, pout, "MT", parts="plot")
  title(main="Metatranscriptomes")
  label.figure("B", font=2, cex=1.7, yfrac=0.955)
  if(pdf) invisible(dev.off())
}

# thermodynamic calculations of relative stabilities along redox gradients
gradox4 <- function(mout, pout, pdf=FALSE) {
  if(pdf) pdf("gradox4.pdf", width=7.5, height=6.5)
  mat <- matrix(0, nrow=4, ncol=15)
  mat[1:2, 1:7] <- 1
  mat[1, 6:8] <- 2
  mat[1, 9:11] <- 3
  mat[2, 6:8] <- 4
  mat[2, 9:11] <- 5
  mat[3, 1:5] <- 6
  mat[3, 6:10] <- 7
  mat[3, 11:15] <- 8
  mat[4, 1:5] <- 9
  mat[4, 6:10] <- 10
  mat[4, 11:15] <- 11
  mat[2, 12:15] <- 12
  layout(mat, heights=c(0.5, 0.5, 0.8, 0.8))
  # add descriptive text in first panel
  plot.new()
  opar <- par(xpd=NA)
  # y-offset is sensitive to letters with dangling parts
  dy <- 0.13; ddy <- 0.015
  yA <- 1.315
  text(-0.28, yA, "A", font=2, cex=1.5)
  text(-0.2, yA-ddy, "Calculate affinity per", adj=0, cex=1.5)
  text(-0.2, yA-ddy-dy, "monomer for overall DNA", adj=0, cex=1.5)
  text(-0.2, yA-ddy-dy-ddy-dy, "and protein composition", adj=0, cex=1.5)
  text(-0.2, yA-ddy-dy-ddy-dy-dy, "in each sample.", adj=0, cex=1.5)
  yB <- 0.315
  text(-0.28, yB, "B", font=2, cex=1.5)
  text(-0.2, yB-ddy, "Subtract sample means", adj=0, cex=1.5)
  text(-0.2, yB-ddy-dy-ddy, quote("to get relative affinity ("*Delta*"A)"*"."), adj=0, cex=1.5)
  text(-0.28, -0.385, "C", font=2, cex=1.5)
  text(-0.2, -0.4, "Compare relative affinity per", adj=0, cex=1.5)
  text(-0.2, -0.53, "monomer in DNA and protein.", adj=0, cex=1.5)
  lines(c(-0.28, -0.28), c(-0.5, -0.85), lwd=2)
  par(opar)
  # plot affinities and relative affinities for proteins and DNA Bison Pool
  aAA_aDNA(mout, pout, datasets="Bison_Pool_IMG_MG", do.mfrow=FALSE)
  # plot relative affinities for proteins vs DNA for multiple datasets
  par(mar=c(3, 3.5, 2.5, 1), mgp=c(1.7, 0.3, 0), xaxs="r", yaxs="r")
  aAA_aDNA(mout, pout, do.mfrow=FALSE)
  # add key for line colors/symbols
  par(mar=c(1.8, 2, 1.2, 1))
  opar <- par(xpd=NA)
  plot.new()
  lines(c(-0.1, 0.4), c(0.9, 0.9), lwd=2, col="blue")
  lines(c(-0.1, 0.4), c(0.7, 0.7), lwd=1, col="blue")
  points(c(0.4, 0.4), c(0.9, 0.7), pch=19, col="blue", cex=2*par("cex"))
  lines(c(-0.1, 0.4), c(0.5, 0.5), lwd=1, col="black")
  lines(c(-0.1, 0.4), c(0.3, 0.3), lwd=1, col="red")
  lines(c(-0.1, 0.4), c(0.1, 0.1), lwd=2, col="red")
  points(c(-0.1, -0.1), c(0.3, 0.1), pch=19, col="red", cex=2*par("cex"))
  text(0.5, 0.5, "ENVIRONMENT", font=3, adj=0)
  text(0.5, 0.9, "oxidizing", adj=0)
  text(0.5, 0.1, "reducing", adj=0)
  text(0.15, 1.3, "MODEL", font=3)
  text(-0.1, 1.1, "reducing")
  text(0.4, 1.1, "oxidizing")
  par(opar)
  if(pdf) invisible(dev.off())
}

# ZC of reads classified to selected abundant species in different datasets 20180529
# revised species plots 20181120
gradox5 <- function(pdf=FALSE, maxdepth=500) {
  #  1 1420917 Marinobacter salarius -- blue
  #  2 28108   Alteromonas macleodii
  #  3 1427364 Candidatus Thioglobus singularis -- red
  #  4 563040  Sulfurimonas autotrophica DSM 16294
  #  5 1977865 Candidatus Pelagibacter sp. RS40 -- orange
  #  6 1410606 Candidatus Nitrosopelagicus brevis -- green
  #  7 1898749 Candidatus Nitrosomarinus catalina -- purple
  #  8 680     Vibrio campbellii
  #  9 146891  Prochlorococcus marinus str. AS9601
  # 10 329     Ralstonia pickettii
  # 11 1747    Cutibacterium acnes
  # 12 374840  Enterobacteria phage phiX174 sensu lato

  # set up figure
  if(pdf) {
    if(!is.null(maxdepth)) pdf("gradox5.pdf", width=8, height=5)
    else pdf("gradoxS3.pdf", width=8, height=5)
  }
  #layout(matrix(c(1, 2, 7, 3, 4, 7, 5, 6, 7), nrow=3, byrow=TRUE))
  layout(matrix(c(1, 2, 6, 3, 4, 5), nrow=2, byrow=TRUE))
  par(mar=c(4, 3.5, 2, 1), mgp=c(2.5, 1, 0))
  labfig <- function(x) label.figure(x, xfrac=0.035, yfrac=0.95, cex=1.6, font=2)

  # now make each plot
  mplot("Diffuse_Vents", "SRA_MG", taxid=c(0, 1420917, 28108, 1427364, 563040, 1977865))
  text(c(2.7, 7.3, 0.7, 0.7, 2.75), c(0.633, 0.604, 0.591, 0.585, 0.58), c(1, 2, 3, 4, 5))
  labfig("A")
  mplot("Menez_Gwen", "SRA_MG", taxid=c(0, 1420917, 28108, 1427364, 1410606, 1898749))
  text(c(-2, -2, -2, -2), c(0.6315, 0.605, 0.596, 0.592, 0.588), c(1, 2, 3, 6, 7))
  labfig("B")
  mplot("ETNP_OMZ", "SRA_MG", taxid=c(0, 28108, 1427364, 1977865))
  text(c(315, 315, 315, 315), c(0.6035, 0.5985, 0.5905), c(2, 3, 5))
  # add line for mean ZC of T. singularis in vents (precomputed value)  20180601
  abline(h=0.5955, lty=3, col="red")
  labfig("C")
  mplot("ETSP_OMZ", "SRA_MG", taxid=c(0, 1427364, 1977865, 1410606, 1898749), maxdepth=maxdepth)
  if(is.null(maxdepth)) xpos <- c(535, 535, 535, 825) else xpos <- rep(525, 4)
  text(xpos, c(0.5935, 0.58, 0.59, 0.577), c(3, 5, 6, 7))
  abline(h=0.5955, lty=3, col="red")
  labfig("D")
  mplot("HOT_ALOHA-2010", "SRA_MG", taxid=c(0, 1427364, 1977865, 1410606, 1898749, 680, 146891), maxdepth=maxdepth)
  if(is.null(maxdepth)) xpos <- rep(1050, 4) else xpos <- rep(520, 4)
  text(xpos, c(0.592, 0.576, 0.588, 0.585, 0.605, 0.571), c(3, 5, 6, 7, 8, 9))
  abline(h=0.5955, lty=3, col="red")
  labfig("E")

  # set up legend
  par(mar=c(0, 0.5, 2, 0), xpd=NA)
  get.legend <- function(abbrvs) {
    # legend text for species
    Ms <- quote(italic("Marinobacter salarius"))
    Am <- quote(italic("Alteromonas macleodii"))
    Ts <- quote(italic("Ca.")~"Thioglobus singularis")
    Sa <- quote(italic("Sulfurimonas autotrophica")~"DSM 16294")
    Ps <- quote(italic("Ca.")~"Pelagibacter sp. RS40")
    Nb <- quote(italic("Ca.")~"Nitrosopelagicus brevis")
    Nc <- quote(italic("Ca.")~"Nitrosomarinus catalina")
    Vc <- quote(italic("Vibrio campbellii"))
    Pm <- quote(italic("Prochlorococcus marinus")~"str. AS9601")
    EM <- quote("Entire metagenome")
    # lines
    lty <- rep(2, length(abbrvs))
    lty[abbrvs=="EM"] <- 1
    lwd <- rep(1.4, length(abbrvs))
    lwd[abbrvs=="EM"] <- 1.9
    # colors
    col <- rep("black", length(abbrvs))
    col[abbrvs=="Ts"] <- "red"
    col[abbrvs=="Ms"] <- "blue"
    col[abbrvs=="Ps"] <- "darkorange"
    col[abbrvs=="Nc"] <- "purple"
    col[abbrvs=="Nb"] <- "darkgreen"
    # get the expressions from the quoted abbreviations
    # use envir= to make this work in the package 20190929
    # https://stackoverflow.com/questions/8016636/using-get-inside-lapply-inside-a-function
    legend <- sapply(abbrvs, get, envir=sys.frame(sys.parent(0)))
    return(list(legend=legend, lty=lty, lwd=lwd, col=col))
  }
  # make legend
  plot.new()
  gl <- get.legend(c("Ms", "Am", "Ts", "Sa", "Ps", "Nb", "Nc", "Vc", "Pm", "EM"))
  leg <- legend("topleft", as.expression(gl$legend), lty=gl$lty, lwd=gl$lwd, col=gl$col, cex=1.1, bty="n")
  text(-0.03, leg$text$y, c(1:9, NA))

  # done!
  if(pdf) invisible(dev.off())
}

############################
### UNEXPORTED FUNCTIONS ###
############################

### next three functions are used to make the different parts of Figure 1 ###

# plot ZC vs nC for nucleobases and nucleosides 20180428
ZC_base <- function(pdf=FALSE) {
  # get ZC for the different compounds
  bases <- c("adenine", "thymine", "uracil", "guanine", "cytosine")
  ZC_bases <- sapply(makeup(info(bases)), ZC)
  RNAsides <- c("adenosine", "thymidine", "uridine", "guanosine", "cytidine")
  ZC_RNAsides <- sapply(makeup(info(RNAsides)), "ZC")
  DNAsides <- c("deoxyadenosine", "deoxythymidine", "deoxyuridine", "deoxyguanosine", "deoxycytidine")
  ZC_DNAsides <- sapply(makeup(info(DNAsides)), "ZC")
  sugars <- c("ribose", "deoxyribose")
  ZC_sugars <- sapply(makeup(info(sugars)), "ZC")
  # set up plot
  if(pdf) pdf("ZC_base.pdf", width=5, height=5)
  par(mar=c(4, 4, 1, 1))
  plot(c(1, 5), c(-0.5, 3.5), xlab="", ylab=NA, type="n", xaxt="n")
  axis(1, at=1:5, labels=c("A", "T", "U", "G", "C"))
  mtext(quote(italic(Z)[C]), side=2, line=2.5, las=0, cex=par("cex"))
  # nucleobases
  points(1:5, ZC_bases, pch=19, cex=1.5)
  # nucleosides: don't plot T for RNA or U for DNA
  iRNA <- c(1, 3:5)
  points(iRNA, ZC_RNAsides[iRNA], pch=15, col="blue", cex=1.5)
  # show A-T and G-C basepairs in DNA
  lines(1:2, ZC_DNAsides[1:2], lty=2)
  lines(4:5, ZC_DNAsides[4:5], lty=2)
  iDNA <- c(1:2, 4:5)
  points(iDNA, ZC_DNAsides[iDNA], pch=17, col="red", cex=1.5)
  # ribose and deoxyribose lines
  abline(h=ZC_sugars, lty=3, lwd=2)
  text(3, ZC_sugars[1]+0.13, "ribose")
  text(3, ZC_sugars[2]+0.11, "deoxyribose")
  legend("topleft", legend=c("nucleobase", "base + ribose (in RNA)", "base + deoxyribose (in DNA)"),
         pch=c(19, 15, 17), col=c("black", "blue", "red"), pt.cex=1.5)
  if(pdf) dev.off()
}

# plot ZC vs G+C content for dsDNA and ssDNA/ssRNA 20180428
ZC_GC <- function(pdf=FALSE) {
  # setup plot
  if(pdf) pdf("ZC_GC.pdf", width=5, height=5)
  par(mar=c(4, 4, 1, 1))
  plot(c(0, 100), c(0.4, 1), xlab="GC content (%)", ylab=NA, type="n")
  mtext(quote(italic(Z)[C]), side=2, line=2.5, las=0, cex=par("cex"))
  # dsDNA
  ZC1 <- ZC(makeup(info(c("deoxyguanosine", "deoxycytidine")), sum=TRUE))
  ZC0 <- ZC(makeup(info(c("deoxyadenosine", "deoxythymidine")), sum=TRUE))
  lines(c(0, 100), c(ZC0, ZC1), col="red")
  # change ribose to deoxyribose
  ZC1 <- ZC(makeup(info(c("guanosine", "cytidine")), sum=TRUE))
  ZC0 <- ZC(makeup(info(c("adenosine", "thymidine")), sum=TRUE))
  lines(c(0, 100), c(ZC0, ZC1), col="blue", lty=2)
  # ssRNA with A=U and G=C
  ZC1 <- ZC(makeup(info(c("guanosine", "cytidine")), sum=TRUE))
  ZC0 <- ZC(makeup(info(c("adenosine", "uridine")), sum=TRUE))
  lines(c(0, 100), c(ZC0, ZC1), col="blue")
  # add legend
  legend("bottomright", legend=c("RNA", "RNA (T instead of U)", "dsDNA"), col=c("blue", "blue", "red"), lty=c(1, 2, 1), lwd=1.5)
  if(pdf) dev.off()
}

# plot ZC of amino acids and DNA codons 20180428
AA_codon <- function(organism=c(146891, 300852), pdf=FALSE) {
  # read aminoacid - codon file
  file <- system.file("extdata/gradox/AA_codon.csv", package = "JMDplots")
  dat <- read.csv(file, as.is=TRUE)
  # remove "End" codons
  dat <- dat[dat$AmAcid!="End", ]
  # get ZC of amino acids
  AAabbrv <- aminoacids(3)
  AAnames <- aminoacids("")
  AA <- AAnames[match(dat$AmAcid, AAabbrv)]
  ZC_AA <- ZC(info(AA))
  # get ZC of codons
  ZC_codon <- numeric()
  for(i in 1:nrow(dat)) {
    nucleosides <- as.character(dat[i, 3:5])
    ZC_codon <- c(ZC_codon, ZC(makeup(info(nucleosides), sum=TRUE)))
  }
  # get point sizes based on frequency (per 1000) of codons in selected organism
  freqs <- dat[, paste0("permil_", organism[1])]
  cex <- sqrt(freqs)/3
  # setup plot
  if(pdf) pdf(paste0("AA_codon_", organism, ".pdf"), width=5, height=5)
  par(mar=c(4, 4, 1, 1))
  plot(c(0.2, 1), c(-1, 1), xlab=quote(DNA~codon~~italic(Z)[C]), ylab=NA, type="n")
  mtext(quote(amino~acid~~italic(Z)[C]), side=2, line=2.5, las=0, cex=par("cex"))
  points(ZC_codon, ZC_AA, pch=19, cex=cex)
  # add weighted regression line
  ZC_lm <- lm(ZC_AA ~ ZC_codon, weights=freqs)
  ZC_codon <- c(0.2, 1.0)
  ZC_AA <- predict(ZC_lm, data.frame(ZC_codon=ZC_codon))
  lines(ZC_codon, ZC_AA, lty=2)
  R2 <- round(summary(ZC_lm)$r.squared, 2)
  legend(0.12, 0.4, legend=substitute(italic(R)^2 == R2, list(R2=R2)), bty="n")
  # add organism name
  if(organism[1]==146891) legend <- c("Prochlorococcus", "  marinus")
  if(organism[1]==300852) legend <- "Thermus thermophilus"
  legend(0.11, 1.15, legend=legend, text.font=3, bty="n")
  if(pdf) dev.off()
}

### next two functions are used to make Figure 4 ###

# protein and DNA affinities 20180505
aAA_aDNA <- function(mout, pout, datasets=c("Bison_Pool_IMG_MG", "Mud_Volcano_SRA_MG", "ETNP_OMZ_SRA_MG",
                                            "BalticSea_Sediment_SRA_MT", "Diffuse_Vents_SRA_MT", "ETNP_OMZ_SRA_MT"),
  do.mfrow=TRUE, do.legend=TRUE, do.abbrev=TRUE, mar1=c(1.8, 2, 1.2, 1), mar2=c(2, 2.5, 1, 1), dataset_in_title=TRUE) {
  # set up figure
  if(is.null(datasets)) {
    datasets <- names(mout)
    par(mfrow=c(4, 5), mar=c(3, 3, 1, 1))
  } else if(length(datasets)==1) {
    if(do.mfrow) par(mfrow=c(2, 2))
  } else {
    if(do.mfrow) par(mfrow=c(2, 3))
  }
  # set up basis species
  #basis("CHNOPS+") # not available in CHNOSZ_1.1.3 (only in later versions)
  basis(c("CO2", "H2O", "NH3", "H3PO4", "H2S", "e-", "H+"))
  basis("H+", -7)
  # use ionized species
  swap.basis("CO2", "HCO3-")
  basis("HCO3-", -3)
  swap.basis("NH3", "NH4+")
  basis("NH4+", -7)
  swap.basis("H2S", "HS-")
  basis("HS-", -9)
  swap.basis("H3PO4", "H2PO4-")
  basis("H2PO4-", -5)
  for(dataset in datasets) {
    # get data
    moutdata <- mout[[dataset]]
    # calculate the thermodynamic properties of the DNA molecule (per nucleotide)
    DNA_obigt <- seqcomp2obigt(moutdata$meancomp, "DNA")
    # add them to the "OBIGT" database in CHNOSZ
    iDNA <- do.call(mod.obigt, DNA_obigt)
    # add species
    species(iDNA)
    # calculate affinity
    aDNA <- affinity(Eh=c(-0.35, -0.15))
    # convert to kJ (2.303RT = 5.71 kJ at 25 degC)
    aDNA$values <- lapply(aDNA$values, "*", 5.71)
    # convert nucleotides to basepairs (multiply by 2)
    aDNA$values <- lapply(aDNA$values, "*", 2)
    # clean up for proteins
    species(delete=TRUE)

    # now do it for proteins
    poutdata <- pout[[paste0(dataset, "P")]]
    AA_obigt <- seqcomp2obigt(poutdata$meancomp, "protein")
    iAA <- do.call(mod.obigt, AA_obigt)
    species(iAA)
    aAA <- affinity(Eh=c(-0.35, -0.15))
    aAA$values <- lapply(aAA$values, "*", 5.71)
    # clean up for next dataset
    species(delete=TRUE)

    # subtract the means
    DNA_affinity <- lapply(aDNA$values, "-", rowMeans(as.data.frame(aDNA$values)))
    AA_affinity <- lapply(aAA$values, "-", rowMeans(as.data.frame(aAA$values)))

    # line colors: 2 red (most reducing samples), black (intermediate), 2 blue (most oxidizing samples)
    n <- length(DNA_affinity)
    col <- rep("black", n)
    col[1:2] <- "red"
    col[c(n-1, n)] <- "blue"
    # line width: 1 bold (most reducing sample), black (intermediate), 1 bold (most oxidizing sample)
    lwd <- rep(1, n)
    lwd[1] <- 2
    lwd[n] <- 2
    # axis labels
    AAlab <- quote(italic(A)~"/ amino acid (kJ)")
    DNAlab <- quote(italic(A)~"/ base pair (kJ)")
    DAAlab <- quote(Delta*italic(A)~"/ amino acid (kJ)")
    DDNAlab <- quote(Delta*italic(A)~"/ base pair (kJ)")
    Ehlab <- "Eh (volt)"

    if(length(datasets) > 1) {
      # set up plot
      xlim <- range(unlist(DNA_affinity))
      ylim <- range(unlist(AA_affinity))
      # extend range (legend obscures lines)
      if(grepl("ETNP_OMZ_SRA_MG", dataset)) xlim[1] <- xlim[1] - 1
      if(grepl("ETNP_OMZ_SRA_MT", dataset)) ylim[2] <- ylim[2] + 0.25
      if(grepl("ETNP_OMZ_SRA_MT", dataset)) xlim[1] <- xlim[1] - 0.2
      if(grepl("Diffuse_Vents", dataset)) ylim[2] <- ylim[2] + 2
      if(grepl("Diffuse_Vents", dataset)) xlim[1] <- xlim[1] - 1
      plot(xlim, ylim, xlab=DDNAlab, ylab=DAAlab, type="n")
      # shade area where affinities of DNA or protein are negative
      rect(0, par("usr")[3], par("usr")[2], 0, col="grey80", lty=0)
      rect(par("usr")[1], 0, 0, par("usr")[4], col="grey80", lty=0)
      rect(par("usr")[1], par("usr")[3], 0, 0, col="grey80", lty=0)
      abline(h=0, v=0, lty=3)
      # redraw plot box and axis tick marks (obscured by rectangles)
      box()
      axis(1, labels=FALSE)
      axis(2, labels=FALSE)
      axis(3, labels=FALSE)
      axis(4, labels=FALSE)
      # draw lines
      for(i in 1:length(DNA_affinity)) lines(DNA_affinity[[i]], AA_affinity[[i]], col=col[i], lwd=lwd[i])
      # draw filled circles at lowest / highest logaH2 (relatively oxidizing / reducing samples)
      for(i in 1:length(DNA_affinity)) {
        if(col[i]=="red") ival <- 1
        if(col[i]=="blue") ival <- length(DNA_affinity[[i]])
        if(col[i]=="black") next  # skip drawing points for intermediate samples
        points(DNA_affinity[[i]][ival], AA_affinity[[i]][ival], pch=19, col=col[i], lwd=lwd[i]*par("cex"), cex=2*par("cex"))
      }
      # add legend
      if(do.legend) {
        ltxt <- rownames(moutdata$meancomp)
        if(dataset=="Bison_Pool_IMG_MG") ltxt <- paste0(studies[["Bison_Pool"]]$xlabels, "m")
        pos <- "topleft"
        if(dataset=="ETNP_OMZ_SRA_MG") pos <- "left"
        legend(pos, legend=rev(ltxt), lwd=rev(lwd), col=rev(col), bg="white", cex=0.8)
      }
      # add title
      main <- dataset2main(dataset, moutdata$abbrev)
      title(main=main, font.main=1)
      # add in-plot label: MG or MT 20180829
      label <- "MG"
      if(grepl("_MT.*", dataset)) label <- "MT"
      # change position for some datasets
      xfrac <- 0.9
      if(grepl("OMZ", dataset)) xfrac <- 0.1
      CHNOSZ::label.plot(label, xfrac=xfrac, yfrac=0.1)
    } else {
      # for a single dataset (i.e. Bison Pool) plot the affinities and relative affinities for protein and DNA
      diagram(aDNA, balance=1, col=col, lwd=lwd, lty=1, names=NULL, mar=mar1, mgp=c(1.3, 0.2, 0), xlab=Ehlab, ylab=DNAlab, yline=2.4)
      if(dataset_in_title) title(main=paste0("DNA (", dataset2main(dataset), ")"), font.main=1)
      else title(main="DNA", font.main=1)
      diagram(aAA, balance=1, col=col, lwd=lwd, lty=1, names=NULL, mar=mar1, mgp=c(1.3, 0.2, 0), xlab=Ehlab, ylab=AAlab, yline=1.9)
      if(dataset_in_title) title(main=paste0("Protein (", dataset2main(dataset), ")"), font.main=1)
      else title(main="protein", font.main=1)
      # relative affinities
      aDNA$values <- DNA_affinity
      aAA$values <- AA_affinity
      diagram(aDNA, balance=1, col=col, lwd=lwd, lty=1, names=NULL, mar=mar2, mgp=c(1.3, 0.2, 0), xlab=Ehlab, ylab=DDNAlab, yline=1.9)
      # shade area where affinity is negative, and re-draw parts of the diagram obscured by the rectangle
      rect(par("usr")[1], par("usr")[3], par("usr")[2], 0, col="grey80", lty=0)
      abline(h=0, lty=3)
      diagram(aDNA, balance=1, col=col, lwd=lwd, lty=1, names=NULL, mar=mar2, mgp=c(1.3, 0.2, 0), xlab="", ylab="", add=TRUE)
      thermo.axis()
      # draw filled circles at lowest / highest logaH2 (relatively oxidizing / reducing samples)
      opar <- par(xpd=TRUE)
      for(i in 1:length(aDNA$values)) {
        if(col[i]=="red") ival <- 1
        if(col[i]=="blue") ival <- length(aDNA$vals[[1]])
        if(col[i]=="black") next  # skip drawing points for intermediate samples
        points(aDNA$vals[[1]][ival], aDNA$values[[i]][ival], pch=19, col=col[i], lwd=lwd[i]*par("cex"), cex=2*par("cex"))
      }
      par(opar)
      diagram(aAA, balance=1, col=col, lwd=lwd, lty=1, names=NULL, mar=mar2, mgp=c(1.3, 0.2, 0), xlab=Ehlab, ylab=DAAlab, yline=1.9)
      rect(par("usr")[1], par("usr")[3], par("usr")[2], 0, col="grey80", lty=0)
      abline(h=0, lty=3)
      diagram(aAA, balance=1, col=col, lwd=lwd, lty=1, names=NULL, mar=mar2, mgp=c(1.3, 0.2, 0), xlab="", ylab="", add=TRUE)
      thermo.axis()
      opar <- par(xpd=TRUE)
      for(i in 1:length(aAA$values)) {
        if(col[i]=="red") ival <- 1
        if(col[i]=="blue") ival <- length(aAA$vals[[1]])
        if(col[i]=="black") next  # skip drawing points for intermediate samples
        points(aAA$vals[[1]][ival], aAA$values[[i]][ival], pch=19, col=col[i], lwd=lwd[i]*par("cex"), cex=2*par("cex"))
      }
      par(opar)
    }
  }
}

# convert data frame of DNA or AA composition to thermodynamic properties 20180505
seqcomp2obigt <- function(seqcomp, type="DNA") {
  if(type=="DNA") {
    # which columns have A, C, G, T
    icol <- match(c("A", "C", "G", "T"), colnames(seqcomp))
    # thermodynamic data for the corresponding nucleotide monophosphates
    iobigt <- info(c("dAMP-2", "dCMP-2", "dGMP-2", "dTMP-2"))
  }
  if(type=="protein") {
    # which columns have the amino acid 3-letter abbreviations
    icol <- match(aminoacids(3), colnames(seqcomp))
    # thermodynamic data for the corresponding amino acids
    iobigt <- info(aminoacids(""))
  }
  # get monomer properties (nucleotides or amino acids)
  # to keep the multipliers, don't use info(iobigt)
  monomer_obigt <- thermo()$obigt[iobigt, ]
  # initialize output
  obigt_out <- monomer_obigt[rep(1, nrow(seqcomp)), ]
  # in CHNOSZ_1.1.3, diagram() (or another function) attempts to interpret names with underscores
  # as proteins, ultimately leading to incorrect diagrams (i.e. only single fields appear)
  obigt_out$name <- gsub("_", "-", rownames(seqcomp))
  obigt_out$date <- today()
  obigt_out$abbrv <- obigt_out$ref1 <- NA
  # loop over sequence compositions (rows)
  for(i in 1:nrow(seqcomp)) {
    # normalize composition to a single monomer unit
    monocomp <- as.numeric(seqcomp[i, icol]) / sum(seqcomp[i, icol])
    # chemical formula
    mkp <- makeup(monomer_obigt$formula, multiplier=monocomp, sum=TRUE)
    obigt_out$formula[i] <- as.chemical.formula(mkp)
    # thermodynamic properties (G, H, S, and HKF parameters)
    obigt_out[i, 8:20] <- colSums(monomer_obigt[, 9:21] * monocomp)
  }
  obigt_out
}

