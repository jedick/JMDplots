# Calculate R-squared of nH2O-ZC and nO2-ZC fits for amino acids,
# using basis species that include different combindation of amino acids 20201011

# works with CHNOSZ 1.3.6
library(CHNOSZ)

# combine different amino acids
aa <- aminoacids("")
aacomb <- combn(20, 3)
# only keep combinations that have sulfur-bearing amino acids
icys <- which(aa=="cysteine")
imet <- which(aa=="methionine")
has.S <- apply(aacomb, 2, function(x) any(x %in% c(icys, imet)))
aacomb <- aacomb[, has.S, drop=FALSE]
# make abbreviations
aa1 <- aminoacids(1)
abbrv1 <- aa1[aacomb[1, ]]
abbrv2 <- aa1[aacomb[2, ]]
abbrv3 <- aa1[aacomb[3, ]]
abbrv <- paste0(abbrv1, abbrv2, abbrv3)
# where to keep the results
O2.R2 <- H2O.R2 <- numeric()
O2.slope <- H2O.slope <- numeric()
# calculate ZC
ZC <- ZC(info(aminoacids(""), "aq"))
# loop over combinations
for(i in 1:ncol(aacomb)) {
  thisbasis <- c(aa[aacomb[, i]], "H2O", "oxygen")
  b <- try(basis(thisbasis), silent = TRUE)
  if(inherits(b, "try-error")) {
    O2.R2 <- c(O2.R2, NA)
    H2O.R2 <- c(H2O.R2, NA)
    O2.slope <- c(O2.slope, NA)
    H2O.slope <- c(H2O.slope, NA)
  } else {
    species(aminoacids(""))
    nO2 <- species()$O2
    nH2O <- species()$H2O
    # calculate R2 and slope for nO2-ZC fit
    O2lm <- lm(nO2 ~ ZC)
    O2.R2 <- c(O2.R2, summary(O2lm)$r.squared)
    O2.slope <- c(O2.slope, coefficients(O2lm)[[2]])
    # calculate R2 and slope for nH2O-ZC fit
    H2Olm <- lm(nH2O ~ ZC)
    H2O.R2 <- c(H2O.R2, summary(H2Olm)$r.squared)
    H2O.slope <- c(H2O.slope, coefficients(H2Olm)[[2]])
  }
}
out <- data.frame(abbrv, O2.R2, H2O.R2, O2.slope, H2O.slope)
write.csv(out, "AAbasis.csv", row.names = FALSE, quote = FALSE)
