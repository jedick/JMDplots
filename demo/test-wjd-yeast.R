## Test that open-system equilibrium distributions reproduce the results of wjd()
## (uses yeast proteins)
# present in CHNOSZ/test/test-wjd.Rd before 20130110
# moved from CHNOSZ to JMDplots 20191019

### FIXME: equil.potentials(w) returns NULL unless we use group additivity parameters from DLH06 20190206
### (issue appears with MKL and OpenBLAS CRAN checks)
reset()
add.obigt("OldAA")
### set up system
# use proteins in the lipid particle (n=19)
y <- yeastgfp("lipid.particle")
# get the amino acid compositions of the proteins
aa <- yeast.aa(y$protein)
# don't use those with NA abundance or sequence (leaves n=17)
ina <- is.na(y$abundance) | is.na(aa$chains)
aa <- aa[!ina, ]
# normalize the proteins to single residues; columns 6:25 are the amino acid counts
aa.625 <- aa[, 6:25]
aa[, 6:25] <- aa.625 / rowSums(aa.625)
# add proteins to thermo$protein
iprotein <- add.protein(aa)
# add proteins to thermo$obigt
iobigt <- info(paste(aa$protein, aa$organism, sep="_"))
### closed system calculation (constant composition of elements)
# use equal initial abundances
Y <- rep(100, length(iobigt))
# run the Gibbs energy minimization (this did not iterate before 20130109,
# due to bug in calculation of free energy derivative)
w <- run.wjd(iobigt, Y=Y, Gfrac=1e-15, nlambda=501)
# the molar abundances
X.closed <- w$X
# get the chemical potentials of the elements
ep <- equil.potentials(w)
# the corresponding logarithms of activities of the basis species
basis("CHNOS")
bl <- basis.logact(ep)
### open system calculation (constant chemical potentials of basis species)
# set the activities of the basis species
basis(names(bl), bl)
# get the affinities of the formation reactions
a <- affinity(iprotein=iprotein)
# then the equilibrium abundances, with total moles of residues as used in wjd
e <- equilibrate(a, loga.balance=log10(sum(Y)))
X.open <- 10^unlist(e$loga.equil)
# the test: abundances calculated both ways are equal
stopifnot(all.equal(X.closed, X.open, tolerance=0.018))
# seems that we could do better than that 1.8% mean difference!
reset()
