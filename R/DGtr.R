# JMDplots/DGtr.R
# Calculate Gibbs energy of transformation
# Code was present in objective.R of CHNOSZ prior to version 2.0.0 20121009
# Moved DGtr() to JMDplots 20221201

# \code{DGtr} calculates the change in Gibbs energy/2.303RT of a system in which species with initial logarithms of activitiy (\code{loga1}) are transformed to the same species with different final logarithms of activity (\code{loga2}) at constant temperature, pressure and chemical potentials of basis species.
# It is calculated as the sum over species of (G2-G1) where G1/RT = -a1*Astar + a1*loga1 - a1 + a constant (where a1 is 10^loga1), likewise for G2.
# \code{Astar} is the starved affinity; that is, the affinity of the reaction to form one mole of the species at unit activity from the basis species in their defined activities.
# The equation is derived by integrating dG = -A/dxi = -A/dn where xi is the reaction progress variable, dn/dxi = 1 is the reaction coefficient on the species, and A = Astar - 2.303RTloga is the chemical affinity (Dick and Shock, 2013).

# NOTE: Astar is provided in the output from CHNOSZ::equilibrate()

DGtr <- function(loga1, loga2, Astar) {
  dgtr <- function(loga1, loga2, Astar) {
    # calculate the Gibbs energy/2.303RT of transformation
    # from an initial to a final assemblage at constant T, P and 
    # chemical potentials of basis species 20120917
    # loga1 - logarithms of activity of species in the initial assemblage
    # loga2 - logarithms of activity of species in the final assemblage
    # Astar - starved (of activity of the species of interest) values of chemical affinity/2.303RT
    # remove the logarithms
    a1 <- 10^loga1
    a2 <- 10^loga2
    # the molal Gibbs energy in the initial and final states
    # (derived from integrating A = Astar - RTln(a) from a1 to a2)
    G1 <- -a1*Astar + a1*loga1 - a1/log(10)
    G2 <- -a2*Astar + a2*loga2 - a2/log(10)
    # calculate the change in molal Gibbs energy for each species
    DG <- G2 - G1
    # return the sum
    return(sum(DG))
  }
  # we need to index both loga1 and Astar
  DGtr <- unlist(lapply(seq(nrow(loga1)), function(i) {
    dgtr(loga1[i, ], loga2, Astar[i, ])
  }))
  return(DGtr)
}
