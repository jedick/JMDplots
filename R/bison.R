# JMDplots/bison.R
# Plots from hot spring (Bison Pool) papers (2011, 2013)
# Code moved from CHNOSZ/hotspring.Rnw 20200712

#add.obigt("OldAA")

bison1 <- function() {
  bison.T <- c(93.3, 79.4, 67.5, 65.3, 57.1)
  bison.pH <- c(7.350, 7.678, 7.933, 7.995, 8.257)

  distance <- c(0, 6, 11, 14, 22)
  par(mfrow = c(1, 2), mar = c(4, 4, 3, 2))
  xpoints <- seq(0, 22, length.out = 128)
  # T plot
  plot(distance, bison.T, xlab = "distance, m", ylab = axis.label("T"))
  Tfun <- splinefun(distance, bison.T, method = "mono")
  lines(xpoints, Tfun(xpoints))

  # pH plot
  plot(distance, bison.pH, xlab = "distance, m", ylab = "pH")
  pHfun <- splinefun(distance, bison.pH, method = "mono")
  lines(xpoints, pHfun(xpoints))
}

### UNEXPORTED OBJECTS ###
