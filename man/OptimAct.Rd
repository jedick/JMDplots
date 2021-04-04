\encoding{UTF-8}
\name{MaximAct}
\alias{MaximAct}
\title{Calculate optimal logaH2O and logfO2 for target proteins}
\description{
Calculate the \logaH2O and \logfO2 that maximize the activities of given target proteins on a background of proteins from the human proteome.
}

\usage{
  MaximAct(aa, seed = 1:100, nbackground = 2000,
    plot.it = TRUE, filebase = NULL, xlab = "sample", names = NULL,
    O2 = c(-72, -67), H2O = c(-2, 6), pH = NULL)
}

\arguments{
  \item{aa}{data frame, amino acid composition in the format of \code{\link[CHNOSZ]{thermo}()$protein}}
  \item{seed}{numeric, seeds to use for subsampling background proteins from human proteome}
  \item{nbackground}{numeric, number of background proteins to sample}
  \item{plot.it}{logical, make a plot?}
  \item{filebase}{character, base name to use for output files}
  \item{xlab}{character or expression, label for \emph{x} axis}
  \item{names}{character, sample names}
  \item{O2}{numeric, range of \logfO2 values}
  \item{H2O}{numeric, range of \logaH2O values}
  \item{pH}{numeric, range of pH values}
}

\details{
\code{MaximAct} is used to compute optimal values of \logaH2O, \logfO2, and optionally pH; that is, those that maximize the equilibrium activity of each of the target proteins (given in \code{aa}) against a background of proteins in the human proteome.
The calculation is described by Dick (2021).

The basis species used for the calculation are \samp{QEC} (glutamine - glutamic acid - cysteine - O2 - H2O) (see Dick et al., 2020), or \samp{QEC+} (the same basis species with the addition of \Hplus) if \code{pH} is not NULL.
The NULL default for \code{pH} means that the system consists of uncharged species, and pH is therefore not defined.
The ranges in \code{O2}, \code{H2O}, and \code{pH} can have an optional third value to specify the resolution; e.g. \code{c(-70, -65, 100)} represents 100 points in the range [-70, -65] (see \code{\link[CHNOSZ]{affinity}}).

Due to memory constraints, a single calculation does not include the entire human proteome, but a random subsample of \code{nbackground} human proteins.
The subsampling is performed for each random seed in \code{seed}.
The results are returned in a list with components \code{O2} and \code{H2O}.

By default, a pair of plots is made that are updated during the calculation; at the end, the mean values of \logfO2 and \logaH2O are plotted as red lines.
Change \code{plot.it} to FALSE to  prevent making the plots.

If \code{filebase} is not NULL, the return value is made \code{\link{invisible}} and the results are saved in CSV files named by concatenating \code{filebase} with \samp{_O2.csv} (for \logfO2) or \samp{_H2O.csv} (for \logaH2O).
The protein names (used in both the plots and CSV files) are taken from \code{aa$protein} or the \code{names} argument if it is not NULL.
}

\seealso{
\code{\link{runMaximAct}} to set up calculations for datasets analyzed in the \code{\link{evdevH2O}} paper.
}

\references{
Dick, J. M., Yu, M. and Tan, J. (2020) Uncovering chemical signatures of salinity gradients through compositional analysis of protein sequences. \emph{Biogeosciences} \bold{17}, 6145--6162. \doi{10.5194/bg-17-6145-2020}

Dick, J. M. (2021) Redox potential linked to water loss from proteins in evolution and developent. bioRxiv. \doi{https://doi.org/10.1101/2021.01.29.428804}
}

\examples{
# An abbreviated version of runMaximAct("transcriptome")
datadir <- system.file("extdata/evdevH2O", package = "JMDplots")
aa <- read.csv(file.path(datadir, "FOK+21_mean_aa.csv"), as.is = TRUE)
aa <- subset(aa, organism == "transcriptome")
OA <- MaximAct(aa, seed = 1:3, nbackground = 200, O2 = c(-72, -65), H2O = c(-2, 5))
## Default arguments used in runMaximAct -- takes much longer
#OA <- MaximAct(aa, seed = 1:100, nbackground = 2000, O2 = c(-72, -65), H2O = c(-2, 5))
# The logaH2O decreases at the end of the time course
H2O <- colMeans(OA$H2O)
stopifnot(H2O["14D"] > H2O["1M"])
stopifnot(H2O["1M"] > H2O["2M"])
# The final value of logaH2O is close to 0
stopifnot(round(H2O["2M"]) == 0)
}
