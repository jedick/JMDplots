\encoding{UTF-8}
\name{chnosz10}
\alias{chnosz10}
\alias{chnosz101}
\alias{chnosz104}
\alias{chnosz105}
\alias{chnosz106}
\alias{chnosz107}
\alias{chnosz10S1}
\alias{chnosz10S2}
\alias{chnosz10S3}
\alias{chnosz10S4}
\alias{chnosz10S5}
\alias{chnosz10S6A}
\alias{chnosz10S6B}
\alias{chnosz10S7}
\alias{chnosz10S8}
\alias{chnosz10S9}
\alias{chnosz10S10}
\title{Plots from CHNOSZ paper}
\description{
Plots from the paper by Dick (2019).
}

\usage{
  chnosz101()
  chnosz104(pdf = FALSE)
  chnosz105(pdf = FALSE)
  chnosz106(pdf = FALSE)
  chnosz107(pdf = FALSE)
  chnosz10S1(pdf = FALSE)
  chnosz10S2(pdf = FALSE)
  chnosz10S3(pdf = FALSE)
  chnosz10S4(pdf = FALSE)
  chnosz10S5(pdf = FALSE)
  chnosz10S6A(pdf = FALSE)
  chnosz10S6B(pdf = FALSE)
  chnosz10S7(pdf = FALSE)
  chnosz10S8(pdf = FALSE)
  chnosz10S9(pdf = FALSE)
  chnosz10S10(pdf = FALSE)
}

\arguments{
  \item{pdf}{logical, make a PDF file?}
}

\details{
This table gives a brief description of each function.

\tabular{ll}{
  \code{chnosz101} \tab Timeline of CHNOSZ development \cr
  \code{chnosz104} \tab Mosaic Eh-pH diagram for the Cu-S-Cl-O-H system \cr
  \code{chnosz105} \tab Corundum solubility diagram \cr
  \code{chnosz106} \tab DEW model (based on Fig. 3 of Sverjensky et al., 2014) \cr
  \code{chnosz107} \tab Thermodynamic properties of reactions involving Al-bearing minerals \cr
  \code{chnosz10S1} \tab Comparison of \logK and maximum affinity methods \cr
  \code{chnosz10S2} \tab Eh-pH diagram like Fig. 5A of Caporuscio et al. (2017) \cr
  \code{chnosz10S3} \tab Gibbs energy of transformation for an assemblage of n-alkanes \cr
  \code{chnosz10S4} \tab \code{\link[CHNOSZ]{findit}} calculations for sulfur species \cr
  \code{chnosz10S5} \tab Debye-Hückel extended term parameter extrapolated from plots of Manning et al. (2013) \cr
  \code{chnosz10S6A} \tab \code{chnosz106} modified to exclude DEW data for acetate \cr
  \code{chnosz10S6B} \tab \code{chnosz10S6A} modified to use default \code{\link[CHNOSZ:nonideal]{bgamma}} equation \cr
  \code{chnosz10S7} \tab \logK of NaCl dissociation \cr
  \code{chnosz10S8} \tab Calcite solubility: comparison with Manning et al. (2013) \cr
  \code{chnosz10S9} \tab Compare gold solubility in HCh and CHNOSZ: hematite-magnetite buffer \cr
  \code{chnosz10S10} \tab Compare gold solubility in HCh and CHNOSZ: pyrite-pyrrhotite-magnetite buffer \cr
}
}

\references{
Caporuscio, F. A., Palaich, S. E. M., Cheshire, M. C. and Jové Colon, C. F. (2017) Corrosion of copper and authigenic sulfide mineral growth in hydrothermal bentonite experiments. \emph{J. Nucl. Mater.} \bold{485}, 137--146. \doi{10.1016/j.jnucmat.2016.12.036}

Dick, J. M., Yu, M., Tan, J. and Lu, A. (2019) Changes in carbon oxidation state of metagenomes along geochemical redox gradients. \emph{Front. Microbiol.} \bold{10}, 120. \doi{10.3389/fmicb.2019.00120}

Manning, C. E., Shock, E. L. and Sverjensky, D. A. (2013) The chemistry of carbon in aqueous fluids at crustal and upper-mantle conditions: Experimental and theoretical constraints. \emph{Rev. Mineral. Geochem.} \bold{75}, 109--148. \doi{10.2138/rmg.2013.75.5}

Sverjensky, D. A., Stagno, V. and Huang, F. (2014) Important role for organic carbon in subduction-zone fluids in the deep carbon cycle. \emph{Nat. Geosci.} \bold{7}, 909--913. \doi{10.1038/ngeo2291}
}

\examples{
# Display timeline of CHNOSZ development in web brower
# Runs Shiny app; use Ctrl-C to return to command line
shinyOK <- require("shiny", quietly = TRUE)
timevisOK <- require("timevis", quietly = TRUE)
if(interactive() & shinyOK & timevisOK) chnosz101()

# Mosaic Eh-pH diagram (creates a PDF file)
if(interactive()) chnosz104()
}
