\encoding{UTF-8}
\name{chnszten}
\alias{chnszten}
\alias{chnszten1}
\alias{chnszten4}
\alias{chnszten5}
\alias{chnszten6}
\alias{chnszten7}
\alias{chnsztenS1}
\alias{chnsztenS2}
\alias{chnsztenS3}
\alias{chnsztenS4}
\alias{chnsztenS5}
\alias{chnsztenS6A}
\alias{chnsztenS6B}
\alias{chnsztenS7}
\alias{chnsztenS8}
\alias{chnsztenS9}
\alias{chnsztenS10}
\title{Plots from CHNOSZ paper}
\description{
Plots from the paper by Dick (2019).
}

\usage{
  chnszten1()
  chnszten4(pdf = FALSE)
  chnszten5(pdf = FALSE)
  chnszten6(pdf = FALSE)
  chnszten7(pdf = FALSE)
  chnsztenS1(pdf = FALSE)
  chnsztenS2(pdf = FALSE)
  chnsztenS3(pdf = FALSE)
  chnsztenS4(pdf = FALSE)
  chnsztenS5(pdf = FALSE)
  chnsztenS6A(pdf = FALSE)
  chnsztenS6B(pdf = FALSE)
  chnsztenS7(pdf = FALSE)
  chnsztenS8(pdf = FALSE)
  chnsztenS9(pdf = FALSE)
  chnsztenS10(pdf = FALSE)
}

\arguments{
  \item{pdf}{logical, make a PDF file?}
}

\details{
This table gives a brief description of each function.

\tabular{ll}{
  \code{chnszten1} \tab Timeline of CHNOSZ development \cr
  \code{chnszten4} \tab Mosaic Eh-pH diagram for the Cu-S-Cl-O-H system \cr
  \code{chnszten5} \tab Corundum solubility diagram \cr
  \code{chnszten6} \tab DEW model (based on Fig. 3 of Sverjensky et al., 2014) \cr
  \code{chnszten7} \tab Thermodynamic properties of reactions involving Al-bearing minerals \cr
  \code{chnsztenS1} \tab Comparison of \logK and maximum affinity methods \cr
  \code{chnsztenS2} \tab Eh-pH diagram like Fig. 5A of Caporuscio et al. (2017) \cr
  \code{chnsztenS3} \tab Gibbs energy of transformation for an assemblage of n-alkanes \cr
  \code{chnsztenS4} \tab \code{\link[CHNOSZ]{findit}} calculations for sulfur species \cr
  \code{chnsztenS5} \tab Debye-Hückel extended term parameter extrapolated from plots of Manning et al. (2013) \cr
  \code{chnsztenS6A} \tab \code{chnszten6} modified to exclude DEW data for acetate \cr
  \code{chnsztenS6B} \tab \code{chnsztenS6A} modified to use default \code{\link[CHNOSZ]{bgamma}} equation \cr
  \code{chnsztenS7} \tab \logK of NaCl dissociation \cr
  \code{chnsztenS8} \tab Calcite solubility: comparison with Manning et al. (2013) \cr
  \code{chnsztenS9} \tab Compare gold solubility in HCh and CHNOSZ: hematite-magnetite buffer \cr
  \code{chnsztenS10} \tab Compare gold solubility in HCh and CHNOSZ: pyrite-pyrrhotite-magnetite buffer \cr
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
if(interactive() & shinyOK & timevisOK) chnszten1()

# Mosaic Eh-pH diagram (creates a PDF file)
if(interactive()) chnszten4()
}
