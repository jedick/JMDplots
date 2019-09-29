\encoding{UTF-8}
\name{gradox}
\alias{gradox}
\alias{gradox1}
\alias{gradox2}
\alias{gradoxS1}
\alias{gradoxS2}
\alias{gradox3}
\alias{gradox4}
\alias{gradox5}
\title{Plots for \dQuote{gradox} paper}
\description{
Plots for the paper by Dick et al. (2019).
}

\usage{
  gradox1(pdf = FALSE)
  gradox2(pdf = FALSE)
  gradoxS1(pdf = FALSE)
  gradoxS2(pdf = FALSE)
  gradox3(mout, pout, pdf = FALSE)
  gradox4(mout, pout, pdf = FALSE)
  gradox5(pdf = FALSE, maxdepth = 500)
}

\arguments{
  \item{pdf}{logical, make a PDF file?}
  \item{mout}{list, value returned by \code{gradoxS1} (DNA \ZC)}
  \item{pout}{list, value returned by \code{gradoxS2} (protein \ZC)}
  \item{maxdepth}{numeric, maximum sample depth (meters)}
}

\details{
This table gives a brief description of each function.

\tabular{ll}{
  \code{gradox1} \tab Some general characteristics of \ZC of DNA, RNA, and proteins \cr
  \code{gradox2} \tab Selected plots of DNA, RNA, and protein \ZC \cr
  \code{gradoxS1} \tab All plots of DNA and RNA \ZC \cr
  \code{gradoxS2} \tab All plots of protein \ZC \cr
  \code{gradox3} \tab \ZC of proteins vs DNA (metagenomes and metatranscriptomes) \cr
  \code{gradox4} \tab Thermodynamic calculations of relative stabilities along redox gradients \cr
  \code{gradox5} \tab \ZC of reads classified to selected abundant species \cr
}
}

\references{
Dick, J. M., Yu, M., Tan, J. and Lu, A. (2019) Changes in carbon oxidation state of metagenomes along geochemical redox gradients. \emph{Front. Microbiol.} \bold{10}, 120. \doi{10.3389/fmicb.2019.00120}
}

\examples{
gradox1()

mout <- gradoxS1()
pout <- gradoxS2()
gradox3(mout, pout)
}
