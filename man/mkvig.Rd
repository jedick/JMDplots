\encoding{UTF-8}
\name{mkvig}
\alias{mkvig}
\title{Compile and View Vignettes from the Command Line}
\description{
Compile the indicated vignette and opens it in the browser.
}

\usage{
  mkvig(vig = NULL)
}

\arguments{
  \item{vig}{character, name of a vignette without \samp{.Rmd} extension}
}

\details{
This function builds and views vignettes for chemical analysis of differential expression datasets.
In order to reduce package building time, these are not included in the pre-built vignettes of the package.


The available vignettes for \code{mkvig} are listed here:
\itemize{
    \item \emph{Cell culture} -- \samp{hypoxia}, \samp{secreted}, \samp{osmotic_bact}, \samp{osmotic_euk}, \samp{osmotic_halo}, \samp{glucose}, \samp{3D}, \samp{osmotic_gene}, \samp{yeast_stress}
    \item \emph{Cancer} -- \samp{breast}, \samp{colorectal}, \samp{liver}, \samp{lung}, \samp{pancreatic}, \samp{prostate}
    \item \emph{Pan-cancer} -- \samp{TCGA}, \samp{HPA}
}

Note that pandoc (including pandoc-citeproc), as a system dependency of \pkg{rmarkdown}, must be installed.
See \pkg{rmarkdown}'s \samp{pandoc} vignette for installation tips.
}

\seealso{
The vignettes can also be run using e.g. \code{\link{demo}("HPA")}, and through the interactive help system (\code{\link{help.start}} > Packages > JMDplots > Code demos).
}

\examples{
\dontrun{
mkvig("osmotic_gene")
}
}
