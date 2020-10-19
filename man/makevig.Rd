\encoding{UTF-8}
\name{makevig}
\alias{makevig}
\title{Compile and View Vignettes from the Command Line}
\description{
This function compiles the indicated vignette and opens it in the browser.
}

\usage{
  makevig(vig = NULL)
}

\arguments{
  \item{vig}{character, name of a vignette without \samp{.Rmd} extension}
}

\details{
This function builds and views the vignettes for compositional analysis of differential expression datasets.
In order to reduce package building time, these are not included in the pre-built vignettes of the package.

The available vignettes are \samp{TCGA}, \samp{HPA}, \samp{osmotic_gene}, and \samp{yeast_stress}.

Note that pandoc (including pandoc-citeproc), as a system dependency of \pkg{rmarkdown}, must be installed.
See \pkg{rmarkdown}'s \samp{pandoc} vignette for installation tips.
}

\seealso{
The vignettes can also be run using e.g. \code{\link{demo}("HPA")}, and through the interactive help system (\code{\link{help.start}} > Packages > JMDplots > Code demos).
The counterpart to this function in \pkg{canprot} is \code{\link[canprot]{mkvig}}.
}

\examples{
\dontrun{
makevig("osmotic_gene")
}
}
