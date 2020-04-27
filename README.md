[![DOI](https://zenodo.org/badge/211601502.svg)](https://zenodo.org/badge/latestdoi/211601502)

# JMDplots

This R package has code and data for plots from papers by [Jeffrey M. Dick](http://chnosz.net/jeff).

## Installation

First install the **remotes** packages from CRAN.

```R
install.packages("remotes")
```

Although it may be possible to build the vignettes without pandoc, having pandoc available will make them look better.
See rmarkdown's [Install Pandoc](https://cran.r-project.org/web/packages/rmarkdown/vignettes/pandoc.html) vignette for tips on installing pandoc.

Install the development version of **canprot**, then install **JMDplots**.

```R
remotes::install_github("jedick/canprot")
remotes::install_github("jedick/JMDplots", build_vignettes = TRUE)
```

Finally, open the vignettes index page to see the plots for the different papers.

```R
browseVignettes("JMDplots")
```

## Online vignettes

The vignettes are available at <http://chnosz.net/JMDplots/doc>.
