[![DOI](https://zenodo.org/badge/211601502.svg)](https://zenodo.org/badge/latestdoi/211601502)

# JMDplots

This R package has code and data for plots from papers by [Jeffrey M. Dick](http://chnosz.net/jeff).
Required R packages are [CHNOSZ](http://chnosz.net) and [canprot](https://github.com/jedick/canprot).
Both of these packages are available on CRAN, but the development version of **canprot** (from GitHub) is required.

## Installation

First install the **remotes** packages from CRAN.

```R
install.packages("remotes")
```

Although it may be possible to build the vignettes without pandoc, having pandoc available will make them look better.
See rmarkdown's [Install Pandoc](https://cran.r-project.org/web/packages/rmarkdown/vignettes/pandoc.html) vignette for tips on installing pandoc.
Then install **JMDplots** from GitHub (including its dependencies).

```R
remotes::install_github("jedick/JMDplots", build_vignettes = TRUE)
```

Finally, open the vignettes index page to see the plots for the different papers.

```R
browseVignettes("JMDplots")
```

## Online vignettes

The vignettes are available at <http://chnosz.net/JMDplots/doc>.
