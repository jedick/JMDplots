[![DOI](https://zenodo.org/badge/211601502.svg)](https://zenodo.org/badge/latestdoi/211601502)

# JMDplots

This R package has code and data for plots from papers by [Jeffrey M. Dick](http://chnosz.net/jeff).
Required R packages are [CHNOSZ](http://chnosz.net), and [canprot](https://github.com/jedick/canprot).
Both of these packages are available on CRAN, but the development version of **canprot** (from GitHub) is required.

## Installation from GitHub

First install the **remotes** packages from CRAN.

```R
install.packages("remotes")
```

Then install **canprot** and **JMDplots** from GitHub, including all dependencies needed to build the vignettes.

```R
remotes::install_github("jedick/canprot", dependencies = TRUE, build_vignettes = TRUE)
remotes::install_github("jedick/JMDplots", dependencies = TRUE, build_vignettes = TRUE)
```

Finally open the vignettes to view the plots.

```R
browseVignettes("JMDplots")
```

## Online vignettes

The HTML vignettes are available at [http://chnosz.net/JMDplots/doc](http://chnosz.net/JMDplots/doc).
