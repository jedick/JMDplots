[![DOI](https://zenodo.org/badge/211601502.svg)](https://zenodo.org/badge/latestdoi/211601502)

# JMDplots

This R package has code and data for plots from papers by [Jeffrey M. Dick](http://chnosz.net/jeff).

Requirements are R (>= 3.1.0), [CHNOSZ](http://chnosz.net) (>= 1.3.2), and [canprot](https://github.com/jedick/canprot) (>= 0.1.3).
knitr and rmarkdown (and a pandoc installation) are needed to build the vignettes.

## Installation from Github

First install some packages from CRAN.
**remotes** provides the `install_github` function.
**knitr** is used to build the vignettes.
**rmarkdown** creates self-contained HTML vignettes.

```R
install.packages(c("remotes", "knitr", "rmarkdown"))
```

Then install **canprot** and **JMDplots** from Github including the vignettes.
This will also install CHNOSZ if it is not installed already.

```R
remotes::install_github("jedick/canprot")
remotes::install_github("jedick/JMDplots", build_vignettes = TRUE)
```

Finally open the vignettes to view the plots.

```R
browseVignettes("JMDplots")
```

## Online vignettes

The HTML vignettes are online at [http://chnosz.net/JMDplots/doc](http://chnosz.net/JMDplots/doc).
