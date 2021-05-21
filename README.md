[![DOI](https://zenodo.org/badge/211601502.svg)](https://zenodo.org/badge/latestdoi/211601502)

# JMDplots

This R package has code and data for plots from papers by [Jeffrey M. Dick](http://chnosz.net/jeff/).

## Installation

First install the **remotes** packages from CRAN.

```R
install.packages("remotes")
```

Then install **JMDplots**.
This command will also install prebuilt vignettes.

```R
remotes::install_github("jedick/JMDplots")
```

To see the plots, use the R help browser to open the vignettes page or open it directly with this command:

```R
browseVignettes("JMDplots")
```

## Building vignettes

If you want to build the vignettes yourself, be aware that it may be possible to build them without [pandoc](https://pandoc.org/), but having pandoc available will make them look better.

```R
remotes::install_github("jedick/JMDplots", build_vignettes = TRUE)
```

## Online vignettes

The vignettes are available at <http://chnosz.net/JMDplots/doc/>.
