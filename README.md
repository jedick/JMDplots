# JMDplots

This R package has code and data for plots from papers by Jeffrey M. Dick.

Requirements are R (>= 3.1.0), [CHNOSZ](http://chnosz.net) (>= 1.3.2), and [canprot](https://github.com/jedick/canprot) (> 0.1.2, i.e. the current development version).
knitr and rmarkdown are suggested for building the vignettes.

## Installation from Github

First install some packages from CRAN.
**remotes** provides the `install_github` function.
**knitr** is needed to build the vignettes.
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
