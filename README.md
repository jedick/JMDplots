[![DOI](https://zenodo.org/badge/211601502.svg)](https://zenodo.org/badge/latestdoi/211601502)

# JMDplots

This R package has code and data for papers by [Jeffrey M. Dick](https://chnosz.net/jeff/).

## Quick links to active papers (pre-publication drafts)

- **orp16S**: Chemical adaptation of estimated bacterial community proteomes to redox potential: [data directory](inst/extdata/orp16S) with sequence processing pipeline,
  [R code](R/orp16S.R), [help page source](man/orp16S.Rd), [vignette source](vignettes/orp16S.Rmd)
  - A compiled HTML version of the vignette is at [chnosz.net](https://chnosz.net/JMDplots/doc/orp16S.html).
  - A preprint is on [bioRxiv](https://doi.org/10.1101/2021.10.12.464155).
- **utogig**: Using thermodynamics to obtain geochemical information from genomes: [data directory](inst/extdata/utogig), [R code](R/utogig.R),
  [help page source](man/utogig.Rd), [vignette source](vignettes/utogig.Rmd)

## Installation

First install the **remotes** packages from CRAN.

```R
install.packages("remotes")
```

Then install [**chem16S**](../chem16S).
This package is used to calculate chemical metrics of estimated community proteomes from 16S rRNA data.

```R
remotes::install_github("jedick/chem16S")

```

Then install **JMDplots**.
This command will also install prebuilt vignettes.

```R
remotes::install_github("jedick/JMDplots")
```

To view the plots, use the R help browser to open the vignettes page or open it directly with this command:

```R
browseVignettes("JMDplots")
```

## Building vignettes

If you want to build the vignettes yourself, be aware that it may be possible to build them without [pandoc](https://pandoc.org/), but having pandoc available will make them look better.

```R
remotes::install_github("jedick/JMDplots", dependencies = TRUE, build_vignettes = TRUE)
```

NOTE (2022-10-24): The development version of CHNOSZ is required for the affinity ranking calculations in the `utogig.Rmd` vignette.
To install the development version, use `install.packages("CHNOSZ", repos = "https://R-Forge.R-project.org")` or `remotes::install_github("jedick/CHNOSZ")`.

## Online vignettes

The vignettes can be viewed at <https://chnosz.net/JMDplots/doc/>.

## Licenses

This package except for the file `inst/extdata/orp16S/metadata/PCL+18.csv` is licensed under the GNU General Public License v3 (GPLv3).

The ORP (mV), DO (mg/L) and Feature (Stream, Spring, Lake, Terrace, or Geyser) data for New Zealand hot springs ([Power et al., 2018](https://doi.org/10.1038/s41467-018-05020-y)) in `PCL+18.csv` were obtained from the [1000 Springs Project](https://1000springs.org.nz) and are licensed under CC-BY-NC-SA.

This package contains a copy of the `dunnTest()` function by Derek H. Ogle from CRAN package [FSA](https://cran.r-project.org/package=FSA), version 0.9.3 (License: GPL (>= 2)), which itself is a wrapper for `dunn.test()` from CRAN package [dunn.test](https://cran.r-project.org/package=dunn.test) by Alexis Dinno.
