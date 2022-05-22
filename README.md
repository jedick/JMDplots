[![DOI](https://zenodo.org/badge/211601502.svg)](https://zenodo.org/badge/latestdoi/211601502)

# JMDplots

This R package has code and data for plots from papers by [Jeffrey M. Dick](http://chnosz.net/jeff/).

## Quick links to active papers (pre-publication drafts)

- **orp16S**: Chemical adaptation of estimated bacterial community proteomes to redox potential: [data directory](inst/extdata/orp16S),
  [R code](R/orp16S.R), [help page source](man/orp16S.Rd), [vignette source](vignettes/orp16S.Rmd)
  - A compiled HTML version of the vignette is at [chnosz.net](https://chnosz.net/JMDplots/doc/orp16S.html) (it might be out of sync with the source).
  - A preprint is on [bioRxiv](https://doi.org/10.1101/2021.10.12.464155).
- **gcbio**: Geochemical biology perspective: [data directory](inst/extdata/gcbio), [R code](R/gcbio.R),
  [help page source](man/gcbio.Rd), [vignette source](vignettes/gcbio.Rmd)

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
remotes::install_github("jedick/JMDplots", build_vignettes = TRUE)
```

## Online vignettes

The vignettes can be viewed at <http://chnosz.net/JMDplots/doc/>.

## Licenses

This package except for the file `inst/extdata/orp16S/metadata/PCL+18.csv` is licensed under the GNU General Public License v3 (GPLv3).

The ORP (mV), DO (mg/L) and Feature (Stream, Spring, Lake, Terrace, or Geyser) data for New Zealand hot springs ([Power et al., 2018](https://doi.org/10.1038/s41467-018-05020-y)) in `PCL+18.csv` were obtained from the [1000 Springs Project](http://1000springs.org.nz) and are licensed under CC-BY-NC-SA.
