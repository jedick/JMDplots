[![DOI](https://zenodo.org/badge/211601502.svg)](https://zenodo.org/badge/latestdoi/211601502)

# JMDplots

This R package has code and data for papers by [Jeffrey M. Dick](https://chnosz.net/jeff/).

## Index to analysis scripts, datafiles, and external links

<details>

<summary>Water and oxygen in the genomic adaptation of human microbiomes</summary>

- [inst/extdata/microhum](inst/extdata/microhum): scripts and processed data files

  - [ARAST](inst/extdata/microhum/ARAST): metagenomes

    - [ARAST.R](inst/extdata/microhum/ARAST/ARAST.R): *script*: metagenome processing pipeline
    - [runARAST.R](inst/extdata/microhum/ARAST/runARAST.R): *script*: run pipeline for particular metagenomes
    - [*_aa.csv](inst/extdata/microhum/ARAST/): *output files*: amino acid composision
    - [*_stats.csv](inst/extdata/microhum/ARAST/): *output files*: processing statistics

  - [KWL22](inst/extdata/microhum/KWL22): metagenome-assembled genomes (MAGs) from [Ke et al. (2022)](https://doi.org/10.1038/s41467-022-32991-w)
  - [metaproteome](inst/extdata/microhum/metaproteome): process metaproteomes

    - [*/mkaa.R](inst/extdata/microhum/metaproteome): *scripts*: metaproteomes processing
    - [*/aa.csv](inst/extdata/microhum/metaproteome): *output files*: amino acid composition

  - [16S](inst/extdata/microhum/16S): 16S rRNA

    - [metadata](inst/extdata/microhum/16S/metadata): *data*: sample metadata for 16S rRNA datasets
    - [pipeline.R](inst/extdata/microhum/16S/pipeline.R): *script*: 16S rRNA processing pipeline
    - [RDP-GTDB](inst/extdata/microhum/16S/RDP-GTDB): *output files*: taxonomic classifications for 16S rRNA datasets

  - [images](inst/extdata/microhum/images): drawings for figure in paper

  - [MR18_Table_S1_modified.csv](inst/extdata/microhum/MR18_Table_S1_modified.csv): *data*: List of Prokaryotes according to their Aerotolerant or Obligate Anaerobic Metabolism, modified from [Million and Raoult (2018)](https://doi.org/10.1016/j.humic.2018.07.002)

- [R/microhum.R](R/microhum.R): code for plots
- [man/microhum.Rd](man/microhum.Rd): manual page
- [vignettes/microhum.Rmd](vignettes/microhum.Rmd): vignette (runs code to make each plot)

  - [microhum.html](https://chnosz.net/JMDplots/vignettes/microhum.html): compiled HTML version of the vignette (external link)

- [bioRxiv](https://doi.org/10.1101/2023.02.12.528246): preprint (external link)

</details>

## Installation

First install the **remotes** packages from CRAN.

```R
install.packages("remotes")
```

Then install other required packages.
- [**canprot**](https://github.com/jedick/canprot) has data and functions for analyzing differential expression of proteins in cancer and cell-culture experiments.
- [**chem16S**](https://github.com/jedick/chem16S) processes 16S rRNA-based taxonomic classifications to calculate chemical metrics of community reference proteomes.

```R
remotes::install_github("jedick/canprot")
remotes::install_github("jedick/chem16S")
```

> **Note**
>
> As of 2023-06-17, **JMDplots** depends on the development version of **canprot** from GitHub, not the released version on CRAN.
> As of 2023-07-31, **JMDplots** depends on the development version of **chem16S** from GitHub, not the released version on CRAN.

Finally, install **JMDplots**.
This command will install prebuilt vignettes; they might not be up-to-date with the code.

```R
remotes::install_github("jedick/JMDplots")
```

To view the plots, use the R help browser to open the vignettes page or open it directly with this command:

```R
browseVignettes("JMDplots")
```

## Building vignettes

If you want to build the vignettes yourself, note that it might be possible to build them without [pandoc](https://pandoc.org/), but having pandoc available will make them look better.

```R
remotes::install_github("jedick/JMDplots", dependencies = TRUE, build_vignettes = TRUE)
```

## Online vignettes

The vignettes can be viewed at <https://chnosz.net/JMDplots/vignettes/>.

## Licenses

This package except for the file `inst/extdata/orp16S/metadata/PCL+18.csv` is licensed under the GNU General Public License v3 (GPLv3).

The ORP (mV), DO (mg/L) and Feature (Stream, Spring, Lake, Terrace, or Geyser) data for New Zealand hot springs ([Power et al., 2018](https://doi.org/10.1038/s41467-018-05020-y)) in `PCL+18.csv` were obtained from the [1000 Springs Project](https://1000springs.org.nz) and are licensed under CC-BY-NC-SA.

This package contains a copy of the `dunnTest()` function by Derek H. Ogle from CRAN package [FSA](https://cran.r-project.org/package=FSA), version 0.9.3 (License: GPL (>= 2)), which itself is a wrapper for `dunn.test()` from CRAN package [dunn.test](https://cran.r-project.org/package=dunn.test) by Alexis Dinno.
