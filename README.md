[![DOI](https://zenodo.org/badge/211601502.svg)](https://zenodo.org/badge/latestdoi/211601502)

# JMDplots

This R package has code and data for papers by [Jeffrey M. Dick](https://chnosz.net/jeff/).
The vignettes in this package have many plots from my papers and can be viewed at <https://chnosz.net/JMDplots/vignettes/>.

## Analysis scripts, data files, and external links for papers

Click on the titles for a detailed list of files; click on the years to go to the published papers.
The manual page for each paper has additional details about scripts, data files, and plotting functions.

<!-- Put a space before <details> to make ghostwriter format the lists correctly -->
 <details>

<summary>Water and oxygen in the genomic adaptation of human microbiomes (*submitted manuscript*)</summary>

- [inst/extdata/microhum](inst/extdata/microhum): scripts and processed data files

  - [ARAST](inst/extdata/microhum/ARAST): analysis of metagenomes

    - [ARAST.R](inst/extdata/microhum/ARAST/ARAST.R): *script*: metagenome processing pipeline
    - [runARAST.R](inst/extdata/microhum/ARAST/runARAST.R): *script*: run pipeline for particular metagenomes
    - [*_aa.csv](inst/extdata/microhum/ARAST/): *output files*: amino acid composision
    - [*_stats.csv](inst/extdata/microhum/ARAST/): *output files*: processing statistics

  - [KWL22](inst/extdata/microhum/KWL22): analysis of metagenome-assembled genomes (MAGs) from [Ke et al. (2022)](https://doi.org/10.1038/s41467-022-32991-w)

    - [mkaa.R](inst/extdata/microhum/KWL22/mkaa.R): *script*: metaproteome processing
    - [KWL22_MAGs_prodigal_aa.csv.xz](inst/extdata/microhum/KWL22/KWL22_MAGs_prodigal_aa.csv.xz): *output file*: amino acid composition
<!--
    - [COVID19_metadata.txt](inst/extdata/microhum/KWL22/COVID19_metadata.txt): *data*: downloaded from <https://github.com/Owenke247/COVID-19/blob/main/Pre-processed_Files/COVID19_metadata.txt>
-->

  - [metaproteome](inst/extdata/microhum/metaproteome): analysis of metaproteomes

    - [*/mkaa.R](inst/extdata/microhum/metaproteome): *scripts*: metaproteome processing
    - [*/aa.csv](inst/extdata/microhum/metaproteome): *output files*: amino acid composition

  - [16S](inst/extdata/microhum/16S): analysis of 16S rRNA gene sequences

    - [metadata](inst/extdata/microhum/16S/metadata): *data*: sample metadata for 16S rRNA datasets
    - [pipeline.R](inst/extdata/microhum/16S/pipeline.R): *script*: 16S rRNA processing pipeline
    - [RDP-GTDB](inst/extdata/microhum/16S/RDP-GTDB): *output files*: taxonomic classifications for 16S rRNA datasets made using the RDP Classifier with a [training set based on GTDB release 207](https://doi.org/10.5281/zenodo.7633100)

  - [images](inst/extdata/microhum/images): drawings for figure in paper

  - [MR18_Table_S1_modified.csv](inst/extdata/microhum/MR18_Table_S1_modified.csv): *data*: List of Prokaryotes according to their Aerotolerant or Obligate Anaerobic Metabolism, modified from [Million and Raoult (2018)](https://doi.org/10.1016/j.humic.2018.07.002)

- [R/microhum.R](R/microhum.R): code for plots
- [man/microhum.Rd](man/microhum.Rd): manual page
- [vignettes/microhum.Rmd](vignettes/microhum.Rmd): vignette including Figures 1--4 and S1--S2

  - [microhum.html](https://chnosz.net/JMDplots/vignettes/microhum.html): compiled HTML version of the vignette (external link)

- [bioRxiv](https://doi.org/10.1101/2023.02.12.528246): preprint (external link)

</details>

 <details>

<summary>Community- and genome-based evidence for a shaping influence of redox potential on bacterial protein evolution (<a href="https://doi.org/10.1128/msystems.00014-23">2023</a>)</summary>

  - [inst/extdata/orp16S](inst/extdata/orp16S): scripts and processed data files

    - [metadata](inst/extdata/orp16S/metadata): *data*: sample metadata for 16S rRNA datasets
    - [pipeline.R](inst/extdata/orp16S/pipeline.R): *script*: 16S rRNA processing pipeline
    - [RDP](inst/extdata/orp16S/RDP): *output files*: taxonomic classifications for 16S rRNA datasets made using the RDP Classifier with its default training set
    - [hydro_p](inst/extdata/orp16S/hydro_p): *data*: shapefiles for the North American Great Lakes, downloaded from [USGS (2010)](https://www.sciencebase.gov/catalog/item/530f8a0ee4b0e7e46bd300dd)
    - [EZdat.csv](inst/extdata/orp16S/EZdat.csv): *output file*: sample data and computed values of Eh7 and \Zc
    - [EZlm.csv](inst/extdata/orp16S/EZlm.csv): *output file*: linear fits between Eh7 and \Zc for each dataset
    - [BKM60.csv](inst/extdata/orp16S/BKM60.csv): *data*: outline of Eh-pH range of natural environments, digitized from Fig. 32 of [Baas Becking et al. (1960)](https://doi.org/10.1086/626659)
    - [MR18_Table_S1.csv](inst/extdata/orp16S/MR18_Table_S1.csv): *data*: list of strictly anaerobic and aerotolerant genera from Table S1 of [Million and Raoult (2018)](https://doi.org/10.1016/j.humic.2018.07.002)

  - [metaproteome](inst/extdata/orp16S/metaproteome): analysis of metaproteomes

    - [*/mkaa.R](inst/extdata/orp16S/metaproteome): *scripts*: metaproteomes processing
    - [*/aa.csv](inst/extdata/orp16S/metaproteome): *output files*: amino acid composition

  - [R/orp16S.R](R/orp16S.R): code for plots
  - [man/orp16S.Rd](man/orp16S.Rd): manual page
  - [vignettes/orp16S.Rmd](vignettes/orp16S.Rmd): vignette including Figures 1--6, S1--S2, and Table 1

    - [orp16S.html](https://chnosz.net/JMDplots/vignettes/orp16S.html): compiled HTML version of the vignette (external link)

  - [bioRxiv](https://doi.org/10.1101/2021.10.12.464155): preprint (external link)

</details>

 <details>

<summary>Using thermodynamics to obtain geochemical information from genomes (<a href="https://doi.org/10.1111/gbi.12532">2023</a>)</summary>

  - [inst/extdata/utogig](inst/extdata/utogig): scripts and processed data files
  - [R/utogig.R](R/utogig.R): code for plots
  - [man/utogig.Rd](man/utogig.Rd): manual page
  - [vignettes/utogig.Rmd](vignettes/utogig.Rmd): vignette including Figures 1--4, S1--S4, Table S6, and conversions between redox scales

    - [utogig.html](https://chnosz.net/JMDplots/vignettes/utogig.html): compiled HTML version of the vignette (external link)

</details>

## Installation

First install the **remotes** packages from CRAN.

```R
install.packages("remotes")
```

Then install other required packages: [**canprot**](https://github.com/jedick/canprot) and [**chem16S**](https://github.com/jedick/chem16S).

```R
remotes::install_github("jedick/canprot")
remotes::install_github("jedick/chem16S")
```

> **Note**
>
> Currently (as of 2023-07-31), **JMDplots** depends on the development versions of
> **canprot** and **chem16S** from GitHub, not the released versions on CRAN.

Finally, install **JMDplots**.
This command will install prebuilt vignettes; they might not be up-to-date with the source code.

```R
remotes::install_github("jedick/JMDplots")
```

To view the plots, use the R help browser or this command to open the vignettes page:

```R
browseVignettes("JMDplots")
```

## Building vignettes

It might be possible to build the vignettes without [pandoc](https://pandoc.org/), but having pandoc available will make them look better.

```R
remotes::install_github("jedick/JMDplots", dependencies = TRUE, build_vignettes = TRUE)
```

## Licenses

This package except for the file `inst/extdata/orp16S/metadata/PCL+18.csv` is licensed under the GNU General Public License v3 (GPLv3).

The ORP (mV), DO (mg/L) and Feature (Stream, Spring, Lake, Terrace, or Geyser) data for New Zealand hot springs ([Power et al., 2018](https://doi.org/10.1038/s41467-018-05020-y)) in `PCL+18.csv` were obtained from the [1000 Springs Project](https://1000springs.org.nz) and are licensed under CC-BY-NC-SA.

This package contains a copy of the `dunnTest()` function by Derek H. Ogle from CRAN package [FSA](https://cran.r-project.org/package=FSA), version 0.9.3 (License: GPL (>= 2)), which itself is a wrapper for `dunn.test()` from CRAN package [dunn.test](https://cran.r-project.org/package=dunn.test) by Alexis Dinno.
