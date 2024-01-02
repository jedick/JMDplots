[![DOI](https://zenodo.org/badge/211601502.svg)](https://zenodo.org/badge/latestdoi/211601502)

# JMDplots

This R package has code and data for papers by [Jeffrey M. Dick](https://chnosz.net/jeff/).
Plots from the papers are reproduced in the vignettes, which are installed with the package and can be viewed at <https://chnosz.net/JMDplots/vignettes/>.

## Analysis scripts and data files

Click on the paper titles for a list of files.
Published papers are indicated by the year with a DOI link.
Links to preprints, if available, are at the end of each list.
See the manual page associated with each paper for additional details about scripts, data files, and plotting functions.

<!-- Put a space before <details> to make ghostwriter format the lists correctly -->
 <details>

<summary>Water and oxygen in the genomic adaptation of human microbiota (<i>submitted manuscript</i>)</summary>

- [inst/extdata/microhum](inst/extdata/microhum): scripts and processed data files

  - [ARAST](inst/extdata/microhum/ARAST): analysis of metagenomes

    - [ARAST.R](inst/extdata/microhum/ARAST/ARAST.R): *script*: metagenome processing pipeline
    - [runARAST.R](inst/extdata/microhum/ARAST/runARAST.R): *script*: run pipeline for particular metagenomes
    - [*_aa.csv](inst/extdata/microhum/ARAST/): *output files*: amino acid composition
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

  - [MR18_Table_S1_modified.csv](inst/extdata/microhum/MR18_Table_S1_modified.csv): *data*: List of Prokaryotes according to their Aerotolerant or Obligate Anaerobic Metabolism, modified from [Million and Raoult (2018)](https://doi.org/10.1016/j.humic.2018.07.002)
  - [Figure_5_genera.txt](inst/extdata/microhum/Figure_5_genera.txt): *data*: List of genera in Figure 5, created from the value invisibly returned by `microhum_5()`.

- [R/microhum.R](R/microhum.R): code for plots
- [man/microhum.Rd](man/microhum.Rd): manual page
- [vignettes/microhum.Rmd](vignettes/microhum.Rmd): vignette including Figures 1&ndash;4 and S1&ndash;S3

  - [microhum.html](https://chnosz.net/JMDplots/vignettes/microhum.html): compiled HTML version of the vignette (external link)

- [bioRxiv](https://doi.org/10.1101/2023.02.12.528246): preprint (external link)

</details>

 <details>

<summary><i>chem16S</i>: community-level chemical metrics for exploring genomic adaptation to environments (<a href="https://doi.org/10.1093/bioinformatics/btad564">2023</a>)</summary>

  - [R/chem16S.R](R/chem16S.R): code for plots
  - [man/chem16S.Rd](man/chem16S.Rd): manual page
  - [vignettes/chem16S.Rmd](vignettes/chem16S.Rmd): vignette including Figure 1

    - [chem16S.html](https://chnosz.net/JMDplots/vignettes/chem16S.html): compiled HTML version of the vignette (external link)

  - [../chem16S/inst/extdata](https://github.com/jedick/chem16S/blob/main/inst/extdata): scripts and processed data files (*NOTE: these files are in the chem16S package; see [chem16S-package.Rd](https://github.com/jedick/chem16S/blob/main/man/chem16S-package.Rd) for details*)

    - [RefSeq](https://github.com/jedick/chem16S/blob/main/inst/extdata/RefSeq): processing scripts and output files of amino acid composition of genus- and higher-level taxa derived from the [RefSeq database](https://www.ncbi.nlm.nih.gov/refseq/)
    - [GTDB](https://github.com/jedick/chem16S/blob/main/inst/extdata/GTDB): processing scripts and output files of amino acid composition of genus- and higher-level taxa derived from the [Genome Taxonomy Database (GTDB)](https://gtdb.ecogenomic.org/)
    - [metadata](https://github.com/jedick/chem16S/blob/main/inst/extdata/metadata): sample metadata for 16S rRNA datasets: Heart Lake Geyser Basin in Yellowstone National Park ([Bowen De León et al., 2012](https://doi.org/10.3389/fmicb.2013.00330)), Baltic Sea ([Herlemann et al., 2016](https://doi.org/10.3389/fmicb.2016.01883)), and Bison Pool in Yellowstone National Park ([Swingley et al., 2012](https://doi.org/10.1371/journal.pone.0038108))
    - [RDP](https://github.com/jedick/chem16S/blob/main/inst/extdata/RDP): output of RDP Classifier for the above 16S rRNA datasets using the default training set
    - [RDP-GTDB](https://github.com/jedick/chem16S/blob/main/inst/extdata/RDP-GTDB): output of RDP Classifier for the above 16S rRNA datasets using a [GTDB-based training set](https://doi.org/10.5281/zenodo.7633100)
    - [DADA2](https://github.com/jedick/chem16S/blob/main/inst/extdata/DADA2): Analysis of two 16S rRNA datasets with [DADA2](https://doi.org/10.18129/B9.bioc.dada2) using a [GTDB-based training set](https://doi.org/10.5281/zenodo.6655692): marine sediment from the Humboldt Sulfuretum ([Fonseca et al., 2022](https://doi.org/10.3389/fmicb.2022.1016418)) and hot springs in the Qinghai-Tibet Plateau ([Zhang et al., 2023](https://doi.org/10.3389/fmicb.2022.994179))

</details>

 <details>

<summary>Community- and genome-based evidence for a shaping influence of redox potential on bacterial protein evolution (<a href="https://doi.org/10.1128/msystems.00014-23">2023</a>)</summary>

  - [inst/extdata/orp16S](inst/extdata/orp16S): scripts and processed data files

    - [metadata](inst/extdata/orp16S/metadata): *data*: sample metadata for 16S rRNA datasets
    - [pipeline.R](inst/extdata/orp16S/pipeline.R): *script*: 16S rRNA processing pipeline
    - [RDP](inst/extdata/orp16S/RDP): *output files*: taxonomic classifications for 16S rRNA datasets made using the RDP Classifier with its default training set
    - [hydro_p](inst/extdata/orp16S/hydro_p): *data*: shapefiles for the North American Great Lakes, downloaded from [USGS (2010)](https://www.sciencebase.gov/catalog/item/530f8a0ee4b0e7e46bd300dd)
    - [EZdat.csv](inst/extdata/orp16S/EZdat.csv): *output file*: sample data and computed values of Eh7 and *Z*<sub>c</sub>
    - [EZlm.csv](inst/extdata/orp16S/EZlm.csv): *output file*: linear fits between Eh7 and *Z*<sub>c</sub> for each dataset
    - [BKM60.csv](inst/extdata/orp16S/BKM60.csv): *data*: outline of Eh-pH range of natural environments, digitized from Fig. 32 of [Baas Becking et al. (1960)](https://doi.org/10.1086/626659)
    - [MR18_Table_S1.csv](inst/extdata/orp16S/MR18_Table_S1.csv): *data*: list of strictly anaerobic and aerotolerant genera from Table S1 of [Million and Raoult (2018)](https://doi.org/10.1016/j.humic.2018.07.002)

  - [metaproteome](inst/extdata/orp16S/metaproteome): analysis of metaproteomes

    - [*/mkaa.R](inst/extdata/orp16S/metaproteome): *scripts*: metaproteomes processing
    - [*/aa.csv](inst/extdata/orp16S/metaproteome): *output files*: amino acid composition

  - [R/orp16S.R](R/orp16S.R): code for plots
  - [man/orp16S.Rd](man/orp16S.Rd): manual page
  - [vignettes/orp16S.Rmd](vignettes/orp16S.Rmd): vignette including Figures 1&ndash;6, S1&ndash;S2, and Table 1

    - [orp16S.html](https://chnosz.net/JMDplots/vignettes/orp16S.html): compiled HTML version of the vignette (external link)

  - [bioRxiv](https://doi.org/10.1101/2021.10.12.464155): preprint (external link)

</details>

 <details>

<summary>Using thermodynamics to obtain geochemical information from genomes (<a href="https://doi.org/10.1111/gbi.12532">2023</a>)</summary>

  - [inst/extdata/utogig](inst/extdata/utogig): scripts and processed data files
  - [R/utogig.R](R/utogig.R): code for plots
  - [man/utogig.Rd](man/utogig.Rd): manual page
  - [vignettes/utogig.Rmd](vignettes/utogig.Rmd): vignette including Figures 1&ndash;4, S1&ndash;S4, Table S6, and conversions between redox scales

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

```R
remotes::install_github("jedick/JMDplots", dependencies = TRUE, build_vignettes = TRUE)
```

> **Note**
> It might be possible to build the vignettes without [pandoc](https://pandoc.org/),
> but having pandoc available will make them look better.

## Licenses

This package except for the file `inst/extdata/orp16S/metadata/PCL+18.csv` is licensed under the GNU General Public License v3 (GPLv3).

The ORP (mV), DO (mg/L) and Feature (Stream, Spring, Lake, Terrace, or Geyser) data for New Zealand hot springs ([Power et al., 2018](https://doi.org/10.1038/s41467-018-05020-y)) in `PCL+18.csv` were obtained from the [1000 Springs Project](https://1000springs.org.nz) and are licensed under CC-BY-NC-SA.

This package contains a copy of the `dunnTest()` function by Derek H. Ogle from CRAN package [FSA](https://cran.r-project.org/package=FSA), version 0.9.3 (License: GPL (>= 2)), which itself is a wrapper for `dunn.test()` from CRAN package [dunn.test](https://cran.r-project.org/package=dunn.test) by Alexis Dinno.
