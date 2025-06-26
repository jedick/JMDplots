[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3544910.svg)](https://doi.org/10.5281/zenodo.3544910)

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

<summary><code>genoGOE</code>: Evolutionary oxidation of proteins (<i>in-preparation manuscript</i>)</summary>

- [inst/extdata/genoGOE](inst/extdata/genoGOE): scripts and processed data files

  - [GTDB](inst/extdata/genoGOE/GTDB): Processing data from GTDB with additional features

    - [methanogen_genomes.csv](inst/extdata/genoGOE/GTDB/methanogen_genomes.csv): genome IDs and taxonomy in GTDB r220 for 19 Class I and 19 Class II methanogen species selected from Fig. 1 of [Lyu and Lu (2018)](https://doi.org/10.1038/ismej.2017.173)
    - [process_GTDB.R](inst/extdata/genoGOE/GTDB/process_GTDB.R): script to obtain DNA and protein sequences for 53 archaeal marker genes in GTDB and amino acid compositions for all proteins in 38 methanogen genomes
    - [ar53_msa_marker_info_r220_XHZ+06.csv](inst/extdata/genoGOE/GTDB/ar53_msa_marker_info_r220_XHZ+06.csv): list of archaeal marker genes from GTDB with protein abundance information for *Methanococcus maripaludis* added from [Xia et al. (2006)](https://doi.org/10.1074/mcp.M500369-MCP200)

  - [methanogen](inst/extdata/genoGOE/methanogen): sequences and amino acid compositions generated using `process_GTDB.R`

    - [marker/fna](inst/extdata/genoGOE/methanogen/marker/fna): Nucleotide sequences of marker genes
    - [marker/faa](inst/extdata/genoGOE/methanogen/marker/faa): Amino acid sequences of marker genes
    - [aa](inst/extdata/genoGOE/methanogen/aa): Amino acid compositions of all proteins

  - [MCK+23](inst/extdata/genoGOE/MCK+23): Data for genomes with sulfur-cycling genes from [Mateos et al. (2023)](https://doi.org/10.1126/sciadv.ade4847)

    - [get_genomes.R](inst/extdata/genoGOE/MCK+23/get_genomes.R): Script to identify genomes from bootstrap files, download genomes from NCBI, and make the following CSV files
    - [genomes.csv](inst/extdata/genoGOE/MCK+23/genomes.csv): Table of genomes with binary labels denoting presence of specific S-cycling gene
    - [genomes_aa.csv](inst/extdata/genoGOE/MCK+23/genomes_aa.csv): Amino acid composition for each genome with availability in NCBI
    - [sulfur_genomes.xlsx](inst/extdata/genoGOE/MCK+23/sulfur_genomes.xlsx): Spreadsheet listing genomes that exclusively contain one S-cycling gene or gene cluster

  - [GMKK20](inst/extdata/genoGOE/GMKK20): Data for extant and reconstructed ancestral nitrogenase sequences taken from [Garcia et al. (2020)](https://doi.org/10.1111/gbi.12381)
  	- [nitrogenase_aa.csv](inst/extdata/genoGOE/GMKK20/nitrogenase_aa.csv): Amino acid compositions computed from the [Extant-MLAnc_Align.fasta](https://github.com/kacarlab/AncientNitrogenase/blob/master/Extant-MLAnc_Align.fasta) file in the [kacarlab](https://github.com/kacarlab) GitHub repo

  - [PIZ+11](inst/extdata/genoGOE/PIZ+11): Data for reconstructed ancestral thioredoxin sequences derived from [Perez-Jimenez et al. (2011)](https://doi.org/10.1038/nsmb.2020)
  	- [thioredoxin.fasta](inst/extdata/genoGOE/PIZ+11/thioredoxin.fasta): Protein sequences obtained from the RCSB PDB using the accessions listed in the next file.
  	- [DAAD19.csv](inst/extdata/genoGOE/PIZ+11/DAAD19.csv): RCSB PDB IDs and ages for proteins listed by [Del Galdo et al. (2019)](https://doi.org/10.1007/s00239-019-09894-4).

  - [CDY+25](inst/extdata/genoGOE/CDY+25): Data for reconstructed ancestral 3-isopropylmalate dehydrogenases sequences taken from [Cui et al. (2025)](https://doi.org/10.1002/pro.70071)
  	- [IPMDH.fasta](inst/extdata/genoGOE/CDY+25/IPMDH.fasta): Protein sequences from the Supporting Information of Cui et al. (2025).

- [inst/extdata/evdevH2O/LMM16](inst/extdata/evdevH2O/LMM16): scripts and processed data files for consensus gene ages from [Liebeskind et al. (2016)](https://doi.org/10.1093/gbe/evw113), modified from the files used by [Dick (2022)](https://doi.org/10.1007/s00239-022-10051-7)

    - [mkaa.R](inst/extdata/evdevH2O/LMM16/mkaa.R): *script*: sum amino acid compositions of proteins in each gene age category
    - [reference_proteomes.csv](inst/extdata/evdevH2O/LMM16/reference_proteomes.csv): *data*: IDs of UniProt reference proteomes for 31 organisms
    - [modeAges_names.csv](inst/extdata/evdevH2O/LMM16/modeAges_names.csv): *output file*: Names of gene age categories for each organism
    - [modeAges_aa.csv](inst/extdata/evdevH2O/LMM16/modeAges_aa.csv): *output file*: Summed amino acid composition for proteins in each gene age category

- [../canprot/inst/extdata/fasta/KHAB17.fasta](https://github.com/jedick/canprot/blob/main/inst/extdata/fasta/KHAB17.fasta): reconstructed ancestral Rubisco sequences taken from [Kaçar et al. (2017)](https://doi.org/10.1111/gbi.12243)

</details>

 <details>

<summary><code>microhum</code>: Adaptations of microbial genomes to human body chemistry (<i>submitted manuscript</i>)</summary>

- [inst/extdata/microhum](inst/extdata/microhum): scripts and processed data files

  - [ARAST](inst/extdata/microhum/ARAST): analysis of metagenomes

    - [ARAST.R](inst/extdata/microhum/ARAST/ARAST.R): *script*: metagenome processing pipeline
    - [runARAST.R](inst/extdata/microhum/ARAST/runARAST.R): *script*: run pipeline for particular metagenomes
    - [*_aa.csv](inst/extdata/microhum/ARAST/): *output files*: amino acid composition
    - [*_stats.csv](inst/extdata/microhum/ARAST/): *output files*: processing statistics

  - [KWL22](inst/extdata/microhum/KWL22): analysis of metagenome-assembled genomes (MAGs) from [Ke et al. (2022)](https://doi.org/10.1038/s41467-022-32991-w)

    - [mkaa.R](inst/extdata/microhum/KWL22/mkaa.R): *script*: metaproteome processing
    - [KWL22_MAGs_prodigal_aa.csv.xz](inst/extdata/microhum/KWL22/KWL22_MAGs_prodigal_aa.csv.xz): *output file*: amino acid composition
    - [BioSample_metadata.txt](inst/extdata/microhum/KWL22/BioSample_metadata.txt): *data*: BioSample metadata for MAGs obtained from NCBI BioProjects [PRJNA624223](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA624223/) and [PRJNA650244](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA650244/).

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
- [vignettes/microhum.Rmd](vignettes/microhum.Rmd): vignette including Figures 1&ndash;6 and figure supplements.

  - [microhum.html](https://chnosz.net/JMDplots/vignettes/microhum.html): compiled HTML version of the vignette (external link)

- [bioRxiv](https://doi.org/10.1101/2023.02.12.528246): preprint (external link)

</details>

 <details>

<summary><code>chem16S</code>: Community-level chemical metrics for exploring genomic adaptation to environments (<a href="https://doi.org/10.1093/bioinformatics/btad564">2023</a>)</summary>

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

<summary><code>orp16S</code>: Community- and genome-based evidence for a shaping influence of redox potential on bacterial protein evolution (<a href="https://doi.org/10.1128/msystems.00014-23">2023</a>)</summary>

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

<summary><code>utogig</code>: Using thermodynamics to obtain geochemical information from genomes (<a href="https://doi.org/10.1111/gbi.12532">2023</a>)</summary>

  - [inst/extdata/utogig](inst/extdata/utogig): scripts and processed data files
  - [R/utogig.R](R/utogig.R): code for plots
  - [man/utogig.Rd](man/utogig.Rd): manual page
  - [vignettes/utogig.Rmd](vignettes/utogig.Rmd): vignette including Figures 1&ndash;4, S1&ndash;S4, Table S6, and conversions between redox scales

    - [utogig.html](https://chnosz.net/JMDplots/vignettes/utogig.html): compiled HTML version of the vignette (external link)

</details>

## Reference databases

Amino acid compositions and taxonomic information have been obtained from the *Saccharomyces* Genome Database (SGD), UniProt, RefSeq, GTDB, and MGnify.
See [man/JMDplots-package.Rd](man/JMDplots-package.Rd) for further details.

<!-- Put a space before <details> to make ghostwriter format the lists correctly -->
 <details>

<summary>Reference databases</summary>

- [inst/extdata/RefDB/organisms](inst/extdata/RefDB/organisms): Data for particular organisms, downloaded from [SGD](https://www.yeastgenome.org/) or [UniProt](https://www.uniprot.org/).

  - [Sce.csv.xz](inst/extdata/RefDB/organisms/Sce.csv.xz): *Saccharomyces cerevisiae* (used in the [scsc](https://chnosz.net/JMDplots/vignettes/scsc.html) and [aoscp](https://chnosz.net/JMDplots/vignettes/aoscp.html) papers)
  - [yeastgfp.csv.xz](inst/extdata/RefDB/organisms/yeastgfp.csv.xz): Subcellular localization and abundance of proteins in *S. cerevisiae* (used in the [scsc](https://chnosz.net/JMDplots/vignettes/scsc.html) paper)
  - [UP000000805_243232.csv.xz](inst/extdata/RefDB/organisms/UP000000805_243232.csv.xz): *Methanocaldococcus jannaschii* (used in the [mjenergy](https://chnosz.net/JMDplots/vignettes/mjenergy.html) paper)
  - [UP000000625_83333.csv.xz](inst/extdata/RefDB/organisms/UP000000625_83333.csv.xz): *Escherichia coli* K12
  - [UP000000803_7227.csv.xz](inst/extdata/RefDB/organisms/UP000000803_7227.csv.xz): *Drosophila melanogaster* (used in the [evdevH2O](https://chnosz.net/JMDplots/vignettes/evdevH2O.html) paper)
  - [UP000001570_224308.csv.xz](inst/extdata/RefDB/organisms/UP000001570_224308.csv.xz): *Bacillus subtilis* strain 168 (used in the [evdevH2O](https://chnosz.net/JMDplots/vignettes/evdevH2O.html) paper)

- [inst/extdata/RefDB/RefSeq](inst/extdata/RefDB/RefSeq): Data files processed from [RefSeq](https://www.ncbi.nlm.nih.gov/refseq/) and used in the [geo16S](https://chnosz.net/JMDplots/vignettes/geo16S.html) and [orp16S](https://chnosz.net/JMDplots/vignettes/orp16S.html) papers
  - [genome_AA.csv.xz](inst/extdata/RefDB/RefSeq/genome_AA.csv.xz): Amino acid compositions of species-level archaeal, bacterial, and viral taxa in the RefSeq database
  - [taxonomy.csv.xz](inst/extdata/RefDB/RefSeq/taxonomy.csv.xz): Taxonomic names for the species
  - Scripts to produce these files are in [**chem16S**](https://github.com/jedick/chem16S)

- [inst/extdata/RefDB/GTDB](inst/extdata/RefDB/GTDB): Data files processed from [GTDB](https://gtdb.ecogenomic.org/) and used in the [microhum](https://chnosz.net/JMDplots/vignettes/microhum.html) manuscript
  - [genome_AA.csv.xz](inst/extdata/RefDB/GTDB/genome_AA.csv.xz): Amino acid compositions of predicted proteins
  - [taxonomy.csv.xz](inst/extdata/RefDB/GTDB/taxonomy.csv.xz): Taxonomic names
  - Scripts to produce these files are in [**chem16S**](https://github.com/jedick/chem16S)

- [inst/extdata/RefDB/UHGG](inst/extdata/RefDB/UHGG): Data files processed from MGnify's [UHGG](https://www.ebi.ac.uk/metagenomics/genome-catalogues/human-gut-v2-0-1) and used in the [microhum](https://chnosz.net/JMDplots/vignettes/microhum.html) manuscript
  - [MGnify_genomes.csv](inst/extdata/RefDB/UHGG/MGnify_genomes.csv): List of 4744 species-level clusters in the Unified Human Gastrointestinal Genome ([UHGG v.2.0.1](https://www.ebi.ac.uk/metagenomics/genome-catalogues/human-gut-v2-0-1))
  - [getMGnify.R](inst/extdata/RefDB/UHGG/getMGnify.R): Commands used to download FASTA files for proteins and to scrape the website for taxonomic information
  - [taxonomy.csv.xz](inst/extdata/RefDB/UHGG/taxonomy.csv.xz): Taxonomy for 2350 selected genomes with contamination < 2% and completeness > 95%
  - [genome_AA.R](inst/extdata/RefDB/UHGG/genome_AA.R): Calculates amino acid compositions of the selected genomes from FASTA files and writes the output file [genome_AA.csv.xz](inst/extdata/RefDB/UHGG/genome_AA.csv.xz)
  - [taxonomy.R](inst/extdata/RefDB/UHGG/taxonomy.R): Combines amino acid compositions of genomes to generate reference proteomes for genera and higher taxonomic levels and writes the output file [taxonomy.csv.xz](inst/extdata/RefDB/UHGG/taxonomy.csv.xz)
  - [fullset](inst/extdata/RefDB/UHGG/fullset): Versions of `taxonomy.csv.xz`, `genome_AA.csv.xz`, and `taxon_AA.csv.xz` for the full set of 4744 genomes

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
