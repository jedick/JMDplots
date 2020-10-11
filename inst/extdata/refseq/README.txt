# This document describes steps for processing the 
# RefSeq database (release 201, 2020-07-13) to produce these files:

protein_refseq.csv: total amino acid composition of all proteins for
  bacteria, archaea, and viral genomes in the RefSeq collection (n = 42787)
taxid_names.csv: taxid, phylum name and species name for 42787 taxa

# These functions/scripts have the following purpose (output files listed in parentheses):
gencat.sh - extract accession, taxid, sequence length from RefSeq release catalog (accession.taxid.txt)
protein.refseq.R - get average amino acid composition for each taxid from gzipped sequence files (protein_refseq.csv)
taxid.names.R - get taxonomic names for each taxid represented (taxid_names.csv)

## Download stuff

1. Download catalog files
wget ftp://ftp.ncbi.nih.gov/refseq/release/release-catalog/RefSeq-release201.catalog.gz (1270585 KB)
wget ftp://ftp.ncbi.nih.gov/refseq/release/release-catalog/release201.MultispeciesAutonomousProtein2taxname.gz (174586 KB)
wget ftp://ftp.ncbi.nih.gov/refseq/release/release-catalog/release201.files.installed (1234 KB)

2. List URLS for the microbial protein sequence files
grep bacteria.*faa release201.files.installed | awk '!($1="")' | sed -e "s/^\ /ftp\:\/\/ftp.ncbi.nih.gov\/refseq\/release\/bacteria\//g" > urllist
grep archaea.*faa release201.files.installed | awk '!($1="")' | sed -e "s/^\ /ftp\:\/\/ftp.ncbi.nih.gov\/refseq\/release\/archaea\//g" >> urllist
grep viral.*faa release201.files.installed | awk '!($1="")' | sed -e "s/^\ /ftp\:\/\/ftp.ncbi.nih.gov\/refseq\/release\/viral\//g" >> urllist

3. Download the sequence files; use -N to skip existing files [1113 files, 25 GB]
wget -N -i urllist

4. Move the *.faa.gz files to a directory named 'protein'

## Protein stuff

# Timings were made on a 2018 Intel i7 laptop
5. gzip -d release201.MultispeciesAutonomousProtein2taxname.gz [2.4 GB, needed for finding all taxids for WP accessions]

6. Use 'gencat.sh' to generate accession.taxid.txt for bacteria, archaea, and viral proteins in the catalog [11 minutes]
   for RefSeq201, 'wc -l accession.taxid.txt' is 148851653

7. Generate protein_refseq.csv
   mkdir csv
   Then, run these R commands in the working directory that contains accession.taxid.txt, protein/, and csv/
   > source("protein.refseq.R")
   > read.allfiles()    # ca. 14 hours (8 cores)
   > protein.refseq()   # 4 minutes

## Taxonomy stuff

8. Edit 'taxid.names.R' so that 'taxdir' points to the directory where the files
    'names.dmp' and 'nodes.dmp' are present. These files can be downloaded from
    ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz (accessed on 2020-07-15)

9. Source 'taxid.names.R' to generate the file 'taxid_names.csv' [~7.5 hours]
