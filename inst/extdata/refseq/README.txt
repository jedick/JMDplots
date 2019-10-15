# the following data files support calculations using the 
# RefSeq database (release 95, 2019-07-12)
protein_refseq.csv: total amino acid composition of all proteins for
  bacteria, archaea, and viral genomes in the RefSeq collection (n = 36425)
taxid_names.csv: taxid, phylum name and species name for 788 microbial taxa

# these functions/scripts have the following purpose (output files listed in parentheses):
gencat.sh - extract accession, taxid, sequence length from RefSeq release catalog (accession.taxid.txt)
protein.refseq.R - get average amino acid composition for each taxid from gzipped sequence files (protein_refseq.csv)
taxid.names.R - get taxonomic names for each taxid represented (taxid_names.csv)

## download stuff
1. Download catalog files
wget ftp://ftp.ncbi.nih.gov/refseq/release/release-catalog/RefSeq-release95.catalog.gz (1129122 KB)
wget ftp://ftp.ncbi.nih.gov/refseq/release/release-catalog/release95.MultispeciesAutonomousProtein2taxname.gz (147264 KB)
wget ftp://ftp.ncbi.nih.gov/refseq/release/release-catalog/release95.files.installed (1132 KB)

2. List URLS for the microbial protein sequence files
grep bacteria.*faa release95.files.installed | awk '!($1="")' | sed -e "s/^\ /ftp\:\/\/ftp.ncbi.nih.gov\/refseq\/release\/bacteria\//g" > urllist
grep archaea.*faa release95.files.installed | awk '!($1="")' | sed -e "s/^\ /ftp\:\/\/ftp.ncbi.nih.gov\/refseq\/release\/archaea\//g" >> urllist
grep viral.*faa release95.files.installed | awk '!($1="")' | sed -e "s/^\ /ftp\:\/\/ftp.ncbi.nih.gov\/refseq\/release\/viral\//g" >> urllist

3. Download the sequence files; use -N to skip existing files [919 files, 25 GB]
wget -N -i urllist

4. move the *.faa.gz files to a directory named 'protein'

## protein stuff
# timings were made on an Intel i7-8550U laptop
5. gzip -d release95.MultispeciesAutonomousProtein2taxname.gz [1.8 GB, needed for finding all taxids for WP accessions]

6. use 'gencat.sh' to generate accession.taxid.txt for bacteria,archaea,viral proteins in the catalog [7 minutes]
   for RefSeq95, 'wc -l accession.taxid.txt' is 121628948

7. generate protein_refseq.csv
   mkdir csv
   Then, run these R commands in the working directory that contains accession.taxid.txt, protein/, and csv/
   > source("protein.refseq.R")
   > read.allfiles()    # ca. 10 hours (8 cores)
   > protein.refseq()   # 2.5 minutes

## taxonomy stuff
8. edit 'taxid.names.R' so that 'taxdir' points to the directory where the files
    'names.dmp' and 'nodes.dmp' are present. these files can be downloaded from
    ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz (accessed on 2013-09-18)

9. source 'taxid.names.R' to generate the file 'taxid_names.csv' [~5.5 hours]
