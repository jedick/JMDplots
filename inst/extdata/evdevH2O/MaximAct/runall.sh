#!/bin/sh

# Script to create all the CSV files in this directory
# Use Rscript to avoid running out of memory in a single R process 20210712
# Add utils to --default-packages 20211228

# NOTE: Generated PNG files are deleted to reduce package size

export NSEED=100

# Trigos phylostrata
Rscript --vanilla --default-packages=JMDplots,utils -e "runMaximAct('TPPG17', bg_organism = 'Hsa', seed=1:$NSEED)"
Rscript --vanilla --default-packages=JMDplots,utils -e "runMaximAct('TPPG17', bg_organism = 'Dme', seed=1:$NSEED)"
Rscript --vanilla --default-packages=JMDplots,utils -e "runMaximAct('TPPG17', bg_organism = 'Bsu', seed=1:$NSEED)"

# Liebeskind gene ages
Rscript --vanilla --default-packages=JMDplots,utils -e "runMaximAct('LMM16', bg_organism = 'Hsa', seed=1:$NSEED)"

# B. subtilis biofilm
Rscript --vanilla --default-packages=JMDplots,utils -e "runMaximAct('transcriptome', bg_organism = 'Bsu', seed=1:$NSEED)"
Rscript --vanilla --default-packages=JMDplots,utils -e "runMaximAct('proteome', bg_organism = 'Bsu', seed=1:$NSEED)"

# Drosophila embryo or adult
Rscript --vanilla --default-packages=JMDplots,utils -e "runMaximAct('fly_embryo', bg_organism = 'Dme', seed=1:$NSEED)"
Rscript --vanilla --default-packages=JMDplots,utils -e "runMaximAct('fly_adult', bg_organism = 'Dme', seed=1:$NSEED)"

# Drosophila developmental time course
Rscript --vanilla --default-packages=JMDplots,utils -e "runMaximAct('fly_development', bg_organism = 'Dme', seed=1:$NSEED)"
