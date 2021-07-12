#!/bin/sh

# Script to create all the CSV files in this directory
# Use Rscript to avoid running out of memory in a single R process 20210712

export NSEED=100

# Trigos phylostrata
Rscript --vanilla --options --default-packages=JMDplots -e "runMaximAct('TPPG17', bg_organism = 'Hsa', seed=1:$NSEED)"
Rscript --vanilla --options --default-packages=JMDplots -e "runMaximAct('TPPG17', bg_organism = 'Sce', seed=1:$NSEED)"
Rscript --vanilla --options --default-packages=JMDplots -e "runMaximAct('TPPG17', bg_organism = 'Eco', seed=1:$NSEED)"
Rscript --vanilla --options --default-packages=JMDplots -e "runMaximAct('TPPG17', bg_organism = 'Mja', seed=1:$NSEED)"

# Liebeskind gene ages
Rscript --vanilla --options --default-packages=JMDplots -e "runMaximAct('LMM16', bg_organism = 'Hsa', seed=1:$NSEED)"

# B. subtilis biofilm
Rscript --vanilla --options --default-packages=JMDplots -e "runMaximAct('transcriptome', bg_organism = 'Hsa', seed=1:$NSEED)"
Rscript --vanilla --options --default-packages=JMDplots -e "runMaximAct('proteome', bg_organism = 'Hsa', seed=1:$NSEED)"

# Drosophila embryo or adult
Rscript --vanilla --options --default-packages=JMDplots -e "runMaximAct('fly_embryo', bg_organism = 'Hsa', seed=1:$NSEED)"
Rscript --vanilla --options --default-packages=JMDplots -e "runMaximAct('fly_adult', bg_organism = 'Hsa', seed=1:$NSEED)"

# Drosophila developmental time course
Rscript --vanilla --options --default-packages=JMDplots -e "runMaximAct('fly_development', bg_organism = 'Hsa', seed=1:$NSEED)"
