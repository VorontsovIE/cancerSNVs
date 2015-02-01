#!/bin/bash

# SITES_FOLDER, FITTING_FOLDER and MOTIF_STATISTICS_FOLDER should be specified in environment variables
cd "$(dirname "$0")"

FITTING_FOLD=1

mkdir -p  $FITTING_FOLDER  ${MOTIF_STATISTICS_FOLDER}/fitting_log  ${MOTIF_STATISTICS_FOLDER}/slices  ${MOTIF_STATISTICS_FOLDER}/full

for CONTEXT in any cpg tpc; do
  mkdir -p  ${FITTING_FOLDER}/${CONTEXT}  ${MOTIF_STATISTICS_FOLDER}/fitting_log/${CONTEXT}/  ${MOTIF_STATISTICS_FOLDER}/slices/${CONTEXT}  ${MOTIF_STATISTICS_FOLDER}/full/${CONTEXT}

  for RANDOM_VARIANT  in  random_genome_13  random_genome_15  random_genome_17  random_shuffle_135  random_shuffle_137  random_shuffle_139; do
    mkdir -p  ${FITTING_FOLDER}/${CONTEXT}  ${MOTIF_STATISTICS_FOLDER}/fitting_log/${CONTEXT}
    ruby fitting_random_sites.rb  ${SITES_FOLDER}/${CONTEXT}/sites_cancer.txt  \
                                  ${SITES_FOLDER}/${CONTEXT}/sites_${RANDOM_VARIANT}.txt  \
                                  --fold $FITTING_FOLD  \
                                  >  ${FITTING_FOLDER}/${CONTEXT}/sites_${RANDOM_VARIANT}.txt  \
                                  2>  ${MOTIF_STATISTICS_FOLDER}/fitting_log/${CONTEXT}/${RANDOM_VARIANT}.log
  done
  ln -f  ${SITES_FOLDER}/${CONTEXT}/sites_cancer.txt  ${FITTING_FOLDER}/${CONTEXT}/sites_cancer.txt
done
