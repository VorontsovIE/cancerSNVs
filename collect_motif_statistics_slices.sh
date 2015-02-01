#!/bin/bash

# FITTING_FOLDER and MOTIF_STATISTICS_FOLDER should be specified in environment variables
cd "$(dirname "$0")"

FOLD_CHANGE=5

for CONTEXT in any cpg tpc; do
  for VARIANT  in  cancer  random_genome_13  random_genome_15  random_genome_17  random_shuffle_135  random_shuffle_137  random_shuffle_139; do
    mkdir -p  ${MOTIF_STATISTICS_FOLDER}/slices/${CONTEXT}/${VARIANT}
    ./generate_all_motif_statistics.sh  ${FITTING_FOLDER}/${CONTEXT}/sites_${VARIANT}.txt  \
                                        ${MOTIF_STATISTICS_FOLDER}/slices/${CONTEXT}/${VARIANT}  \
                                        0.0005  $FOLD_CHANGE

  done
done
