#!/bin/bash

# FITTING_FOLDER and MOTIF_STATISTICS_FOLDER should be specified in environment variables
cd "$(dirname "$0")"

FOLD_CHANGE=4

for CONTEXT in ${CONTEXTS}; do
  for VARIANT  in  cancer  ${RANDOM_VARIANTS}; do
    mkdir -p  ${MOTIF_STATISTICS_FOLDER}/slices/${CONTEXT}/${VARIANT}
    ./generate_all_motif_statistics.sh  ${FITTING_FOLDER}/${CONTEXT}/sites_${VARIANT}.txt  \
                                        ${MOTIF_STATISTICS_FOLDER}/slices/${CONTEXT}/${VARIANT}  \
                                        0.0005  $FOLD_CHANGE

  done
done
