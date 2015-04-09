#!/bin/bash

# FITTING_FOLDER and MOTIF_STATISTICS_FOLDER should be specified in environment variables
cd "$(dirname "$0")"

for VARIANT  in  cancer  ${RANDOM_VARIANTS}; do
  for CONTEXT in ${CONTEXTS}; do
    mkdir -p  ${MOTIF_STATISTICS_FOLDER}/slices/${CONTEXT}/${VARIANT}
    ./generate_all_motif_statistics.sh  ${FITTING_FOLDER}/${CONTEXT}/sites_${VARIANT}.txt  \
                                        ${MOTIF_STATISTICS_FOLDER}/slices/${CONTEXT}/${VARIANT}
  done
done
