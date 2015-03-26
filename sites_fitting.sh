#!/bin/bash

# SITES_FOLDER, FITTING_FOLDER and MOTIF_STATISTICS_FOLDER should be specified in environment variables
cd "$(dirname "$0")"

mkdir -p  $FITTING_FOLDER  ${MOTIF_STATISTICS_FOLDER}/fitting_log  ${MOTIF_STATISTICS_FOLDER}/slices  ${MOTIF_STATISTICS_FOLDER}/full

for CONTEXT in ${CONTEXTS}; do
  mkdir -p  ${FITTING_FOLDER}/${CONTEXT}  ${MOTIF_STATISTICS_FOLDER}/fitting_log/${CONTEXT}/  ${MOTIF_STATISTICS_FOLDER}/slices/${CONTEXT}  ${MOTIF_STATISTICS_FOLDER}/full/${CONTEXT}

  for RANDOM_VARIANT  in  ${RANDOM_VARIANTS}; do
    mkdir -p  ${FITTING_FOLDER}/${CONTEXT}  ${MOTIF_STATISTICS_FOLDER}/fitting_log/${CONTEXT}
    ruby fitting_random_sites.rb  ${SITES_FOLDER}/${CONTEXT}/sites_cancer.txt  \
                                  ${SITES_FOLDER}/${CONTEXT}/sites_${RANDOM_VARIANT}.txt  \
                                  ${SNV_FOLDER}/SNV_infos_cancer.txt  \
                                  ${SNV_FOLDER}/SNV_infos_${RANDOM_VARIANT}.txt  \
                                  --fold $FITTING_FOLD  \
                                  >  ${FITTING_FOLDER}/${CONTEXT}/sites_${RANDOM_VARIANT}.txt  \
                                  2>  ${MOTIF_STATISTICS_FOLDER}/fitting_log/${CONTEXT}/${RANDOM_VARIANT}.log
  done
  ln -f  ${SITES_FOLDER}/${CONTEXT}/sites_cancer.txt  ${FITTING_FOLDER}/${CONTEXT}/sites_cancer.txt
done
