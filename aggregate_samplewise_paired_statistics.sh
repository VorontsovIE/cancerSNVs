#!/bin/bash

cd "$(dirname "$0")"

CANCER_TYPE_1=$1
CANCER_TYPE_2=$2

for CONTEXT in ${CONTEXTS}; do
  mkdir -p  ${MOTIF_STATISTICS_FOLDER}/full/${CONTEXT}/samples_vs_samples/${CANCER_TYPE_1}/

  ruby summary.rb   ${MOTIF_STATISTICS_FOLDER}/slices/${CONTEXT}/samples/${CANCER_TYPE_1}/cancer  \
                    ${MOTIF_STATISTICS_FOLDER}/slices/${CONTEXT}/samples/${CANCER_TYPE_2}/cancer  \
                    ./source_data/motif_names.txt  ./source_data/hocomoco_genes_infos.csv  \
                    --correction ${CORRECTION_METHOD}  \
                    >  ${MOTIF_STATISTICS_FOLDER}/full/${CONTEXT}/samples_vs_samples/${CANCER_TYPE_1}/${CANCER_TYPE_2}.csv

  for SUBJECTED_OR_PROTECTED  in  subjected  protected; do
    for DISRUPTION_OR_EMERGENCE  in  disruption  emergence; do
      mkdir -p ${MOTIF_STATISTICS_FOLDER}/filtered/${SUBJECTED_OR_PROTECTED}/${DISRUPTION_OR_EMERGENCE}/${CONTEXT}/samples_vs_samples/${CANCER_TYPE_1}
      ruby filter_summary.rb  ${MOTIF_STATISTICS_FOLDER}/full/${CONTEXT}/samples_vs_samples/${CANCER_TYPE_1}/${CANCER_TYPE_2}.csv  \
                              --motif-qualities A,B,C,D  --significance 0.05  \
                              --${DISRUPTION_OR_EMERGENCE}  --${SUBJECTED_OR_PROTECTED}  \
                              >  ${MOTIF_STATISTICS_FOLDER}/filtered/${SUBJECTED_OR_PROTECTED}/${DISRUPTION_OR_EMERGENCE}/${CONTEXT}/samples_vs_samples/${CANCER_TYPE_1}/${CANCER_TYPE_2}.csv
    done
  done
done