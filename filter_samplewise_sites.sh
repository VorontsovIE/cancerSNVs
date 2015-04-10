#!/bin/bash

# SITES_FOLDER, FITTING_FOLDER and MOTIF_STATISTICS_FOLDER and so on should be specified in environment variables
# Also CANCER_SAMPLES_BY_TYPE, FITTING_FOLD, FOLD_CHANGE, MIN_FITTING_RATE should be specified

cd "$(dirname "$0")"

CANCER_TYPE=$1 # ER_plus_ve_HER2_plus_ve
CANCER_SAMPLES=$2 # PD4194a,PD4198a
VARIANT=$3 # random_genome_123

declare -A CONTEXT_DECLARATIONS # All except any
CONTEXT_DECLARATIONS=( \
                        [tpc]="--contexts TCN  --mutation-types promoter,intronic" \
                        [cpg]="--contexts NCG  --mutation-types promoter,intronic" \
                        [non_tpc]="--invert-context-request  --contexts TCN  --mutation-types promoter,intronic" \
                      )

for CONTEXT  in  ${CONTEXTS}; do
  mkdir -p  ${SITES_FOLDER}/${CONTEXT}/samples
done

ruby bin/preparations/filter_mutations.rb ${SNV_FOLDER}/SNV_infos_${VARIANT}.txt \
                                          ${SITES_FOLDER}/any/sites_${VARIANT}.txt \
                                          --cancer-samples ${CANCER_SAMPLES} \
                                          >  ${SITES_FOLDER}/any/samples/${CANCER_TYPE}/sites_${VARIANT}.txt

for CONTEXT in ${CONTEXTS}; do
  if [[ "$CONTEXT" != "any" ]]; then
    ruby bin/preparations/filter_mutations.rb   ${SNV_FOLDER}/SNV_infos_${VARIANT}.txt  \
                                                ${SITES_FOLDER}/any/samples/${CANCER_TYPE}/sites_${VARIANT}.txt  \
                                                ${CONTEXT_DECLARATIONS[$CONTEXT]} \
                                                >  ${SITES_FOLDER}/${CONTEXT}/samples/${CANCER_TYPE}/sites_${VARIANT}.txt
  fi
done
