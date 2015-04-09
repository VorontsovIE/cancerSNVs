#!/bin/bash

# SITES_FOLDER, FITTING_FOLDER and MOTIF_STATISTICS_FOLDER and so on should be specified in environment variables
# Also CANCER_SAMPLES_BY_TYPE, FITTING_FOLD, FOLD_CHANGE, MIN_FITTING_RATE should be specified

cd "$(dirname "$0")"

CANCER_TYPE=$1 # ER_plus_ve_HER2_plus_ve
CANCER_SAMPLES=$2 # PD4194a,PD4198a

declare -A CONTEXT_DECLARATIONS # All except any
CONTEXT_DECLARATIONS=( \
                        [tpc]="--contexts TCN  --mutation-types promoter,intronic" \
                        [cpg]="--contexts NCG  --mutation-types promoter,intronic" \
                        [non_tpc]="--invert-context-request  --contexts TCN  --mutation-types promoter,intronic" \
                      )

for CONTEXT  in  ${CONTEXTS}; do
  mkdir -p  ${SITES_FOLDER}/${CONTEXT}/samples
done

ruby bin/preparations/filter_mutations.rb ${SNV_FOLDER}/SNV_infos_cancer.txt \
                                          ${SITES_FOLDER}/any/sites_cancer.txt \
                                          --cancer-samples ${CANCER_SAMPLES} \
                                          >  ${SITES_FOLDER}/any/samples/sites_cancer_${CANCER_TYPE}.txt

for CONTEXT in ${CONTEXTS}; do
  if [[ "$CONTEXT" != "any" ]]; then
    ruby bin/preparations/filter_mutations.rb   ${SNV_FOLDER}/SNV_infos_cancer.txt  \
                                                ${SITES_FOLDER}/any/samples/sites_cancer_${CANCER_TYPE}.txt  \
                                                ${CONTEXT_DECLARATIONS[$CONTEXT]} \
                                                >  ${SITES_FOLDER}/${CONTEXT}/samples/sites_cancer_${CANCER_TYPE}.txt
  fi
done

#####
mkdir -p  $FITTING_FOLDER  ${MOTIF_STATISTICS_FOLDER}/fitting_log/  ${MOTIF_STATISTICS_FOLDER}/slices/samples  ${MOTIF_STATISTICS_FOLDER}/full/samples

for CONTEXT in ${CONTEXTS}; do
  mkdir -p  ${FITTING_FOLDER}/${CONTEXT}/samples/${CANCER_TYPE}  ${MOTIF_STATISTICS_FOLDER}/fitting_log/${CONTEXT}/samples/${CANCER_TYPE}  ${MOTIF_STATISTICS_FOLDER}/slices/${CONTEXT}  ${MOTIF_STATISTICS_FOLDER}/full/${CONTEXT}/samples

  for RANDOM_VARIANT  in  ${RANDOM_VARIANTS}; do
    mkdir -p  ${MOTIF_STATISTICS_FOLDER}/fitting_log/${CONTEXT}/samples

    ruby fitting_random_sites.rb  ${SITES_FOLDER}/${CONTEXT}/samples/sites_cancer_${CANCER_TYPE}.txt  \
                                  ${SITES_FOLDER}/${CONTEXT}/sites_${RANDOM_VARIANT}.txt  \
                                  ${SNV_FOLDER}/SNV_infos_cancer.txt  \
                                  ${SNV_FOLDER}/SNV_infos_${RANDOM_VARIANT}.txt  \
                                  --fold $FITTING_FOLD  \
                                  >  ${FITTING_FOLDER}/${CONTEXT}/samples/${CANCER_TYPE}/sites_${RANDOM_VARIANT}.txt  \
                                  2>  ${MOTIF_STATISTICS_FOLDER}/fitting_log/${CONTEXT}/samples/${CANCER_TYPE}/${RANDOM_VARIANT}.log

    ln -f  ${SITES_FOLDER}/${CONTEXT}/samples/sites_cancer_${CANCER_TYPE}.txt  ${FITTING_FOLDER}/${CONTEXT}/samples/${CANCER_TYPE}/sites_cancer.txt
  done
done

#####
for CONTEXT in ${CONTEXTS}; do
  for VARIANT  in  cancer  ${RANDOM_VARIANTS}; do
    mkdir -p  ${MOTIF_STATISTICS_FOLDER}/slices/${CONTEXT}/samples/${CANCER_TYPE}/${VARIANT}
    ./generate_all_motif_statistics.sh  ${FITTING_FOLDER}/${CONTEXT}/samples/${CANCER_TYPE}/sites_${VARIANT}.txt  \
                                        ${MOTIF_STATISTICS_FOLDER}/slices/${CONTEXT}/samples/${CANCER_TYPE}/${VARIANT}  \
                                        0.0005  $FOLD_CHANGE

  done
done
