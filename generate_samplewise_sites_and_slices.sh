#!/bin/bash

# SITES_FOLDER, FITTING_FOLDER and MOTIF_STATISTICS_FOLDER and so on should be specified in environment variables
# Also CANCER_SAMPLES_BY_TYPE, FITTING_FOLD, FOLD_CHANGE, MIN_FITTING_RATE should be specified

cd "$(dirname "$0")"

CANCER_TYPE=$1 # ER_plus_ve_HER2_plus_ve

#####
mkdir -p  $FITTING_FOLDER  ${MOTIF_STATISTICS_FOLDER}/fitting_log/samples/${CANCER_TYPE}  ${MOTIF_STATISTICS_FOLDER}/slices  ${MOTIF_STATISTICS_FOLDER}/full

for CONTEXT in ${CONTEXTS}; do
  mkdir -p  ${FITTING_FOLDER}/${CONTEXT}/samples/${CANCER_TYPE}  ${MOTIF_STATISTICS_FOLDER}/fitting_log/${CONTEXT}/samples/${CANCER_TYPE}  ${MOTIF_STATISTICS_FOLDER}/slices/${CONTEXT}/samples/${CANCER_TYPE}  ${MOTIF_STATISTICS_FOLDER}/full/${CONTEXT}/samples/${CANCER_TYPE}

  for RANDOM_VARIANT  in  ${RANDOM_GENOME_VARIANTS}; do
    mkdir -p  ${FITTING_FOLDER}/${CONTEXT}/samples/${CANCER_TYPE}  ${MOTIF_STATISTICS_FOLDER}/fitting_log/${CONTEXT}/samples/${CANCER_TYPE}
    ruby fitting_random_sites.rb  ${SITES_FOLDER}/${CONTEXT}/samples/${CANCER_TYPE}/sites_cancer.txt  \
                                  ${SITES_FOLDER}/${CONTEXT}/sites_${RANDOM_VARIANT}.txt  \
                                  ${SNV_FOLDER}/SNV_infos_cancer.txt  \
                                  ${SNV_FOLDER}/SNV_infos_${RANDOM_VARIANT}.txt  \
                                  --fold $FITTING_FOLD_GENOME  \
                                  >  ${FITTING_FOLDER}/${CONTEXT}/samples/${CANCER_TYPE}/sites_${RANDOM_VARIANT}.txt  \
                                  2>  ${MOTIF_STATISTICS_FOLDER}/fitting_log/${CONTEXT}/samples/${CANCER_TYPE}/${RANDOM_VARIANT}.log
  done

  # Use shuffle variants generated from the same sample
  for RANDOM_VARIANT  in  ${RANDOM_SHUFFLE_VARIANTS}; do
    mkdir -p  ${FITTING_FOLDER}/${CONTEXT}/samples/${CANCER_TYPE}  ${MOTIF_STATISTICS_FOLDER}/fitting_log/${CONTEXT}/samples/${CANCER_TYPE}
    ruby fitting_random_sites.rb  ${SITES_FOLDER}/${CONTEXT}/samples/${CANCER_TYPE}/sites_cancer.txt  \
                                  ${SITES_FOLDER}/${CONTEXT}/samples/${CANCER_TYPE}/sites_${RANDOM_VARIANT}.txt  \
                                  ${SNV_FOLDER}/SNV_infos_cancer.txt  \
                                  ${SNV_FOLDER}/SNV_infos_${RANDOM_VARIANT}.txt  \
                                  --fold $FITTING_FOLD_SHUFFLE  \
                                  >  ${FITTING_FOLDER}/${CONTEXT}/samples/${CANCER_TYPE}/sites_${RANDOM_VARIANT}.txt  \
                                  2>  ${MOTIF_STATISTICS_FOLDER}/fitting_log/${CONTEXT}/samples/${CANCER_TYPE}/${RANDOM_VARIANT}.log
  done
  ln -f  ${SITES_FOLDER}/${CONTEXT}/samples/${CANCER_TYPE}/sites_cancer.txt  ${FITTING_FOLDER}/${CONTEXT}/samples/${CANCER_TYPE}/sites_cancer.txt
done

#####
for RANDOM_VARIANT  in  ${RANDOM_VARIANTS}; do
  for CONTEXT in ${CONTEXTS}; do
    mkdir -p  ${MOTIF_STATISTICS_FOLDER}/slices/${CONTEXT}/samples/${CANCER_TYPE}/${RANDOM_VARIANT}
    ./generate_all_motif_statistics.sh  ${FITTING_FOLDER}/${CONTEXT}/samples/${CANCER_TYPE}/sites_${RANDOM_VARIANT}.txt  \
                                        ${MOTIF_STATISTICS_FOLDER}/slices/${CONTEXT}/samples/${CANCER_TYPE}/${RANDOM_VARIANT}
  done
done
