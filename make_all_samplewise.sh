#!/bin/bash

# SITES_FOLDER, FITTING_FOLDER and MOTIF_STATISTICS_FOLDER and so on should be specified in environment variables

cd "$(dirname "$0")"

# Choose only subset of cancer types
declare -A CANCER_SAMPLES_BY_TYPE
CANCER_SAMPLES_BY_TYPE=( \
                        [ER_plus_ve_HER2_plus_ve]=PD4194a,PD4198a \
                        [ER_plus_ve_HER2_minus_ve]=PD4120a,PD4103a,PD4088a,PD4085a,PD3851a \
                        [ER_minus_ve_HER2_plus_ve]=PD4199a,PD4192a \
                        [triple_negative]=PD4248a,PD4086a,PD4109a \
                        [BRCA1]=PD4107a,PD3890a,PD3905a,PD4005a,PD4006a \
                        [BRCA2]=PD3904a,PD3945a,PD4115a,PD4116a \
                        [PD3851a]=PD3851a \
                        [PD3890a]=PD3890a \
                        [PD3904a]=PD3904a \
                        [PD3905a]=PD3905a \
                        [PD3945a]=PD3945a \
                        [PD4005a]=PD4005a \
                        [PD4006a]=PD4006a \
                        [PD4085a]=PD4085a \
                        [PD4086a]=PD4086a \
                        [PD4088a]=PD4088a \
                        [PD4103a]=PD4103a \
                        [PD4107a]=PD4107a \
                        [PD4109a]=PD4109a \
                        [PD4115a]=PD4115a \
                        [PD4116a]=PD4116a \
                        [PD4120a]=PD4120a \
                        [PD4192a]=PD4192a \
                        [PD4194a]=PD4194a \
                        [PD4198a]=PD4198a \
                        [PD4199a]=PD4199a \
                        [PD4248a]=PD4248a \
                        [ER_plus_ve_HER2_minus_ve_without_PD4120a]=PD4103a,PD4088a,PD4085a,PD3851a \
                        [all_except_PD4120a]=PD3851a,PD3890a,PD3904a,PD3905a,PD3945a,PD4005a,PD4006a,PD4085a,PD4086a,PD4088a,PD4103a,PD4107a,PD4109a,PD4115a,PD4116a,PD4192a,PD4194a,PD4198a,PD4199a,PD4248a \
                      )

export FITTING_FOLD=1
export FOLD_CHANGE=5
export MIN_FITTING_RATE=0.99

for CANCER_TYPE  in  ${!CANCER_SAMPLES_BY_TYPE[@]}; do
  ./generate_samplewise_sites_and_slices.sh ${CANCER_TYPE} ${CANCER_SAMPLES_BY_TYPE[$CANCER_TYPE]}
end

for CANCER_TYPE  in  ${!CANCER_SAMPLES_BY_TYPE[@]}; do
  ./aggregate_samplewise_statistics.sh ${CANCER_TYPE}
done

for CANCER_TYPE_1  in  ${!CANCER_SAMPLES_BY_TYPE[@]}; do
  for CANCER_TYPE_2  in  ${!CANCER_SAMPLES_BY_TYPE[@]}; do
    ./aggregate_samplewise_paired_statistics.sh ${CANCER_TYPE_1} ${CANCER_TYPE_2}
  done
done
