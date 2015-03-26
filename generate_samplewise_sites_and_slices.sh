#!/bin/bash

# SITES_FOLDER, FITTING_FOLDER and MOTIF_STATISTICS_FOLDER and so on should be specified in environment variables
# Also CANCER_SAMPLES_BY_TYPE, FITTING_FOLD, FOLD_CHANGE, MIN_FITTING_RATE should be specified

cd "$(dirname "$0")"

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

declare -A CONTEXT_DECLARATIONS # All except any
CONTEXT_DECLARATIONS=( \
                        [tpc]="--contexts TCN  --mutation-types promoter,intronic" \
                        [cpg]="--contexts NCG  --mutation-types promoter,intronic" \
                        [non_tpc]="--invert-context-request  --contexts TCN  --mutation-types promoter,intronic" \
                      )

for CONTEXT  in  ${CONTEXTS}; do
  mkdir -p  ${SITES_FOLDER}/${CONTEXT}/samples
done

for CANCER_TYPE  in  ${!CANCER_SAMPLES_BY_TYPE[@]}; do
  ruby bin/preparations/filter_mutations.rb  --cancer-samples ${CANCER_SAMPLES_BY_TYPE[$CANCER_TYPE]}  ${SNV_FOLDER}/SNV_infos_cancer.txt  ${SITES_FOLDER}/any/sites_cancer.txt  >  ${SITES_FOLDER}/any/samples/sites_cancer_${CANCER_TYPE}.txt

  for CONTEXT in ${CONTEXTS}; do
    if [[ "$CONTEXT" != "any" ]]; then
      ruby bin/preparations/filter_mutations.rb   ${SNV_FOLDER}/SNV_infos_cancer.txt  \
                                                  ${SITES_FOLDER}/any/samples/sites_cancer_${CANCER_TYPE}.txt  \
                                                  ${CONTEXT_DECLARATIONS[$CONTEXT]} \
                                                  >  ${SITES_FOLDER}/${CONTEXT}/samples/sites_cancer_${CANCER_TYPE}.txt
    fi
  done
done


for CANCER_TYPE  in  ${!CANCER_SAMPLES_BY_TYPE[@]}; do
  mkdir -p  $FITTING_FOLDER  ${MOTIF_STATISTICS_FOLDER}/fitting_log/  ${MOTIF_STATISTICS_FOLDER}/slices/samples  ${MOTIF_STATISTICS_FOLDER}/full/samples

  for CONTEXT in ${CONTEXTS}; do
    mkdir -p  ${FITTING_FOLDER}/${CONTEXT}/samples/${CANCER_TYPE}  ${MOTIF_STATISTICS_FOLDER}/fitting_log/${CONTEXT}/samples/${CANCER_TYPE}  ${MOTIF_STATISTICS_FOLDER}/slices/${CONTEXT}  ${MOTIF_STATISTICS_FOLDER}/full/${CONTEXT}/samples

    for RANDOM_VARIANT  in  ${RANDOM_VARIANTS}; do
      mkdir -p  ${MOTIF_STATISTICS_FOLDER}/fitting_log/${CONTEXT}/samples

      # echo "ruby fitting_random_sites.rb  ${SITES_FOLDER}/${CONTEXT}/samples/sites_cancer_${CANCER_TYPE}.txt  ${SITES_FOLDER}/${CONTEXT}/sites_${RANDOM_VARIANT}.txt  --fold $FITTING_FOLD  >  ${FITTING_FOLDER}/${CONTEXT}/samples/${CANCER_TYPE}/sites_${RANDOM_VARIANT}.txt  2>  ${MOTIF_STATISTICS_FOLDER}/fitting_log/${CONTEXT}/samples/${CANCER_TYPE}/${RANDOM_VARIANT}.log"

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
done

for CANCER_TYPE  in  ${!CANCER_SAMPLES_BY_TYPE[@]}; do
  for CONTEXT in ${CONTEXTS}; do
    for VARIANT  in  cancer  ${RANDOM_VARIANTS}; do
      mkdir -p  ${MOTIF_STATISTICS_FOLDER}/slices/${CONTEXT}/samples/${CANCER_TYPE}/${VARIANT}
      ./generate_all_motif_statistics.sh  ${FITTING_FOLDER}/${CONTEXT}/samples/${CANCER_TYPE}/sites_${VARIANT}.txt  \
                                          ${MOTIF_STATISTICS_FOLDER}/slices/${CONTEXT}/samples/${CANCER_TYPE}/${VARIANT}  \
                                          0.0005  $FOLD_CHANGE

    done
  done
done
