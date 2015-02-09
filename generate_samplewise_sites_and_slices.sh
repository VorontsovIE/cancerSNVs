#!/bin/bash

# SITES_FOLDER, FITTING_FOLDER and MOTIF_STATISTICS_FOLDER and so on should be specified in environment variables
# Also CANCER_SAMPLES_BY_TYPE, FITTING_FOLD, FOLD_CHANGE, MIN_FITTING_RATE should be specified

cd "$(dirname "$0")"

for CONTEXT  in  any cpg tpc; do
  mkdir -p  ${SITES_FOLDER}/${CONTEXT}/samples
done

for CANCER_TYPE  in  ${!CANCER_SAMPLES_BY_TYPE[@]};  do
  ruby bin/preparations/filter_mutations.rb  --cancer-samples ${CANCER_SAMPLES_BY_TYPE[$CANCER_TYPE]}  ${SNV_FOLDER}/SNV_infos_cancer.txt  ${SITES_FOLDER}/any/sites_cancer.txt  >  ${SITES_FOLDER}/any/samples/sites_cancer_${CANCER_TYPE}.txt

  ruby bin/preparations/filter_mutations.rb   ${SNV_FOLDER}/SNV_infos_cancer.txt  \
                                              ${SITES_FOLDER}/any/samples/sites_cancer_${CANCER_TYPE}.txt  \
                                              --contexts TCN  --mutation-types promoter,intronic  \
                                              >  ${SITES_FOLDER}/tpc/samples/sites_cancer_${CANCER_TYPE}.txt

  ruby bin/preparations/filter_mutations.rb   ${SNV_FOLDER}/SNV_infos_cancer.txt  \
                                              ${SITES_FOLDER}/any/samples/sites_cancer_${CANCER_TYPE}.txt  \
                                              --contexts NCG  --mutation-types promoter,intronic  \
                                              >  ${SITES_FOLDER}/cpg/samples/sites_cancer_${CANCER_TYPE}.txt
done


for CANCER_TYPE  in  ${!CANCER_SAMPLES_BY_TYPE[@]};  do
  mkdir -p  $FITTING_FOLDER  ${MOTIF_STATISTICS_FOLDER}/fitting_log/  ${MOTIF_STATISTICS_FOLDER}/slices/samples  ${MOTIF_STATISTICS_FOLDER}/full/samples

  for CONTEXT in any cpg tpc; do
    mkdir -p  ${FITTING_FOLDER}/${CONTEXT}/samples/${CANCER_TYPE}  ${MOTIF_STATISTICS_FOLDER}/fitting_log/${CONTEXT}/samples/${CANCER_TYPE}  ${MOTIF_STATISTICS_FOLDER}/slices/${CONTEXT}  ${MOTIF_STATISTICS_FOLDER}/full/${CONTEXT}/samples

    for RANDOM_VARIANT  in  random_genome_13  random_genome_15  random_genome_17  random_shuffle_135  random_shuffle_137  random_shuffle_139; do
      mkdir -p  ${MOTIF_STATISTICS_FOLDER}/fitting_log/${CONTEXT}/samples

      # echo "ruby fitting_random_sites.rb  ${SITES_FOLDER}/${CONTEXT}/samples/sites_cancer_${CANCER_TYPE}.txt  ${SITES_FOLDER}/${CONTEXT}/sites_${RANDOM_VARIANT}.txt  --fold $FITTING_FOLD  >  ${FITTING_FOLDER}/${CONTEXT}/samples/${CANCER_TYPE}/sites_${RANDOM_VARIANT}.txt  2>  ${MOTIF_STATISTICS_FOLDER}/fitting_log/${CONTEXT}/samples/${CANCER_TYPE}/${RANDOM_VARIANT}.log"

      ruby fitting_random_sites.rb  ${SITES_FOLDER}/${CONTEXT}/samples/sites_cancer_${CANCER_TYPE}.txt  \
                                    ${SITES_FOLDER}/${CONTEXT}/sites_${RANDOM_VARIANT}.txt  \
                                    --fold $FITTING_FOLD  \
                                    >  ${FITTING_FOLDER}/${CONTEXT}/samples/${CANCER_TYPE}/sites_${RANDOM_VARIANT}.txt  \
                                    2>  ${MOTIF_STATISTICS_FOLDER}/fitting_log/${CONTEXT}/samples/${CANCER_TYPE}/${RANDOM_VARIANT}.log

      ln -f  ${SITES_FOLDER}/${CONTEXT}/samples/sites_cancer_${CANCER_TYPE}.txt  ${FITTING_FOLDER}/${CONTEXT}/samples/${CANCER_TYPE}/sites_cancer.txt
    done
  done
done

for CANCER_TYPE  in  ${!CANCER_SAMPLES_BY_TYPE[@]};  do
  for CONTEXT in any cpg tpc; do
    for VARIANT  in  cancer random_genome_13  random_genome_15  random_genome_17  random_shuffle_135  random_shuffle_137  random_shuffle_139; do
      mkdir -p  ${MOTIF_STATISTICS_FOLDER}/slices/${CONTEXT}/samples/${CANCER_TYPE}/${VARIANT}
      ./generate_all_motif_statistics.sh  ${FITTING_FOLDER}/${CONTEXT}/samples/${CANCER_TYPE}/sites_${VARIANT}.txt  \
                                          ${MOTIF_STATISTICS_FOLDER}/slices/${CONTEXT}/samples/${CANCER_TYPE}/${VARIANT}  \
                                          0.0005  $FOLD_CHANGE

    done
  done
done
