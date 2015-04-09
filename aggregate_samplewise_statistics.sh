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

for CANCER_TYPE  in  ${!CANCER_SAMPLES_BY_TYPE[@]};  do
  for CONTEXT in ${CONTEXTS}; do
    mkdir -p  ${MOTIF_STATISTICS_FOLDER}/full/${CONTEXT}/samples/${CANCER_TYPE}

    for RANDOM_VARIANT  in  ${RANDOM_VARIANTS}; do
      ruby summary.rb   ${MOTIF_STATISTICS_FOLDER}/slices/${CONTEXT}/samples/${CANCER_TYPE}/cancer  \
                        ${MOTIF_STATISTICS_FOLDER}/slices/${CONTEXT}/samples/${CANCER_TYPE}/${RANDOM_VARIANT}  \
                        ./source_data/motif_names.txt  ./source_data/hocomoco_genes_infos.csv  \
                        ${MOTIF_STATISTICS_FOLDER}/fitting_log/${CONTEXT}/samples/${CANCER_TYPE}/${RANDOM_VARIANT}.log  \
                        --correction ${CORRECTION_METHOD}  \
                        >  ${MOTIF_STATISTICS_FOLDER}/full/${CONTEXT}/samples/${CANCER_TYPE}/${RANDOM_VARIANT}.csv

      for SUBJECTED_OR_PROTECTED  in  subjected  protected; do
        for DISRUPTION_OR_EMERGENCE  in  disruption  emergence; do
          mkdir -p ${MOTIF_STATISTICS_FOLDER}/filtered/${SUBJECTED_OR_PROTECTED}/${DISRUPTION_OR_EMERGENCE}/${CONTEXT}/samples/${CANCER_TYPE}
          ruby filter_summary.rb  ${MOTIF_STATISTICS_FOLDER}/full/${CONTEXT}/samples/${CANCER_TYPE}/${RANDOM_VARIANT}.csv  \
                                  --motif-qualities A,B,C,D  --significance 0.05  \
                                  --${DISRUPTION_OR_EMERGENCE}  --${SUBJECTED_OR_PROTECTED}  \
                                  >  ${MOTIF_STATISTICS_FOLDER}/filtered/${SUBJECTED_OR_PROTECTED}/${DISRUPTION_OR_EMERGENCE}/${CONTEXT}/samples/${CANCER_TYPE}/${RANDOM_VARIANT}.csv
        done
      done
    done

    mkdir -p  ${MOTIF_STATISTICS_FOLDER}/pvalue_statitics/${CONTEXT}/samples/${CANCER_TYPE}
    ls ${MOTIF_STATISTICS_FOLDER}/full/${CONTEXT}/samples/${CANCER_TYPE}/random_genome_*.csv | ruby motif_pvalue_stats.rb  >  ${MOTIF_STATISTICS_FOLDER}/pvalue_statitics/${CONTEXT}/samples/${CANCER_TYPE}/compared_to_each_genome.csv
    ls ${MOTIF_STATISTICS_FOLDER}/full/${CONTEXT}/samples/${CANCER_TYPE}/random_shuffle_*.csv | ruby motif_pvalue_stats.rb  >  ${MOTIF_STATISTICS_FOLDER}/pvalue_statitics/${CONTEXT}/samples/${CANCER_TYPE}/compared_to_each_shuffle.csv
    ls ${MOTIF_STATISTICS_FOLDER}/full/${CONTEXT}/samples/${CANCER_TYPE}/random_*.csv | ruby motif_pvalue_stats.rb  >  ${MOTIF_STATISTICS_FOLDER}/pvalue_statitics/${CONTEXT}/samples/${CANCER_TYPE}/compared_to_each.csv
  done

  for CONTEXT in ${CONTEXTS}; do
    for SUBJECTED_OR_PROTECTED  in  subjected  protected; do
      for DISRUPTION_OR_EMERGENCE  in  disruption  emergence; do

        FILTERED_FOLDER=${MOTIF_STATISTICS_FOLDER}/filtered/${SUBJECTED_OR_PROTECTED}/${DISRUPTION_OR_EMERGENCE}/${CONTEXT}/samples/${CANCER_TYPE}
        COMMON_MOTIFS_FOLDER=${MOTIF_STATISTICS_FOLDER}/common_motifs/${SUBJECTED_OR_PROTECTED}/${DISRUPTION_OR_EMERGENCE}/${CONTEXT}/samples/${CANCER_TYPE}

        mkdir -p  ${COMMON_MOTIFS_FOLDER}

        ls  ${FILTERED_FOLDER}/random_genome_*.csv | ruby common_motifs.rb  >  ${COMMON_MOTIFS_FOLDER}/compared_to_each_genome.txt
        ls  ${FILTERED_FOLDER}/random_shuffle_*.csv | ruby common_motifs.rb  >  ${COMMON_MOTIFS_FOLDER}/compared_to_each_shuffle.txt
        ls  ${FILTERED_FOLDER}/random_*.csv | ruby common_motifs.rb  >  ${COMMON_MOTIFS_FOLDER}/compared_to_each.txt
      done
    done
  done
done


for CANCER_TYPE_1  in  ${!CANCER_SAMPLES_BY_TYPE[@]};  do
  for CONTEXT in ${CONTEXTS}; do
    mkdir -p  ${MOTIF_STATISTICS_FOLDER}/full/${CONTEXT}/samples_vs_samples/${CANCER_TYPE_1}/

    for CANCER_TYPE_2  in  ${!CANCER_SAMPLES_BY_TYPE[@]};  do
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
  done
done
