#!/bin/bash

# SITES_FOLDER, FITTING_FOLDER and MOTIF_STATISTICS_FOLDER and so on should be specified in environment variables
# Also CANCER_SAMPLES_BY_TYPE, FITTING_FOLD, FOLD_CHANGE, MIN_FITTING_RATE should be specified

cd "$(dirname "$0")"

CANCER_TYPE=$1

for CONTEXT in ${CONTEXTS}; do
  mkdir -p  ${MOTIF_STATISTICS_FOLDER}/full/${CONTEXT}/samples/${CANCER_TYPE}

  for RANDOM_VARIANT  in  ${RANDOM_VARIANTS}; do
    ruby summary.rb   ${MOTIF_STATISTICS_FOLDER}/slices/${CONTEXT}/samples/${CANCER_TYPE}/cancer  \
                      ${MOTIF_STATISTICS_FOLDER}/slices/${CONTEXT}/samples/${CANCER_TYPE}/${RANDOM_VARIANT}  \
                      ./source_data/motif_names.txt  ./source_data/hocomoco_genes_infos.csv  \
                      ${MOTIF_STATISTICS_FOLDER}/fitting_log/${CONTEXT}/samples/${CANCER_TYPE}/${RANDOM_VARIANT}.log  \
                      --correction ${CORRECTION_METHOD}  \
                      --expand-control-set ${EXPAND_CONTROL_SET_FOLD} \
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
