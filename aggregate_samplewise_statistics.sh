#!/bin/bash

# SITES_FOLDER, FITTING_FOLDER and MOTIF_STATISTICS_FOLDER and so on should be specified in environment variables
# Also CANCER_SAMPLES_BY_TYPE, FITTING_FOLD, FOLD_CHANGE, MIN_FITTING_RATE should be specified

cd "$(dirname "$0")"

CANCER_TYPE=$1
CONTEXT=$2

mkdir -p  ${MOTIF_STATISTICS_FOLDER}/full/${CONTEXT}/samples/${CANCER_TYPE}

for RANDOM_VARIANT  in  ${RANDOM_VARIANTS}; do
  ruby summary.rb   ${MOTIF_STATISTICS_FOLDER}/slices/${CONTEXT}/samples/${CANCER_TYPE}/cancer  \
                    ${MOTIF_STATISTICS_FOLDER}/slices/${CONTEXT}/samples/${CANCER_TYPE}/${RANDOM_VARIANT}  \
                    --fitting-log ${MOTIF_STATISTICS_FOLDER}/fitting_log/${CONTEXT}/samples/${CANCER_TYPE}/${RANDOM_VARIANT}.log  \
                    --correction ${CORRECTION_METHOD}  \
                    --expand-control-set ${EXPAND_CONTROL_SET_FOLD} \
                    >  ${MOTIF_STATISTICS_FOLDER}/full/${CONTEXT}/samples/${CANCER_TYPE}/${RANDOM_VARIANT}.csv

  for SUBJECTED_OR_PROTECTED  in  subjected  protected; do
    for CHARACTERISTIC  in  disruption  emergence 'substitution-in-core'; do
      mkdir -p ${MOTIF_STATISTICS_FOLDER}/filtered/${SUBJECTED_OR_PROTECTED}/${CHARACTERISTIC}/${CONTEXT}/samples/${CANCER_TYPE}
      ruby filter_summary.rb  ${MOTIF_STATISTICS_FOLDER}/full/${CONTEXT}/samples/${CANCER_TYPE}/${RANDOM_VARIANT}.csv  \
                              --motif-qualities A,B,C,D  --significance 0.05  \
                              --characteristic ${CHARACTERISTIC}  --${SUBJECTED_OR_PROTECTED}  \
                              >  ${MOTIF_STATISTICS_FOLDER}/filtered/${SUBJECTED_OR_PROTECTED}/${CHARACTERISTIC}/${CONTEXT}/samples/${CANCER_TYPE}/${RANDOM_VARIANT}.csv
    done
  done
done

# mkdir -p  ${MOTIF_STATISTICS_FOLDER}/pvalue_statitics/${CONTEXT}/samples/${CANCER_TYPE}
# ls ${MOTIF_STATISTICS_FOLDER}/full/${CONTEXT}/samples/${CANCER_TYPE}/random_genome_*.csv | ruby motif_pvalue_stats.rb  >  ${MOTIF_STATISTICS_FOLDER}/pvalue_statitics/${CONTEXT}/samples/${CANCER_TYPE}/compared_to_each_genome.csv
# ls ${MOTIF_STATISTICS_FOLDER}/full/${CONTEXT}/samples/${CANCER_TYPE}/random_shuffle_*.csv | ruby motif_pvalue_stats.rb  >  ${MOTIF_STATISTICS_FOLDER}/pvalue_statitics/${CONTEXT}/samples/${CANCER_TYPE}/compared_to_each_shuffle.csv
# ls ${MOTIF_STATISTICS_FOLDER}/full/${CONTEXT}/samples/${CANCER_TYPE}/random_*.csv | ruby motif_pvalue_stats.rb  >  ${MOTIF_STATISTICS_FOLDER}/pvalue_statitics/${CONTEXT}/samples/${CANCER_TYPE}/compared_to_each.csv

for SUBJECTED_OR_PROTECTED  in  subjected  protected; do
  for CHARACTERISTIC  in  disruption  emergence 'substitution-in-core'; do

    FILTERED_FOLDER=${MOTIF_STATISTICS_FOLDER}/filtered/${SUBJECTED_OR_PROTECTED}/${CHARACTERISTIC}/${CONTEXT}/samples/${CANCER_TYPE}
    COMMON_MOTIFS_FOLDER=${MOTIF_STATISTICS_FOLDER}/common_motifs/${SUBJECTED_OR_PROTECTED}/${CHARACTERISTIC}/${CONTEXT}/samples/${CANCER_TYPE}

    mkdir -p  ${COMMON_MOTIFS_FOLDER}

    RANDOM_GENOME_VARIANTS_FILTERED=""
    for SEED in ${RANDOM_GENOME_SEEDS}; do
      RANDOM_GENOME_VARIANTS_FILTERED+=" ${FILTERED_FOLDER}/random_genome_${SEED}.csv"
    done

    RANDOM_SHUFFLE_VARIANTS_FILTERED=""
    for SEED in ${RANDOM_SHUFFLE_SEEDS}; do
      RANDOM_SHUFFLE_VARIANTS_FILTERED+=" ${FILTERED_FOLDER}/random_shuffle_${SEED}.csv"
    done

    echo ${RANDOM_GENOME_VARIANTS_FILTERED} | ruby common_motifs.rb  >  ${COMMON_MOTIFS_FOLDER}/compared_to_each_genome.txt
    echo ${RANDOM_SHUFFLE_VARIANTS_FILTERED} | ruby common_motifs.rb  >  ${COMMON_MOTIFS_FOLDER}/compared_to_each_shuffle.txt
    echo "${RANDOM_GENOME_VARIANTS_FILTERED} ${RANDOM_SHUFFLE_VARIANTS_FILTERED}" | ruby common_motifs.rb  >  ${COMMON_MOTIFS_FOLDER}/compared_to_each.txt
  done
done
