#!/bin/bash

# MOTIF_STATISTICS_FOLDER should be specified in environment variables

cd "$(dirname "$0")"

EXPAND_CONTROL_SET_FOLD=1

mkdir -p $MOTIF_STATISTICS_FOLDER/filtered

for CONTEXT in ${CONTEXTS}; do
  mkdir -p  ${MOTIF_STATISTICS_FOLDER}/full/${CONTEXT}

  for RANDOM_VARIANT  in  ${RANDOM_VARIANTS}; do
    # --ignore-underfitting \
    ruby summary.rb   ${MOTIF_STATISTICS_FOLDER}/slices/${CONTEXT}/cancer  \
                      ${MOTIF_STATISTICS_FOLDER}/slices/${CONTEXT}/${RANDOM_VARIANT}  \
                      ./source_data/motif_names.txt  ./source_data/hocomoco_genes_infos.csv  \
                      ${MOTIF_STATISTICS_FOLDER}/fitting_log/${CONTEXT}/${RANDOM_VARIANT}.log  \
                      --correction ${CORRECTION_METHOD}  \
                      --expand-control-set ${EXPAND_CONTROL_SET_FOLD} \
                      >  ${MOTIF_STATISTICS_FOLDER}/full/${CONTEXT}/${RANDOM_VARIANT}.csv

    for SUBJECTED_OR_PROTECTED  in  subjected  protected; do
      for CHARACTERISTIC  in  disruption  emergence 'substitution-in-core'; do
        mkdir -p ${MOTIF_STATISTICS_FOLDER}/filtered/${SUBJECTED_OR_PROTECTED}/${CHARACTERISTIC}/${CONTEXT}
        ruby filter_summary.rb  ${MOTIF_STATISTICS_FOLDER}/full/${CONTEXT}/${RANDOM_VARIANT}.csv  \
                                --motif-qualities A,B,C,D  --significance 0.05  \
                                --characteristic ${CHARACTERISTIC}  --${SUBJECTED_OR_PROTECTED}  \
                                >  ${MOTIF_STATISTICS_FOLDER}/filtered/${SUBJECTED_OR_PROTECTED}/${CHARACTERISTIC}/${CONTEXT}/${RANDOM_VARIANT}.csv
      done
    done
  done

  for SUBJECTED_OR_PROTECTED  in  subjected  protected; do
    for CHARACTERISTIC  in  disruption  emergence 'substitution-in-core'; do

      FILTERED_FOLDER=${MOTIF_STATISTICS_FOLDER}/filtered/${SUBJECTED_OR_PROTECTED}/${CHARACTERISTIC}/${CONTEXT}
      COMMON_MOTIFS_FOLDER=${MOTIF_STATISTICS_FOLDER}/common_motifs/${SUBJECTED_OR_PROTECTED}/${CHARACTERISTIC}/${CONTEXT}

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
done
