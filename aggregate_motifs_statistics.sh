#!/bin/bash

# MOTIF_STATISTICS_FOLDER should be specified in environment variables

cd "$(dirname "$0")"

mkdir -p $MOTIF_STATISTICS_FOLDER/filtered

for CONTEXT in ${CONTEXTS}; do
  mkdir -p  ${MOTIF_STATISTICS_FOLDER}/full/${CONTEXT}

  for RANDOM_VARIANT  in  ${RANDOM_VARIANTS}; do
    ruby summary.rb   ${MOTIF_STATISTICS_FOLDER}/slices/${CONTEXT}/cancer  \
                      ${MOTIF_STATISTICS_FOLDER}/slices/${CONTEXT}/${RANDOM_VARIANT}  \
                      ./source_data/motif_names.txt  ./source_data/hocomoco_genes_infos.csv  \
                      ${MOTIF_STATISTICS_FOLDER}/fitting_log/${CONTEXT}/${RANDOM_VARIANT}.log  \
                      --correction fdr  \
                      >  ${MOTIF_STATISTICS_FOLDER}/full/${CONTEXT}/${RANDOM_VARIANT}.csv

    for SUBJECTED_OR_PROTECTED  in  subjected  protected; do
      for DISRUPTION_OR_EMERGENCE  in  disruption  emergence; do
        mkdir -p ${MOTIF_STATISTICS_FOLDER}/filtered/${SUBJECTED_OR_PROTECTED}/${DISRUPTION_OR_EMERGENCE}/${CONTEXT}
        ruby filter_summary.rb  ${MOTIF_STATISTICS_FOLDER}/full/${CONTEXT}/${RANDOM_VARIANT}.csv  \
                                --motif-qualities A,B,C,D  --significance 0.05  \
                                --${DISRUPTION_OR_EMERGENCE}  --${SUBJECTED_OR_PROTECTED}  \
                                >  ${MOTIF_STATISTICS_FOLDER}/filtered/${SUBJECTED_OR_PROTECTED}/${DISRUPTION_OR_EMERGENCE}/${CONTEXT}/${RANDOM_VARIANT}.csv
      done
    done
  done

  for SUBJECTED_OR_PROTECTED  in  subjected  protected; do
    for DISRUPTION_OR_EMERGENCE  in  disruption  emergence; do

      FILTERED_FOLDER=${MOTIF_STATISTICS_FOLDER}/filtered/${SUBJECTED_OR_PROTECTED}/${DISRUPTION_OR_EMERGENCE}/${CONTEXT}
      COMMON_MOTIFS_FOLDER=${MOTIF_STATISTICS_FOLDER}/common_motifs/${SUBJECTED_OR_PROTECTED}/${DISRUPTION_OR_EMERGENCE}/${CONTEXT}

      mkdir -p  ${COMMON_MOTIFS_FOLDER}

      ls  ${FILTERED_FOLDER}/random_genome_*.csv | ruby common_motifs.rb  >  ${COMMON_MOTIFS_FOLDER}/compared_to_each_genome.txt
      ls  ${FILTERED_FOLDER}/random_shuffle_*.csv | ruby common_motifs.rb  >  ${COMMON_MOTIFS_FOLDER}/compared_to_each_shuffle.txt
      ls  ${FILTERED_FOLDER}/random_*.csv | ruby common_motifs.rb  >  ${COMMON_MOTIFS_FOLDER}/compared_to_each.txt
    done
  done
done
