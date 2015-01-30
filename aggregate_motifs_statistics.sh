#!/bin/bash
cd "$(dirname "$0")"

MOTIF_STATISTICS_FOLDER=$1

MOTIF_NAMES=./source_data/motif_names.txt
MOTIF_INFOS=./source_data/hocomoco_genes_infos.csv

mkdir -p $MOTIF_STATISTICS_FOLDER/filtered

for CONTEXT in any cpg tpc; do
  mkdir -p  ${MOTIF_STATISTICS_FOLDER}/full/${CONTEXT}  ${MOTIF_STATISTICS_FOLDER}/filtered/${CONTEXT}
  
  for RANDOM_VARIANT  in  random_shuffle_135  random_shuffle_137  random_shuffle_139  random_genome_13  random_genome_15  random_genome_17; do
    ruby summary.rb   ${MOTIF_STATISTICS_FOLDER}/slices/${CONTEXT}/cancer  \
                      ${MOTIF_STATISTICS_FOLDER}/slices/${CONTEXT}/${RANDOM_VARIANT}  \
                      $MOTIF_NAMES  $MOTIF_INFOS --correction fdr  \
                      >  ${MOTIF_STATISTICS_FOLDER}/full/${CONTEXT}/${RANDOM_VARIANT}.csv

    ruby filter_summary.rb  ${MOTIF_STATISTICS_FOLDER}/full/${CONTEXT}/${RANDOM_VARIANT}.csv  \
                            --motif-qualities A,B,C,D  --significance 0.05  \
                            >  ${MOTIF_STATISTICS_FOLDER}/filtered/${CONTEXT}/${RANDOM_VARIANT}.csv
  done
done

for CONTEXT in any cpg tpc; do
  mkdir -p  ${MOTIF_STATISTICS_FOLDER}/common_motifs/${CONTEXT}
  
  SHUFFLE_FILTERED_RESULTS_FILES=""
  for RANDOM_VARIANT  in  random_shuffle_135  random_shuffle_137  random_shuffle_139; do
    SHUFFLE_FILTERED_RESULTS_FILES="${SHUFFLE_FILTERED_RESULTS_FILES}  ${MOTIF_STATISTICS_FOLDER}/filtered/${CONTEXT}/${RANDOM_VARIANT}.csv"
  done

  GENOME_FILTERED_RESULTS_FILES=""
  for RANDOM_VARIANT  in  random_genome_13  random_genome_15  random_genome_17; do
    GENOME_FILTERED_RESULTS_FILES="${GENOME_FILTERED_RESULTS_FILES}  ${MOTIF_STATISTICS_FOLDER}/filtered/${CONTEXT}/${RANDOM_VARIANT}.csv"
  done

  ALL_FILTERED_RESULTS="${GENOME_FILTERED_RESULTS_FILES}  ${SHUFFLE_FILTERED_RESULTS_FILES}"

  ruby common_motifs.rb ${GENOME_FILTERED_RESULTS_FILES}  >  ${MOTIF_STATISTICS_FOLDER}/common_motifs/${CONTEXT}/compared_to_each_genome.txt
  ruby common_motifs.rb ${SHUFFLE_FILTERED_RESULTS_FILES}  >  ${MOTIF_STATISTICS_FOLDER}/common_motifs/${CONTEXT}/compared_to_each_shuffle.txt
  ruby common_motifs.rb ${ALL_FILTERED_RESULTS}  >  ${MOTIF_STATISTICS_FOLDER}/common_motifs/${CONTEXT}/compared_to_each.txt
done
