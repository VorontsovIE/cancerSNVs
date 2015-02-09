#!/bin/bash

# MOTIF_STATISTICS_FOLDER should be specified in environment variables

cd "$(dirname "$0")"

MIN_FITTING_RATE=0.99

mkdir -p $MOTIF_STATISTICS_FOLDER/filtered

for CONTEXT in any cpg tpc; do
  mkdir -p  ${MOTIF_STATISTICS_FOLDER}/full/${CONTEXT}

  for RANDOM_VARIANT  in  random_shuffle_135  random_shuffle_137  random_shuffle_139  random_genome_13  random_genome_15  random_genome_17; do
    ruby summary.rb   ${MOTIF_STATISTICS_FOLDER}/slices/${CONTEXT}/cancer  \
                      ${MOTIF_STATISTICS_FOLDER}/slices/${CONTEXT}/${RANDOM_VARIANT}  \
                      ./source_data/motif_names.txt  ./source_data/hocomoco_genes_infos.csv  \
                      --correction fdr  \
                      >  ${MOTIF_STATISTICS_FOLDER}/full/${CONTEXT}/${RANDOM_VARIANT}.csv

    for SUBJECTED_OR_PROTECTED  in  subjected  protected; do
      for DISRUPTION_OR_EMERGENCE  in  disruption  emergence; do
        mkdir -p ${MOTIF_STATISTICS_FOLDER}/filtered/${SUBJECTED_OR_PROTECTED}/${DISRUPTION_OR_EMERGENCE}/${CONTEXT}
        ruby filter_summary.rb  ${MOTIF_STATISTICS_FOLDER}/full/${CONTEXT}/${RANDOM_VARIANT}.csv  \
                                --motif-qualities A,B,C,D  --significance 0.05  \
                                --${DISRUPTION_OR_EMERGENCE}  --${SUBJECTED_OR_PROTECTED}  \
                                --fitting-log ${MOTIF_STATISTICS_FOLDER}/fitting_log/${CONTEXT}/${RANDOM_VARIANT}.log  \
                                --min-fitting-rate $MIN_FITTING_RATE  \
                                >  ${MOTIF_STATISTICS_FOLDER}/filtered/${SUBJECTED_OR_PROTECTED}/${DISRUPTION_OR_EMERGENCE}/${CONTEXT}/${RANDOM_VARIANT}.csv
      done
    done
  done

  mkdir -p  ${MOTIF_STATISTICS_FOLDER}/pvalue_statitics/${CONTEXT}
  ls ${MOTIF_STATISTICS_FOLDER}/full/${CONTEXT}/random_genome_*.csv | ruby motif_pvalue_stats.rb  >  ${MOTIF_STATISTICS_FOLDER}/pvalue_statitics/${CONTEXT}/compared_to_each_genome.csv
  ls ${MOTIF_STATISTICS_FOLDER}/full/${CONTEXT}/random_shuffle_*.csv | ruby motif_pvalue_stats.rb  >  ${MOTIF_STATISTICS_FOLDER}/pvalue_statitics/${CONTEXT}/compared_to_each_shuffle.csv
  ls ${MOTIF_STATISTICS_FOLDER}/full/${CONTEXT}/random_*.csv | ruby motif_pvalue_stats.rb  >  ${MOTIF_STATISTICS_FOLDER}/pvalue_statitics/${CONTEXT}/compared_to_each.csv
done

for CONTEXT in any cpg tpc; do
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
