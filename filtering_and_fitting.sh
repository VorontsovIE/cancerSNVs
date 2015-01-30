#!/bin/bash

# SNV_FOLDER, CHUNK_FOLDER, SITES_FOLDER, FITTING_FOLDER and MOTIF_STATISTICS_FOLDER should be specified in environment variables

cd "$(dirname "$0")"

FITTING_FOLD=1
FOLD_CHANGE=5


for CONTEXT in any cpg tpc
do
  mkdir -p ${SITES_FOLDER}/${CONTEXT}
done

for VARIANT  in  cancer  random_genome_13  random_genome_15  random_genome_17  random_shuffle_135  random_shuffle_137  random_shuffle_139; do
  ln  ${CHUNK_FOLDER}/sites_${VARIANT}.txt  ${SITES_FOLDER}/any/sites_${VARIANT}.txt
done

# generate subsets of sites in TpC/CpG mutration contexts
for RANDOM_VARIANT  in  cancer  random_genome_13  random_genome_15  random_genome_17  random_shuffle_135  random_shuffle_137  random_shuffle_139; do

  ruby bin/preparations/filter_mutations.rb   ${SNV_FOLDER}/SNV_infos_${RANDOM_VARIANT}.txt  \
                                              ${SITES_FOLDER}/any/sites_${RANDOM_VARIANT}.txt  \
                                              --contexts TCN  --mutation-types promoter,intronic  \
                                              >  ${SITES_FOLDER}/tpc/sites_${RANDOM_VARIANT}.txt

  ruby bin/preparations/filter_mutations.rb   ${SNV_FOLDER}/SNV_infos_${RANDOM_VARIANT}.txt  \
                                              ${SITES_FOLDER}/any/sites_${RANDOM_VARIANT}.txt  \
                                              --contexts NCG  --mutation-types promoter,intronic  \
                                              >  ${SITES_FOLDER}/cpg/sites_${RANDOM_VARIANT}.txt
done

mkdir -p  $FITTING_FOLDER  ${MOTIF_STATISTICS_FOLDER}/fitting_log  ${MOTIF_STATISTICS_FOLDER}/slices  ${MOTIF_STATISTICS_FOLDER}/full

for CONTEXT in any cpg tpc; do
  mkdir -p  ${FITTING_FOLDER}/${CONTEXT}  ${MOTIF_STATISTICS_FOLDER}/fitting_log/${CONTEXT}/  ${MOTIF_STATISTICS_FOLDER}/slices/${CONTEXT}  ${MOTIF_STATISTICS_FOLDER}/full/${CONTEXT}

  for RANDOM_VARIANT  in  random_genome_13  random_genome_15  random_genome_17  random_shuffle_135  random_shuffle_137  random_shuffle_139; do

    ruby fitting_random_sites.rb  ${SITES_FOLDER}/${CONTEXT}/sites_cancer.txt  \
                                  ${SITES_FOLDER}/${CONTEXT}/sites_${RANDOM_VARIANT}.txt  \
                                  --fold $FITTING_FOLD  \
                                  >  ${FITTING_FOLDER}/${CONTEXT}/sites_${RANDOM_VARIANT}.txt  \
                                  2>  ${MOTIF_STATISTICS_FOLDER}/fitting_log/${CONTEXT}/${RANDOM_VARIANT}.log
  done
  ln  ${SITES_FOLDER}/${CONTEXT}/sites_cancer.txt  ${FITTING_FOLDER}/${CONTEXT}/sites_cancer.txt

  for VARIANT  in  cancer  random_genome_13  random_genome_15  random_genome_17  random_shuffle_135  random_shuffle_137  random_shuffle_139; do
    mkdir -p  ${MOTIF_STATISTICS_FOLDER}/slices/${CONTEXT}/${VARIANT}
    ./generate_all_motif_statistics.sh  ${FITTING_FOLDER}/${CONTEXT}/sites_${VARIANT}.txt  \
                                        ${MOTIF_STATISTICS_FOLDER}/slices/${CONTEXT}/${VARIANT}  \
                                        0.0005  $FOLD_CHANGE

  done
done
