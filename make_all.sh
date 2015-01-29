#!/bin/bash

SEQ_FOLDER=./results/sequences
SNV_FOLDER=./results/SNVs
CHUNK_FOLDER=./results/sequence_chunks

FITTING_FOLD=1
FOLD_CHANGE=5

SITES_FOLDER=./results/sites
FITTING_FOLDER=./results/fitted_sites
MOTIF_STATISTICS_FOLDER=./results/motif_statistics



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

mkdir -p $FITTING_FOLDER
mkdir -p ${MOTIF_STATISTICS_FOLDER}/fitting_log
mkdir -p ${MOTIF_STATISTICS_FOLDER}/slices
mkdir -p ${MOTIF_STATISTICS_FOLDER}/full

for CONTEXT in any cpg tpc; do
  for RANDOM_VARIANT  in  random_genome_13  random_genome_15  random_genome_17  random_shuffle_135  random_shuffle_137  random_shuffle_139; do

    ruby fitting_random_sites.rb  ${SITES_FOLDER}/${CONTEXT}/cancer_sites.txt  \
                                  ${SITES_FOLDER}/${CONTEXT}/sites_${RANDOM_VARIANT}.txt  \
                                  --fold $FITTING_FOLD  \
                                  >  ${FITTING_FOLDER}/${CONTEXT}/sites_${RANDOM_VARIANT}.txt  \
                                  2>  ${MOTIF_STATISTICS_FOLDER}/fitting_log/${CONTEXT}/${RANDOM_VARIANT}.log
  done
  ln  ${SITES_FOLDER}/${CONTEXT}/cancer_sites.txt  ${FITTING_FOLDER}/${CONTEXT}/sites_cancer.txt

  for VARIANT  in  cancer  random_genome_13  random_genome_15  random_genome_17  random_shuffle_135  random_shuffle_137  random_shuffle_139; do
    ./generate_all_motif_statistics.sh  ${FITTING_FOLDER}/${CONTEXT}/sites_${VARIANT}.txt  \
                                        ${MOTIF_STATISTICS_FOLDER}/slices/${CONTEXT}/${VARIANT}  \
                                        0.0005  $FOLD_CHANGE

  done
done

./agregate_motifs_statistics.sh  $MOTIF_STATISTICS_FOLDER
