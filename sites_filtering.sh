#!/bin/bash

# SNV_FOLDER, CHUNK_FOLDER and SITES_FOLDER should be specified in environment variables
cd "$(dirname "$0")"

for CONTEXT in any cpg tpc
do
  mkdir -p ${SITES_FOLDER}/${CONTEXT}
done

for VARIANT  in  cancer  random_genome_13  random_genome_15  random_genome_17  random_shuffle_135  random_shuffle_137  random_shuffle_139; do
  ln -f  ${CHUNK_FOLDER}/sites_${VARIANT}.txt  ${SITES_FOLDER}/any/sites_${VARIANT}.txt
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
