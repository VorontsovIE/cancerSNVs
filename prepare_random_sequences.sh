#!/bin/bash

# SEQ_FOLDER and SNV_FOLDER should be specified in environment variables

cd "$(dirname "$0")"

# Generate random sequences from genome
for SEED in ${RANDOM_GENOME_SEEDS}; do
  # generate random sequences
  ruby bin/preparations/generate_random_genome_sequences.rb  ${SNV_FOLDER}/SNV_infos_cancer.txt  ./source_data/exons.txt  ./source_data/genome  --fold 100  --random-seed $SEED --flank-length 50  >  ${SEQ_FOLDER}/sequences_random_genome_${SEED}.txt

  # generate SNV infos
  ruby bin/preparations/create_snv_infos_by_genome_sequences.rb  ${SEQ_FOLDER}/sequences_random_genome_${SEED}.txt  >  ${SNV_FOLDER}/SNV_infos_random_genome_${SEED}_non_marked_up.txt
  ruby bin/preparations/snv_markup.rb ${SNV_FOLDER}/SNV_infos_random_genome_${SEED}_non_marked_up.txt  ./source_data/exons.txt  >  ${SNV_FOLDER}/SNV_infos_random_genome_${SEED}.txt
  rm ${SNV_FOLDER}/SNV_infos_random_genome_${SEED}_non_marked_up.txt
done


# Generate random sequences by shuffling contexts
for SEED in ${RANDOM_SHUFFLE_SEEDS}; do
  # generate random sequences
  ruby ./bin/preparations/shuffle_SNVs.rb  ${SEQ_FOLDER}/sequences_cancer.txt  --random-seed $SEED  --fold 100  >  ${SEQ_FOLDER}/sequences_random_shuffle_${SEED}.txt
  # generate SNV infos
  ln -f  ${SNV_FOLDER}/SNV_infos_cancer.txt  ${SNV_FOLDER}/SNV_infos_random_shuffle_${SEED}.txt
done
