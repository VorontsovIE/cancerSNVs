#!/bin/bash
cd "$(dirname "$0")"

SEQ_FOLDER=$1
SNV_FOLDER=$2

# Generate random sequences from genome
for SEED in 13 15 17; do
  # generate random sequences
  ruby bin/preparations/generate_random_genome_sequences.rb  ${SNV_FOLDER}/SNV_infos_cancer.txt  ./source_data/exons.txt  ./source_data/genome  --fold 10  --random-seed $SEED  >  ${SEQ_FOLDER}/sequences_random_genome_${SEED}.txt

  # generate SNV infos
  ruby bin/preparations/create_snv_infos_by_genome_sequences.rb  ${SEQ_FOLDER}/sequences_random_genome_${SEED}.txt  >  ${SNV_FOLDER}/SNV_infos_random_genome_${SEED}_non_marked_up.txt
  ruby bin/preparations/snv_markup.rb ${SNV_FOLDER}/SNV_infos_random_genome_${SEED}_non_marked_up.txt  ./source_data/cage_peaks.txt  ./source_data/exons.txt  >  ${SNV_FOLDER}/SNV_infos_random_genome_${SEED}.txt
  rm ${SNV_FOLDER}/SNV_infos_random_genome_${SEED}_non_marked_up.txt
done


# Generate random sequences by shuffling contexts
for SEED in 135 137 139; do
  # generate random sequences
  ruby ./bin/preparations/shuffle_SNVs.rb  ${SEQ_FOLDER}/sequences_cancer.txt  --random-seed $SEED  --fold 10  >  ${SEQ_FOLDER}/sequences_random_shuffle_${SEED}.txt
  # generate SNV infos
  ln  ${SNV_FOLDER}/SNV_infos_cancer.txt  ${SNV_FOLDER}/SNV_infos_random_shuffle_${SEED}.txt
done
