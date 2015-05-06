#!/bin/bash

# RESULTS_FOLDER, SEQ_FOLDER, SNV_FOLDER and CHUNK_FOLDER should be specified in environment variables

cd "$(dirname "$0")"

mkdir -p  $SEQ_FOLDER  $SNV_FOLDER  $CHUNK_FOLDER

ruby bin/preparations/load_cancer_mutations_sequences.rb --promoter-upstream 5000 --promoter-downstream 500 --kataegis-expansion 1000

# Markup SNVs and filter regulatory only, non-duplicated SNVs
ruby bin/preparations/snv_markup.rb  ./source_data/SNV_infos_original.txt  >  ${RESULTS_FOLDER}/SNV_infos_marked_up.txt
ruby bin/preparations/filter_snv_infos.rb  ${RESULTS_FOLDER}/SNV_infos_marked_up.txt  >  ${RESULTS_FOLDER}/SNV_infos_regulatory.txt


#  Later we work only with regulatory non-duplicated SNVs
ln -f  ${RESULTS_FOLDER}/SNV_infos_regulatory.txt  ${SNV_FOLDER}/SNV_infos_cancer.txt

##################################

# Extract sequences around SNVs (not more necessary as SNV infos and sequences are joined in a single file now)
# ruby bin/preparations/extract_snv_sequences.rb  ${SNV_FOLDER}/SNV_infos_cancer.txt  ./source_data/genome  --flank-length 50  >  ${SEQ_FOLDER}/sequences_cancer.txt

# Generate random sequences
./prepare_random_sequences.sh

# Split sequences into equal chunks in order to run chunks in parallel
NUMBER_OF_CORES=8  ./prepare_sequences_for_perfectosape_run.sh
