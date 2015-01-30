#!/bin/bash
cd "$(dirname "$0")"

export SEQ_FOLDER=./results/sequences
export SNV_FOLDER=./results/SNVs
export CHUNK_FOLDER=./results/sequence_chunks
export SITES_FOLDER=./results/sites
export FITTING_FOLDER=./results/fitted_sites
export MOTIF_STATISTICS_FOLDER=./results/motif_statistics

# Prepare marked up SNVs and sequences.
# Generate random sequences.
# Split files into chunks for running computations in parallel.
./preparations.sh

# ATTENTION! This command will run lots of background java workers.
# This step can take several days.
# If possible, put chunks on different machines so that workers run on several # cores and fix script so that it run only chunks on a current machine.
# In that case one also need to defer concatenating results and do it manually
# after joining processed chunks from all machines.
${CHUNK_FOLDER}/run_perfectosape_multithread.sh

# Separate sites by contexts. Perform random sites fitting. Extract motif statistics slices.
./filtering_and_fitting

# Aggregate motif statistics slices and find motifs of interest
./aggregate_motifs_statistics.sh
