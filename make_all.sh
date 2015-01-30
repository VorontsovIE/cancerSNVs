#!/bin/bash
cd "$(dirname "$0")"

# Prepare marked up SNVs and sequences.
# Generate random sequences.
# Split files into chunks for running computations in parallel.
./preparations.sh

# ATTENTION! This command will run lots of background java workers.
# This step can take several days.
# If possible, put chunks on different machines so that workers run on several # cores and fix script so that it run only chunks on a current machine.
# In that case one also need to defer concatenating results and do it manually
# after joining processed chunks from all machines.
./results/sequence_chunks/run_perfectosape_multithread.sh



# Separate sites by contexts. Perform random sites fitting. Extract motif statistics slices.
./filtering_and_fitting

# Aggregate motif statistics slices and find motifs of interest
./aggregate_motifs_statistics.sh  $MOTIF_STATISTICS_FOLDER
