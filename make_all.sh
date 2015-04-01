#!/bin/bash
cd "$(dirname "$0")"

export RESULTS_FOLDER=./results
export SEQ_FOLDER=${RESULTS_FOLDER}/sequences
export SNV_FOLDER=${RESULTS_FOLDER}/SNVs
export CHUNK_FOLDER=${RESULTS_FOLDER}/sequence_chunks
export SITES_FOLDER=${RESULTS_FOLDER}/sites
export FITTING_FOLDER=${RESULTS_FOLDER}/fitted_sites
export MOTIF_STATISTICS_FOLDER=${RESULTS_FOLDER}/motif_statistics


FITTING_FOLD_GENOME=35
FITTING_FOLD_SHUFFLE=25
FITTING_FOLD=1

CONTEXTS="any cpg tpc non_tpc"

RANDOM_GENOME_SEEDS="123" #"13-15-17" #"13 15 17"
RANDOM_SHUFFLE_SEEDS="124" #"135-137-139" #"135 137 139"

RANDOM_GENOME_VARIANTS=""
for SEED in ${RANDOM_GENOME_SEEDS}; do
  RANDOM_GENOME_VARIANTS+=" random_genome_${SEED}"
done

RANDOM_SHUFFLE_VARIANTS=""
for SEED in ${RANDOM_SHUFFLE_SEEDS}; do
  RANDOM_SHUFFLE_VARIANTS+=" random_shuffle_${SEED}"
done

RANDOM_VARIANTS=""
RANDOM_VARIANTS+=" ${RANDOM_GENOME_VARIANTS}"
RANDOM_VARIANTS+=" ${RANDOM_SHUFFLE_VARIANTS}"

export RANDOM_GENOME_SEEDS
export RANDOM_GENOME_VARIANTS
export RANDOM_SHUFFLE_SEEDS
export RANDOM_SHUFFLE_VARIANTS
export RANDOM_VARIANTS
export CONTEXTS
export FITTING_FOLD
export FITTING_FOLD_GENOME
export FITTING_FOLD_SHUFFLE

# Prepare marked up SNVs and sequences.
# Generate random sequences.
# Split files into chunks for running computations in parallel.
# Estimated time is about 40 min to complete preparations
./preparations.sh

# ATTENTION! This command will run lots of background java workers.
# This step can take several days (it took about 1.5 days using 8 cores).
# If possible, put chunks on different machines so that workers run on several # cores and fix script so that it run only chunks on a current machine.
# In that case one also need to defer concatenating results and do it manually
# after joining processed chunks from all machines.
# By default we use 4 workers (for 4 cores), but it's easily configurable.
# If necessary one also can provide JVM options for PerfectosAPE run in a file ./prepare_sequences_for_perfectosape_run.sh
# For example it is worth to allow JVM to take up to 1-2Gb RAM (note, each worker will consume that much)
${CHUNK_FOLDER}/run_perfectosape_multithread.sh

# Separate sites by contexts.
time ./sites_filtering.sh
# Perform random sites fitting.
time ./sites_fitting.sh

# Extract motif statistics slices.
time ./collect_motif_statistics_slices.sh
# Aggregate motif statistics slices and find motifs of interest
time ./aggregate_motifs_statistics.sh
time ./aggregate_additional_motifs_statistics.sh
