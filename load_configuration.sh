#!/bin/bash

# ATTENTION! `./load_configuration.sh` will have no effect. It'd be used as `source ./load_configuration.sh`

export RESULTS_FOLDER=./results
export SEQ_FOLDER=${RESULTS_FOLDER}/sequences
export SNV_FOLDER=${RESULTS_FOLDER}/SNVs
export CHUNK_FOLDER=${RESULTS_FOLDER}/sequence_chunks
export SITES_FOLDER=${RESULTS_FOLDER}/sites
export FITTING_FOLDER=${RESULTS_FOLDER}/fitted_sites
export MOTIF_STATISTICS_FOLDER=${RESULTS_FOLDER}/motif_statistics

export EXPAND_FLANKS_LENGTH=11
export CORRECTION_METHOD="fdr"

export FITTING_FOLD_GENOME=35
export FITTING_FOLD_SHUFFLE=25
export FOLD_CHANGE_CUTOFF=4
export SITE_PVALUE_CUTOFF=0.0005
export MOTIF_NAMES=./source_data/motif_names.txt
export GENE_INFOS=./source_data/hocomoco_genes_infos.csv

export EXPAND_CONTROL_SET_FOLD=1

export CONTEXTS="any cpg tpc non_tpc"

export RANDOM_GENOME_SEEDS="123" #"13-15-17" #"13 15 17"
export RANDOM_SHUFFLE_SEEDS="124" #"135-137-139" #"135 137 139"

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

export RANDOM_GENOME_VARIANTS
export RANDOM_SHUFFLE_VARIANTS
export RANDOM_VARIANTS
