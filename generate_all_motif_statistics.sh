#!/bin/bash
cd "$(dirname "$0")"

MOTIF_NAMES=./source_data/motif_names.txt
SITES_FILE=$1
OUTPUT_DIR=$2
PVALUE_CUTOFF=$3
FOLD_CHANGE_CUTOFF=$4
# ruby calculate_motif_statistics.rb  $SITES_FILE  $MOTIF_NAMES  >  ${OUTPUT_DIR}/all.txt
ruby calculate_motif_statistics.rb  $SITES_FILE  $MOTIF_NAMES  --site-before $PVALUE_CUTOFF  >  ${OUTPUT_DIR}/sites_before.txt
ruby calculate_motif_statistics.rb  $SITES_FILE  $MOTIF_NAMES  --site-after $PVALUE_CUTOFF  >  ${OUTPUT_DIR}/sites_after.txt
ruby calculate_motif_statistics.rb  $SITES_FILE  $MOTIF_NAMES  --site-before $PVALUE_CUTOFF  --disrupted $FOLD_CHANGE_CUTOFF  >  ${OUTPUT_DIR}/sites_disrupted.txt
ruby calculate_motif_statistics.rb  $SITES_FILE  $MOTIF_NAMES  --site-after $PVALUE_CUTOFF  --emerged $FOLD_CHANGE_CUTOFF  >  ${OUTPUT_DIR}/sites_emerged.txt
