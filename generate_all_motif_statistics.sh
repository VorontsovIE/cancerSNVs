#!/bin/bash
cd "$(dirname "$0")"

SITES_FILE=$1
OUTPUT_DIR=$2

mkdir -p  ${OUTPUT_DIR}
# ruby calculate_motif_statistics.rb  ${SITES_FILE}  ${MOTIF_NAMES}  >  ${OUTPUT_DIR}/all.txt
ruby calculate_motif_statistics.rb  ${SITES_FILE}  ${MOTIF_NAMES}  --site-before ${SITE_PVALUE_CUTOFF}  >  ${OUTPUT_DIR}/sites_before.txt &
ruby calculate_motif_statistics.rb  ${SITES_FILE}  ${MOTIF_NAMES}  --site-after ${SITE_PVALUE_CUTOFF}  >  ${OUTPUT_DIR}/sites_after.txt &
ruby calculate_motif_statistics.rb  ${SITES_FILE}  ${MOTIF_NAMES}  --site-before ${SITE_PVALUE_CUTOFF}  --disrupted ${FOLD_CHANGE_CUTOFF}  >  ${OUTPUT_DIR}/sites_disrupted.txt &
ruby calculate_motif_statistics.rb  ${SITES_FILE}  ${MOTIF_NAMES}  --site-after ${SITE_PVALUE_CUTOFF}  --emerged ${FOLD_CHANGE_CUTOFF}  >  ${OUTPUT_DIR}/sites_emerged.txt &

ruby calculate_motif_statistics.rb  ${SITES_FILE}  ${MOTIF_NAMES}  --substitution-in-core  >  ${OUTPUT_DIR}/substitutions_in_core.txt &
ruby calculate_motif_statistics.rb  ${SITES_FILE}  ${MOTIF_NAMES}  --substitution-in-flank  >  ${OUTPUT_DIR}/substitutions_in_flank.txt &
wait
