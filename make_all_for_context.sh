#!/bin/bash

CONTEXT      = $1 # any
SITES_CANCER = $2 # ./results/intermediate/site_subsets/cancer_any.txt
SITES_RANDOM = $3 # ./results/intermediate/site_subsets/random_any.txt

FITTING_FOLD = 1
FOLD_CHANGE = 5

MOTIF_NAMES = ./source_data/motif_names.txt
MOTIF_INFOS = ./source_data/hocomoco_genes_infos.csv
SNV_INFOS = ./source_data/SNV_infos_regulatory.txt

FITTING_FOLDER = ./results/intermediate/site_subsets/fitted
SITES_FITTED = $FITTING_FOLDER/random_$CONTEXT.txt

LOG_FOLDER = $FITTING_FOLDER/log
FITTING_LOG_FILE = $LOG_FOLDER/random_$CONTEXT.log

MOTIF_STATISTICS_FOLDER = ./results/motif_statistics
## not files, but prefices for many files
MOTIF_STATISTICS_CANCER = $MOTIF_STATISTICS_FOLDER/$CONTEXT/cancer.txt
MOTIF_STATISTICS_RANDOM = $MOTIF_STATISTICS_FOLDER/$CONTEXT/random.txt

MOTIF_STATISTICS_ALL = $MOTIF_STATISTICS_FOLDER/$CONTEXT.csv
MOTIF_STATISTICS_FILTERED = $MOTIF_STATISTICS_FOLDER/$CONTEXT_filtered.csv

mkdir -p FITTING_FOLDER
mkdir -p LOG_FOLDER


ruby fitting_random_sites.rb $SITES_CANCER $SITES_RANDOM --fold $FITTING_FOLD > $SITES_FITTED 2> $FITTING_LOG_FILE

mkdir -p $MOTIF_STATISTICS_FOLDER

ruby generate_all_motif_statistics.rb  $SITES_CANCER  $MOTIF_STATISTICS_CANCER  $MOTIF_NAMES  --pvalue 0.0005  --fold-change $FOLD_CHANGE
ruby generate_all_motif_statistics.rb  $SITES_FITTED  $MOTIF_STATISTICS_RANDOM  $MOTIF_NAMES  --pvalue 0.0005  --fold-change $FOLD_CHANGE
ruby summary.rb  $MOTIF_STATISTICS_CANCER  $MOTIF_STATISTICS_RANDOM  $MOTIF_NAMES  $MOTIF_INFOS  >  $MOTIF_STATISTICS_ALL
ruby filter_summary.rb  $MOTIF_STATISTICS_ALL  >  $MOTIF_STATISTICS_FILTERED

cp -r $LOG_FOLDER  $MOTIF_STATISTICS_FOLDER
