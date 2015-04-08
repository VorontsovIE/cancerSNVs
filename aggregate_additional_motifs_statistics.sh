#!/bin/bash

cd "$(dirname "$0")"

for CONTEXT in ${CONTEXTS}; do
  mkdir -p  ${MOTIF_STATISTICS_FOLDER}/pvalue_statitics/${CONTEXT}
  ls ${MOTIF_STATISTICS_FOLDER}/full/${CONTEXT}/random_genome_*.csv | ruby motif_pvalue_stats.rb  >  ${MOTIF_STATISTICS_FOLDER}/pvalue_statitics/${CONTEXT}/compared_to_each_genome.csv
  ls ${MOTIF_STATISTICS_FOLDER}/full/${CONTEXT}/random_shuffle_*.csv | ruby motif_pvalue_stats.rb  >  ${MOTIF_STATISTICS_FOLDER}/pvalue_statitics/${CONTEXT}/compared_to_each_shuffle.csv
  ls ${MOTIF_STATISTICS_FOLDER}/full/${CONTEXT}/random_*.csv | ruby motif_pvalue_stats.rb  >  ${MOTIF_STATISTICS_FOLDER}/pvalue_statitics/${CONTEXT}/compared_to_each.csv
done

mkdir -p ./results/core_flank_content/
for CONTEXT in ${CONTEXTS}; do
  mkdir -p ./results/core_flank_content/${CONTEXT}/all
  mkdir -p ./results/core_flank_content/${CONTEXT}/filtered
  for RANDOM_VARIANT in ${RANDOM_VARIANTS}; do
    ruby core_flank_content.rb ./results/fitted_sites/${CONTEXT}/sites_cancer.txt ./results/fitted_sites/${CONTEXT}/sites_${RANDOM_VARIANT}.txt > ./results/core_flank_content/${CONTEXT}/all/${RANDOM_VARIANT}.txt
    ruby core_flank_content.rb ./results/fitted_sites/${CONTEXT}/sites_cancer.txt ./results/fitted_sites/${CONTEXT}/sites_${RANDOM_VARIANT}.txt --significance 0.05 > ./results/core_flank_content/${CONTEXT}/filtered/${RANDOM_VARIANT}.txt
  done
done

for CONTEXT in ${CONTEXTS}; do
  ls ./results/fitted_sites/${CONTEXT}/sites_cancer.txt | ruby bin/supplementary/snv_positions_in_motif_site.rb  --expand-region ${EXPAND_FLANKS_LENGTH}  --folder ./results/disruption_position_profile/${CONTEXT}/cancer
  ls ./results/fitted_sites/${CONTEXT}/sites_random_genome_*.txt | ruby bin/supplementary/snv_positions_in_motif_site.rb  --expand-region ${EXPAND_FLANKS_LENGTH}  --folder ./results/disruption_position_profile/${CONTEXT}/genome/
  ls ./results/fitted_sites/${CONTEXT}/sites_random_shuffle_*.txt | ruby bin/supplementary/snv_positions_in_motif_site.rb  --expand-region ${EXPAND_FLANKS_LENGTH}  --folder ./results/disruption_position_profile/${CONTEXT}/shuffle/
done
for CONTEXT in ${CONTEXTS}; do
  mkdir -p ./results/disruption_position_profile/${CONTEXT}/combined
  for MOTIF in `cat ./source_data/motif_names.txt`; do

    rm ./results/disruption_position_profile/${CONTEXT}/combined/${MOTIF}.txt
    echo -en "cancer\t" >> ./results/disruption_position_profile/${CONTEXT}/combined/${MOTIF}.txt
    cat ./results/disruption_position_profile/${CONTEXT}/cancer/${MOTIF}.txt >> ./results/disruption_position_profile/${CONTEXT}/combined/${MOTIF}.txt

    echo -en "genome\t" >> ./results/disruption_position_profile/${CONTEXT}/combined/${MOTIF}.txt
    cat ./results/disruption_position_profile/${CONTEXT}/genome/${MOTIF}.txt >> ./results/disruption_position_profile/${CONTEXT}/combined/${MOTIF}.txt

    echo -en "shuffle\t" >> ./results/disruption_position_profile/${CONTEXT}/combined/${MOTIF}.txt
    cat ./results/disruption_position_profile/${CONTEXT}/shuffle/${MOTIF}.txt >> ./results/disruption_position_profile/${CONTEXT}/combined/${MOTIF}.txt

    ruby -e 'File.write(ARGV[0], File.readlines(ARGV[0]).map{|l| l.chomp.split("\t")}.transpose.map{|l| l.join("\t") }.join("\n"))' -- ./results/disruption_position_profile/${CONTEXT}/combined/${MOTIF}.txt

  done
done
