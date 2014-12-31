# TODO: rename SUBSTITUTIONS_13Apr2012_snz_promoter_markup2.txt to smth more plausible (SNV_infos.txt);
# TODO: add download links (wget tasks)
# TODO: cleanup source_data folder; make folders for intermediates
# TODO: change order of actions to perform random filter and regulatory&context filtering on SNV stage, not sites stage

# 15 sec
ruby bin/preparations/markup_mutations_in_promoters.rb ./source_data/gene_tss.txt /home/ilya/iogen/cages/hg19/freeze1/hg19.cage_peak_tpm_ann.osc.txt  ./source_data/SUBSTITUTIONS_13Apr2012_snz.txt > ./source_data/SUBSTITUTIONS_13Apr2012_snz_promoter_markup3.txt

# 1 min
ruby bin/preparations/extract_snv_sequences.rb ./source_data/SUBSTITUTIONS_13Apr2012_snz_promoter_markup3.txt /home/ilya/iogen/genome/hg19/ --flank-length 25 > ./source_data/SNV_sequences.txt


# 5 min
ruby bin/preparations/filter_repeats.rb ./source_data/sites_cancer.txt > ./source_data/sites_cancer_wo_repeats.txt
# 40 min
ruby bin/preparations/filter_repeats.rb ./source_data/sites_random.txt > ./source_data/sites_random_wo_repeats.txt


# 40 sec each
ruby bin/preparations/filter_mutations.rb ./source_data/sites_cancer_wo_repeats.txt --contexts TCN --mutation-types promoter,intronic > ./source_data/sites_cancer_tpc.txt
ruby bin/preparations/filter_mutations.rb ./source_data/sites_cancer_wo_repeats.txt --contexts NCG --mutation-types promoter,intronic > ./source_data/sites_cancer_cpg.txt
ruby bin/preparations/filter_mutations.rb ./source_data/sites_cancer_wo_repeats.txt --mutation-types promoter,intronic > ./source_data/sites_cancer_any.txt
# 15 min each
ruby bin/preparations/filter_mutations.rb ./source_data/sites_random_wo_repeats.txt --contexts TCN --mutation-types promoter,intronic > ./source_data/sites_random_tpc.txt
ruby bin/preparations/filter_mutations.rb ./source_data/sites_random_wo_repeats.txt --contexts NCG --mutation-types promoter,intronic > ./source_data/sites_random_cpg.txt
ruby bin/preparations/filter_mutations.rb ./source_data/sites_random_wo_repeats.txt --mutation-types promoter,intronic > ./source_data/sites_random_any.txt

# 10 min
ruby fitting_random_sites.rb ./source_data/sites_cancer_cpg.txt ./source_data/sites_random_cpg.txt --fold 2 > ./source_data/sites_random_fitted_cpg.txt 2> >(tee ./source_data/sites_random_fitted_cpg.log >&2)
ruby fitting_random_sites.rb ./source_data/sites_cancer_tpc.txt ./source_data/sites_random_tpc.txt --fold 2 > ./source_data/sites_random_fitted_tpc.txt 2> >(tee ./source_data/sites_random_fitted_tpc.log >&2)
ruby fitting_random_sites.rb ./source_data/sites_cancer_any.txt ./source_data/sites_random_any.txt --fold 2 > ./source_data/sites_random_fitted_any.txt 2> >(tee ./source_data/sites_random_fitted_any.log >&2)

# ~1 minute
mkdir -p ./results/motif_statistics

ruby generate_all_motif_statistics.rb ./source_data/sites_cancer_cpg.txt ./results/motif_statistics/cpg/cancer.txt --pvalue 0.0005 --fold-change 5
ruby generate_all_motif_statistics.rb ./source_data/sites_random_cpg.txt ./results/motif_statistics/cpg/random.txt --pvalue 0.0005 --fold-change 5

ruby generate_all_motif_statistics.rb ./source_data/sites_cancer_tpc.txt ./results/motif_statistics/tpc/cancer.txt --pvalue 0.0005 --fold-change 5
ruby generate_all_motif_statistics.rb ./source_data/sites_random_tpc.txt ./results/motif_statistics/tpc/random.txt --pvalue 0.0005 --fold-change 5

ruby generate_all_motif_statistics.rb ./source_data/sites_cancer_any.txt ./results/motif_statistics/any/cancer.txt --pvalue 0.0005 --fold-change 5
ruby generate_all_motif_statistics.rb ./source_data/sites_random_any.txt ./results/motif_statistics/any/random.txt --pvalue 0.0005 --fold-change 5

ruby summary.rb ./results/motif_statistics/cpg/cancer.txt ./results/motif_statistics/cpg/random.txt > ./results/motif_statistics/cpg.csv
ruby summary.rb ./results/motif_statistics/tpc/cancer.txt ./results/motif_statistics/tpc/random.txt > ./results/motif_statistics/tpc.csv
ruby summary.rb ./results/motif_statistics/any/cancer.txt ./results/motif_statistics/any/random.txt > ./results/motif_statistics/any.csv

ruby filter_summary.rb ./results/motif_statistics/cpg.csv > ./results/motif_statistics/cpg_filtered.csv
ruby filter_summary.rb ./results/motif_statistics/tpc.csv > ./results/motif_statistics/tpc_filtered.csv
ruby filter_summary.rb ./results/motif_statistics/any.csv > ./results/motif_statistics/any_filtered.csv
