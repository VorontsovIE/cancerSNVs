# TODO: rename SUBSTITUTIONS_13Apr2012_snz_promoter_markup2.txt to smth more plausible (SNV_infos.txt); SNPs.txt --> SNV_sequences.txt
# TODO: add download links (wget tasks)
# TODO: cleanup source_data folder; make folders for intermediates
# TODO: change order of actions to perform random filter and regulatory&context filtering on SNV stage, not sites stage

# 15 sec
ruby bin/preparations/markup_mutations_in_promoters.rb source_data/gene_tss.txt source_data/SUBSTITUTIONS_13Apr2012_snz.txt > source_data/SUBSTITUTIONS_13Apr2012_snz_promoter_markup2.txt

# 1 min
ruby bin/preparations/extract_snv_sequences.rb source_data/SUBSTITUTIONS_13Apr2012_snz_promoter_markup2.txt /home/ilya/iogen/genome/hg19/ --flank-length 25 > source_data/SNPs.txt


# 5 min
ruby bin/preparations/filter_repeats.rb source_data/sites_cancer.txt > source_data/sites_cancer_wo_repeats.txt
# 40 min
ruby bin/preparations/filter_repeats.rb source_data/sites_random.txt > source_data/sites_random_wo_repeats.txt


# 40 sec each
ruby bin/preparations/filter_mutations.rb source_data/sites_cancer_wo_repeats.txt --contexts TCN --mutation-types promoter,intronic > source_data/sites_cancer_tpc.txt
ruby bin/preparations/filter_mutations.rb source_data/sites_cancer_wo_repeats.txt --contexts NCG --mutation-types promoter,intronic > source_data/sites_cancer_cpg.txt
# 15 min each
ruby bin/preparations/filter_mutations.rb source_data/sites_random_wo_repeats.txt --contexts TCN --mutation-types promoter,intronic > source_data/sites_random_tpc.txt
ruby bin/preparations/filter_mutations.rb source_data/sites_random_wo_repeats.txt --contexts NCG --mutation-types promoter,intronic > source_data/sites_random_cpg.txt


# 2 min
ruby fitting_random_sites.rb source_data/sites_cancer_cpg.txt source_data/sites_random_cpg.txt > source_data/sites_random_fitted_cpg.txt 2> >(tee source_data/sites_random_fitted_cpg.log >&2)
ruby fitting_random_sites.rb source_data/sites_cancer_tpc.txt source_data/sites_random_tpc.txt > source_data/sites_random_fitted_tpc.txt 2> >(tee source_data/sites_random_fitted_tpc.log >&2)

# ~1 minute
mkdir -p results/motif_statistics

ruby calculate_motif_statistics.rb source_data/sites_cancer_cpg.txt > results/motif_statistics/cancer_cpg_all.txt
ruby calculate_motif_statistics.rb source_data/sites_cancer_cpg.txt --site-after 0.0005 --emerged 5 > results/motif_statistics/cancer_cpg_sites_emerged.txt
ruby calculate_motif_statistics.rb source_data/sites_cancer_cpg.txt --site-before 0.0005 --disrupted 5 > results/motif_statistics/cancer_cpg_sites_disrupted.txt
ruby calculate_motif_statistics.rb source_data/sites_cancer_cpg.txt --site-after 0.0005 > results/motif_statistics/cancer_cpg_sites_after.txt
ruby calculate_motif_statistics.rb source_data/sites_cancer_cpg.txt --site-before 0.0005 > results/motif_statistics/cancer_cpg_sites_before.txt

ruby calculate_motif_statistics.rb source_data/sites_cancer_tpc.txt > results/motif_statistics/cancer_tpc_all.txt
ruby calculate_motif_statistics.rb source_data/sites_cancer_tpc.txt --site-after 0.0005 --emerged 5 > results/motif_statistics/cancer_tpc_sites_emerged.txt
ruby calculate_motif_statistics.rb source_data/sites_cancer_tpc.txt --site-before 0.0005 --disrupted 5 > results/motif_statistics/cancer_tpc_sites_disrupted.txt
ruby calculate_motif_statistics.rb source_data/sites_cancer_tpc.txt --site-after 0.0005 > results/motif_statistics/cancer_tpc_sites_after.txt
ruby calculate_motif_statistics.rb source_data/sites_cancer_tpc.txt --site-before 0.0005 > results/motif_statistics/cancer_tpc_sites_before.txt

ruby calculate_motif_statistics.rb source_data/sites_random_fitted_cpg.txt > results/motif_statistics/random_fitted_cpg_all.txt
ruby calculate_motif_statistics.rb source_data/sites_random_fitted_cpg.txt --site-after 0.0005 --emerged 5 > results/motif_statistics/random_fitted_cpg_sites_emerged.txt
ruby calculate_motif_statistics.rb source_data/sites_random_fitted_cpg.txt --site-before 0.0005 --disrupted 5 > results/motif_statistics/random_fitted_cpg_sites_disrupted.txt
ruby calculate_motif_statistics.rb source_data/sites_random_fitted_cpg.txt --site-after 0.0005 > results/motif_statistics/random_fitted_cpg_sites_after.txt
ruby calculate_motif_statistics.rb source_data/sites_random_fitted_cpg.txt --site-before 0.0005 > results/motif_statistics/random_fitted_cpg_sites_before.txt

ruby calculate_motif_statistics.rb source_data/sites_random_fitted_tpc.txt > results/motif_statistics/random_fitted_tpc_all.txt
ruby calculate_motif_statistics.rb source_data/sites_random_fitted_tpc.txt --site-after 0.0005 --emerged 5 > results/motif_statistics/random_fitted_tpc_sites_emerged.txt
ruby calculate_motif_statistics.rb source_data/sites_random_fitted_tpc.txt --site-before 0.0005 --disrupted 5 > results/motif_statistics/random_fitted_tpc_sites_disrupted.txt
ruby calculate_motif_statistics.rb source_data/sites_random_fitted_tpc.txt --site-after 0.0005 > results/motif_statistics/random_fitted_tpc_sites_after.txt
ruby calculate_motif_statistics.rb source_data/sites_random_fitted_tpc.txt --site-before 0.0005 > results/motif_statistics/random_fitted_tpc_sites_before.txt
