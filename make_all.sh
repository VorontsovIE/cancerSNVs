# 1 min
ruby bin/preparations/extract_snv_sequences.rb source_data/SUBSTITUTIONS_13Apr2012_snz_promoter_markup2.txt /home/ilya/iogen/genome/hg19/ --flank-length 25 > source_data/SNPs.txt



# 5 min
ruby bin/preparations/filter_repeats.rb source_data/sites_cancer.txt > source_data/sites_cancer_wo_repeats.txt 
# 40 min
ruby bin/preparations/filter_repeats.rb source_data/sites_random.txt > source_data/sites_random_wo_repeats.txt 

# 15 min each
ruby bin/preparations/filter_mutations.rb source_data/sites_random_wo_repeats.txt --contexts TCN --mutation-types promoter,intronic > source_data/sites_random_tpc.txt
ruby bin/preparations/filter_mutations.rb source_data/sites_random_wo_repeats.txt --contexts NCG --mutation-types promoter,intronic > source_data/sites_random_cpg.txt
# 40 sec each
ruby bin/preparations/filter_mutations.rb source_data/sites_cancer_wo_repeats.txt --contexts TCN --mutation-types promoter,intronic > source_data/sites_cancer_tpc.txt 
ruby bin/preparations/filter_mutations.rb source_data/sites_cancer_wo_repeats.txt --contexts NCG --mutation-types promoter,intronic > source_data/sites_cancer_cpg.txt 

# 35 sec
ruby fitting_random_sites.rb source_data/sites_cancer_cpg.txt source_data/sites_random_cpg.txt > source_data/fitted_shuffle_cpg.txt 2> >(tee source_data/fitted_shuffle_cpg.log >&2)

