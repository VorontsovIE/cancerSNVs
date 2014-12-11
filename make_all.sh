# 1 min
ruby bin/preparations/extract_snv_sequences.rb source_data/SUBSTITUTIONS_13Apr2012_snz_promoter_markup2.txt /home/ilya/iogen/genome/hg19/ --flank-length 25 > source_data/SNPs.txt

# 15 min each
ruby bin/preparations/filter_mutations.rb source_data/mutated_sites_shuffled.txt --contexts TCN --mutation-types promoter,intronic > source_data/mutated_sites_shuffled_tpc.txt
ruby bin/preparations/filter_mutations.rb source_data/mutated_sites_shuffled.txt --contexts NCG --mutation-types promoter,intronic > source_data/mutated_sites_shuffled_cpg.txt
# 40 sec each
ruby bin/preparations/filter_mutations.rb source_data/cancer_SNPs.txt --contexts TCN --mutation-types promoter,intronic > source_data/cancer_SNPs_tpc.txt 
ruby bin/preparations/filter_mutations.rb source_data/cancer_SNPs.txt --contexts NCG --mutation-types promoter,intronic > source_data/cancer_SNPs_cpg.txt 

# 35 sec
ruby fitting_random_sites.rb source_data/cancer_SNPs_cpg.txt source_data/mutated_sites_shuffled_cpg.txt > source_data/fitted_shuffle_cpg.txt 2> >(tee source_data/fitted_shuffle_cpg.log >&2)

