CONSIDERED_MOTIFS=AP2D_f1,BATF_si,FOSB_f1,FOSL2_f1,HXD13_f1,JUNB_f1,JUND_f1,KLF15_f1,KLF1_f1,PBX2_f1,PKNX1_si,SMRC1_f1,SP3_f1,TFCP2_f1,UBIP1_f1

ls ./results/fitted_sites/any/sites_cancer.txt | ruby bin/supplementary/snv_positions_in_motif_site.rb --motif ${CONSIDERED_MOTIFS} --folder ./results/disruption_position_profile/cancer

ls ./results/fitted_sites/any/sites_random_shuffle_*.txt | ruby bin/supplementary/snv_positions_in_motif_site.rb --motif ${CONSIDERED_MOTIFS} --folder ./results/disruption_position_profile/shuffle/ 

ls ./results/fitted_sites/any/sites_random_genome_*.txt | ruby bin/supplementary/snv_positions_in_motif_site.rb --motif ${CONSIDERED_MOTIFS} --folder ./results/disruption_position_profile/genome/

ls ./results/fitted_sites/any/sites_random_*.txt | ruby bin/supplementary/snv_positions_in_motif_site.rb --motif ${CONSIDERED_MOTIFS} --folder ./results/disruption_position_profile/shuffle_and_genome/


ruby significances_by_samples.rb > results/significances_by_sample_all.csv
ruby significances_by_samples.rb > results/significances_by_sample_genome_only.csv
