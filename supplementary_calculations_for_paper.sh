
# # CONSIDERED_MOTIFS=AP2D_f1,BATF_si,FOSB_f1,FOSL2_f1,HXD13_f1,JUNB_f1,JUND_f1,KLF15_f1,KLF1_f1,PBX2_f1,PKNX1_si,SMRC1_f1,SP3_f1,TFCP2_f1,UBIP1_f1
# CONSIDERED_MOTIFS=AP2A_f2,AP2B_f1,AP2D_f1,ATF5_si,CEBPA_do,CEBPB_f1,EGR1_f2,KLF15_f1,KLF1_f1,MAZ_f1,PAX5_si,PURA_f1,SP1_f2,SP2_si,SP3_f1,SP4_f1,ZBT7B_si,ZN148_si,

# ls ./results/fitted_sites/any/sites_cancer.txt | ruby bin/supplementary/snv_positions_in_motif_site.rb --motif ${CONSIDERED_MOTIFS} --folder ./results/disruption_position_profile/any/cancer
# ls ./results/fitted_sites/any/sites_random_genome_*.txt | ruby bin/supplementary/snv_positions_in_motif_site.rb --motif ${CONSIDERED_MOTIFS} --folder ./results/disruption_position_profile/any/genome/
# ls ./results/fitted_sites/any/sites_random_shuffle_*.txt | ruby bin/supplementary/snv_positions_in_motif_site.rb --motif ${CONSIDERED_MOTIFS} --folder ./results/disruption_position_profile/any/shuffle/ 

# ls ./results/fitted_sites/any/sites_random_*.txt | ruby bin/supplementary/snv_positions_in_motif_site.rb --motif ${CONSIDERED_MOTIFS} --folder ./results/disruption_position_profile/any/shuffle_and_genome/


ruby significances_by_samples.rb ${RANDOM_VARIANTS} > results/significances_by_sample_all.csv
ruby significances_by_samples.rb ${RANDOM_VARIANTS} > results/significances_by_sample_genome_only.csv
