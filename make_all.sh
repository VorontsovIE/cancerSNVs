# TODO: add download links (wget tasks): SNV_infos_original.txt from paper; genome hg19; repeat masker track; ensembl gene markup; fantom peaks;
# TODO: cleanup source_data folder; make folders for intermediates
# TODO: change order of actions to perform random filter and regulatory&context filtering on SNV stage, not sites stage (?)

# http://feb2014.archive.ensembl.org/biomart/martview/b139ef98cf27cbd5649f7c5f6d3e2c0c?VIRTUALSCHEMANAME=default&ATTRIBUTES=hsapiens_gene_ensembl.default.structure.ensembl_gene_id|hsapiens_gene_ensembl.default.structure.ensembl_transcript_id|hsapiens_gene_ensembl.default.structure.exon_chrom_start|hsapiens_gene_ensembl.default.structure.exon_chrom_end|hsapiens_gene_ensembl.default.structure.is_constitutive|hsapiens_gene_ensembl.default.structure.rank|hsapiens_gene_ensembl.default.structure.phase|hsapiens_gene_ensembl.default.structure.cdna_coding_start|hsapiens_gene_ensembl.default.structure.cdna_coding_end|hsapiens_gene_ensembl.default.structure.genomic_coding_start|hsapiens_gene_ensembl.default.structure.genomic_coding_end|hsapiens_gene_ensembl.default.structure.ensembl_exon_id|hsapiens_gene_ensembl.default.structure.cds_start|hsapiens_gene_ensembl.default.structure.cds_end|hsapiens_gene_ensembl.default.structure.ensembl_peptide_id|hsapiens_gene_ensembl.default.structure.chromosome_name|hsapiens_gene_ensembl.default.structure.start_position|hsapiens_gene_ensembl.default.structure.end_position|hsapiens_gene_ensembl.default.structure.transcript_start|hsapiens_gene_ensembl.default.structure.transcript_end|hsapiens_gene_ensembl.default.structure.strand|hsapiens_gene_ensembl.default.structure.external_gene_id|hsapiens_gene_ensembl.default.structure.external_gene_db|hsapiens_gene_ensembl.default.structure.5_utr_start|hsapiens_gene_ensembl.default.structure.5_utr_end|hsapiens_gene_ensembl.default.structure.3_utr_start|hsapiens_gene_ensembl.default.structure.3_utr_end|hsapiens_gene_ensembl.default.structure.cds_length|hsapiens_gene_ensembl.default.structure.transcript_count|hsapiens_gene_ensembl.default.structure.description|hsapiens_gene_ensembl.default.structure.gene_biotype&FILTERS=

# 7 min
ruby bin/preparations/markup_mutations_in_promoters.rb ./source_data/gene_tss.txt /home/ilya/iogen/cages/hg19/freeze1/hg19.cage_peak_tpm_ann.osc.txt ./source_data/SNV_infos_original.txt /home/ilya/iogen/genome/hg19_repeatMasker "/home/ilya/iogen/genome/hg19_exons(ensembl,GRCh37.p13).txt" > ./source_data/SNV_infos.txt

mkdir -p results/intermediate
mkdir -p results/intermediate/site_subsets

# 1 min
ruby bin/preparations/extract_snv_sequences.rb ./source_data/SNV_infos.txt /home/ilya/iogen/genome/hg19/ --flank-length 25 > ./results/intermediate/SNV_sequences.txt

# 5 min
ruby bin/preparations/filter_repeats.rb ./source_data/sites_cancer.txt ./source_data/SNV_infos.txt /home/ilya/iogen/genome/hg19_repeatMasker > ./results/intermediate/site_subsets/cancer_wo_repeats.txt
# 40 min
ruby bin/preparations/filter_repeats.rb ./source_data/sites_random.txt ./source_data/SNV_infos.txt /home/ilya/iogen/genome/hg19_repeatMasker > ./results/intermediate/site_subsets/random_wo_repeats.txt


# 40 sec each
ruby bin/preparations/filter_mutations.rb ./source_data/SNV_infos.txt ./results/intermediate/site_subsets/cancer_wo_repeats.txt --contexts TCN --mutation-types promoter,intronic > ./results/intermediate/site_subsets/cancer_tpc.txt
ruby bin/preparations/filter_mutations.rb ./source_data/SNV_infos.txt ./results/intermediate/site_subsets/cancer_wo_repeats.txt --contexts NCG --mutation-types promoter,intronic > ./results/intermediate/site_subsets/cancer_cpg.txt
ruby bin/preparations/filter_mutations.rb ./source_data/SNV_infos.txt ./results/intermediate/site_subsets/cancer_wo_repeats.txt --mutation-types promoter,intronic > ./results/intermediate/site_subsets/cancer_any.txt
# 15 min each
ruby bin/preparations/filter_mutations.rb ./source_data/SNV_infos.txt ./results/intermediate/site_subsets/random_wo_repeats.txt --contexts TCN --mutation-types promoter,intronic > ./results/intermediate/site_subsets/random_tpc.txt
ruby bin/preparations/filter_mutations.rb ./source_data/SNV_infos.txt ./results/intermediate/site_subsets/random_wo_repeats.txt --contexts NCG --mutation-types promoter,intronic > ./results/intermediate/site_subsets/random_cpg.txt
ruby bin/preparations/filter_mutations.rb ./source_data/SNV_infos.txt ./results/intermediate/site_subsets/random_wo_repeats.txt --mutation-types promoter,intronic > ./results/intermediate/site_subsets/random_any.txt

# 10 min
ruby fitting_random_sites.rb ./results/intermediate/site_subsets/cancer_cpg.txt ./results/intermediate/site_subsets/random_cpg.txt --fold 2 > ./results/intermediate/site_subsets/fitted/random_cpg.txt 2> >(tee ./results/intermediate/site_subsets/fitted/log/random_cpg.log >&2)
ruby fitting_random_sites.rb ./results/intermediate/site_subsets/cancer_tpc.txt ./results/intermediate/site_subsets/random_tpc.txt --fold 2 > ./results/intermediate/site_subsets/fitted/random_tpc.txt 2> >(tee ./results/intermediate/site_subsets/fitted/log/random_tpc.log >&2)
ruby fitting_random_sites.rb ./results/intermediate/site_subsets/cancer_any.txt ./results/intermediate/site_subsets/random_any.txt --fold 2 > ./results/intermediate/site_subsets/fitted/random_any.txt 2> >(tee ./results/intermediate/site_subsets/fitted/log/random_any.log >&2)

# ~1 minute
mkdir -p ./results/motif_statistics

ruby generate_all_motif_statistics.rb ./results/intermediate/site_subsets/cancer_cpg.txt ./results/motif_statistics/cpg/cancer.txt ./source_data/motif_names.txt --pvalue 0.0005 --fold-change 5
ruby generate_all_motif_statistics.rb ./results/intermediate/site_subsets/fitted/random_cpg.txt ./results/motif_statistics/cpg/random.txt ./source_data/motif_names.txt --pvalue 0.0005 --fold-change 5

ruby generate_all_motif_statistics.rb ./results/intermediate/site_subsets/cancer_tpc.txt ./results/motif_statistics/tpc/cancer.txt ./source_data/motif_names.txt --pvalue 0.0005 --fold-change 5
ruby generate_all_motif_statistics.rb ./results/intermediate/site_subsets/fitted/random_tpc.txt ./results/motif_statistics/tpc/random.txt ./source_data/motif_names.txt --pvalue 0.0005 --fold-change 5

ruby generate_all_motif_statistics.rb ./results/intermediate/site_subsets/cancer_any.txt ./results/motif_statistics/any/cancer.txt ./source_data/motif_names.txt --pvalue 0.0005 --fold-change 5
ruby generate_all_motif_statistics.rb ./results/intermediate/site_subsets/fitted/random_any.txt ./results/motif_statistics/any/random.txt ./source_data/motif_names.txt --pvalue 0.0005 --fold-change 5

ruby summary.rb ./results/motif_statistics/cpg/cancer.txt ./results/motif_statistics/cpg/random.txt ./source_data/motif_names.txt ./source_data/hocomoco_genes_infos.csv > ./results/motif_statistics/cpg.csv
ruby summary.rb ./results/motif_statistics/tpc/cancer.txt ./results/motif_statistics/tpc/random.txt ./source_data/motif_names.txt ./source_data/hocomoco_genes_infos.csv > ./results/motif_statistics/tpc.csv
ruby summary.rb ./results/motif_statistics/any/cancer.txt ./results/motif_statistics/any/random.txt ./source_data/motif_names.txt ./source_data/hocomoco_genes_infos.csv > ./results/motif_statistics/any.csv

ruby filter_summary.rb ./results/motif_statistics/cpg.csv > ./results/motif_statistics/cpg_filtered.csv
ruby filter_summary.rb ./results/motif_statistics/tpc.csv > ./results/motif_statistics/tpc_filtered.csv
ruby filter_summary.rb ./results/motif_statistics/any.csv > ./results/motif_statistics/any_filtered.csv
