#!/bin/bash

# # TODO: add download links (wget tasks): SNV_infos_original.txt from paper; genome hg19; repeat masker track; ensembl gene markup; fantom peaks;
# # TODO: cleanup source_data folder; make folders for intermediates
# # TODO: change order of actions to perform random filter and regulatory&context filtering on SNV stage, not sites stage (?)

# # http://feb2014.archive.ensembl.org/biomart/martview/b139ef98cf27cbd5649f7c5f6d3e2c0c?VIRTUALSCHEMANAME=default&ATTRIBUTES=hsapiens_gene_ensembl.default.structure.ensembl_gene_id|hsapiens_gene_ensembl.default.structure.ensembl_transcript_id|hsapiens_gene_ensembl.default.structure.exon_chrom_start|hsapiens_gene_ensembl.default.structure.exon_chrom_end|hsapiens_gene_ensembl.default.structure.is_constitutive|hsapiens_gene_ensembl.default.structure.rank|hsapiens_gene_ensembl.default.structure.phase|hsapiens_gene_ensembl.default.structure.cdna_coding_start|hsapiens_gene_ensembl.default.structure.cdna_coding_end|hsapiens_gene_ensembl.default.structure.genomic_coding_start|hsapiens_gene_ensembl.default.structure.genomic_coding_end|hsapiens_gene_ensembl.default.structure.ensembl_exon_id|hsapiens_gene_ensembl.default.structure.cds_start|hsapiens_gene_ensembl.default.structure.cds_end|hsapiens_gene_ensembl.default.structure.ensembl_peptide_id|hsapiens_gene_ensembl.default.structure.chromosome_name|hsapiens_gene_ensembl.default.structure.start_position|hsapiens_gene_ensembl.default.structure.end_position|hsapiens_gene_ensembl.default.structure.transcript_start|hsapiens_gene_ensembl.default.structure.transcript_end|hsapiens_gene_ensembl.default.structure.strand|hsapiens_gene_ensembl.default.structure.external_gene_id|hsapiens_gene_ensembl.default.structure.external_gene_db|hsapiens_gene_ensembl.default.structure.5_utr_start|hsapiens_gene_ensembl.default.structure.5_utr_end|hsapiens_gene_ensembl.default.structure.3_utr_start|hsapiens_gene_ensembl.default.structure.3_utr_end|hsapiens_gene_ensembl.default.structure.cds_length|hsapiens_gene_ensembl.default.structure.transcript_count|hsapiens_gene_ensembl.default.structure.description|hsapiens_gene_ensembl.default.structure.gene_biotype&FILTERS=

GENOME_FOLDER = /home/ilya/iogen/genome/hg19
CAGE_PEAKS = /home/ilya/iogen/cages/hg19/freeze1/hg19.cage_peak_tpm_ann.osc.txt
REPEAT_MASKER_FOLDER = /home/ilya/iogen/genome/hg19_repeatMasker
ENSEMBL_EXONS = /home/ilya/iogen/genome/hg19_exons\(ensembl\,GRCh37.p13\).txt

# Markup SNVs, filter regulatory only and extract sequences
# mkdir -p results/intermediate
# ruby bin/preparations/snv_markup.rb $CAGE_PEAKS ./source_data/SNV_infos_original.txt $REPEAT_MASKER_FOLDER $ENSEMBL_EXONS > ./source_data/SNV_infos_full.txt
# ruby bin/preparations/filter_snv_infos.rb source_data/SNV_infos_full.txt  $GENOME_FOLDER  >  source_data/SNV_infos_regulatory.txt
# ruby bin/preparations/extract_snv_sequences.rb ./source_data/SNV_infos_regulatory.txt  $GENOME_FOLDER --flank-length 25 > ./results/intermediate/SNV_sequences.txt


# # Generate random sequences from genome
# ruby bin/preparations/generate_random_genome_sequences.rb source_data/SNV_infos_regulatory.txt  $ENSEMBL_EXONS  --fold 10  --random-seed 13  >  ./results/intermediate/genome_seqs_seed_13.txt
# ruby bin/preparations/generate_random_genome_sequences.rb source_data/SNV_infos_regulatory.txt  $ENSEMBL_EXONS  --fold 10  --random-seed 15  >  ./results/intermediate/genome_seqs_seed_15.txt
# ruby bin/preparations/generate_random_genome_sequences.rb source_data/SNV_infos_regulatory.txt  $ENSEMBL_EXONS  --fold 10  --random-seed 17  >  ./results/intermediate/genome_seqs_seed_17.txt
#
# ruby bin/preparations/create_snv_infos_by_genome_sequences.rb ./results/intermediate/genome_seqs_seed_13.txt > ./results/intermediate/SNV_infos_genome_seqs_seed_13_non_marked_up.txt
# ruby bin/preparations/create_snv_infos_by_genome_sequences.rb ./results/intermediate/genome_seqs_seed_15.txt > ./results/intermediate/SNV_infos_genome_seqs_seed_15_non_marked_up.txt
# ruby bin/preparations/create_snv_infos_by_genome_sequences.rb ./results/intermediate/genome_seqs_seed_17.txt > ./results/intermediate/SNV_infos_genome_seqs_seed_17_non_marked_up.txt
#
# ruby bin/preparations/snv_markup.rb /home/ilya/iogen/cages/hg19/freeze1/hg19.cage_peak_tpm_ann.osc.txt ./results/intermediate/SNV_infos_genome_seqs_seed_13_non_marked_up.txt $REPEAT_MASKER_FOLDER "/home/ilya/iogen/genome/hg19_exons(ensembl,GRCh37.p13).txt" > ./results/intermediate/SNV_infos_genome_seqs_seed_13.txt
# ruby bin/preparations/snv_markup.rb /home/ilya/iogen/cages/hg19/freeze1/hg19.cage_peak_tpm_ann.osc.txt ./results/intermediate/SNV_infos_genome_seqs_seed_15_non_marked_up.txt $REPEAT_MASKER_FOLDER "/home/ilya/iogen/genome/hg19_exons(ensembl,GRCh37.p13).txt" > ./results/intermediate/SNV_infos_genome_seqs_seed_15.txt
# ruby bin/preparations/snv_markup.rb /home/ilya/iogen/cages/hg19/freeze1/hg19.cage_peak_tpm_ann.osc.txt ./results/intermediate/SNV_infos_genome_seqs_seed_17_non_marked_up.txt $REPEAT_MASKER_FOLDER "/home/ilya/iogen/genome/hg19_exons(ensembl,GRCh37.p13).txt" > ./results/intermediate/SNV_infos_genome_seqs_seed_17.txt

# # Generate random sequences by shuffling contexts
# ruby ./bin/preparations/shuffle_SNVs.rb ./results/intermediate/SNV_sequences.txt --random-seed 135 --fold 10 > ./results/intermediate/sequences_shuffled_seed_135.txt
# ruby ./bin/preparations/shuffle_SNVs.rb ./results/intermediate/SNV_sequences.txt --random-seed 137 --fold 10 > ./results/intermediate/sequences_shuffled_seed_137.txt
# ruby ./bin/preparations/shuffle_SNVs.rb ./results/intermediate/SNV_sequences.txt --random-seed 139 --fold 10 > ./results/intermediate/sequences_shuffled_seed_139.txt

# mkdir -p ./results/intermediate/chunks/
# split --additional-suffix .txt --suffix-length 1 --number l/4  ./results/intermediate/sequences_shuffled_seed_135.txt ./results/intermediate/chunks/sequences_shuffled_seed_135_chunk_
# split --additional-suffix .txt --suffix-length 1 --number l/4  ./results/intermediate/sequences_shuffled_seed_137.txt ./results/intermediate/chunks/sequences_shuffled_seed_137_chunk_
# split --additional-suffix .txt --suffix-length 1 --number l/4  ./results/intermediate/sequences_shuffled_seed_139.txt ./results/intermediate/chunks/sequences_shuffled_seed_139_chunk_
# split --additional-suffix .txt --suffix-length 1 --number l/4  ./results/intermediate/genome_seqs_seed_13.txt ./results/intermediate/chunks/genome_seqs_seed_13_chunk_
# split --additional-suffix .txt --suffix-length 1 --number l/4  ./results/intermediate/genome_seqs_seed_15.txt ./results/intermediate/chunks/genome_seqs_seed_15_chunk_
# split --additional-suffix .txt --suffix-length 1 --number l/4  ./results/intermediate/genome_seqs_seed_17.txt ./results/intermediate/chunks/genome_seqs_seed_17_chunk_


# mkdir -p ./chunks/sites/
# java -cp ape.jar ru.autosome.perfectosape.SNPScan ./hocomoco/ ./chunks/sequences_shuffled_seed_135_chunk_a.txt --fold-change-cutoff 1 --precalc ./thresholds/> ./chunks/sites/sequences_shuffled_seed_135_chunk_a.txt
# java -cp ape.jar ru.autosome.perfectosape.SNPScan ./hocomoco/ ./chunks/sequences_shuffled_seed_137_chunk_a.txt --fold-change-cutoff 1 --precalc ./thresholds/> ./chunks/sites/sequences_shuffled_seed_137_chunk_a.txt
# java -cp ape.jar ru.autosome.perfectosape.SNPScan ./hocomoco/ ./chunks/sequences_shuffled_seed_139_chunk_a.txt --fold-change-cutoff 1 --precalc ./thresholds/> ./chunks/sites/sequences_shuffled_seed_139_chunk_a.txt

# java -cp ape.jar ru.autosome.perfectosape.SNPScan ./hocomoco/ ./chunks/sequences_shuffled_seed_135_chunk_b.txt --fold-change-cutoff 1 --precalc ./thresholds/> ./chunks/sites/sequences_shuffled_seed_135_chunk_b.txt
# java -cp ape.jar ru.autosome.perfectosape.SNPScan ./hocomoco/ ./chunks/sequences_shuffled_seed_137_chunk_b.txt --fold-change-cutoff 1 --precalc ./thresholds/> ./chunks/sites/sequences_shuffled_seed_137_chunk_b.txt
# java -cp ape.jar ru.autosome.perfectosape.SNPScan ./hocomoco/ ./chunks/sequences_shuffled_seed_139_chunk_b.txt --fold-change-cutoff 1 --precalc ./thresholds/> ./chunks/sites/sequences_shuffled_seed_139_chunk_b.txt

# java -cp ape.jar ru.autosome.perfectosape.SNPScan ./hocomoco/ ./chunks/sequences_shuffled_seed_135_chunk_c.txt --fold-change-cutoff 1 --precalc ./thresholds/> ./chunks/sites/sequences_shuffled_seed_135_chunk_c.txt
# java -cp ape.jar ru.autosome.perfectosape.SNPScan ./hocomoco/ ./chunks/sequences_shuffled_seed_137_chunk_c.txt --fold-change-cutoff 1 --precalc ./thresholds/> ./chunks/sites/sequences_shuffled_seed_137_chunk_c.txt
# java -cp ape.jar ru.autosome.perfectosape.SNPScan ./hocomoco/ ./chunks/sequences_shuffled_seed_139_chunk_c.txt --fold-change-cutoff 1 --precalc ./thresholds/> ./chunks/sites/sequences_shuffled_seed_139_chunk_c.txt

# java -cp ape.jar ru.autosome.perfectosape.SNPScan ./hocomoco/ ./chunks/sequences_shuffled_seed_135_chunk_d.txt --fold-change-cutoff 1 --precalc ./thresholds/> ./chunks/sites/sequences_shuffled_seed_135_chunk_d.txt
# java -cp ape.jar ru.autosome.perfectosape.SNPScan ./hocomoco/ ./chunks/sequences_shuffled_seed_137_chunk_d.txt --fold-change-cutoff 1 --precalc ./thresholds/> ./chunks/sites/sequences_shuffled_seed_137_chunk_d.txt
# java -cp ape.jar ru.autosome.perfectosape.SNPScan ./hocomoco/ ./chunks/sequences_shuffled_seed_139_chunk_d.txt --fold-change-cutoff 1 --precalc ./thresholds/> ./chunks/sites/sequences_shuffled_seed_139_chunk_d.txt

# java -cp ape.jar ru.autosome.perfectosape.SNPScan ./hocomoco/ ./chunks/SNV_sequences.txt --fold-change-cutoff 1 --precalc ./thresholds/> ./chunks/sites/SNV_sequences.txt

