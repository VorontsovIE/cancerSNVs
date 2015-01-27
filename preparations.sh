#!/bin/bash

# # TODO: add download links (wget tasks): SNV_infos_original.txt from paper; genome hg19; repeat masker track; ensembl gene markup; fantom peaks;
# # TODO: cleanup source_data folder; make folders for intermediates
# # TODO: change order of actions to perform random filter and regulatory&context filtering on SNV stage, not sites stage (?)

# # http://feb2014.archive.ensembl.org/biomart/martview/b139ef98cf27cbd5649f7c5f6d3e2c0c?VIRTUALSCHEMANAME=default&ATTRIBUTES=hsapiens_gene_ensembl.default.structure.ensembl_gene_id|hsapiens_gene_ensembl.default.structure.ensembl_transcript_id|hsapiens_gene_ensembl.default.structure.exon_chrom_start|hsapiens_gene_ensembl.default.structure.exon_chrom_end|hsapiens_gene_ensembl.default.structure.is_constitutive|hsapiens_gene_ensembl.default.structure.rank|hsapiens_gene_ensembl.default.structure.phase|hsapiens_gene_ensembl.default.structure.cdna_coding_start|hsapiens_gene_ensembl.default.structure.cdna_coding_end|hsapiens_gene_ensembl.default.structure.genomic_coding_start|hsapiens_gene_ensembl.default.structure.genomic_coding_end|hsapiens_gene_ensembl.default.structure.ensembl_exon_id|hsapiens_gene_ensembl.default.structure.cds_start|hsapiens_gene_ensembl.default.structure.cds_end|hsapiens_gene_ensembl.default.structure.ensembl_peptide_id|hsapiens_gene_ensembl.default.structure.chromosome_name|hsapiens_gene_ensembl.default.structure.start_position|hsapiens_gene_ensembl.default.structure.end_position|hsapiens_gene_ensembl.default.structure.transcript_start|hsapiens_gene_ensembl.default.structure.transcript_end|hsapiens_gene_ensembl.default.structure.strand|hsapiens_gene_ensembl.default.structure.external_gene_id|hsapiens_gene_ensembl.default.structure.external_gene_db|hsapiens_gene_ensembl.default.structure.5_utr_start|hsapiens_gene_ensembl.default.structure.5_utr_end|hsapiens_gene_ensembl.default.structure.3_utr_start|hsapiens_gene_ensembl.default.structure.3_utr_end|hsapiens_gene_ensembl.default.structure.cds_length|hsapiens_gene_ensembl.default.structure.transcript_count|hsapiens_gene_ensembl.default.structure.description|hsapiens_gene_ensembl.default.structure.gene_biotype&FILTERS=

GENOME_FOLDER = /home/ilya/iogen/genome/hg19
CAGE_PEAKS = /home/ilya/iogen/cages/hg19/freeze1/hg19.cage_peak_tpm_ann.osc.txt
REPEAT_MASKER_FOLDER = /home/ilya/iogen/genome/hg19_repeatMasker
ENSEMBL_EXONS = /home/ilya/iogen/genome/hg19_exons\(ensembl\,GRCh37.p13\).txt

SNV_INFOS = ./source_data/SNV_infos_regulatory.txt
INTERMEDIATE_FOLDER = ./results/intermediate
SNV_SEQUENCES = $INTERMEDIATE_FOLDER/SNV_sequences.txt

MOTIF_COLLECTION = ./hocomoco
MOTIF_COLLECTION_THRESHOLDS = ./thresholds


mkdir -p $INTERMEDIATE_FOLDER

# Markup SNVs, filter regulatory only and extract sequences
ruby bin/preparations/snv_markup.rb  ./source_data/SNV_infos_original.txt  $CAGE_PEAKS  $REPEAT_MASKER_FOLDER  $ENSEMBL_EXONS  >  $INTERMEDIATE_FOLDER/SNV_infos_marked_up.txt
ruby bin/preparations/filter_snv_infos.rb  $INTERMEDIATE_FOLDER/SNV_infos_marked_up.txt  $GENOME_FOLDER  >  $SNV_INFOS
ruby bin/preparations/extract_snv_sequences.rb  $SNV_INFOS  $GENOME_FOLDER --flank-length 25  >  $SNV_SEQUENCES


CHUNK_FOLDER = $INTERMEDIATE_FOLDER/chunks
mkdir -p  $CHUNK_FOLDER

# java -cp ape.jar ru.autosome.perfectosape.SNPScan  $MOTIF_COLLECTION  $CHUNK_FOLDER/SNV_sequences.txt  --fold-change-cutoff 1  --precalc $MOTIF_COLLECTION_THRESHOLDS  >  $CHUNK_FOLDER/sites/SNV_sequences.txt

# Generate random sequences from genome
for SEED in 13 15 17
do
  RANDOM_SEQUENCES = $INTERMEDIATE_FOLDER/background_genome_sequences_${SEED}.txt
  CHUNK_PREFIX = $CHUNK_FOLDER/background_genome_sequences_${SEED}_chunk_

  ruby bin/preparations/generate_random_genome_sequences.rb  $SNV_INFOS  $ENSEMBL_EXONS  --fold 10  --random-seed $SEED  >  $RANDOM_SEQUENCES
  ruby bin/preparations/create_snv_infos_by_genome_sequences.rb $RANDOM_SEQUENCES > $INTERMEDIATE_FOLDER/background_genome_SNV_${SEED}_non_marked_up.txt
  ruby bin/preparations/snv_markup.rb $INTERMEDIATE_FOLDER/background_genome_SNV_${SEED}_non_marked_up.txt  $CAGE_PEAKS  $REPEAT_MASKER_FOLDER  $ENSEMBL_EXONS  >  $INTERMEDIATE_FOLDER/background_genome_SNV_${SEED}.txt

  split  --additional-suffix .txt  --suffix-length 1  --number l/4  $RANDOM_SEQUENCES  $CHUNK_PREFIX

  # for SUFFIX in a b c d
  # do
  #   java -cp ape.jar ru.autosome.perfectosape.SNPScan  $MOTIF_COLLECTION  ${CHUNK_PREFIX}${SUFFIX}.txt  --fold-change-cutoff 1  --precalc $MOTIF_COLLECTION_THRESHOLDS  >  $CHUNK_FOLDER/sites/background_genome_sites_${SEED}_chunk_#{SUFFIX}.txt
  # done

  # rm $INTERMEDIATE_FOLDER/random_genome_sites_${SEED}.txt
  # for SUFFIX in a b c d
  # do
  #   cat $CHUNK_FOLDER/sites/background_genome_sites_${SEED}_chunk_#{SUFFIX}.txt >> $INTERMEDIATE_FOLDER/random_genome_sites_${SEED}.txt
  # done
  # for SUFFIX in a b c d
  # do
  #   rm $CHUNK_FOLDER/sites/background_genome_sites_${SEED}_chunk_#{SUFFIX}.txt
  # done
done

# Generate random sequences by shuffling contexts
for SEED in 135 137 139
do
  RANDOM_SEQUENCES = $INTERMEDIATE_FOLDER/sequences_shuffled_${SEED}.txt
  CHUNK_PREFIX = $CHUNK_FOLDER/sequences_shuffled_${SEED}_chunk_

  ruby ./bin/preparations/shuffle_SNVs.rb  $SNV_SEQUENCES  --random-seed $SEED  --fold 10  >  $RANDOM_SEQUENCES

  split  --additional-suffix .txt  --suffix-length 1  --number l/4  $RANDOM_SEQUENCES  $CHUNK_PREFIX
  # for SUFFIX in a b c d
  # do
  #   java -cp ape.jar ru.autosome.perfectosape.SNPScan  $MOTIF_COLLECTION  ${CHUNK_PREFIX}${SUFFIX}.txt  --fold-change-cutoff 1  --precalc $MOTIF_COLLECTION_THRESHOLDS  >  $CHUNK_FOLDER/sites/sites_shuffled_${SEED}_chunk_#{SUFFIX}.txt
  # done
  # rm $INTERMEDIATE_FOLDER/random_shuffled_sites_${SEED}.txt
  # for SUFFIX in a b c d
  # do
  #   cat $CHUNK_FOLDER/sites/sites_shuffled_${SEED}_chunk_#{SUFFIX}.txt >> $INTERMEDIATE_FOLDER/random_shuffled_sites_${SEED}.txt
  # done
  # for SUFFIX in a b c d
  # do
  #   rm $CHUNK_FOLDER/sites/sites_shuffled_${SEED}_chunk_#{SUFFIX}.txt
  # done
done
