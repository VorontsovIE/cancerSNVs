#!/bin/bash

# # TODO: add download links (wget tasks): SNV_infos_original.txt from paper; genome hg19; repeat masker track; ensembl gene markup; fantom peaks;
# # TODO: cleanup source_data folder; make folders for intermediates
# # TODO: change order of actions to perform random filter and regulatory&context filtering on SNV stage, not sites stage (?)

# # http://feb2014.archive.ensembl.org/biomart/martview/b139ef98cf27cbd5649f7c5f6d3e2c0c?VIRTUALSCHEMANAME=default&ATTRIBUTES=hsapiens_gene_ensembl.default.structure.ensembl_gene_id|hsapiens_gene_ensembl.default.structure.ensembl_transcript_id|hsapiens_gene_ensembl.default.structure.exon_chrom_start|hsapiens_gene_ensembl.default.structure.exon_chrom_end|hsapiens_gene_ensembl.default.structure.is_constitutive|hsapiens_gene_ensembl.default.structure.rank|hsapiens_gene_ensembl.default.structure.phase|hsapiens_gene_ensembl.default.structure.cdna_coding_start|hsapiens_gene_ensembl.default.structure.cdna_coding_end|hsapiens_gene_ensembl.default.structure.genomic_coding_start|hsapiens_gene_ensembl.default.structure.genomic_coding_end|hsapiens_gene_ensembl.default.structure.ensembl_exon_id|hsapiens_gene_ensembl.default.structure.cds_start|hsapiens_gene_ensembl.default.structure.cds_end|hsapiens_gene_ensembl.default.structure.ensembl_peptide_id|hsapiens_gene_ensembl.default.structure.chromosome_name|hsapiens_gene_ensembl.default.structure.start_position|hsapiens_gene_ensembl.default.structure.end_position|hsapiens_gene_ensembl.default.structure.transcript_start|hsapiens_gene_ensembl.default.structure.transcript_end|hsapiens_gene_ensembl.default.structure.strand|hsapiens_gene_ensembl.default.structure.external_gene_id|hsapiens_gene_ensembl.default.structure.external_gene_db|hsapiens_gene_ensembl.default.structure.5_utr_start|hsapiens_gene_ensembl.default.structure.5_utr_end|hsapiens_gene_ensembl.default.structure.3_utr_start|hsapiens_gene_ensembl.default.structure.3_utr_end|hsapiens_gene_ensembl.default.structure.cds_length|hsapiens_gene_ensembl.default.structure.transcript_count|hsapiens_gene_ensembl.default.structure.description|hsapiens_gene_ensembl.default.structure.gene_biotype&FILTERS=

ln -s  /home/ilya/iogen/genome/hg19  ./source_data/genome
ln -s  /home/ilya/iogen/cages/hg19/freeze1/hg19.cage_peak_tpm_ann.osc.txt  ./source_data/cage_peaks.txt
ln -s  /home/ilya/iogen/genome/hg19_repeatMasker  ./source_data/genome_repeats
ln -s  /home/ilya/iogen/genome/hg19_exons\(ensembl\,GRCh37.p13\).txt  ./source_data/exons.txt
ln  /home/ilya/iogen/hocomoco  ./source_data/motif_collection
ln  /home/ilya/iogen/hocomoco-thresholds  ./source_data/motif_thresholds


SEQ_FOLDER=./results/sequences
SNV_FOLDER=./results/SNVs
CHUNK_FOLDER=./results/sequence_chunks

MOTIF_COLLECTION=./source_data/motif_collection
MOTIF_COLLECTION_THRESHOLDS=./source_data/motif_thresholds

NUMBER_OF_CORES=4


# Markup SNVs and filter regulatory only, non-duplicated SNVs
ruby bin/preparations/snv_markup.rb  ./source_data/SNV_infos_original.txt  ./source_data/cage_peaks.txt  ./source_data/genome_repeats  ./source_data/exons.txt  >  ./source_data/SNV_infos_marked_up.txt
ruby bin/preparations/filter_snv_infos.rb  ./source_data/SNV_infos_marked_up.txt  ./source_data/genome  >  ./source_data/SNV_infos_regulatory.txt


#  Later we work only with regulatory non-duplicated SNVs
ln  ./source_data/SNV_infos_regulatory.txt  ${SNV_FOLDER}/SNV_infos_cancer.txt

##################################

# Extract sequences around SNVs
ruby bin/preparations/extract_snv_sequences.rb  ${SNV_FOLDER}/SNV_infos_cancer.txt  ./source_data/genome  --flank-length 25  >  ${SEQ_FOLDER}/sequences_cancer.txt


# Generate random sequences from genome
for SEED in 13 15 17; do
  RANDOM_SEQUENCES=${SEQ_FOLDER}/sequences_random_genome_${SEED}.txt

  # generate random sequences
  ######## TODO: include genome folder as a script parameter !!!!!!!!!!
  ruby bin/preparations/generate_random_genome_sequences.rb  ${SNV_FOLDER}/SNV_infos_cancer.txt  ./source_data/exons.txt  --fold 10  --random-seed $SEED  >  $RANDOM_SEQUENCES

  # generate SNV infos
  ruby bin/preparations/create_snv_infos_by_genome_sequences.rb  $RANDOM_SEQUENCES  >  ${SNV_FOLDER}/SNV_infos_random_genome_${SEED}_non_marked_up.txt
  ruby bin/preparations/snv_markup.rb ${SNV_FOLDER}/SNV_infos_random_genome_${SEED}_non_marked_up.txt  ./source_data/cage_peaks.txt  ./source_data/genome_repeats  ./source_data/exons.txt  >  ${SNV_FOLDER}/SNV_infos_random_genome_${SEED}.txt
  rm ${SNV_FOLDER}/SNV_infos_random_genome_${SEED}_non_marked_up.txt
done


# Generate random sequences by shuffling contexts
for SEED in 135 137 139; do
  RANDOM_SEQUENCES=${SEQ_FOLDER}/sequences_random_shuffle_${SEED}.txt

  # generate random sequences
  ruby ./bin/preparations/shuffle_SNVs.rb  ${SEQ_FOLDER}/sequences_cancer.txt  --random-seed $SEED  --fold 10  >  $RANDOM_SEQUENCES
  # generate SNV infos
  ln  ${SNV_FOLDER}/SNV_infos_cancer.txt  ${SNV_FOLDER}/SNV_infos_random_shuffle_${SEED}.txt
done


##################################

./prepare_sequences_for_perfectosape_run.sh  $SEQ_FOLDER  $CHUNK_FOLDER  4
