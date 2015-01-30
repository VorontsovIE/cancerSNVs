#!/bin/bash

# SEQ_FOLDER, SNV_FOLDER and CHUNK_FOLDER should be specified in environment variables

cd "$(dirname "$0")"

# PerfectosAPE package
# wget http://opera.autosome.ru/downloads/ape.jar -O ape.jar

# Genome assembly (strictly hg19)
# wget --timestamping 'ftp://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz' -O chromFa.tar.gz
# tar -zxf chromFa.tar.gz
# TODO: rename folder
# TODO: decode FASTA into plain sequences (remove possible headers and newline breaks; rename *.fa --> *.plain)

# Ensembl exons from (also hg19 version)
# # http://feb2014.archive.ensembl.org/biomart/martview/b139ef98cf27cbd5649f7c5f6d3e2c0c?VIRTUALSCHEMANAME=default&ATTRIBUTES=hsapiens_gene_ensembl.default.structure.ensembl_gene_id|hsapiens_gene_ensembl.default.structure.ensembl_transcript_id|hsapiens_gene_ensembl.default.structure.exon_chrom_start|hsapiens_gene_ensembl.default.structure.exon_chrom_end|hsapiens_gene_ensembl.default.structure.is_constitutive|hsapiens_gene_ensembl.default.structure.rank|hsapiens_gene_ensembl.default.structure.phase|hsapiens_gene_ensembl.default.structure.cdna_coding_start|hsapiens_gene_ensembl.default.structure.cdna_coding_end|hsapiens_gene_ensembl.default.structure.genomic_coding_start|hsapiens_gene_ensembl.default.structure.genomic_coding_end|hsapiens_gene_ensembl.default.structure.ensembl_exon_id|hsapiens_gene_ensembl.default.structure.cds_start|hsapiens_gene_ensembl.default.structure.cds_end|hsapiens_gene_ensembl.default.structure.ensembl_peptide_id|hsapiens_gene_ensembl.default.structure.chromosome_name|hsapiens_gene_ensembl.default.structure.start_position|hsapiens_gene_ensembl.default.structure.end_position|hsapiens_gene_ensembl.default.structure.transcript_start|hsapiens_gene_ensembl.default.structure.transcript_end|hsapiens_gene_ensembl.default.structure.strand|hsapiens_gene_ensembl.default.structure.external_gene_id|hsapiens_gene_ensembl.default.structure.external_gene_db|hsapiens_gene_ensembl.default.structure.5_utr_start|hsapiens_gene_ensembl.default.structure.5_utr_end|hsapiens_gene_ensembl.default.structure.3_utr_start|hsapiens_gene_ensembl.default.structure.3_utr_end|hsapiens_gene_ensembl.default.structure.cds_length|hsapiens_gene_ensembl.default.structure.transcript_count|hsapiens_gene_ensembl.default.structure.description|hsapiens_gene_ensembl.default.structure.gene_biotype&FILTERS=

ln -sf  /home/ilya/iogen/genome/hg19  ./source_data/genome
ln -sf  /home/ilya/iogen/genome/hg19_exons\(ensembl\,GRCh37.p13\).txt  ./source_data/exons.txt
ln -f  /home/ilya/iogen/hocomoco  ./source_data/motif_collection
ln -f  /home/ilya/iogen/hocomoco-thresholds  ./source_data/motif_thresholds

mkdir -p  $SEQ_FOLDER  $SNV_FOLDER  $CHUNK_FOLDER

# Markup SNVs and filter regulatory only, non-duplicated SNVs
ruby bin/preparations/snv_markup.rb  ./source_data/SNV_infos_original.txt  ./source_data/exons.txt  >  ./source_data/SNV_infos_marked_up.txt
ruby bin/preparations/filter_snv_infos.rb  ./source_data/SNV_infos_marked_up.txt  ./source_data/genome  >  ./source_data/SNV_infos_regulatory.txt


#  Later we work only with regulatory non-duplicated SNVs
ln  ./source_data/SNV_infos_regulatory.txt  ${SNV_FOLDER}/SNV_infos_cancer.txt

##################################

# Extract sequences around SNVs
ruby bin/preparations/extract_snv_sequences.rb  ${SNV_FOLDER}/SNV_infos_cancer.txt  ./source_data/genome  --flank-length 25  >  ${SEQ_FOLDER}/sequences_cancer.txt

# Generate random sequences
./prepare_random_sequences.sh

# Split sequences into equal chunks in order to run chunks in parallel
NUMBER_OF_CORES=4  ./prepare_sequences_for_perfectosape_run.sh
