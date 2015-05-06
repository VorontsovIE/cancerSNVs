#!/bin/bash

# RESULTS_FOLDER, SEQ_FOLDER, SNV_FOLDER and CHUNK_FOLDER should be specified in environment variables

cd "$(dirname "$0")"

# Cancer somatic SNVs from "Mutational Processes Molding the Genomes of 21 Breast Cancers. Nik-Zainal et.al."
wget ftp://ftp.sanger.ac.uk/pub/cancer/Nik-ZainalEtAl/SUBSTITUTIONS_13Apr2012_snz.txt -O ./source_data/SNV_infos_original.txt

# Cancer somatic SNVs from "Signatures of mutational processes in human cancer. Alexandrov et al."
wget --recursive --directory-prefix='./source_data/AlexandrovEtAl' ftp://ftp.sanger.ac.uk/pub/cancer/AlexandrovEtAl/
mv ./source_data/AlexandrovEtAl/ftp.sanger.ac.uk/pub/cancer/AlexandrovEtAl/* ./source_data/AlexandrovEtAl/
rm -r ./source_data/AlexandrovEtAl/ftp.sanger.ac.uk/

# PerfectosAPE package
wget http://opera.autosome.ru/downloads/ape.jar -O ape.jar

# Genome assembly (strictly hg19)
GENOME_FOLDER='/home/ilya/genomes/Ensembl-GRCh37.p13'
mkdir -p ${GENOME_FOLDER}
wget --directory-prefix=${GENOME_FOLDER} ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/dna/README
wget --directory-prefix=${GENOME_FOLDER} ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.75.dna_sm.chromosome.*.fa.gz
gzip -d ${GENOME_FOLDER}/*.fa.gz
ruby bin/preparations/convert_fasta_to_plain.rb ${GENOME_FOLDER}
rm ${GENOME_FOLDER}/*.fa

# Download ensembl exons table
wget -O /home/ilya/iogen/genome/hg19_exons\(ensembl\,GRCh37.p13\).txt 'http://feb2014.archive.ensembl.org/biomart/martservice?query=<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query><Query  virtualSchemaName = "default" formatter = "TSV" header = "1" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" ><Dataset name = "hsapiens_gene_ensembl" interface = "default" ><Attribute name = "ensembl_gene_id" /><Attribute name = "ensembl_transcript_id" /><Attribute name = "exon_chrom_start" /><Attribute name = "exon_chrom_end" /><Attribute name = "is_constitutive" /><Attribute name = "rank" /><Attribute name = "phase" /><Attribute name = "cdna_coding_start" /><Attribute name = "cdna_coding_end" /><Attribute name = "genomic_coding_start" /><Attribute name = "genomic_coding_end" /><Attribute name = "ensembl_exon_id" /><Attribute name = "cds_start" /><Attribute name = "cds_end" /><Attribute name = "ensembl_peptide_id" /><Attribute name = "chromosome_name" /><Attribute name = "start_position" /><Attribute name = "end_position" /><Attribute name = "transcript_start" /><Attribute name = "transcript_end" /><Attribute name = "strand" /><Attribute name = "external_gene_id" /><Attribute name = "external_gene_db" /><Attribute name = "5_utr_start" /><Attribute name = "5_utr_end" /><Attribute name = "3_utr_start" /><Attribute name = "3_utr_end" /><Attribute name = "cds_length" /><Attribute name = "transcript_count" /><Attribute name = "description" /><Attribute name = "gene_biotype" /></Dataset></Query>'

ln -sf  ${GENOME_FOLDER}  ./source_data/genome
ln -sf  /home/ilya/iogen/genome/hg19_exons\(ensembl\,GRCh37.p13\).txt  ./source_data/exons.txt
ln -sf  /home/ilya/iogen/hocomoco/  ./source_data/motif_collection
ln -f  /home/ilya/iogen/hocomoco-thresholds  ./source_data/motif_thresholds

mkdir -p  $SEQ_FOLDER  $SNV_FOLDER  $CHUNK_FOLDER

ruby bin/preparations/load_cancer_mutations_sequences.rb --promoter-upstream 5000 --promoter-downstream 500 --kataegis-expansion 1000

# Markup SNVs and filter regulatory only, non-duplicated SNVs
ruby bin/preparations/snv_markup.rb  ./source_data/SNV_infos_original.txt  >  ${RESULTS_FOLDER}/SNV_infos_marked_up.txt
ruby bin/preparations/filter_snv_infos.rb  ${RESULTS_FOLDER}/SNV_infos_marked_up.txt  >  ${RESULTS_FOLDER}/SNV_infos_regulatory.txt


#  Later we work only with regulatory non-duplicated SNVs
ln -f  ${RESULTS_FOLDER}/SNV_infos_regulatory.txt  ${SNV_FOLDER}/SNV_infos_cancer.txt

##################################

# Extract sequences around SNVs (not more necessary as SNV infos and sequences are joined in a single file now)
# ruby bin/preparations/extract_snv_sequences.rb  ${SNV_FOLDER}/SNV_infos_cancer.txt  ./source_data/genome  --flank-length 50  >  ${SEQ_FOLDER}/sequences_cancer.txt

# Generate random sequences
./prepare_random_sequences.sh

# Split sequences into equal chunks in order to run chunks in parallel
NUMBER_OF_CORES=8  ./prepare_sequences_for_perfectosape_run.sh
