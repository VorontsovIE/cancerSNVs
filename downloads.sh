#!/bin/bash
cd "$(dirname "$0")"

# Genome assembly and list of all exons are huge and can be downloaded into any 
# common place to be reusable across different projects. By default they are downloaded right here
# In project we will use filesystem links to them

# Ensembl genome assembly (GRCh37)
GENOME_FOLDER='./source_data/Ensembl-GRCh37.p13'
# Ensembl exons table
EXONS_FILE='./source_data/hg19_exons\(ensembl\,GRCh37.p13\).txt'

MOTIF_COLLECTION='./source_data/hocomoco'
MOTIF_COLLECTION_THRESHOLDS='./source_data/hocomoco-thresholds'

##### Source data:
# Cancer somatic SNVs from "Mutational Processes Molding the Genomes of 21 Breast Cancers. Nik-Zainal et.al."
wget ftp://ftp.sanger.ac.uk/pub/cancer/Nik-ZainalEtAl/SUBSTITUTIONS_13Apr2012_snz.txt -O ./source_data/SNV_infos_original.txt

# Cancer somatic SNVs from "Signatures of mutational processes in human cancer. Alexandrov et al."
wget --recursive --directory-prefix='./source_data/AlexandrovEtAl' ftp://ftp.sanger.ac.uk/pub/cancer/AlexandrovEtAl/
mv ./source_data/AlexandrovEtAl/ftp.sanger.ac.uk/pub/cancer/AlexandrovEtAl/* ./source_data/AlexandrovEtAl/
rm -r ./source_data/AlexandrovEtAl/ftp.sanger.ac.uk/

# Kataegis regions
wget -O ./source_data/coordinates_of_kataegis.xls  http://www.nature.com/nature/journal/v500/n7463/extref/nature12477-s2.xls
ruby ./bin/preparations/extract_kataegis_regions.rb  ./source_data/coordinates_of_kataegis.xls > ./source_data/AlexandrovEtAl/coordinates_of_kataegis.csv
rm ./source_data/coordinates_of_kataegis.xls

##### Supplementary data:
# PerfectosAPE package
wget http://opera.autosome.ru/downloads/ape.jar -O ape.jar

# Download genome and convert FASTA files into plain text 
# (i.e. no gaps, no unnecessary symbols, so that we can seek files by chromosome position)
mkdir -p ${GENOME_FOLDER}
wget --directory-prefix=${GENOME_FOLDER} ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/dna/README
wget --directory-prefix=${GENOME_FOLDER} ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.75.dna_sm.chromosome.*.fa.gz
gzip -d ${GENOME_FOLDER}/*.fa.gz
ruby bin/preparations/convert_fasta_to_plain.rb ${GENOME_FOLDER}
rm ${GENOME_FOLDER}/*.fa

# Download exonic Ensembl markup
wget -O ${EXONS_FILE} 'http://feb2014.archive.ensembl.org/biomart/martservice?query=<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query><Query  virtualSchemaName = "default" formatter = "TSV" header = "1" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" ><Dataset name = "hsapiens_gene_ensembl" interface = "default" ><Attribute name = "ensembl_gene_id" /><Attribute name = "ensembl_transcript_id" /><Attribute name = "exon_chrom_start" /><Attribute name = "exon_chrom_end" /><Attribute name = "is_constitutive" /><Attribute name = "rank" /><Attribute name = "phase" /><Attribute name = "cdna_coding_start" /><Attribute name = "cdna_coding_end" /><Attribute name = "genomic_coding_start" /><Attribute name = "genomic_coding_end" /><Attribute name = "ensembl_exon_id" /><Attribute name = "cds_start" /><Attribute name = "cds_end" /><Attribute name = "ensembl_peptide_id" /><Attribute name = "chromosome_name" /><Attribute name = "start_position" /><Attribute name = "end_position" /><Attribute name = "transcript_start" /><Attribute name = "transcript_end" /><Attribute name = "strand" /><Attribute name = "external_gene_id" /><Attribute name = "external_gene_db" /><Attribute name = "5_utr_start" /><Attribute name = "5_utr_end" /><Attribute name = "3_utr_start" /><Attribute name = "3_utr_end" /><Attribute name = "cds_length" /><Attribute name = "transcript_count" /><Attribute name = "description" /><Attribute name = "gene_biotype" /></Dataset></Query>'

# Download hocomoco motif collection and calculate thresholds
wget -O ./source_data/hocomoco_v9_pwm.tar.gz  http://opera.autosome.ru/downloads/motif_collections/hocomoco_v9_pwm.tar.gz
mkdir -p ${MOTIF_COLLECTION}
tar --directory=${MOTIF_COLLECTION} -zxf ./source_data/hocomoco_v9_pwm.tar.gz
rm ./source_data/hocomoco_v9_pwm.tar.gz
java -cp ape.jar ru.autosome.ape.PrecalculateThresholds ${MOTIF_COLLECTION} ${MOTIF_COLLECTION_THRESHOLDS} --silent

# Links to all supplementary resources:
ln -nsf  `realpath ${GENOME_FOLDER}`  ./source_data/genome
ln -sf   `realpath ${EXONS_FILE}`  ./source_data/exons.txt
ln -nsf  `realpath ${MOTIF_COLLECTION}`  ./source_data/motif_collection
ln -nsf  `realpath ${MOTIF_COLLECTION_THRESHOLDS}`  ./source_data/motif_thresholds
