## Workflow
Workflow is following:
* Get cancer SNVs
* Remove and recreate location type markup: promoter and intronic regions
* Filter SNVs to obtain regulatory SNVs only (located either in promoter or in intronic regions). This SNV collection is used for the rest of computations. We don't use non-regulatory SNVs anymore.
* Generate sequences with substitutions for site extraction:
    * Extract sequences around cancer SNVs (regulatory, you remember) from genome.
    * Generate several copies (e.g. 10; more helps in further fitting, but linearly slows down site extraction, which is the slowest step in workflow) of shuffled sequences. Shuffling preserves 1-bp context around SNV, so distribution of mutation contexts is the same as in cancer.
    * Take random sequences from genome (not equal to cancer positions) and put random mutations into them, so that mutation types distribution is the same as in cancer. Better to have several times more mutations than in cancer (e.g. 10 times more)
* (parallelizable task, use several cores on several servers because it's a bottleneck. It can take about 1.5 days using 8 CPU cores for multiplication factor 10 as we used). Run PerfectosAPE in order to obtain sites in all sequences collected in a previous step (cancer, genomic random and shuffling random).
* If necessary, take subset of sites for a certain mutation type, such as TpCpN or NpCpG (on either strand). Because random sequences preserve mutation context distribution, this step is safe. If one need the only mutation context, now and further (not our case), he can make this filtering before sequence extraction and generation.
* Make fitting of random sites to real cancer sites. Sites around somatic mutations and around randomly created mutations can present different P-value distributions, we sample sites of random mutations so that each site's P-value distribution for each different mutation context was proportional to a corresponding in-cancer distribution.
* Calculate statistics of sites total and of disrupted sites in cancer substitutions and in choosen random substitutions. 
* Evaluate Fischer exact test for disrupted sites (across all sites) in cancer vs the same value in random mutations.
* Now we can get rate of site disruption in cancer vs rate of site disruption  in random mutations and calculate their ratio. Significance of this difference is given by a P-value from Fisher exact test (we should make multiple comparison correction such as Bonferony/Holms/FDR).
* Filter only motifs whose sites are disrupted significantly more often than in random, choose only A-,B,-C- quality motifs.
* Repeat the same computations for several random mutation sets. Intersect results in order to get reliable answer.

## Reproduce an experiment
In order to reproduce an experiment one should at first resolve dependencies and download source data. Look through `preparations.sh` and fix links to downloaded data. Normally you'd have a genome in a separate folder and have a symlink pointing to it: `./source_data/genome/`. Same is true for almost all other external data. Fix `ln` commands in `preparations.sh` according to your machine data location.

Workflow is carried out by `make_all.sh` script in a folder root. Before running it, read comments inside, it'll save you lots of CPU-time.

### Dependencies:
    * bunch of GNU tools (bash, echo, cat, split, etc...) included in any Linux distributive
    * ruby (version 2.1 or greater)
    * java (version 1.6 or higher)
    * PerfectosAPE [http://opera.autosome.ru/perfectosape/description] package. Download *.jar file here[http://opera.autosome.ru/downloads/ape.jar] and put it into a root folder.
    * ruby gem bundler
    * Some other ruby gems (run `bundle install` to install them)

### Source files:
    * ./source_data/SNV_infos_original.txt -- Breast cancer SNVs from Mutational processes molding the genomes of 21 breast cancers. Cell, 149(5):979–993. Nik-Zainal, et al. (2012).
    These data use hg19 assembly, so one need to download other data (such as genome itself or exon markup) related to hg19.
    * Human genome hg19 assembly from UCSC[http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz]. Chromosome sequences in plain text format. One should preprocess files so that sequences had no headers and line-breaks (in order to obtain sequences using file seek). Files should be named like `genome_folder/chr1.plain`
    * Ensembl exon markup[http://feb2014.archive.ensembl.org/biomart/martview/b139ef98cf27cbd5649f7c5f6d3e2c0c?VIRTUALSCHEMANAME=default&ATTRIBUTES=hsapiens_gene_ensembl.default.structure.ensembl_gene_id|hsapiens_gene_ensembl.default.structure.ensembl_transcript_id|hsapiens_gene_ensembl.default.structure.exon_chrom_start|hsapiens_gene_ensembl.default.structure.exon_chrom_end|hsapiens_gene_ensembl.default.structure.is_constitutive|hsapiens_gene_ensembl.default.structure.rank|hsapiens_gene_ensembl.default.structure.phase|hsapiens_gene_ensembl.default.structure.cdna_coding_start|hsapiens_gene_ensembl.default.structure.cdna_coding_end|hsapiens_gene_ensembl.default.structure.genomic_coding_start|hsapiens_gene_ensembl.default.structure.genomic_coding_end|hsapiens_gene_ensembl.default.structure.ensembl_exon_id|hsapiens_gene_ensembl.default.structure.cds_start|hsapiens_gene_ensembl.default.structure.cds_end|hsapiens_gene_ensembl.default.structure.ensembl_peptide_id|hsapiens_gene_ensembl.default.structure.chromosome_name|hsapiens_gene_ensembl.default.structure.start_position|hsapiens_gene_ensembl.default.structure.end_position|hsapiens_gene_ensembl.default.structure.transcript_start|hsapiens_gene_ensembl.default.structure.transcript_end|hsapiens_gene_ensembl.default.structure.strand|hsapiens_gene_ensembl.default.structure.external_gene_id|hsapiens_gene_ensembl.default.structure.external_gene_db|hsapiens_gene_ensembl.default.structure.5_utr_start|hsapiens_gene_ensembl.default.structure.5_utr_end|hsapiens_gene_ensembl.default.structure.3_utr_start|hsapiens_gene_ensembl.default.structure.3_utr_end|hsapiens_gene_ensembl.default.structure.cds_length|hsapiens_gene_ensembl.default.structure.transcript_count|hsapiens_gene_ensembl.default.structure.description|hsapiens_gene_ensembl.default.structure.gene_biotype&FILTERS=].
    * Fantom 5 TSS-peaks (deprecated dependency, not used and soon will be removed).
