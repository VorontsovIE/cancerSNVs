query_xml = <<-EOS
<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query  virtualSchemaName = "default" formatter = "TSV" header = "1" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >
  <Dataset name = "hsapiens_gene_ensembl" interface = "default" >
    <Attribute name = "ensembl_gene_id" />
    <Attribute name = "ensembl_transcript_id" />
    <Attribute name = "exon_chrom_start" />
    <Attribute name = "exon_chrom_end" />
    <Attribute name = "is_constitutive" />
    <Attribute name = "rank" />
    <Attribute name = "phase" />
    <Attribute name = "cdna_coding_start" />
    <Attribute name = "cdna_coding_end" />
    <Attribute name = "genomic_coding_start" />
    <Attribute name = "genomic_coding_end" />
    <Attribute name = "ensembl_exon_id" />
    <Attribute name = "cds_start" />
    <Attribute name = "cds_end" />
    <Attribute name = "ensembl_peptide_id" />
    <Attribute name = "chromosome_name" />
    <Attribute name = "start_position" />
    <Attribute name = "end_position" />
    <Attribute name = "transcript_start" />
    <Attribute name = "transcript_end" />
    <Attribute name = "strand" />
    <Attribute name = "external_gene_id" />
    <Attribute name = "external_gene_db" />
    <Attribute name = "5_utr_start" />
    <Attribute name = "5_utr_end" />
    <Attribute name = "3_utr_start" />
    <Attribute name = "3_utr_end" />
    <Attribute name = "cds_length" />
    <Attribute name = "transcript_count" />
    <Attribute name = "description" />
    <Attribute name = "gene_biotype" />
  </Dataset>
</Query>
EOS

namespace 'source_data' do
  namespace 'exonic_markup' do
    desc 'Download and prepare exonic markup (Ensembl, GRCh37.p13)'
    task :prepare => [LocalPaths::ExonicMarkup]
    file LocalPaths::ExonicMarkup => [SystemPaths::ExonicMarkup] do
      ln_sf SystemPaths::ExonicMarkup, LocalPaths::ExonicMarkup
    end

    file  SystemPaths::ExonicMarkup do
      sh 'wget', '-O', SystemPaths::ExonicMarkup, "http://feb2014.archive.ensembl.org/biomart/martservice?query=#{query_xml}"
    end
  end
end
