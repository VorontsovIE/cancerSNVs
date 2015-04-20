$:.unshift File.absolute_path('../../lib', __dir__)
require 'breast_cancer_snv'

raise 'Specify file with sequences'  unless sequences_filename = ARGV[0] # './source_data/sequences_background_genome.txt'

puts BreastCancerSNV::FILE_HEADER
File.readlines(sequences_filename).each.lazy.map {|line|
  full_pos, seq = line.chomp.split("\t")
  chromosome, pos = full_pos.split(':')
  pos = pos.to_i
  chromosome = chromosome.sub(/\Achr/, '')

  seq_w_snv = SequenceWithSNP.from_string(seq)
  seq_w_snv_pyrimidine = seq_w_snv.in_pyrimidine_context

  short_seq_w_snv = seq_w_snv.subsequence(before: 10, after: 10)
  short_seq_w_snv_pyrimidine = seq_w_snv_pyrimidine.subsequence(before: 10, after: 10)

  BreastCancerSNV.new(full_pos,
                      'random genome mutations', chromosome, pos, 'GRCh37',

                      seq_w_snv.allele_variants[0], seq_w_snv.allele_variants[1],

                      short_seq_w_snv_pyrimidine.left,
                      seq_w_snv_pyrimidine.allele_variants[0],
                      seq_w_snv_pyrimidine.allele_variants[1],
                      short_seq_w_snv_pyrimidine.right,

                      seq_w_snv.pyrimidine_context? ? '+' : '-',
                      '', '', '', '', #:gene, :gene_id, :ccds_id, :transcript_id,

                      '', RegionType.new, #:gene_type, :mutation_region_types,

                      '', '', '', # :mRNA_mut_syntax, :cds_mut_syntax, :aa_mut_syntax,
                      '', '' # :current_conf_status, :validation_platform
                    )
}.each {|snv_info|
  puts snv_info.to_s
}
