$:.unshift File.absolute_path('../../lib', __dir__)
require 'data_import/breast_cancer_snv'
require 'load_genome_structure'
require 'optparse'

flank_length = 25

OptionParser.new do |opts|
  opts.banner = 'Usage: #{program_name} <SNVs file> <genome folder> [options]'
  opts.separator 'Options:'
  opts.on('--flank-length LENGTH') {|value| flank_length = value.to_i }
end.parse!(ARGV)

raise 'Specify file with SNV infos'  unless snvs_filename = ARGV[0] # './source_data/SNV_infos.txt'
raise 'Specify genome folder'  unless genome_folder = ARGV[1] # './source_data/SNV_infos.txt'

puts BreastCancerSNV::FILE_HEADER
snvs = BreastCancerSNV.each_in_file(snvs_filename).select{|snv|
  snv.promoter? || snv.intronic?
}.to_a.uniq{|snv| # remove duplicates
  snv.snp_sequence_from_genome(genome_folder, flank_length, flank_length, name: snv.variant_id).in_pyrimidine_context
}.each{|snv|
  puts snv
}
