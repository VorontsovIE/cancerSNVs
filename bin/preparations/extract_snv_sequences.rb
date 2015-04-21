$:.unshift File.absolute_path('../../lib', __dir__)
require 'data_import/breast_cancer_snv'
require 'optparse'

flank_length = 25

OptionParser.new do |opts|
  opts.banner = 'Usage: #{program_name} <SNVs file> <genome folder> [options]'
  opts.separator 'Options:'
  opts.on('--flank-length LENGTH') {|value| flank_length = value.to_i }
end.parse!(ARGV)

raise 'Specify a file with substitutions infos'  unless substitutions_filename = ARGV[0] # './source_data/SNV_infos.txt'
raise 'Specify genome folder'  unless genome_folder = ARGV[1] # '/home/ilya/iogen/genome/hg19/'

BreastCancerSNV.each_in_file(substitutions_filename).each do |snv|
  sequence = snv.snp_sequence_from_genome(genome_folder, flank_length, flank_length, name: snv.variant_id)
  puts sequence
end
