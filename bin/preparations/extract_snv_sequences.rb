$:.unshift File.absolute_path('../../lib', __dir__)
require 'breast_cancer_snv'
require 'optparse'

flank_length = 25

OptionParser.new do |opts|
  opts.banner = 'Usage: #{program_name} <SNVs file> <genome folder> [options]'
  opts.separator 'Options:'
  opts.on('--flank-length LENGTH') {|value| flank_length = value.to_i }
end.parse!(ARGV)

substitutions_filename = ARGV[0] # 'source_data/SUBSTITUTIONS_13Apr2012_snz_promoter_markup2.txt'
genome_folder = ARGV[1] # '/home/ilya/iogen/genome/hg19/'
raise 'Specify a file with substitutions infos and a genome folder'  unless substitutions_filename && genome_folder

BreastCancerSNV.each_substitution_in_file(substitutions_filename).each do |snv|
  sequence = snv.snp_sequence_from_genome(genome_folder, flank_length, flank_length)
  puts "#{snv.variant_id}\t#{sequence}"
end
