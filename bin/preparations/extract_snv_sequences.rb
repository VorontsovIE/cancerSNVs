$:.unshift File.absolute_path('../../lib', __dir__)
require 'snv_info'
require 'optparse'

flank_length = 25

OptionParser.new do |opts|
  opts.banner = 'Usage: #{program_name} <SNVs file> <genome folder> [options]'
  opts.separator 'Options:'
  opts.on('--flank-length LENGTH') {|value| flank_length = value.to_i }
end.parse!(ARGV)

raise 'Specify a file with substitutions infos'  unless substitutions_filename = ARGV[0] # './source_data/SNV_infos.txt'
raise 'Specify genome folder'  unless genome_folder = ARGV[1] # '/home/ilya/iogen/genome/hg19/'

SNVInfo.each_in_file(substitutions_filename) do |snv|
  puts [snv.variant_id, snv.snv_sequence].join("\t")
end
