$:.unshift File.absolute_path('../../lib', __dir__)
require 'snv_info'
require 'load_genome_structure'
require 'optparse'

OptionParser.new do |opts|
  opts.banner = 'Usage: #{program_name} <SNVs file>'
end.parse!(ARGV)

raise 'Specify file with SNV infos'  unless snvs_filename = ARGV[0] # './source_data/SNV_infos.txt'

puts SNVInfo::HEADER
SNVInfo.each_in_file(snvs_filename).select{|snv|
  snv.regulatory?
}.to_a.uniq{|snv| # remove duplicates
  snv.in_pyrimidine_context.snv_sequence
}.each{|snv|
  puts snv
}
