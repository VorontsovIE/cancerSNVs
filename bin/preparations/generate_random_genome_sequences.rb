$:.unshift File.absolute_path('../../lib', __dir__)
require_relative '../../experiment_configuration'
require 'random_genome_distribution'
require 'snv_info'
require 'sequence_encoding'
require 'load_genome_structure'
require 'set'
require 'optparse'

####################################################
flank_length = 25
fold = 10 # how many times we should multiply original distribution
seed = nil

OptionParser.new do |opts|
  opts.banner = 'Usage: #{program_name} <SNVs file> [options]'
  opts.separator 'Options:'
  opts.on('--flank-length LENGTH', 'Length of substitution sequence flanks') {|value| flank_length = value.to_i }
  opts.on('--fold FOLD', 'Multiply original context distribution FOLD times') {|value| fold = value.to_i }
  opts.on('--random-seed SEED', 'Seed for random generator') {|value| seed = Integer(value) }
end.parse!(ARGV)

raise 'Specify SNV infos'  unless snvs_filename = ARGV[0] # 'source_data/SNV_infos.txt'

GENOME_MARKUP = GENOME_MARKUP_LOADER.load_markup

GENOME_READER_IN_MEMORY = GenomeReader::MemoryReader.load_from_disk(
  LocalPaths::Genome,
  chromosome_name_matcher: /^Homo_sapiens\.GRCh37.75\.dna_sm\.chromosome\.(?<chromosome>\w+)\.plain$/
)

genomic_content = calculate_genomic_context_distribution(GENOME_READER_IN_MEMORY,
                                                        exclude_N: true,
                                                        exclude_chromosome: ->(chr){ chr == :MT })

generate_random_genome_according_to_snvs(snvs_filename, genomic_content: genomic_content, seed: , stream: $stdout)
