$:.unshift File.absolute_path('../../lib', __dir__)
require_relative '../../experiment_configuration'
require 'snv_info'
require 'sequence_encoding'
require 'load_genome_structure'
require 'set'
require 'optparse'

class Range
  def random_step(from, to)
    return enum_for(:random_step, from, to)  unless block_given?
    pos = self.begin
    delta = to-from
    if exclude_end?
      while pos < self.end
        yield pos
        pos += (from + (rand * delta).round)
      end
    else
      while pos <= self.end
        yield pos
        pos += (from + (rand * delta).round)
      end
    end
  end
end

# chromosome should be encoded in numeric code with range of 0...alphabet_length
def calculate_context_distribution(chromosome, context_length:, alphabet_length:, initial_content: nil)
  raise 'Context length should be a positive number'  unless context_length >= 1
  raise 'Sequence is shorter than context length'  unless chromosome.length >= context_length
  content = initial_content || Hash.new(0)

  number_of_contexts = alphabet_length ** context_length

  # initial window
  window_code = 0
  chromosome.take(context_length){|letter_code|
    window_code = (window_code * alphabet_length + letter_code)
  }
  content[window_code] += 1

  # slide window
  chromosome.drop(context_length).each do |letter_code|
    window_code = (window_code * alphabet_length + letter_code) % number_of_contexts
    content[window_code] += 1
  end
  content
end

flank_length = 25
fold = 10 # how many times we should multiply original distribution

OptionParser.new do |opts|
  opts.banner = 'Usage: #{program_name} <SNVs file> [options]'
  opts.separator 'Options:'
  opts.on('--flank-length LENGTH', 'Length of substitution sequence flanks') {|value| flank_length = value.to_i }
  opts.on('--fold FOLD', 'Multiply original context distribution FOLD times') {|value| fold = value.to_i }
  opts.on('--random-seed SEED', 'Seed for random generator') {|value| srand(Integer(value)) }
end.parse!(ARGV)

raise 'Specify SNV infos'  unless site_infos_filename = ARGV[0] # 'source_data/SNV_infos.txt'

GENOME_MARKUP = GENOME_MARKUP_LOADER.load_markup

def calculate_genomic_context_distribution(genome_reader, exclude_N: true, exclude_chromosome: ->(chr){ false })
  genomic_content = nil
  encoder = SequenceEncoder.default_encoder
  genome_reader.chromosome_names.reject(&exclude_chromosome).each do |chromosome|
    sequence = genome_reader.read_sequence(chromosome, ZERO_BASED_EXCLUSIVE, 0, Float::INFINITY).upcase
    sequence_code = encoder.encode_sequence(sequence)
    genomic_content = calculate_context_distribution(sequence_code,
                                                context_length: 3,
                                                alphabet_length: encoder.alphabet_length,
                                                initial_content: genomic_content)
  end
  # genomic_content = {0=>110995028, 1=>42125960, 8=>46557445, 43=>57819574, 91=>57133135, 80=>56717307, 27=>58762946, 12=>51523761, 61=>34547919, 58=>40536627, 41=>48861274, 83=>64119652, 42=>58800571, 88=>58497976, 66=>27417189, 85=>56740626, 53=>38691998, 18=>72131745, 40=>37300780, 77=>37346116, 13=>46615697, 67=>43640460, 86=>41753576, 92=>54944630, 68=>42283909, 90=>60157436, 78=>59611666, 93=>111360191, 81=>44797566, 30=>53461290, 26=>43541886, 6=>33722474, 28=>53159342, 15=>59562478, 16=>38646516, 55=>41724811, 76=>32820831, 33=>51568473, 5=>58303460, 10=>64039825, 50=>57071784, 2=>57725267, 87=>53554990, 63=>33755618, 31=>38200529, 52=>48868606, 62=>38247266, 60=>44823835, 3=>72039236, 75=>60075904, 25=>54742690, 17=>53161506, 56=>34556602, 57=>6940766, 35=>6411783, 11=>40523652, 32=>8046785, 38=>7306090, 51=>27383861, 65=>32833817, 7=>7289162, 37=>8048543, 82=>6423366, 36=>6933663, 4=>33, 24=>46, 124=>239849850, 123=>33, 115=>12, 29=>13, 122=>283, 110=>183, 84=>177, 49=>285, 120=>47, 100=>25, 14=>14, 74=>41, 89=>12, 44=>14, 99=>47, 94=>32, 64=>19, 121=>33, 108=>18, 107=>5, 34=>102, 112=>98, 79=>6, 118=>5, 39=>5, 102=>14, 69=>10, 98=>6, 117=>14, 106=>13, 19=>7, 95=>5, 71=>4, 105=>10, 116=>11, 103=>11, 101=>10, 45=>5, 9=>7, 59=>6, 72=>3, 20=>2, 113=>5, 21=>4, 73=>2, 119=>2, 22=>1, 70=>1, 111=>3, 54=>3, 96=>5, 23=>2, 97=>2, 114=>1, 109=>1, 47=>1, 46=>1, 48=>1}
  genomic_content = genomic_content.map{|k,v| [encoder.decode_sequence(k, code_length: 3), v] }.to_h
  genomic_content = genomic_content.reject{|k,v| k.match(/N/i) }  if exclude_N
  genomic_content
end

# fills `known_snv_positions_by_chromosome` and returns extended fold times SNV context distribution
def calculate_SNV_context_distribution_from_stream(snv_stream, known_snv_positions_by_chromosome:)
  snv_context_distribution = Hash.new{|hsh, context| hsh[context] = Hash.new(0) }
  snv_stream.each{|snv|
    known_snv_positions_by_chromosome[snv.chromosome] << snv.position

    context = snv.context_before
    mutation_to = snv.mutant_base
    snv_context_distribution[context][mutation_to] += 1
  }
  snv_context_distribution
end

def multiply_context_distribution(context_distribution, fold)
  extended_context_distribution = Hash.new{|hsh, context| hsh[context] = Hash.new(0) }
  context_distribution.each_key do |context|
    context_distribution[context].each do |mutation_to, goal|
      extended_context_distribution[context][mutation_to] = fold * goal
    end
  end
  extended_context_distribution
end

genomic_content = calculate_genomic_context_distribution(GENOME_READER, exclude_N: true, exclude_chromosome: ->(chr){ chr == :MT })

known_snv_positions_by_chromosome = Hash.new {|hsh, key| hsh[key] = Set.new }
snv_context_distribution = calculate_SNV_context_distribution_from_stream(
  SNVInfo.each_in_file(site_infos_filename),
  known_snv_positions_by_chromosome: known_snv_positions_by_chromosome
)
necessary_context_distribution = multiply_context_distribution(snv_context_distribution, fold)

is_known_snv = ->(chr, pos) { known_snv_positions_by_chromosome[chr].include?(pos) }

rates = {}
necessary_context_distribution.each_key {|context|
  rates[context] = genomic_content[context] / necessary_context_distribution[context].inject(0, &:+)
}
contexts = rates.keys

sequence_hashes = Set.new

miss = 0

puts SNVInfo::HEADER

marked_up_chromosomes = GENOME_READER.chromosome_names.sort.select{|chromosome|
  GENOME_MARKUP.chromosome_marked_up?(chromosome)
}.reject{|chr| chr == :MT }

marked_up_chromosomes.each do |chromosome|
  sequence = GENOME_READER.read_sequence(chromosome, ZERO_BASED_EXCLUSIVE, 0, Float::INFINITY).upcase

  contexts.each do |context|
    rate = rates[context]
    step = (rate / 1.95).to_i # heuristic
    # $stderr.puts "context: #{context}; step: #{step}"
    start_pos = flank_length + rand(step) # chromosome start (with padding)
    end_pos = sequence.length - flank_length  # chromosome end (with padding)
    (start_pos...end_pos).random_step(1, 2*step - 1) ### .select{ rand <= 1.0/rate  } shouldn't be used instead of .step: it is too slow
      .select{|pos| GENOME_MARKUP.regulatory?(chromosome, pos) }
      .reject{|pos| is_known_snv.(chromosome, pos) }
      .select{|pos| sequence[pos-1, 3] == context }
      .map{|pos| [pos, sequence[pos - flank_length, 2*flank_length + 1]] }
      .reject{|pos,seq| seq.match(/N/) }
      .each do |pos,seq|

      mut = necessary_context_distribution[context].select{|k,v| v > 0 }.keys.sample
      if !mut
        miss += 1
        next
      end

      position = pos + 1
      synthetic_snv_name = "#{chromosome}:#{position}/#{mut}"
      seq_w_snv = SequenceWithSNV.new(seq[0, flank_length], [seq[flank_length], mut], seq[flank_length + 1, flank_length])


      snv_info = SNVInfo.new(synthetic_snv_name, seq_w_snv,
                            'random genome mutations', 'random genome mutations',
                            chromosome, position, :+,
                            GENOME_MARKUP.get_region_type(chromosome, position)
              ).in_pyrimidine_context # BUG: why resulting strands are all :+ ?!?!?!

      hash = snv_info.snv_sequence.hash
      next  if sequence_hashes.include?(hash) # possibly duplicate
      sequence_hashes << hash

      puts snv_info
      necessary_context_distribution[context][mut] -= 1
    end
  end
end

$stderr.puts "\n================================\n"
$stderr.puts "Missed (skipped due to overfill): #{miss}"
$stderr.puts "\n--------------------------------\n"
$stderr.puts necessary_context_distribution.select{|k,hsh| hsh.any?{|k2,v| v > 0} }
$stderr.puts "\n================================\n"
if necessary_context_distribution.select{|k,hsh| hsh.any?{|k2,v| v > 0} }.empty?
  $stderr.puts "OK"
else
  $stderr.puts "Not enough sequences"
  exit 1
end
