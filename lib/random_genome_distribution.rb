require_relative '../experiment_configuration'
require_relative 'snv_info'
require_relative 'sequence_encoding'
require_relative 'load_genome_structure'
require 'set'

class Range
  # go through range by steps of size from `from` to `to`
  def random_step(from, to, random_generator: Random::DEFAULT)
    return enum_for(:random_step, from, to, random_generator: random_generator)  unless block_given?
    pos = self.begin
    if exclude_end?
      while pos < self.end
        yield pos
        pos += random_generator.rand(from..to)
      end
    else
      while pos <= self.end
        yield pos
        pos += random_generator.rand(from..to)
      end
    end
  end
end

class RandomGenomeGenerator
  attr_reader :necessary_context_distribution, :flank_length, :sequence_hashes, :miss, :is_known_snv
  def initialize(necessary_context_distribution:, flank_length:, is_known_snv:)
    @necessary_context_distribution = necessary_context_distribution
    @flank_length = flank_length
    @sequence_hashes = Set.new
    @miss = 0
    @is_known_snv = is_known_snv
  end

  def random_regulatory_positions(sequence, position_generator:, chromosome_name:)
    unless block_given?
      return enum_for(:random_regulatory_positions, sequence,
                                                    position_generator: position_generator,
                                                    chromosome_name: chromosome_name)
    end

    position_generator
      .select{|pos| GENOME_MARKUP.regulatory?(chromosome_name, pos) }
      .reject{|pos| is_known_snv.(chromosome_name, pos) }
      .each{|pos|
        yield pos
      }
  end

  # choose random mutation
  def choose_mutation(context)
    necessary_context_distribution[context].select{|k,v| v > 0 }.keys.sample
  end

  # zero-based position
  def snv_around_position(chromosome_sequence:, chromosome_name:, position:, mutation_to:)
    sequence = chromosome_sequence[position - flank_length, 2*flank_length + 1]
    position_one_based = position + 1
    synthetic_snv_name = "#{chromosome_name}:#{position_one_based}"
    seq_w_snv = SequenceWithSNV.new(chromosome_sequence[position - flank_length, flank_length],
                                    [chromosome_sequence[position], mutation_to],
                                    chromosome_sequence[position + 1, flank_length])

    context = seq_w_snv.subsequence(before: 1, after: 1)
    SNVInfo.new("#{synthetic_snv_name}@#{context}", seq_w_snv,
                '', '', # Random genome (but it's too large and expands file to enormous size)
                chromosome_name, position_one_based, :+,
                GENOME_MARKUP.get_region_type(chromosome_name, position_one_based)
    )
  end

  def yield_uniq_mutations(sequence, position_generator:, chromosome_name:, context:)
    total_goal = necessary_context_distribution[context].each_value.inject(0, &:+)
    random_regulatory_positions(sequence, position_generator: position_generator, chromosome_name: chromosome_name) do |pos|
      next  unless sequence[pos-1, 3] == context
      break  if total_goal.zero?

      mutation_to = choose_mutation(context)
      if !mutation_to
        @miss += 1
        next
      end

      snv_info = snv_around_position(
        chromosome_sequence: sequence,
        chromosome_name: chromosome_name,
        position: pos,
        mutation_to: mutation_to
      )

      next  if snv_info.snv_sequence.sequence_variant(0).match(/N/)

      hash = snv_info.snv_sequence.in_pyrimidine_context.hash
      next  if sequence_hashes.include?(hash) # possibly duplicate
      sequence_hashes << hash

      necessary_context_distribution[context][mutation_to] -= 1
      total_goal -= 1
      yield snv_info
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

def calculate_genomic_context_distribution(genome_reader, exclude_N: true, exclude_chromosome: ->(chr){ false })
  genomic_content = nil
  encoder = SequenceEncoder.default_encoder
  genome_reader.chromosome_names.reject(&exclude_chromosome).each do |chromosome|
    $stderr.puts "Start loading chromosome #{chromosome}"

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
def calculate_SNV_context_distribution(snv_stream, known_snv_positions_by_chromosome:, exclude_N: true)
  snv_context_distribution = Hash.new{|hsh, context| hsh[context] = Hash.new(0) }
  snv_stream.each{|snv|
    known_snv_positions_by_chromosome[snv.chromosome] << snv.position

    context = snv.context_before
    mutation_to = snv.mutant_base
    snv_context_distribution[context][mutation_to] += 1
  }

  snv_context_distribution = snv_context_distribution.reject{|k,v| k.match(/N/i) }  if exclude_N
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

def generate_random_genome_according_to_snvs(from_filename:, genome_reader:, genomic_content:, fold:, random_generator: Random::DEFAULT, flank_length:, output_stream: $stdout)
  known_snv_positions_by_chromosome = Hash.new {|hsh, key| hsh[key] = Set.new }
  snv_context_distribution = calculate_SNV_context_distribution(SNVInfo.each_in_file(from_filename),
                                                                known_snv_positions_by_chromosome: known_snv_positions_by_chromosome,
                                                                exclude_N: true)
  is_known_snv = ->(chr, pos) { known_snv_positions_by_chromosome[chr].include?(pos) }
  necessary_context_distribution = multiply_context_distribution(snv_context_distribution, fold)

  rates = {}
  necessary_context_distribution.each_key {|context|
    rates[context] = genomic_content[context] / necessary_context_distribution[context].each_value.inject(0, &:+)
  }
  contexts = rates.keys

  marked_up_chromosomes = genome_reader.chromosome_names.sort.select{|chromosome|
    GENOME_MARKUP.chromosome_marked_up?(chromosome)
  }.reject{|chr| chr == :MT }

  random_genome_generator = RandomGenomeGenerator.new(necessary_context_distribution: necessary_context_distribution,
                                                      flank_length: flank_length,
                                                      is_known_snv: is_known_snv)

  contexts.each do |context|
    rate = rates[context]
    step = (rate / 1.95).to_i # heuristic
    $stderr.puts("#{context}: step #{step}")

    output_stream.puts SNVInfo::HEADER
    marked_up_chromosomes.each do |chromosome|
      sequence = genome_reader.read_sequence(chromosome, ZERO_BASED_EXCLUSIVE, 0, Float::INFINITY).upcase

      start_pos = flank_length + random_generator.rand(step) # chromosome start (with padding)
      end_pos = sequence.length - flank_length  # chromosome end (with padding)
      random_positions = (start_pos...end_pos).random_step(1, 2*step - 1, random_generator: random_generator)

      random_genome_generator.yield_uniq_mutations(sequence, position_generator: random_positions, chromosome_name: chromosome, context: context) do |snv_info|
        output_stream.puts snv_info
      end
    end
  end

  necessary_context_distribution = random_genome_generator.necessary_context_distribution
  $stderr.puts "\n================================\n"
  $stderr.puts "Missed (skipped due to overfill): #{random_genome_generator.miss}"
  $stderr.puts "\n--------------------------------\n"
  $stderr.puts necessary_context_distribution.select{|k,hsh| hsh.any?{|k2,v| v > 0} }
  $stderr.puts "\n================================\n"
  if necessary_context_distribution.select{|k,hsh| hsh.any?{|k2,v| v > 0} }.empty?
    $stderr.puts "OK"
  else
    $stderr.puts "Not enough sequences"
    raise 'Not enough sequences to choose'
  end
end
