require 'interval_notation'
require_relative '../../experiment_configuration'
require_relative '../../lib/snv_info'
require_relative '../../lib/sequence'

def calculate_SNV_context_distribution(snv_stream)
  snv_context_distribution = Hash.new{|hsh, context| hsh[context] = Hash.new(0) }
  snv_stream.each{|snv|
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

def main_chromosomes(genome_reader)
  genome_reader.chromosome_names.select{|chromosome|
    chromosome.match(/^(\d+|X|Y)$/i)
  }.sort
end

def regulatory_fasta(sequence:, markup:)
  markup.intervals.map{|basic_interval|
    sequence[basic_interval.integer_points]
  }.join('N'.b) # sentinel
end

# Reduce sequence in such a way that regions between N-s can contain a word of size `window_size`.
# Consequent N-s are also removed
def remove_small_windows!(sequence, window_size:)
  sequence.gsub!(/N[ACGT]{1,#{window_size - 1}}N/ni, 'N'.b)
  sequence.gsub!(/^[ACGT]{1,#{window_size - 1}}N/ni, 'N'.b)
  sequence.gsub!(/N[ACGT]{1,#{window_size - 1}}$/ni, 'N'.b)
  sequence.gsub!(/N{2,}/ni, 'N'.b) # We don't need polyN-sequences longer than window size
  sequence
end

# sentinel delimits intervals, markup is a hash: {chr_name => interval_set}
def regulatory_fasta_all_chromosomes(genome_reader:, markup_by_chromosome:)
  main_chromosomes(genome_reader).map{|chr_name|
    regulatory_fasta(sequence: genome_reader.read_chromosome(chr_name),
                    markup: markup_by_chromosome[chr_name])
  }.join('N'.b) # sentinel
end

def snv_neighborhood_by_chromosome(snv_stream, flank_length:)
  positions_by_chromosome = Hash.new{|h,k| h[k] = [] }
  if flank_length > 0
    snv_stream.each{|snv|
      positions_by_chromosome[snv.chromosome] << IntervalNotation::Syntax::Long.closed_closed(snv.position - flank_length, snv.position + flank_length)
    }
  elsif flank_length == 0
    snv_stream.each{|snv|
      positions_by_chromosome[snv.chromosome] << IntervalNotation::Syntax::Long.point(snv.position)
    }
  else
    raise ArgumentError, 'Negative flank length'
  end

  result = positions_by_chromosome.map{|chr, mutation_points|
    [chr, IntervalNotation::Operations.union(mutation_points)]
  }.to_h
  Hash.new(IntervalNotation::Syntax::Long::Empty).merge(result)
end

# 3nt-context such that (2*flank_length + 1)-window has no N-s
def calculate_context_distribution(sequence, flank_length:)
  result = Hash.new(0)
  sequence.split('N').each do |chunk|
    (flank_length ... (chunk.length - flank_length)).each do |ind|
      context = chunk[ind-1, 3]
      result[context] += 1
    end
  end
  result
end

def choose_random_element(element_counts, random_generator:)
  partial_sums = element_counts.each_with_object([]) do |(elem, count), arr|
    arr << (arr.empty? ? [elem, count] : [elem, arr.last.last + count])
  end

  r = random_generator.rand(1 .. partial_sums.last.last)
  partial_sums.bsearch{|elem, cumulative_count| r <= cumulative_count }.first
end

# zero-based position
def snv_around_position(chromosome_sequence:, chromosome_name:, position:, mutation_to:, flank_length:)
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
              'Regulatory' # don't try to etimate real type (we mixed all regulatory sequences together)
  )
end


PYRIMIDINES = ['C', 'T'].map(&:b)
def pyrimidine_context_hash(sequence, flank_length:)
  (PYRIMIDINES.include?(sequence[flank_length]) ? sequence : Sequence.revcomp(sequence)).hash
end

##################

fold = 1
flank_length = 50
random_generator = Random.new

OptionParser.new{|opts|
  opts.on('--seed SEED', 'Random generator seed'){|value| random_generator = Random.new(Integer(value)) }
  opts.on('--flank-length LENGTH', 'Flank length'){|value| flank_length = Integer(value) }
  opts.on('--fold TIMES', 'Random set size is that times more than a cancer one'){|value| fold = Integer(value) }
}.parse!(ARGV)

cancer_snvs_filename = ARGV[0]

# Markup to cut regulatory FASTA (with actual mutations cutted out)
genome_markup = GENOME_MARKUP_LOADER.load_markup
ensembl_markup_by_chromosome = genome_markup.regulatory_by_chromosome
cancer_mutations = snv_neighborhood_by_chromosome(SNVInfo.each_in_file(cancer_snvs_filename), flank_length: flank_length)

markup_by_chromosome = Hash.new(IntervalNotation::Syntax::Long::Empty)
main_chromosomes(GENOME_READER).each{|chr|
  markup_by_chromosome[chr] = ensembl_markup_by_chromosome[chr] - cancer_mutations[chr]
}

# Cut FASTA
sequence = regulatory_fasta_all_chromosomes(
  genome_reader: GENOME_READER,
  markup_by_chromosome: markup_by_chromosome
)
sequence.upcase!
remove_small_windows!(sequence, window_size: 2*flank_length + 1)

# Calculate 1bp-context (trinucleotides around substitution) distributions
# We don't expect N-s anywhere in sequence, such sequences should be removed upstream for cancer SNVs.
genomic_context_distribution = calculate_context_distribution(sequence, flank_length: flank_length)

snv_context_distribution = calculate_SNV_context_distribution(SNVInfo.each_in_file(cancer_snvs_filename))
snv_context_goal = multiply_context_distribution(snv_context_distribution, fold)
snv_context_goal_mutations_unaware = snv_context_goal.map{|k, vals| [k, vals.values.inject(0,:+)] }.to_h


sequence_hashes = Set.new
singleN = 'N'.b
window_length = 2*flank_length + 1
polyN = ('N' * window_length).b # This should be ascii
range = flank_length ... (sequence.length - flank_length) # We don't change sequence length (polyN-chunk have the same size as replaced chunk)
contexts = snv_context_goal_mutations_unaware.keys \
  .sort_by{|context|
    # first take context which cover most part of genome contexts of the same type
    snv_context_goal_mutations_unaware[context].to_f / genomic_context_distribution[context]
  }.reverse
  .map(&:b) # use ASCII

contexts.each do |context|

  num_substitutions = snv_context_goal_mutations_unaware[context]

  num_substitutions.times do
    pos = random_generator.rand(range)
    seq = sequence[pos - flank_length, window_length]
    hash = pyrimidine_context_hash(seq, flank_length: flank_length)

    while seq[flank_length - 1, 3] != context || seq.index(singleN) || sequence_hashes.include?(hash)
      pos = random_generator.rand(range)
      seq = sequence[pos - flank_length, window_length]
      hash = pyrimidine_context_hash(seq, flank_length: flank_length)
    end

    sequence_hashes << hash

    mutation_to = choose_random_element(snv_context_goal[context], random_generator: random_generator)
    snv_context_goal[context][mutation_to] -= 1

    snv_info = snv_around_position(
      chromosome_sequence: sequence,
      chromosome_name: 'Glued',
      position: pos,
      flank_length: flank_length,
      mutation_to: mutation_to
    )
    puts snv_info
    sequence[(pos - flank_length), window_length] = polyN
  end
end
