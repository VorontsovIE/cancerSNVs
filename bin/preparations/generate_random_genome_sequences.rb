$:.unshift File.absolute_path('../../lib', __dir__)
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
def calculate_context_content(chromosome, context_length:, alphabet_length:, initial_content: nil)
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
  opts.banner = 'Usage: #{program_name} <SNVs file> <ensembl genome markup> [options]'
  opts.separator 'Options:'
  opts.on('--flank-length LENGTH', 'Length of substitution sequence flanks') {|value| flank_length = value.to_i }
  opts.on('--fold FOLD', 'Multiply original context distribution FOLD times') {|value| fold = value.to_i }
  opts.on('--random-seed SEED', 'Seed for random generator') {|value| srand(Integer(value)) }
end.parse!(ARGV)

raise 'Specify SNV infos'  unless site_infos_filename = ARGV[0] # 'source_data/SNV_infos.txt'
raise 'Specify ensembl exons markup'  unless exons_filename = ARGV[1] # '/home/ilya/iogen/genome/hg19_exons(ensembl,GRCh37.p13).txt'

raise 'Specify genome folder'  unless genome_folder = ARGV[2] # '/home/ilya/iogen/genome/hg19'

# TODO:make promoter expansion configurable!
promoters_by_chromosome = load_promoters_by_chromosome(exons_filename, length_5_prime: 2000, length_3_prime: 500)
introns_by_chromosome = read_introns_by_chromosome(exons_filename)

is_promoter = ->(chr, pos) { promoters_by_chromosome[chr] && promoters_by_chromosome[chr].include_position?(pos) }
is_intron = ->(chr, pos) { introns_by_chromosome[chr] && introns_by_chromosome[chr].include_position?(pos) }

encoder = SequenceEncoder.default_encoder

# genomic_content = nil
# Dir.glob(File.join(genome_folder, '*.plain')).sort.each do |chromosome_filename|
#   $stderr.puts chromosome_filename
#   sequence_code = encoder.encode_sequence(File.read(chromosome_filename))
#   genomic_content = calculate_context_content(sequence_code, context_length: 3, alphabet_length: encoder.alphabet_length, initial_content: genomic_content)
# end
genomic_content = {0=>110995028, 1=>42125960, 8=>46557445, 43=>57819574, 91=>57133135, 80=>56717307, 27=>58762946, 12=>51523761, 61=>34547919, 58=>40536627, 41=>48861274, 83=>64119652, 42=>58800571, 88=>58497976, 66=>27417189, 85=>56740626, 53=>38691998, 18=>72131745, 40=>37300780, 77=>37346116, 13=>46615697, 67=>43640460, 86=>41753576, 92=>54944630, 68=>42283909, 90=>60157436, 78=>59611666, 93=>111360191, 81=>44797566, 30=>53461290, 26=>43541886, 6=>33722474, 28=>53159342, 15=>59562478, 16=>38646516, 55=>41724811, 76=>32820831, 33=>51568473, 5=>58303460, 10=>64039825, 50=>57071784, 2=>57725267, 87=>53554990, 63=>33755618, 31=>38200529, 52=>48868606, 62=>38247266, 60=>44823835, 3=>72039236, 75=>60075904, 25=>54742690, 17=>53161506, 56=>34556602, 57=>6940766, 35=>6411783, 11=>40523652, 32=>8046785, 38=>7306090, 51=>27383861, 65=>32833817, 7=>7289162, 37=>8048543, 82=>6423366, 36=>6933663, 4=>33, 24=>46, 124=>239849850, 123=>33, 115=>12, 29=>13, 122=>283, 110=>183, 84=>177, 49=>285, 120=>47, 100=>25, 14=>14, 74=>41, 89=>12, 44=>14, 99=>47, 94=>32, 64=>19, 121=>33, 108=>18, 107=>5, 34=>102, 112=>98, 79=>6, 118=>5, 39=>5, 102=>14, 69=>10, 98=>6, 117=>14, 106=>13, 19=>7, 95=>5, 71=>4, 105=>10, 116=>11, 103=>11, 101=>10, 45=>5, 9=>7, 59=>6, 72=>3, 20=>2, 113=>5, 21=>4, 73=>2, 119=>2, 22=>1, 70=>1, 111=>3, 54=>3, 96=>5, 23=>2, 97=>2, 114=>1, 109=>1, 47=>1, 46=>1, 48=>1}
genomic_content = genomic_content.map{|k,v| [encoder.decode_sequence(k, code_length: 3), v] }.to_h
genomic_content = genomic_content.reject{|k,v| k.match(/N/i) }

##################
snv_context_content = Hash.new(0)
snv_context_content_mut = Hash.new{|hsh, mutation_to| hsh[mutation_to] = Hash.new(0) }

SNVInfo.each_in_file(site_infos_filename).map{|snv|
  [snv.context_before, snv.mutant_base]
}.each{|context, mut|
  snv_context_content[context] += 1
  snv_context_content_mut[context][mut] += 1
}

# snv_context_content = {"AGG"=>2615, "AGC"=>1782, "ACG"=>1718, "GTC"=>687, "AAT"=>2015, "CCG"=>1247, "CGT"=>1813, "GGC"=>1480, "AAA"=>2013, "CGG"=>1248, "CCT"=>2629, "GCT"=>1804, "TGC"=>1858, "GAT"=>1052, "ACT"=>2362, "AAG"=>1642, "GGG"=>1814, "TCA"=>21019, "TGT"=>2930, "ATC"=>963, "ACA"=>2960, "ATG"=>1510, "TGG"=>2677, "CTC"=>1271, "CTT"=>1740, "ATA"=>1699, "TCC"=>6312, "GGT"=>1859, "TTT"=>2002, "CAT"=>1500, "GGA"=>6274, "TAG"=>1003, "AGA"=>17869, "GAG"=>1264, "GCA"=>1894, "CCA"=>2559, "CAA"=>930, "TAC"=>919, "GCC"=>1416, "GCG"=>1180, "AAC"=>1207, "GTT"=>1229, "CCC"=>1757, "CGC"=>1243, "TAT"=>1745, "TTA"=>1544, "AGT"=>2392, "TCG"=>2236, "CGA"=>2281, "ACC"=>1839, "TCT"=>17944, "GTA"=>908, "CAG"=>1305, "TGA"=>21125, "ATT"=>1988, "GTG"=>991, "CAC"=>1013, "CTA"=>970, "TAA"=>1487, "TTG"=>1010, "CTG"=>1316, "GAC"=>724, "TTC"=>1042, "GAA"=>1091}


snv_positions = Hash.new {|hash, key| hash[key] = [] }
SNVInfo.each_in_file(site_infos_filename).each{|snv|
  snv_positions[snv.chromosome] << snv.position
}
snv_positions = snv_positions.map{|k,v| ["chr#{k}".to_sym, v.to_set] }.to_h


rates = {}
snv_context_content.each_key {|context|
  rates[context] = genomic_content[context] / snv_context_content[context]
}
contexts = rates.keys

necessary_seqs = {}
snv_context_content.each do |context, goal|
  necessary_seqs[context] = fold * goal
end

necessary_seqs_mut = Hash.new { |hash, key| hash[key] = {} }
snv_context_content_mut.each_key do |context|
  snv_context_content_mut[context].each do |mut,goal|
    necessary_seqs_mut[context][mut] = fold * goal
  end
end


sequence_hashes = Set.new

miss = 0
Dir.glob(File.join(genome_folder, '*.plain')).sort.select{|chromosome_filename|
  chr_name = File.basename(chromosome_filename, File.extname(chromosome_filename)).to_sym
  promoters_by_chromosome.has_key?(chr_name) || introns_by_chromosome.has_key?(chr_name)  
}.each do |chromosome_filename|
  sequence = File.read(chromosome_filename).upcase
  chr_name = File.basename(chromosome_filename, File.extname(chromosome_filename)).to_sym
  contexts.each do |context|
    rate = rates[context]
    step = (rate / (fold * 1.95)).to_i
    # $stderr.puts "context: #{context}; step: #{step}"
    start_pos = flank_length + rand(step) # chromosome start (with padding)
    end_pos = sequence.length - flank_length  # chromosome end (with padding)
    (start_pos...end_pos).random_step(1, 2*step - 1) ### .select{ rand <= 1.0/rate  } instead of using step is too slow
      .select{|pos| sequence[pos-1, 3] == context }
      .select{|pos| is_intron.call(chr_name, pos) || is_promoter.call(chr_name, pos) }
      .map{|pos| [pos, sequence[pos - flank_length, 2*flank_length + 1]] }
      .reject{|pos,seq| seq.match(/N/) }
      .reject{|pos,seq| snv_positions[chr_name] && snv_positions[chr_name].include?(pos) }
      .each do |pos,seq|

      mut = necessary_seqs_mut[context].select{|k,v| v > 0 }.keys.sample
      if !mut
        miss += 1
        next
      end

      synthetic_snv_name = "#{chr_name}:#{pos + 1}/#{mut}"
      seq_w_snv = SequenceWithSNV.new(seq[0, flank_length], [seq[flank_length], mut], seq[flank_length + 1, flank_length])

      next  if sequence_hashes.include?(seq_w_snv.in_pyrimidine_context.hash) # possibly duplicate
      sequence_hashes << seq_w_snv.in_pyrimidine_context.hash

      puts [synthetic_snv_name, seq_w_snv].join("\t")
      necessary_seqs[context] -= 1
      necessary_seqs_mut[context][mut] -= 1
    end
  end
end

$stderr.puts "\n================================\n"
$stderr.puts "Missed (skipped due to overfill): #{miss}"
$stderr.puts "\n--------------------------------\n"
$stderr.puts necessary_seqs_mut.select{|k,hsh| hsh.any?{|k2,v| v > 0} }
$stderr.puts "\n================================\n"
if necessary_seqs_mut.select{|k,hsh| hsh.any?{|k2,v| v > 0} }.empty?
  $stderr.puts "OK"
else
  $stderr.puts "Not enough sequences"
  exit 1
end
