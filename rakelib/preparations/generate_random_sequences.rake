$:.unshift File.join(LocalPaths::Root, 'lib')
require 'set'
require 'snv_info'
require 'random_genome_distribution'

def shuffle_snvs(from_filename:, output_stream:, random_generator: Random::DEFAULT, fold: 1)
  sequence_hashes = Set.new

  output_stream.puts SNVInfo::HEADER
  SNVInfo.each_in_file(from_filename) do |snv_info|
    fold.times do |suffix|
      name_of_shuffled = "#{snv_info.variant_id.split('@').first}_#{suffix}"

      seq_w_snv = snv_info.snv_sequence
      shuffled_seq = seq_w_snv.with_flanks_shuffled(random_generator: random_generator)
      shuffled_seq = seq_w_snv.with_flanks_shuffled(random_generator: random_generator)  while sequence_hashes.include?(shuffled_seq.hash) # ignore possible duplicates
      sequence_hashes << shuffled_seq.hash

      context = shuffled_seq.subsequence(before: 1, after: 1)
      shuffled_snv = SNVInfo.new("#{name_of_shuffled}@#{context}", shuffled_seq,
          '', '', # appropriate sample/cancer type names will be too long (smth like "Shuffled Lung Adeno" which is repeated each line)
          snv_info.chromosome, snv_info.position, snv_info.strand,
          snv_info.mutation_region_types)
      output_stream.puts shuffled_snv
    end
  end

end

def generate_random_shuffle_task(output_filename:, task_name:, cancer_filename:, random_generator:, fold:)
  file cancer_filename
  file output_filename => [cancer_filename] do
    File.open(output_filename, 'w') do |fw|
      shuffle_snvs(from_filename: cancer_filename,
                  output_stream: fw,
                  fold: fold,
                  random_generator: random_generator)
    end
  end
  task task_name => output_filename
end

def generate_random_genome_task(output_filename:, task_name:, cancer_filename:, random_generator:, fold:)
  file cancer_filename
  file output_filename => [cancer_filename] do
    File.open(output_filename, 'w') do |fw|
      generate_random_genome_according_to_snvs(from_filename: cancer_filename,
                                              genome_markup: GENOME_MARKUP_LOADER.load_markup,
                                              output_stream: fw,
                                              fold: fold,
                                              flank_length: 50,
                                              genome_reader: GENOME_READER,
                                              genomic_content: get_genomic_content,
                                              random_generator: random_generator)
    end
  end
  task task_name => output_filename
end

def get_genomic_content
  # return $genomic_content  if $genomic_content
  # $genomic_content = calculate_genomic_context_distribution(
  #                       GENOME_READER,
  #                       exclude_N: true,
  #                       exclude_chromosome: ->(chr){
  #                         chr_name = chr.to_s
  #                         !chr_name.match(/^(\d+|X|Y)$/i)
  #                       })
  # $stderr.puts "Genomic content loaded"
  # $stderr.puts $genomic_content
  # File.write('./genomic_content_distribution.txt', $genomic_content.to_s)
  {"AAA"=>109671348, "GAA"=>56334225, "AAT"=>71230656, "ATT"=>71328720, "TTC"=>56404541, "TCT"=>63269496, "CTA"=>36843040, "TAC"=>32424799, "ACA"=>57551081, "CAT"=>52501415, "TTA"=>59519778, "TAG"=>36890662, "AGA"=>63171516, "ATA"=>58916160, "TAA"=>59429582, "AAC"=>41595681, "ACC"=>33217755, "CCA"=>52665758, "AGC"=>39955536, "GCC"=>34011113, "CCT"=>50784312, "CTC"=>48107749, "TCA"=>55986236, "ATC"=>38142449, "CAC"=>42903625, "CAG"=>57891673, "AGG"=>50728463, "GGC"=>33995097, "GCA"=>41147623, "ACT"=>45964533, "CTT"=>57103064, "CTG"=>57934337, "TGA"=>56003090, "AAG"=>56991469, "GCT"=>39969847, "TGC"=>41183604, "CAA"=>54053983, "TAT"=>58974661, "TCG"=>6310520, "CGT"=>7186417, "GTT"=>41761912, "TTT"=>110084438, "TTG"=>54269349, "ACG"=>7168287, "GTA"=>32440719, "TCC"=>44104934, "TGG"=>52767918, "GGG"=>37585305, "GGA"=>44128179, "GAG"=>48106270, "ATG"=>52502831, "GAT"=>38183425, "GGT"=>33256365, "GTG"=>43009093, "TGT"=>57760996, "GTC"=>27016608, "AGT"=>46024558, "GAC"=>26977547, "CGA"=>6298507, "CCC"=>37553264, "GCG"=>6799994, "CGC"=>6794339, "CCG"=>7883731, "CGG"=>7883270}
end

namespace 'preparations' do
  desc 'Generate random SNVs as a control group: shuffled and from genome.'
  task generate_random_SNVs: ['preparations:generate_random_SNVs:shuffle', 'preparations:generate_random_SNVs:genome']
  namespace 'generate_random_SNVs' do
    desc 'Generate random SNVs with shuffled flanks.'
    task :shuffle => ['preparations:generate_random_SNVs:NikZainal:shuffle', 'preparations:generate_random_SNVs:Alexandrov:shuffle']

    desc 'Generate random SNVs from genome, mimic context distribution of original SNVs.'
    task :genome => ['preparations:generate_random_SNVs:NikZainal:genome', 'preparations:generate_random_SNVs:Alexandrov:genome']
  end
end

# Alexandrov
AlexandrovWholeGenomeCancers.each do |cancer_type|
  cancer_filename = File.join(LocalPaths::Secondary::SNVs, 'Alexandrov', cancer_type.to_s, 'cancer.txt')

  generate_random_genome_task(output_filename: File.join(LocalPaths::Secondary::SNVs, 'Alexandrov', cancer_type.to_s, 'random_genome.txt'),
                              task_name: 'preparations:generate_random_SNVs:Alexandrov:genome',
                              cancer_filename: cancer_filename,
                              fold: Configuration::Alexandrov::RandomGenomeFold,
                              random_generator: Random.new("#{Configuration::AlexandrovRandomGenomeSeeds}_#{cancer_type}".hash))

  generate_random_shuffle_task(output_filename: File.join(LocalPaths::Secondary::SNVs, 'Alexandrov', cancer_type.to_s, 'random_shuffle.txt'),
                              task_name: 'preparations:generate_random_SNVs:Alexandrov:shuffle',
                              cancer_filename: cancer_filename,
                              fold: Configuration::Alexandrov::RandomShuffleFold,
                              random_generator: Random.new("{Configuration::AlexandrovRandomShuffleSeeds}_#{cancer_type}".hash))
end


# Nik-Zainal
Configuration::RandomGenomeSeeds.each do |seed|
  generate_random_genome_task(output_filename: File.join(LocalPaths::Secondary::SNVs, 'NikZainal', "random_genome_#{seed}.txt"),
                              task_name: 'preparations:generate_random_SNVs:NikZainal:genome',
                              cancer_filename: File.join(LocalPaths::Secondary::SNVs, 'NikZainal', 'cancer.txt'),
                              fold: Configuration::Alexandrov::RandomGenomeFold,
                              random_generator: Random.new(seed))
end

Configuration::RandomShuffleSeeds.each do |seed|
  generate_random_shuffle_task(output_filename: File.join(LocalPaths::Secondary::SNVs, 'NikZainal', "random_shuffle_#{seed}.txt"),
                              task_name: 'preparations:generate_random_SNVs:NikZainal:shuffle',
                              cancer_filename: File.join(LocalPaths::Secondary::SNVs, 'NikZainal', 'cancer.txt'),
                              fold: Configuration::Alexandrov::RandomShuffleFold,
                              random_generator: Random.new(seed))
end
