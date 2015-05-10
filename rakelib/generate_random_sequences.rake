$:.unshift File.join(LocalPaths::Root, 'lib')
require 'set'
require 'snv_info'

def shuffle_snvs(from_filename:, to_filename:, seed: nil, fold: 1)
  srand(seed)  if seed
  sequence_hashes = Set.new
  File.open(to_filename, 'w') do |fw|
    fw.puts SNVInfo::HEADER
    SNVInfo.each_in_file(from_filename) do |snv_info|
      fold.times do |suffix|
        name_of_shuffled = "#{snv_info.variant_id}_#{suffix}"

        seq_w_snv = snv_info.snv_sequence
        shuffled_seq = seq_w_snv.with_flanks_shuffled
        shuffled_seq = seq_w_snv.with_flanks_shuffled  while sequence_hashes.include?(shuffled_seq.hash) # ignore possible duplicates
        sequence_hashes << shuffled_seq.hash

        shuffled_snv = SNVInfo.new(name_of_shuffled, shuffled_seq,
            "Shuffled #{snv_info.cancer_type}", "Shuffled #{snv_info.sample_id}",
            snv_info.chromosome, snv_info.position, snv_info.strand,
            snv_info.mutation_region_types)
        fw.puts shuffled_snv
      end
    end
  end
end

namespace 'preparations' do
  desc 'Generate random SNVs as a control group: shuffled and from genome.'
  task generate_random_SNVs: ['generate_random_SNVs:shuffle', 'generate_random_SNVs:genome']
  namespace 'generate_random_SNVs' do
    desc 'Generate random SNVs with shuffled flanks.'
    task :shuffle

    desc 'Generate random SNVs from genome, mimic context distribution of original SNVs.'
    task :genome
    AlexandrovCancerTypes.each do |cancer_type|
      shuffle_snvs_filename = File.join(LocalPaths::Secondary::SNVs, 'Alexandrov', cancer_type.to_s, 'random_shuffle.txt')
      random_genome_snvs_filename = File.join(LocalPaths::Secondary::SNVs, 'Alexandrov', cancer_type.to_s, 'random_genome.txt')
      
      task shuffle: [shuffle_snvs_filename]
      file shuffle_snvs_filename => File.join(LocalPaths::Secondary::SNVs, 'Alexandrov', cancer_type.to_s, "#{cancer_type}.txt") do |t|
        Configuration::RandomShuffleSeeds.each do |seed|
          shuffle_snvs(from_filename: t.prerequisites.first, to_filename: t.name, seed: seed, fold: Configuration::RandomShuffleFold)
        end
      end

      task genome: [random_genome_snvs_filename]
      file random_genome_snvs_filename => File.join(LocalPaths::Secondary::SNVs, 'Alexandrov', cancer_type.to_s, "#{cancer_type}.txt") do |t|
        Configuration::RandomGenomeSeeds.each do |seed|
          ruby 'bin/preparations/generate_random_genome_sequences.rb',
                t.prerequisites.first,
                "--fold=#{Configuration::RandomGenomeFold}",
                "--random-seed=#{seed}",
                '--flank-length=50',
                {out: t.name}, {}
        end
      end
    end

  end
end
