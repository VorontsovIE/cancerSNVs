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
            '', '', # appropriate sample/cancer type names will be too long (smth like "Shuffled Lung Adeno" which is repeated each line)
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
    AlexandrovWholeGenomeCancers.each do |cancer_type|
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
          # ruby 'bin/preparations/generate_random_genome_sequences.rb',
          #       t.prerequisites.first,
          #       "--fold=#{Configuration::RandomGenomeFold}",
          #       "--random-seed=#{seed}",
          #       '--flank-length=50',
          #       {out: t.name}, {}

          require 'random_genome_distribution'
          GENOME_MARKUP ||= GENOME_MARKUP_LOADER.load_markup

          # GENOME_READER_IN_MEMORY ||= GenomeReader::MemoryReader.load_from_disk(
          #   LocalPaths::Genome,
          #   chromosome_name_matcher: /^Homo_sapiens\.GRCh37.75\.dna_sm\.chromosome\.(?<chromosome>\w+)\.plain$/
          # )
          $stderr.puts "Genome markup loaded"

          # GENOMIC_CONTENT ||= calculate_genomic_context_distribution(
          #                       GENOME_READER,
          #                       exclude_N: true,
          #                       exclude_chromosome: ->(chr){
          #                         chr_name = chr.to_s
          #                         chr_name == 'MT' || chr_name.start_with?('HG') || chr_name.start_with?('HS')
          #                       })
          # $stderr.puts "Genomic content loaded"
          # $stderr.puts GENOMIC_CONTENT
          # File.write('./genomic_content_distribution.txt', GENOMIC_CONTENT.to_s)

          GENOMIC_CONTENT = {"AAA"=>109671348, "GAA"=>56334225, "AAT"=>71230656, "ATT"=>71328720, "TTC"=>56404541, "TCT"=>63269496, "CTA"=>36843040, "TAC"=>32424799, "ACA"=>57551081, "CAT"=>52501415, "TTA"=>59519778, "TAG"=>36890662, "AGA"=>63171516, "ATA"=>58916160, "TAA"=>59429582, "AAC"=>41595681, "ACC"=>33217755, "CCA"=>52665758, "AGC"=>39955536, "GCC"=>34011113, "CCT"=>50784312, "CTC"=>48107749, "TCA"=>55986236, "ATC"=>38142449, "CAC"=>42903625, "CAG"=>57891673, "AGG"=>50728463, "GGC"=>33995097, "GCA"=>41147623, "ACT"=>45964533, "CTT"=>57103064, "CTG"=>57934337, "TGA"=>56003090, "AAG"=>56991469, "GCT"=>39969847, "TGC"=>41183604, "CAA"=>54053983, "TAT"=>58974661, "TCG"=>6310520, "CGT"=>7186417, "GTT"=>41761912, "TTT"=>110084438, "TTG"=>54269349, "ACG"=>7168287, "GTA"=>32440719, "TCC"=>44104934, "TGG"=>52767918, "GGG"=>37585305, "GGA"=>44128179, "GAG"=>48106270, "ATG"=>52502831, "GAT"=>38183425, "GGT"=>33256365, "GTG"=>43009093, "TGT"=>57760996, "GTC"=>27016608, "AGT"=>46024558, "GAC"=>26977547, "CGA"=>6298507, "CCC"=>37553264, "GCG"=>6799994, "CGC"=>6794339, "CCG"=>7883731, "CGG"=>7883270}

          File.open(t.name, 'w') do |fw|
            generate_random_genome_according_to_snvs(t.prerequisites.first, genome_reader: GENOME_READER, genomic_content: GENOMIC_CONTENT, fold: Configuration::RandomGenomeFold, seed: seed, stream: fw)
          end
        end
      end
    end

  end
end
