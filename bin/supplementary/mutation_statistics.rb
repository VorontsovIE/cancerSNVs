# Statistics over cancer samples from Alexandrov et al.

$:.unshift File.absolute_path('../../lib', __dir__)
require 'breast_cancer_snv'
require 'snv_info'
require 'sample_info'
require 'load_genome_structure'
require 'set'
require 'interval_notation'

GENOME_FOLDER = './source_data/genome'
EXONS_FILENAME = './source_data/exons.txt'

ALEXANDROV_ET_AL_FOLDER = './source_data/AlexandrovEtAl/'
SOMATIC_MUTATIONS_FOLDER = File.join(ALEXANDROV_ET_AL_FOLDER, 'somatic_mutation_data')
KATAEGIS_COORDINATES_FILENAME = File.join(ALEXANDROV_ET_AL_FOLDER, 'coordinates_of_kataegis.csv')
SAMPLE_INFOS_FILENAME = File.join(ALEXANDROV_ET_AL_FOLDER, 'samples_summary.txt')

def whole_genome_samples_by_cancer(sample_infos_filename)
  samples = SampleInfo.each_in_file(sample_infos_filename).to_a

  whole_genome_samples = samples.group_by(&:cancer_type).map{|cancer_type, cancer_samples|
    [cancer_type, cancer_samples.select(&:whole_genome?).map(&:sample_name).to_set]
  }.to_h
end

# only whole-genome, non-mitochondrial
def load_cancer_mutations_by_cancer_type(somatic_mutations_folder, whole_genome_samples)
  cancer_types = Dir.glob(File.join(somatic_mutations_folder, '*')).map{|dirname| File.basename(dirname).to_sym }

  cancer_types.reject{|cancer_type|
    whole_genome_samples[cancer_type].empty?
  }.map{|cancer_type|
    filename = File.join(somatic_mutations_folder, cancer_type.to_s, "#{cancer_type}_clean_somatic_mutations_for_signature_analysis.txt")
    filtered_mutations = MutationInfo.each_in_file(filename).reject{|mutation|
      mutation.chromosome == :MT
    }.select{|mutation|
      whole_genome_samples[cancer_type].include?(mutation.sample_id)
    }

    [cancer_type, filtered_mutations]
  }.to_h
end

promoters_by_chromosome = load_promoters_by_chromosome(EXONS_FILENAME, length_5_prime: 5000, length_3_prime: 500, convert_chromosome_names: false)
introns_by_chromosome = read_introns_by_chromosome(EXONS_FILENAME, convert_chromosome_names: false)
kataegis_regions_by_chromosome = load_kataegis_regions_by_chromosome(KATAEGIS_COORDINATES_FILENAME, expanded_intervals: 1000)

whole_genome_samples = whole_genome_samples_by_cancer(SAMPLE_INFOS_FILENAME)
mutations_by_cancer = load_cancer_mutations_by_cancer_type(SOMATIC_MUTATIONS_FOLDER, whole_genome_samples)


File.open('./results/alexandrov_somatic_mutations_promoter_intronic_kataegis.txt', 'w') do |fw|
  mutations_by_cancer.sort.each do |cancer_type, mutations|
    fw.puts "> #{cancer_type}"
    [true, false].each do |should_be_kataegis|
      [true, false].each do |should_be_intronic|
        [true, false].each do |should_be_promoter|
          regulatory_mutations = mutations.select{|mutation|
            is_promoter = promoters_by_chromosome[mutation.chromosome].intersect?( mutation.interval )
            is_intronic = introns_by_chromosome[mutation.chromosome].intersect?( mutation.interval )
            is_kataegis = kataegis_regions_by_chromosome[mutation.chromosome].intersect?( mutation.interval )
            (is_promoter == should_be_promoter) && (is_intronic == should_be_intronic) && (is_kataegis == should_be_kataegis)
          }
          mutation_types_counts = regulatory_mutations.group_by(&:mutation_type).map{|mut_type, muts| [mut_type, muts.size] }.sort.to_h
          fw.puts ["promoter:#{should_be_promoter} intronic:#{should_be_intronic} kataegis:#{should_be_kataegis}", mutation_types_counts].join("\t")
        end
      end
    end
  end
end


File.open('./results/alexandrov_somatic_mutations_contexts.txt', 'w') do |fw|
  fw.puts '(promoter || intronic) && !kataegis'
  mutations_by_cancer.sort.each do |cancer_type, mutations|
    fw.puts "> #{cancer_type}"

    regulatory_mutations = mutations.select(&:snv?).select{|mutation|
      is_promoter = promoters_by_chromosome[mutation.chromosome].intersect?( mutation.interval )
      is_intronic = introns_by_chromosome[mutation.chromosome].intersect?( mutation.interval )
      is_kataegis = kataegis_regions_by_chromosome[mutation.chromosome].intersect?( mutation.interval )
      (is_promoter || is_intronic) && !is_kataegis
    }

    context_counts = regulatory_mutations.map{|mutation|
      context = mutation.to_snv_info.load_sequence(GENOME_FOLDER, five_prime_flank_length: 1, three_prime_flank_length: 1).upcase
      ['C', 'T'].include?(context[1]) ? context : context.reverse.tr('ACGTN', 'TGCAN')
    }.group_by(&:itself).map{|context, mutations| [context, mutations.size] }.sort_by{|k,v| v }.reverse.to_h

    fw.puts context_counts
  end
end
