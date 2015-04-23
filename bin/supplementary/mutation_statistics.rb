# Statistics over cancer samples from Alexandrov et al.

$:.unshift File.absolute_path('../../lib', __dir__)
require 'load_genome_structure'
require 'data_import/cancer_mutations_loading'
require 'region_type'
require 'data_import/cancer_mutations_from_alexandrov_et_al'
require 'experiment_configuration'
require 'optparse'

promoter_length_5_prime = 5000
promoter_length_3_prime = 500
kataegis_expansion_length = 1000

OptionParser.new do |opts|
  opts.banner = "Usage: #{opts.program_name} [options]"
  opts.on('--promoter-upstream LENGTH', "Promoter's length upstream of TSS") {|value|
    promoter_length_5_prime = Integer(value)
  }
  opts.on('--promoter-downstream LENGTH', "Promoter's length downstream of TSS") {|value|
    promoter_length_3_prime = Integer(value)
  }
  opts.on('--kataegis-expansion LENGTH', "Kataegis region expansion radius") {|value|
    kataegis_expansion_length = Integer(value)
  }
end.parse!(ARGV)

ALEXANDROV_ET_AL_FOLDER = './source_data/AlexandrovEtAl/'
SOMATIC_MUTATIONS_FOLDER = File.join(ALEXANDROV_ET_AL_FOLDER, 'somatic_mutation_data')
KATAEGIS_COORDINATES_FILENAME = File.join(ALEXANDROV_ET_AL_FOLDER, 'coordinates_of_kataegis.csv')
SAMPLE_INFOS_FILENAME = File.join(ALEXANDROV_ET_AL_FOLDER, 'samples_summary.txt')

introns_by_chromosome = read_introns_by_chromosome(EXONS_FILENAME, convert_chromosome_names: false)
promoters_by_chromosome = load_promoters_by_chromosome(EXONS_FILENAME,
                                                      length_5_prime: promoter_length_5_prime,
                                                      length_3_prime: promoter_length_3_prime,
                                                      convert_chromosome_names: false)
kataegis_regions_by_chromosome = load_kataegis_regions_by_chromosome(KATAEGIS_COORDINATES_FILENAME,
                                                                    expansion_length: kataegis_expansion_length,
                                                                    convert_chromosome_names: false)

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
      context = mutation.to_snv_info(GENOME_FOLDER, cancer_type: cancer_type).context_before
      ['C', 'T'].include?(context[1]) ? context : context.reverse.tr('ACGTN', 'TGCAN')
    }.group_by(&:itself).map{|context, mutations| [context, mutations.size] }.sort_by{|k,v| v }.reverse.to_h

    fw.puts context_counts
  end
end
