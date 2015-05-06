$:.unshift File.absolute_path('../../lib', __dir__)
require 'load_genome_structure'
require 'data_import/cancer_mutations_loading'
require 'experiment_configuration'
require 'fileutils'
require 'optparse'

output_folder = './results/snv_infos/'

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

introns_by_chromosome = read_introns_by_chromosome(EXONS_FILENAME)
promoters_by_chromosome = load_promoters_by_chromosome(EXONS_FILENAME,
                                                      length_5_prime: promoter_length_5_prime,
                                                      length_3_prime: promoter_length_3_prime)
kataegis_regions_by_chromosome = load_kataegis_regions_by_chromosome(KATAEGIS_COORDINATES_FILENAME,
                                                                    expansion_length: kataegis_expansion_length)

is_promoter = ->(chr, pos) { promoters_by_chromosome[chr].include_position?(pos) }
is_intronic = ->(chr, pos) { introns_by_chromosome[chr].include_position?(pos) }
is_kataegis = ->(chr, pos) { kataegis_regions_by_chromosome[chr].include_position?(pos) }


whole_genome_samples = whole_genome_samples_by_cancer(SAMPLE_INFOS_FILENAME)
mutations_by_cancer = load_cancer_mutations_by_cancer_type(SOMATIC_MUTATIONS_FOLDER, whole_genome_samples)

  
FileUtils.mkdir_p(output_folder)  unless Dir.exist?(output_folder)

mutations_by_cancer.each{|cancer_type, mutations|
  File.open(File.join(output_folder,"#{cancer_type}.txt"), 'w') do |fw|
    fw.puts SNVInfo::HEADER
    mutations
      .select(&:snv?)
      .map{|mutation|
        snv_info = mutation.to_snv_info(GENOME_READER,
          cancer_type: cancer_type,
          variant_id: "#{mutation.sample_id}_chr#{mutation.chromosome}:#{mutation.position_start}/#{mutation.after_substitution}",
          flank_length: 50)
        snv_info.mutation_region_types << :promoter  if is_promoter.call(snv_info.chromosome, snv_info.position)
        snv_info.mutation_region_types << :intronic  if is_intronic.call(snv_info.chromosome, snv_info.position)
        snv_info.mutation_region_types << :kataegis  if is_kataegis.call(snv_info.chromosome, snv_info.position)
        snv_info
      }
      .select(&:regulatory?)
      .each {|snv_info|
        fw.puts snv_info
      }
  end
}
