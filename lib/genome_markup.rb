require_relative 'load_genome_structure'
require_relative 'region_type'

class GenomeMarkup
  attr_reader :introns_by_chromosome
  attr_reader :promoters_by_chromosome
  attr_reader :kataegis_regions_by_chromosome

  def initialize(introns_by_chromosome:, promoters_by_chromosome:, kataegis_regions_by_chromosome:)
    @introns_by_chromosome = introns_by_chromosome
    @promoters_by_chromosome = promoters_by_chromosome
    @kataegis_regions_by_chromosome = kataegis_regions_by_chromosome
  end

  def promoter?(chromosome, position)
    promoters_by_chromosome[chromosome].include_position?(position)
  end

  def intronic?(chromosome, position)
    introns_by_chromosome[chromosome].include_position?(position)
  end

  def kataegis?(chromosome, position)
    kataegis_regions_by_chromosome[chromosome].include_position?(position)
  end

  def get_region_type(chromosome, position)
    result = RegionType.new
    result << :intronic  if genome_markup.intronic?(chromosome, position)
    result << :promoter  if genome_markup.promoter?(chromosome, position)
    result << :kataegis  if genome_markup.kataegis?(chromosome, position)
    result
  end

  def self.load_structure(exonic_markup, kataegis_coordinates,
                          promoter_length_5_prime: 5000,
                          promoter_length_3_prime: 500,
                          kataegis_expansion_length: 1000)
    introns = read_introns_by_chromosome(exonic_markup)
    promoters = load_promoters_by_chromosome(exonic_markup, length_5_prime: promoter_length_5_prime,
                                                            length_3_prime: promoter_length_3_prime)
    kataegis = load_kataegis_regions_by_chromosome(kataegis_coordinates, expansion_length: kataegis_expansion_length)
    self.new(introns_by_chromosome: introns,
             promoters_by_chromosome: promoters,
             kataegis_regions_by_chromosome: kataegis)
  end
end
