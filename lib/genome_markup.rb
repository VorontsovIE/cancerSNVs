require 'interval_notation'
require_relative 'load_genome_structure'
require_relative 'region_type'

class GenomeMarkup
  attr_reader :introns_by_chromosome
  attr_reader :promoters_by_chromosome
  attr_reader :kataegis_regions_by_chromosome
  attr_reader :regulatory_by_chromosome

  # To load genome markup from source files, use GenomeMarkupLoader
  def initialize(introns_by_chromosome:, promoters_by_chromosome:, kataegis_regions_by_chromosome:)
    @introns_by_chromosome = Hash.new(IntervalNotation::Syntax::Long::Empty).merge(introns_by_chromosome)
    @promoters_by_chromosome = Hash.new(IntervalNotation::Syntax::Long::Empty).merge(promoters_by_chromosome)
    @kataegis_regions_by_chromosome = Hash.new(IntervalNotation::Syntax::Long::Empty).merge(kataegis_regions_by_chromosome)

    @regulatory_by_chromosome = Hash.new(IntervalNotation::Syntax::Long::Empty)
    (@introns_by_chromosome.keys + @promoters_by_chromosome.keys + kataegis_regions_by_chromosome.keys).uniq.each do |chr|
      @regulatory_by_chromosome[chr] = (@promoters_by_chromosome[chr] | @introns_by_chromosome[chr]) - kataegis_regions_by_chromosome[chr]
    end
  end

  def promoter?(chromosome, position)
    case position
    when Integer
      promoters_by_chromosome[chromosome].include_position?(position)
    when IntervalNotation::IntervalSet # position is an interval
      promoters_by_chromosome[chromosome].intersect?(position)
    else
      raise TypeError, 'Position should be either integer or interval set'
    end
  end

  def intronic?(chromosome, position)
    case position
    when Integer
      introns_by_chromosome[chromosome].include_position?(position)
    when IntervalNotation::IntervalSet # position is an interval
      introns_by_chromosome[chromosome].intersect?(position)
    else
      raise TypeError, 'Position should be either integer or interval set'
    end
  end

  def kataegis?(chromosome, position)
    case position
    when Integer
      kataegis_regions_by_chromosome[chromosome].include_position?(position)
    when IntervalNotation::IntervalSet # position is an interval
      kataegis_regions_by_chromosome[chromosome].intersect?(position)
    else
      raise TypeError, 'Position should be either integer or interval set'
    end
  end

  # only singular positions
  def regulatory?(chromosome, position)
    # get_region_type(chromosome, position).regulatory?

    regulatory_by_chromosome[chromosome].include_position?(position)
  end

  def chromosome_marked_up?(chromosome)
    promoters_by_chromosome.has_key?(chromosome) || \
    introns_by_chromosome.has_key?(chromosome) || \
    kataegis_regions_by_chromosome.has_key?(chromosome)
  end

  def get_region_type(chromosome, position)
    result = RegionType.new
    result << :intronic  if intronic?(chromosome, position)
    result << :promoter  if promoter?(chromosome, position)
    result << :kataegis  if kataegis?(chromosome, position)
    result
  end
end

# This class can store data necessary to load markup (but do it lazily and caches the result)
GenomeMarkupLoader = Struct.new(:exonic_markup_filename, :kataegis_coordinates_filename,
                                :promoter_length_5_prime,
                                :promoter_length_3_prime,
                                :kataegis_expansion_length) do

  # Just a constructor with human-readable params
  def self.create(exonic_markup_filename:, kataegis_coordinates_filename:,
                  promoter_length_5_prime: 5000,
                  promoter_length_3_prime: 500,
                  kataegis_expansion_length: 1000)
    self.new(exonic_markup_filename, kataegis_coordinates_filename,
            promoter_length_5_prime,
            promoter_length_3_prime,
            kataegis_expansion_length)
  end

  # It can last very long time
  def load_markup
    if !@cache
      introns = read_introns_by_chromosome(exonic_markup_filename)
      promoters = load_promoters_by_chromosome(exonic_markup_filename,
                                              length_5_prime: promoter_length_5_prime,
                                              length_3_prime: promoter_length_3_prime)
      kataegis = load_kataegis_regions_by_chromosome(kataegis_coordinates_filename,
                                                    expansion_length: kataegis_expansion_length)
      @cache = GenomeMarkup.new(introns_by_chromosome: introns,
                               promoters_by_chromosome: promoters,
                               kataegis_regions_by_chromosome: kataegis)
    end
    @cache
  end
end
