require 'interval_notation'
require_relative 'sequence'

class GenomeRegionReader
  attr_reader :genome_folder
  attr_reader :mode

  def initialize(genome_folder, mode:)
    raise 'Mode can be either :zero_based or :one_based'  unless [:zero_based, :one_based].include? mode
    @genome_folder = genome_folder
    @mode = mode
  end

  def chromosome_filename(chromosome)
    chromosome = chromosome.to_s.sub(/^chr/i, '')
    File.join(genome_folder, "chr#{chromosome}.plain")
  end

  def has_chromosome?(chromosome)
    File.exist?(chromosome_filename(chromosome))
  end

  # start zero-based
  def load_sequence(chromosome, start, length)
    case mode
    when :zero_based
      seek_position = start
    when :one_based
      seek_position = start - 1
    else
      raise NotImplementedError
    end
    File.open(chromosome_filename(chromosome)) do |f|
      f.seek(seek_position)
      f.read(length)
    end
  end

  def load_interval(chromosome, interval)
    if !interval.from.is_a?(Fixnum) && interval.from_finite?  &&  !interval.to.is_a?(Fixnum) && interval.to_finite?
      raise NotImplementedError, 'Interval should have integer boundaries or be infinity'
    end
    start = interval.include_from? ? interval.from : (interval.from + 1)
    length = interval.length + 1 - (interval.include_from? ? 0 : 1) - (interval.include_to? ? 0 : 1)
    if length > 0
      load_sequence(chromosome, start, length)
    elsif length == 0
      ''
    else
      raise ArgumentError, 'Negative interval length'
    end
  end

  def load_region(chromosome, region)
    region.intervals.map{|interval| load_interval(chromosome, interval) }.join
  end

  def load_region_respecting_strand(chromosome, region, strand)
    seq = load_region(chromosome, region)
    case strand
    when :+
      seq
    when :-
      Sequence.revcomp(seq)
    else
      raise ArgumentError, "Unknown strand `#{strand.inspect}`"
    end
  end
end
