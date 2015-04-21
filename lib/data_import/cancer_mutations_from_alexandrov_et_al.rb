require 'interval_notation'
require_relative '../snv_info'
require_relative '../sequence_with_snv'

# SNV info from Alexandrov et al.
MutationInfo = Struct.new(:sample_id, :mutation_type, :chromosome, :position_start, :position_end, :before_substitution, :after_substitution, :quality) do
  # PD3952a subs  X 41206199  41206199  C T SOMATIC-FOR-SIGNATURE-ANALYSIS
  def self.from_string(str)
    sample_id, mutation_type, \
      chromosome, position_start, position_end, \
      before_substitution, after_substitution, \
      quality = str.chomp.split("\t")

    self.new( sample_id.to_sym, mutation_type.to_sym, chromosome.to_sym,
              position_start.to_i, position_end.to_i,
              before_substitution.upcase.to_sym, after_substitution.upcase.to_sym,
              quality.to_sym )
  end

  def self.each_in_stream(stream, &block)
    stream.each_line.map{|line|
      self.from_string(line)
    }.each(&block)
  end

  def self.each_in_file(filename, &block)
    File.open(filename){|f|
      self.each_in_stream(f, &block)
    }
  end

  def to_s
    [ sample_id, mutation_type,
      chromosome, position_start, position_end,
      before_substitution, after_substitution,
      quality ].join("\t")
  end

  def to_snv_info(genome_folder, variant_id: nil, mutation_region_types: RegionType.new, flank_length: 25)
    raise 'Can\'t convert non-SNV into SNVInfo'  unless snv?
    position = position_start
    seq = File.open( File.join(genome_folder, "chr#{chromosome}.plain") ){|f|
      f.seek(position - flank_length - 1)
      f.read(2*flank_length + 1).upcase
    }
    unless seq[flank_length] == before_substitution.to_s
      raise "Base before substitution `#{before_substitution}` is not consistent with reference genome `#{seq[flank_length]}`"
    end
    five_prime_flank = seq[0, flank_length]
    three_prime_flank = seq[flank_length + 1, flank_length]
    allele_variants = [before_substitution, after_substitution]
    snv_seq = SequenceWithSNV.new(five_prime_flank, allele_variants, three_prime_flank)
    SNVInfo.new(variant_id, sample_id, chromosome, position, :+, snv_seq, mutation_region_types)
  end

  def snv? # single nucleotide substitution
    mutation_type == :subs
  end

  def indel?
    mutation_type == :indel
  end

  def interval
    if position_start == position_end
      @interval ||= IntervalNotation::Syntax::Long.point(position_start)
    elsif position_start < position_end
      @interval ||= IntervalNotation::Syntax::Long.closed_closed(position_start, position_end)
    else
      @interval ||= IntervalNotation::Syntax::Long.closed_closed(position_end, position_start)
    end
  end
end
