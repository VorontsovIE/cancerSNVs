require 'interval_notation'

# SNV info from Alexandrov et al.
MutationInfo = Struct.new(:sample_id, :mutation_type, :chromosome, :position_start, :position_end, :ref_base_plus_strand, :mut_base_plus_strand, :quality) do
  # PD3952a subs  X 41206199  41206199  C T SOMATIC-FOR-SIGNATURE-ANALYSIS
  def self.from_string(str)
    sample_id, mutation_type, \
      chromosome, position_start, position_end, \
      ref_base_plus_strand, mut_base_plus_strand, \
      quality = str.chomp.split("\t")

    self.new( sample_id.to_sym, mutation_type.to_sym, chromosome.to_sym,
                      position_start.to_i, position_end.to_i,
                      ref_base_plus_strand.upcase.to_sym, mut_base_plus_strand.upcase.to_sym,
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
      ref_base_plus_strand, mut_base_plus_strand,
      quality ].join("\t")
  end

  def to_snv_info
    raise 'Can\'t convert non-SNV into SNVInfo'  unless snv?
    SNVInfo.new(sample_id, chromosome, position_start, ref_base_plus_strand, mut_base_plus_strand)
  end

  def snv? # single nucleotide substitution
    mutation_type == :subs
  end

  def indel?
    mutation_type == :indel
  end

  # def load_sequence(genome_folder, five_prime_flank_length: 0, three_prime_flank_length: 0)
  #   left, right = [position_start, position_end].sort
  #   File.open (File.join(genome_folder, "chr#{chromosome}.plain")) do |f|
  #     f.seek(left - five_prime_flank_length - 1)
  #     f.read(five_prime_flank_length + three_prime_flank_length + right - left + 1).upcase
  #   end
  # end

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

SNVInfo = Struct.new(:sample_id, :chromosome, :position, :ref_base_plus_strand, :mut_base_plus_strand) do
  def to_s
    [ sample_id, :subs,
      chromosome, position, position,
      ref_base_plus_strand, mut_base_plus_strand,
    ].join("\t")
  end

  # def self.each_in_stream(stream, &block)
  #   MutationInfo.each_in_stream(stream).select(&:snv?).map(&:to_snv_info).each(&block)
  # end

  # def self.each_in_file(filename, &block)
  #   MutationInfo.each_in_file(filename).select(&:snv?).map(&:to_snv_info).each(&block)
  # end

  def load_sequence(genome_folder, five_prime_flank_length: 0, three_prime_flank_length: 0)
    File.open (File.join(genome_folder, "chr#{chromosome}.plain")) do |f|
      f.seek(position - five_prime_flank_length - 1)
      f.read(five_prime_flank_length + three_prime_flank_length + 1).upcase
    end
  end
end
