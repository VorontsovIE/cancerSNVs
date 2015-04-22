require_relative 'region_type'

SNVInfo = Struct.new(:variant_id, :snv_sequence,
                     :sample_id, :chromosome, :position, :strand,
                     :mutation_region_types) do

  def self.from_string(str)
    variant_id, snv_sequence, \
      sample_id,  chromosome, position, strand, \
      mutation_region_types = str.chomp.split("\t", 7)
    new(variant_id,
        SequenceWithSNV.from_string(snv_sequence),
        sample_id,  chromosome.to_sym, position.to_i, strand.to_sym,
        RegionType.from_string(mutation_region_types))
  end

  def self.each_in_stream(stream, &block)
    stream.each_line.lazy.map{|line| from_string(line) }.each(&block)
  end

  def self.each_in_file(filename, &block)
    return enum_for(:each_in_file, filename).lazy  unless block_given?
    File.open(filename) do |f|
      f.readline # skip header
      each_in_stream(f, &block)
    end
  end

  def to_s
    SNVInfo::COLUMN_ORDER.map{|column| send(column) }.join("\t")
  end

  def context_before
    snv_sequence.context(before: 1, after: 1, allele_variant_number: 0)
  end

  def context_after
    snv_sequence.context(before: 1, after: 1, allele_variant_number: 1)
  end

  def reference_base
    snv_sequence.allele_variants[0]
  end

  def mutant_base
    snv_sequence.allele_variants[1]
  end

  def context_full
    snv_sequence.subsequence(before: 1, after: 1).to_s
  end

  def self.reverse_strand(original_strand)
    case original_strand
    when :+
      :-
    when :-
      :+
    else
      raise 'Unknown strand'
    end
  end

  def revcomp
    SNVInfo.new(variant_id, snv_sequence.revcomp,
                sample_id, chromosome, position, SNVInfo.reverse_strand(strand),
                mutation_region_types)
  end

  def in_pyrimidine_context?
    snv_sequence.in_pyrimidine_context?
  end

  def in_pyrimidine_context
    in_pyrimidine_context? ? self : revcomp
  end


  # Deprecated
  def load_sequence(genome_folder, five_prime_flank_length, three_prime_flank_length)
    File.open( File.join(genome_folder, "chr#{chromosome}.plain") ) do |f|
      f.seek(position - five_prime_flank_length - 1)
      f.read(five_prime_flank_length + three_prime_flank_length + 1).upcase
    end
  end

  # Deprecated
  def load_site_sequence(genome_folder, site, flank_length = 0)
    load_sequence(genome_folder,
                  flank_length + site.seq_1_five_flank_length,
                  flank_length + site.seq_1_three_flank_length)
  end
end

SNVInfo::COLUMN_ORDER = [:variant_id, :snv_sequence,
                         :sample_id, :chromosome, :position,
                         :strand, :mutation_region_types]

SNVInfo::COLUMN_TITLES = {
  variant_id: 'Variant id', sample_id: 'Sample',
  chromosome: 'Chromosome', position: 'Position',
  strand: 'Strand',
  snv_sequence: 'SNV sequence',
  mutation_region_types: 'Mutation region type'
}

SNVInfo::HEADER = SNVInfo::COLUMN_ORDER.map{|column| SNVInfo::COLUMN_TITLES[column] }.join("\t")
