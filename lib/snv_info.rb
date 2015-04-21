require_relative 'region_type'

SNVInfo = Struct.new(:variant_id, :sample_id, :chromosome, :position, :strand,
                     :snv_sequence,
                     :mutation_region_types) do

  def self.from_string(str)
    variant_id, sample_id,  chromosome, position, strand, \
                snv_sequence, \
                mutation_region_types = str.chomp.split("\t", 7)
    new(variant_id, sample_id,  chromosome, position, strand.to_sym,
        SequenceWithSNV.from_string(snv_sequence),
        RegionType.from_string(mutation_region_types))
  end

  def self.each_in_stream(stream, &block)
    stream.each_line.lazy.map{|line| from_string(line) }.each(&block)
  end

  def self.each_in_file(filename, &block)
    return enum_for(:each_in_file, filename).lazy  unless block_given?
    File.open(filename) do |f|
      # f.readline # skip header
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

  def revcomp
    case strand
    when :+
      reverse_strand = :-
    when :-
      reverse_strand = :+
    else
      raise 'Unknown strand'
    end

    SNVInfo.new(variant_id, sample_id, chromosome, position, reverse_strand,
                snv_sequence.revcomp, mutation_region_types)
  end

  def in_pyrimidine_context?
    snv_sequence.pyrimidine_context?
  end

  def in_pyrimidine_context
    in_pyrimidine_context? ? self : revcomp
  end
end

SNVInfo::COLUMN_ORDER = [:variant_id, :sample_id, :chromosome, :position, :strand,
                         :snv_sequence, :mutation_region_types]

SNVInfo::COLUMN_TITLES = {
  variant_id: 'Variant id', sample_id: 'Sample',
  chromosome: 'Chromosome', position: 'Position',
  strand: 'Strand',
  snv_sequence: 'SNV sequence',
  mutation_region_types: 'Mutation region type'
}

SNVInfo::HEADER = SNVInfo::COLUMN_ORDER.map{|column| SNVInfo::COLUMN_TITLES[column] }.join("\t")
