require_relative 'region_type'

SNVInfo = Struct.new(:variant_id, :sample_id, :chromosome, :position,
                     :snv_sequence,
                     :mutation_region_types) do

  def self.from_string(str)
    variant_id, sample_id,  chromosome, position, \
                snv_sequence, \
                mutation_region_types = str.chomp.split("\t", 6)
    SNVInfo.new(variant_id, sample_id,  chromosome, position,
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
end

SNVInfo::COLUMN_ORDER = [:variant_id, :sample_id, :chromosome, :position,
                         :snv_sequence, :mutation_region_types]

SNVInfo::COLUMN_TITLES = {
  variant_id: 'Variant id', sample_id: 'Sample',
  chromosome: 'Chromosome', position: 'Position',
  snv_sequence: 'SNV sequence',
  mutation_region_types: 'Mutation region type'
}

SNVInfo::HEADER = SNVInfo::COLUMN_ORDER.map{|column| SNVInfo::COLUMN_TITLES[column] }.join("\t")
