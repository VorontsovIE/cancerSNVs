SNVInfo = Struct.new(:variant_id, :sample_id, :chromosome, :position,
                     :snv_sequence, :context_before, :context_after,
                     :mutation_region_types) do
  COLUMN_ORDER = [:variant_id, :sample_id, :chromosome, :position,
                  :snv_sequence_wo_name, :context_before, :context_after,
                  :mutation_region_types]
  COLUMN_TITLES = { variant_id: 'Variant id', sample_id: 'Sample',
                    chromosome: 'Chromosome', position: 'Position',
                    snv_sequence_wo_name: 'SNV sequence',
                    context_before: 'Context before substitution (1bp around SNV)',
                    context_after: 'Context after substitution (1bp around SNV)',
                    mutation_region_types: 'Mutation region type' }
  HEADER = COLUMN_ORDER.map{|column| COLUMN_TITLES[column] }

  def self.from_string(str)
    variant_id, sample_id,  chromosome, position, \
                snv_sequence, context_before, context_after, \
                mutation_region_types = str.chomp.split("\t")
    SNVInfo.new(variant_id, sample_id,  chromosome, position,
                SequenceWithSNP.from_string(snv_sequence),
                context_before, context_after,
                mutation_region_types)
  end

  def self.each_in_stream(stream)
    stream.each_line.lazy.map{|line| from_string(line) }.each(&block)
  end

  def self.each_in_file(filename, &block)
    return enum_for(:each_in_file, filename).lazy  unless block_given?
    File.open(filename) do |f|
      f.readline # skip header
      each_in_stream(f, &block)
    end
  end

  def snv_sequence_wo_name
    snv_sequence.to_s_without_name
  end

  def to_s
    COLUMN_ORDER.map{|column| send(column) }.join("\t")
  end
end
