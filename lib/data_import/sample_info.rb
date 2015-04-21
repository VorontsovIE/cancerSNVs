# Cancer samples from Alexandrov et al.
SampleInfo = Struct.new(:cancer_type, :sample_name, :sequencing_type, :data_source) do
  def exome?
    sequencing_type == :'Exome'
  end

  def whole_genome?
    sequencing_type == :'Whole genome'
  end

  def to_s
    [cancer_type, sample_name, sequence_type, data_source].join("\t")
  end

  # Colorectum  TCGA-AA-3552-01A-01W-0831-10  Exome TCGA data portal
  def self.from_string(line)
    cancer_type, sample_name, sequence_type, data_source = line.chomp.split("\t").map(&:to_sym)
    self.new(cancer_type, sample_name, sequence_type, data_source)
  end

  def self.each_in_file(filename, &block)
    File.open(filename){|f|
      f.readline
      f.each_line.map{|line| self.from_string(line) }.each(&block)
    }
  end
end
