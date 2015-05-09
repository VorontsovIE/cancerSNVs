require 'set'
require_relative 'sample_info'
require_relative 'cancer_mutations_from_alexandrov_et_al'
require_relative '../snv_info'

def whole_genome_samples_by_cancer(sample_infos_filename)
  samples = SampleInfo.each_in_file(sample_infos_filename).to_a

  whole_genome_samples = samples.group_by(&:cancer_type).map{|cancer_type, cancer_samples|
    [cancer_type, cancer_samples.select(&:whole_genome?).map(&:sample_name).to_set]
  }.to_h
  whole_genome_samples.default = Set.new
  whole_genome_samples
end

# non-mitochondrial, from specified samples only
def mutations_filered_by_sample(filename, samples)
  MutationInfo.each_in_file(filename).reject{|mutation|
    mutation.chromosome == :MT
  }.select{|mutation|
    samples.include?(mutation.sample_id)
  }
end

# only whole-genome, non-mitochondrial
def load_cancer_mutations_by_cancer_type(somatic_mutations_folder, sample_infos_filename)
  whole_genome_samples_by_cancer(sample_infos_filename)
  .reject{|cancer_type, samples|
    samples.empty?
  }.map{|cancer_type, samples|
    filename = File.join(somatic_mutations_folder, 
                        cancer_type.to_s,
                        "#{cancer_type}_clean_somatic_mutations_for_signature_analysis.txt")
    [cancer_type, mutations_filered_by_sample(filename, samples)]
  }.sort_by{|cancer_type, mutations|
    cancer_type
  }.to_h
end

MutationsByCancerTypeLoader = Struct.new(:mutations_folder, :samples_summary_filename) do
  def self.create(mutations_folder:, samples_summary_filename:)
    self.new(mutations_folder, samples_summary_filename)
  end

  def load
    if !@cache
      @cache = load_cancer_mutations_by_cancer_type(mutations_folder, samples_summary_filename)
    end
    @cache
  end
end
