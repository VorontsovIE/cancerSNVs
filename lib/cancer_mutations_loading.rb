require 'set'
require_relative 'sample_info'
require_relative 'snv_info'

def whole_genome_samples_by_cancer(sample_infos_filename)
  samples = SampleInfo.each_in_file(sample_infos_filename).to_a

  whole_genome_samples = samples.group_by(&:cancer_type).map{|cancer_type, cancer_samples|
    [cancer_type, cancer_samples.select(&:whole_genome?).map(&:sample_name).to_set]
  }.to_h
end

# only whole-genome, non-mitochondrial
def load_cancer_mutations_by_cancer_type(somatic_mutations_folder, whole_genome_samples)
  cancer_types = Dir.glob(File.join(somatic_mutations_folder, '*')).map{|dirname| File.basename(dirname).to_sym }

  cancer_types.reject{|cancer_type|
    !whole_genome_samples[cancer_type] || whole_genome_samples[cancer_type].empty?
  }.map{|cancer_type|
    filename = File.join(somatic_mutations_folder, cancer_type.to_s, "#{cancer_type}_clean_somatic_mutations_for_signature_analysis.txt")
    filtered_mutations = MutationInfo.each_in_file(filename).reject{|mutation|
      mutation.chromosome == :MT
    }.select{|mutation|
      whole_genome_samples[cancer_type].include?(mutation.sample_id)
    }

    [cancer_type, filtered_mutations]
  }.to_h
end
