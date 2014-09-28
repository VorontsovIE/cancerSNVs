def count_each_element(values)
  result = {}
  values.each do |value|
    if result.has_key?(value)
        result[value] += 1
    else
        result[value] = 1
    end
  end
  result
end

module BreastCancerSNPs
  class OutputConfigurator
    attr_reader :filename_result, :motif_names
    def initialize(filename_result, motif_names)
      @filename_result, @motif_names = filename_result, motif_names
      Dir.mkdir dirname unless Dir.exist?(dirname)
    end

    def dirname
      File.dirname(filename_result)
    end

    def basename
      File.basename(filename_result)
    end

    def filename_w_prefix(prefix)
      File.join(dirname, prefix + basename)
    end

    def output_counts_for_each_motif(output_filename, motif_counts)
      File.open(output_filename, 'w') do |fw|
        motif_names.each do |motif_name|
          fw.puts "#{motif_name}\t#{motif_counts.fetch(motif_name, 0)}"
        end
      end
    end

    def motif_statistics(file_prefix, motif_counts)
      output_counts_for_each_motif(filename_w_prefix(file_prefix), motif_counts)
    end
  end
end


def combine_conditions(*conditions)
  ->(*args, &block) do
    conditions.all?{|condition| condition.call(*args, &block) }
  end
end

def mutation_in_set_checker(mutations_subset)
  ->(name_snp, motif_name, fold_change, pvalue_1, pvalue_2) { mutations_subset.include?(name_snp.split("_")[0]) }
end

def disrupted_site_checker(fold_change_cutoff)
  ->(name_snp, motif_name, fold_change, pvalue_1, pvalue_2) { fold_change <= (1.0 / fold_change_cutoff) }
end

def created_site_checker(fold_change_cutoff)
  ->(name_snp, motif_name, fold_change, pvalue_1, pvalue_2) { fold_change >= fold_change_cutoff }
end
