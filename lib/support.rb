require 'fileutils'

module BreastCancerSNPs
  class OutputConfigurator
    attr_reader :filename_result, :motif_names
    def initialize(filename_result, motif_names)
      @filename_result, @motif_names = filename_result, motif_names
      FileUtils.mkdir_p(dirname) unless Dir.exist?(dirname)
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

def count_each_motif_mutations(all_mutations_filename)
  File.open(all_mutations_filename) do |f|
    # "27610826_3 MAZ_f1  -7  direct  cggctgaGgaggaggag -7  direct  cggctgaCgaggaggag G/C 1.1218764110455249E-4 9.602413003842941E-4  0.11683275970285215"
    mutated_sites = f.each_line.lazy.drop(1).select do |line|
      name_snp, motif_name,
                pos_1,orientation_1,seq_1,
                pos_2,orientation_2,seq_2,
                variants,
                pvalue_1, pvalue_2, fold_change = line.split("\t")

      pos_1 = pos_1.to_i
      pos_2 = pos_2.to_i
      orientation_1 = orientation_1.to_sym
      orientation_2 = orientation_2.to_sym
      pvalue_1 = pvalue_1.to_f
      pvalue_2 = pvalue_2.to_f
      fold_change = fold_change.to_f

      yield name_snp, motif_name, fold_change, pvalue_1, pvalue_2
    end
    mutated_motifs = mutated_sites.map{|el| el.split("\t")[1] }
    count_each_element(mutated_motifs)
  end
end
