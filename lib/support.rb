require 'fileutils'

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

# Procs instead of lambdas because &->(...) being passed as a block obtains 1 array argument instead of five actual arguments
def combine_conditions(*conditions)
  Proc.new do |*args, &block|
    conditions.all?{|condition| condition.call(*args, &block) }
  end
end

def mutation_in_set_checker(mutations_subset)
  Proc.new do |line, name_snp, motif_name, fold_change, pvalue_1, pvalue_2| mutations_subset.include?(name_snp.split("_")[0]); end
end

def disrupted_site_checker(fold_change_cutoff)
  Proc.new do |line, name_snp, motif_name, fold_change, pvalue_1, pvalue_2| fold_change <= (1.0 / fold_change_cutoff); end
end

def created_site_checker(fold_change_cutoff)
  Proc.new do |line, name_snp, motif_name, fold_change, pvalue_1, pvalue_2| fold_change >= fold_change_cutoff; end
end

def is_site_checker(pvalue_threshold)
  Proc.new do |line, name_snp, motif_name, fold_change, pvalue_1, pvalue_2| pvalue_1 <= pvalue_threshold; end
end


def each_mutation_infos(all_mutations_filename)
  return enum_for(:each_mutation_infos, all_mutations_filename)  unless block_given?
  File.open(all_mutations_filename) do |f|
    f.readline #skip header
    # "27610826_3 MAZ_f1  -7  direct  cggctgaGgaggaggag -7  direct  cggctgaCgaggaggag G/C 1.1218764110455249E-4 9.602413003842941E-4  0.11683275970285215"
    mutated_sites = f.each_line do |line|
      name_snp, motif_name,
                pos_1,orientation_1,seq_1,
                pos_2,orientation_2,seq_2,
                variants,
                pvalue_1, pvalue_2, fold_change = line.split("\t")

      pos_1 = pos_1.to_i
      pos_2 = pos_2.to_i
      # orientation_1 = orientation_1.to_sym
      # orientation_2 = orientation_2.to_sym
      pvalue_1 = pvalue_1.to_f
      pvalue_2 = pvalue_2.to_f
      fold_change = fold_change.to_f

      yield line, name_snp, motif_name, fold_change, pvalue_1, pvalue_2
    end
  end
end

def count_each_motif_mutations(all_mutations_filename, &block)
  mutated_sites = each_mutation_infos(all_mutations_filename).select(&block)
  mutated_motifs = mutated_sites.map{|line, name_snp, motif_name, fold_change, pvalue_1, pvalue_2| motif_name }
  count_each_element(mutated_motifs)
end
