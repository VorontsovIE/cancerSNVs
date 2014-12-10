require 'fileutils'
require 'set'

def complement(str)
  str.tr('acgtnACGTN'.freeze, 'tgcanTGCAN'.freeze)
end

def revcomp(str)
  complement(str.reverse)
end


MutatatedSiteInfo = Struct.new( :line,
                                :variant_id, :motif_name,
                                :fold_change, :pvalue_1, :pvalue_2,
                                :pos_1, :orientation_1, :seq_1,
                                :pos_2, :orientation_2, :seq_2,
                                :variants ) do

  # "27610826_3 MAZ_f1  -7  direct  cggctgaGgaggaggag -7  direct  cggctgaCgaggaggag G/C 1.1218764110455249E-4 9.602413003842941E-4  0.11683275970285215"
  def self.from_string(line)
    variant_id, motif_name,
              pos_1, orientation_1, seq_1,
              pos_2, orientation_2, seq_2,
              variants,
              pvalue_1, pvalue_2, fold_change = line.split("\t")
    MutatatedSiteInfo.new(line,
                          variant_id, motif_name,
                          fold_change.to_f, pvalue_1.to_f, pvalue_2.to_f,
                          pos_1.to_i, orientation_1.to_sym, seq_1,
                          pos_2.to_i, orientation_2.to_sym, seq_2,
                          variants)
  end

  def normalized_snp_name
    variant_id.split("_")[0]
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

# Procs instead of lambdas because &->(...) being passed as a block obtains 1 array argument instead of five actual arguments
def combine_conditions(*conditions)
  Proc.new do |*args, &block|
    conditions.all?{|condition| condition.call(*args, &block) }
  end
end

#############
def mutation_in_set_checker(mutations_subset)
  Proc.new do |mutated_site| mutations_subset.include?(mutated_site.normalized_snp_name); end
end

def disrupted_site_checker(fold_change_cutoff)
  Proc.new do |mutated_site| mutated_site.fold_change <= (1.0 / fold_change_cutoff); end
end

def created_site_checker(fold_change_cutoff)
  Proc.new do |mutated_site| mutated_site.fold_change >= fold_change_cutoff; end
end

def was_site_checker(pvalue_threshold)
  Proc.new do |mutated_site| mutated_site.pvalue_1 <= pvalue_threshold; end
end

def becomes_site_checker(pvalue_threshold)
  Proc.new do |mutated_site| mutated_site.pvalue_2 <= pvalue_threshold; end
end


##############
def disrupted_and_in_set_checker(mutations_subset)
  combine_conditions(mutation_in_set_checker(mutations_subset), was_site_checker(0.0005), disrupted_site_checker(5))
end

def created_and_in_set_checker(mutations_subset)
  combine_conditions(mutation_in_set_checker(mutations_subset), becomes_site_checker(0.0005), created_site_checker(5))
end

def improved_and_in_set_checker(mutations_subset)
  combine_conditions(mutation_in_set_checker(mutations_subset), becomes_site_checker(0.0005), created_site_checker(1))
end

def was_site_and_in_set_checker(mutations_subset)
  combine_conditions(mutation_in_set_checker(mutations_subset), was_site_checker(0.0005))
end

def becomes_site_and_in_set_checker(mutations_subset)
  combine_conditions(mutation_in_set_checker(mutations_subset), becomes_site_checker(0.0005))
end

#####

def disrupted_motifs_in_set(mutated_site_infos_filename, mutations_subset)
  count_each_motif_mutations(mutated_site_infos_filename, &disrupted_and_in_set_checker(mutations_subset))
end

def created_motifs_in_set(mutated_site_infos_filename, mutations_subset)
  count_each_motif_mutations(mutated_site_infos_filename, &created_and_in_set_checker(mutations_subset))
end

def improved_motifs_in_set(mutated_site_infos_filename, mutations_subset)
  count_each_motif_mutations(mutated_site_infos_filename, &improved_and_in_set_checker(mutations_subset))
end

def was_site_motifs_in_set(mutated_site_infos_filename, mutations_subset)
  count_each_motif_mutations(mutated_site_infos_filename, &was_site_and_in_set_checker(mutations_subset))
end

def becomes_site_motifs_in_set(mutated_site_infos_filename, mutations_subset)
  count_each_motif_mutations(mutated_site_infos_filename, &becomes_site_and_in_set_checker(mutations_subset))
end

##############

def each_mutated_site_info_in_stream(stream, &block)
  stream.each_line.lazy.reject{|line|
    line.start_with?('#')
  }.map{|line|
    MutatatedSiteInfo.from_string(line)
  }.each(&block)
end

def each_mutated_site_info(all_mutations_filename, &block)
  return enum_for(:each_mutated_site_info, all_mutations_filename).lazy  unless block_given?
  File.open(all_mutations_filename) do |f|
    each_mutated_site_info_in_stream(f, &block)
  end
end

def count_each_motif_mutations(all_mutations_filename, &block)
  mutated_sites = each_mutated_site_info(all_mutations_filename).select(&block)
  mutated_motifs = mutated_sites.map{|mutated_site| mutated_site.motif_name }
  count_each_element(mutated_motifs)
end
