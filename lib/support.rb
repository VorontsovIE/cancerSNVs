require 'fileutils'
require 'set'
require_relative 'site_info'

def complement(str)
  str.tr('acgtnACGTN'.freeze, 'tgcanTGCAN'.freeze)
end

def revcomp(str)
  complement(str.reverse)
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

def count_each_motif_mutations(all_mutations_filename, &block)
  mutated_sites = MutatatedSiteInfo.each_site(all_mutations_filename).select(&block)
  mutated_motifs = mutated_sites.map{|mutated_site| mutated_site.motif_name }
  count_each_element(mutated_motifs)
end
