$:.unshift File.absolute_path('lib', __dir__)
require 'set'
require 'statistical_significance'
require 'mutation_features'
require 'support'
require 'histogram'

mutated_site_infos_filename = ARGV[0] # 'source_data/cancer_SNPs.txt'
raise "Specify file with mutations"  unless mutated_site_infos_filename

mut_types = File.readlines('source_data/SUBSTITUTIONS_13Apr2012_snz_promoter_markup2.txt').drop(1).map{|el| data = el.split("\t"); [data[0], data[17]] };
regulatory_mutation_names = mutation_names_by_mutation_type(mut_types){|mut_name, mut_type| intronic_mutation?(mut_type) || promoter_mutation?(mut_type) }

each_mutated_site_info(mutated_site_infos_filename) do |mutated_site_info|
  puts mutated_site_info.line  if regulatory_mutation_names.include?(mutated_site_info.normalized_snp_name) && mutated_site_info.pvalue_1 <= 0.0005 # && mutated_site_info.fold_change <= 1.0 &&
end
