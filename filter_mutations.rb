$:.unshift File.absolute_path('lib', __dir__)
require 'set'
require 'statistical_significance'
require 'mutation_features'
require 'support'
require 'histogram'

mutations_filename = ARGV[0] # 'source_data/cancer_SNPs.txt'
raise "Specify file with mutations"  unless mutations_filename

mut_types = File.readlines('source_data/SUBSTITUTIONS_13Apr2012_snz_promoter_markup2.txt').drop(1).map{|el| data = el.split("\t"); [data[0], data[17]] };
regulatory_mutation_names = mutation_names_by_mutation_type(mut_types){|mut_name, mut_type| intronic_mutation?(mut_type) || promoter_mutation?(mut_type) }

each_mutation_infos(mutations_filename) do |line, name_snp, motif_name, fold_change, pvalue_1, pvalue_2|
  puts line  if regulatory_mutation_names.include?(name_snp.split("_")[0]) && pvalue_1 <= 0.0005 # && fold_change <= 1.0 &&
end
