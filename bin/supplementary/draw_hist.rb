$:.unshift File.absolute_path('lib', __dir__)
require 'set'
require 'statistical_significance'
require 'mutation_features'
require 'support'
require 'histogram'

motif = ARGV[0]
mutated_site_infos_filename = ARGV[1] # 'source_data/sites_cancer.txt'

raise "Specify motif name and file with mutation infos "  unless motif && mutated_site_infos_filename
# motifs = ['AP2A_f2', 'ESR1_do', 'NFKB1_f1']

create_histogram = ->{ Histogram.new(1e-7, 0.0005, 0.25){|pvalue| - Math.log2(pvalue) } }

motif_names = File.readlines('source_data/motif_names.txt').map(&:strip)

mut_types = File.readlines('source_data/SUBSTITUTIONS_13Apr2012_snz_promoter_markup2.txt').drop(1).map{|el| data = el.split("\t"); [data[0], data[17]] };
intronic_mutation_names = mutation_names_by_mutation_type(mut_types){|mut_name, mut_type| intronic_mutation?(mut_type) }
promoter_mutation_names = mutation_names_by_mutation_type(mut_types){|mut_name, mut_type| promoter_mutation?(mut_type) }
regulatory_mutation_names = mutation_names_by_mutation_type(mut_types){|mut_name, mut_type| intronic_mutation?(mut_type) || promoter_mutation?(mut_type) }

histogram = create_histogram.call

total = 0
MutatatedSiteInfo.each_site(mutated_site_infos_filename) do |mutated_site_info|
  next  unless motif == mutated_site_info.motif_name
  next  unless mutated_site_info.pvalue_1 <= 0.0005
  next  unless regulatory_mutation_names.include?(mutated_site_info.normalized_snp_name)
  histogram.add_element( mutated_site_info.pvalue_1 )
end

$stderr.puts "Loaded #{total}"

# $stderr.puts "#{motif} -- #{random_histogram[motif].elements_total}; #{histogram[motif].elements_total}"
histogram.each_bin do |bin_range, bin_count|
  fraction = bin_count.to_f / histogram.elements_total_in_range
  puts "#{bin_range.begin}\t#{bin_range.end}\t#{fraction}"
end
