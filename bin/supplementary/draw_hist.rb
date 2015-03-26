$:.unshift File.absolute_path('../../lib', __dir__)
require 'set'
require 'histogram'
require 'breast_cancer_snv'
require 'site_info'

motif = ARGV[0].to_sym # AHR_si
mutated_site_infos_filename = ARGV[1] # './source_data/sites_cancer.txt'
snv_infos_filename = ARGV[2] #'./source_data/SNV_infos.txt'

raise "Specify motif name, file with sites infos and file with SNV infos"  unless motif && mutated_site_infos_filename && snv_infos_filename

regulatory_mutation_names = BreastCancerSNV.each_substitution_in_file(snv_infos_filename).select(&:regulatory?).map(&:variant_id).to_set

histogram = Histogram.new(1e-7, 0.0005, 1.0){|pvalue| - Math.log2(pvalue) }

MutatatedSiteInfo.each_site(mutated_site_infos_filename).select{|site_info|
  site_info.motif_name == motif &&
  site_info.site_before_substitution?(pvalue_cutoff: 0.0005) &&
  regulatory_mutation_names.include?(site_info.normalized_snp_name)
}.each{|site_info|
  histogram.add_element( site_info.pvalue_1 )
}

puts histogram.bin_counts_info(ignore_flanks: false)
