# Sampling random sites to make distribution of random site pvalues to be the same as an actual pvalue distribution

$:.unshift File.absolute_path('lib', __dir__)
require 'support'
require 'histogram'
require 'motifwise_histogram_fitting'

motif_names = File.readlines('source_data/motif_names.txt').map(&:strip)

mutated_site_infos_cancer_filename = ARGV[0] # 'source_data/sites_cancer_cpg.txt'
mutated_site_infos_random_filename = ARGV[1] # 'source_data/sites_random_cpg.txt'

fitters = MotifHistogramFitter.new(motif_names, Histogram.new(1e-7, 0.0005, 1.0){|pvalue| - Math.log2(pvalue) })

MutatatedSiteInfo.each_site(mutated_site_infos_cancer_filename).select(&:site_before_substitution?).each do |mutated_site_info|
  fitters.add_element(mutated_site_info.motif_name, mutated_site_info.pvalue_1)
end

$stderr.puts "Loaded #{fitters.goal_total} original sites"

MutatatedSiteInfo.each_site(mutated_site_infos_random_filename).select(&:site_before_substitution?).each do |mutated_site_info|
  fitters.fit_element(mutated_site_info.motif_name, mutated_site_info.pvalue_1) do
    puts mutated_site_info.line
  end
end

$stderr.puts "Loaded #{fitters.current_total} fitted sites (#{100.0 * fitters.current_total / fitters.goal_total}%)"
if fitters.current_total >= fitters.goal_total
  $stderr.puts "Mutations fitted"
else
  fitters.print_discrepancies(msg: 'Mutations not fitted', output_stream: $stderr)
end
