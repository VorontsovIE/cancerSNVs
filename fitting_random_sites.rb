# Sampling random sites to make distribution of random site pvalues to be the same as an actual pvalue distribution

$:.unshift File.absolute_path('lib', __dir__)
require 'histogram'
require 'motifwise_histogram_fitting'
require 'optparse'

fitting_fold = 1
OptionParser.new do |opts|
  opts.banner = "Usage: #{opts.program_name} <original sites (cancer)> <sites to be fitted (random)> [options]\n" +
                "Calculates P-value distribution for each site in cancer and takes a subset of random sites\n" +
                "so that their P-value distributions fit original distributions."
  opts.on('--fold N', 'Take N times more sites (dilate histogram of distribution N times)') {|value|
    fitting_fold = value.to_f
  }
end.parse!(ARGV)

mutated_site_infos_cancer_filename = ARGV[0] # './results/intermediate/site_subsets/cancer_cpg.txt'
mutated_site_infos_random_filename = ARGV[1] # './results/intermediate/site_subsets/random_cpg.txt'

raise 'Specify cancer and random site files'  unless mutated_site_infos_cancer_filename && mutated_site_infos_random_filename

histograms_for_motifs = Hash.new{|hsh, motif_name|
  hsh[motif_name] = Histogram.new(1e-7, 1, 1.0){|pvalue| - Math.log2(pvalue)}
}

MutatatedSiteInfo.each_site(mutated_site_infos_cancer_filename).each do |site|
  histograms_for_motifs[site.motif_name].add_element(site.pvalue_1)
end


fitters_for_motifs = histograms_for_motifs.map{|motif, histogram|
  goal_histogram = histogram.multiply(fitting_fold)
  [motif, HistogramFitting.new(goal_histogram)]
}.to_h

fitters = MotifHistogramFitter.new(fitters_for_motifs, raise_on_missing: false)

$stderr.puts "Loaded #{fitters.goal_total} original sites"

MutatatedSiteInfo.each_site(mutated_site_infos_random_filename).each_with_index do |site, index|
  $stderr.puts "#{index} sites processed"  if index % 100000 == 0
  fitters.fit_element(site.motif_name, site.pvalue_1) do
    puts site.line
  end
end
$stderr.puts "Loaded #{fitters.current_total} fitted sites"

fitters.print_discrepancies(output_stream: $stderr)
