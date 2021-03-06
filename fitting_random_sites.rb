# Sampling random sites to make distribution of random site pvalues to be the same as an actual pvalue distribution

$:.unshift File.absolute_path('lib', __dir__)
require 'fitting/histogram'
require 'perfectosape/results_short'
require 'snv_info'
require 'fitting/multi_histogram_fitter'
require 'optparse'

def context_by_snv_name(variant_id)
  expanded_context = variant_id.split('@').last
  SequenceWithSNV.from_string(expanded_context).sequence_variant(0)
end

fitting_fold = 1
pvalue_cutoff = Configuration::PvalueCutoff
OptionParser.new do |opts|
  opts.banner = "Usage: #{opts.program_name} <original sites (cancer)> <sites to be fitted (random)> <SNV infos cancer> <SNV infos random> [options]\n" +
                "Calculates P-value distribution for each site in cancer and takes a subset of random sites\n" +
                "so that their P-value distributions fit original distributions."
  opts.on('--fold N', 'Take N times more sites (dilate histogram of distribution N times)') {|value|
    fitting_fold = Integer(value)
  }
  opts.on('--pvalue-cutoff CUTOFF', 'P-value of an actual site') {|value|
    pvalue_cutoff = Float(value)
  }
end.parse!(ARGV)

raise 'Specify file with mutated sites in cancer'  unless mutated_site_infos_cancer_filename = ARGV[0] # './results/intermediate/site_subsets/cancer_cpg.txt'
raise 'Specify file with mutated sites in control group'  unless mutated_site_infos_random_filename = ARGV[1] # './results/intermediate/site_subsets/random_cpg.txt'

histograms = MultiHistogram.new{
  Histogram.new(-1, 1, 2)
}

PerfectosAPE::ResultShort.each_in_file(mutated_site_infos_cancer_filename).each do |site|
  motif = site.motif_name
  context = context_by_snv_name(site.variant_id)
  histograms.add_element([motif, context], 0)
end

fitters = histograms.multiply(fitting_fold).fitter(raise_on_missing: false)

$stderr.puts "# Loaded original #{fitters.goal_total / fitting_fold} sites. Now need #{fitters.goal_total} sites"

PerfectosAPE::ResultShort.each_in_file(mutated_site_infos_random_filename).each_with_index do |site, index|
  motif = site.motif_name
  context = context_by_snv_name(site.variant_id)
  if fitters.element_can_be_added?([motif, context], 0)
    fitters.fit_element([motif, context], 0)
    puts site.line
  end
end
$stderr.puts "# Loaded #{fitters.current_total} fitted sites"
fitters.print_discrepancies(output_stream: $stderr)
