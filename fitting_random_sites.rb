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

histograms_affinity = MultiHistogram.new{
  Histogram.new(1.0 / 2**40, 1.0, 0.5){|pvalue_1| - Math.log2(pvalue_1) } # Take actual site with any P-value into the same bin
}

histograms_context = MultiHistogram.new{
  Histogram.new(-1, 1, 2)
}

PerfectosAPE::ResultShort.each_in_file(mutated_site_infos_cancer_filename).each do |site|
  motif = site.motif_name
  context = context_by_snv_name(site.variant_id)
  histograms_context.add_element([motif, context], 0)
  histograms_affinity.add_element([motif], site.pvalue_1)
end

fitters_context = histograms_context.multiply(fitting_fold).fitter(raise_on_missing: false)
fitters_affinity = histograms_affinity.multiply(fitting_fold).fitter(raise_on_missing: false)

$stderr.puts "# Loaded original #{fitters_context.goal_total / fitting_fold} sites (fitting by context). Now need #{fitters_context.goal_total} sites."
$stderr.puts "# Loaded original #{fitters_affinity.goal_total / fitting_fold} sites (fitting by affinity). Now need #{fitters_affinity.goal_total} sites."

PerfectosAPE::ResultShort.each_in_file(mutated_site_infos_random_filename).each_with_index do |site, index|
  motif = site.motif_name
  context = context_by_snv_name(site.variant_id)
  if fitters_context.element_can_be_added?([motif, context], 0) && fitters_affinity.element_can_be_added?([motif], site.pvalue_1)
    fitters_context.fit_element([motif, context], 0){ }
    fitters_affinity.fit_element([motif], site.pvalue_1){ }
    puts site.line
  end
end

$stderr.puts "# Loaded #{fitters_context.current_total} fitted sites (by context)."
$stderr.puts "# Loaded #{fitters.current_total} fitted sites (by affinity)"

fitters_context.print_discrepancies(output_stream: $stderr)
$stderr.puts '#######################' # after this line fitting log won't be processed as total underfitting should be the same for both fitters
fitters_affinity.print_discrepancies(output_stream: $stderr)
