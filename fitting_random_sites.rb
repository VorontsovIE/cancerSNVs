# Sampling random sites to make distribution of random site pvalues to be the same as an actual pvalue distribution

$:.unshift File.absolute_path('lib', __dir__)
require 'fitting/histogram'
require 'perfectosape/results_short'
require 'snv_info'
require 'fitting/multi_histogram_fitter'
require 'optparse'

def context_by_snv_name(snv_infos_filename)
  results = {}
  SNVInfo
    .each_in_file(snv_infos_filename)
    .map{|snv|
      [snv.variant_id, snv.in_pyrimidine_context.snv_sequence.context(before: 1, after: 1, allele_variant_number: 0)]
    }.to_h
end

fitting_fold = 1
OptionParser.new do |opts|
  opts.banner = "Usage: #{opts.program_name} <original sites (cancer)> <sites to be fitted (random)> <SNV infos cancer> <SNV infos random> [options]\n" +
                "Calculates P-value distribution for each site in cancer and takes a subset of random sites\n" +
                "so that their P-value distributions fit original distributions."
  opts.on('--fold N', 'Take N times more sites (dilate histogram of distribution N times)') {|value|
    fitting_fold = value.to_i
  }
end.parse!(ARGV)

raise 'Specify file with mutated sites in cancer'  unless mutated_site_infos_cancer_filename = ARGV[0] # './results/intermediate/site_subsets/cancer_cpg.txt'
raise 'Specify file with mutated sites in control group'  unless mutated_site_infos_random_filename = ARGV[1] # './results/intermediate/site_subsets/random_cpg.txt'
raise 'Specify file with SNV infos for cancer'  unless cancer_snv_sequences_file = ARGV[2] # ./results/SNVs/SNV_infos_cancer.txt
raise 'Specify file with SNV infos for control group'  unless random_snv_sequences_file = ARGV[3] # ./results/SNVs/SNV_infos_random_genome_13.txt

cancer_snv_contexts = context_by_snv_name(cancer_snv_sequences_file)
random_snv_contexts = context_by_snv_name(random_snv_sequences_file)

histograms = MultiHistogram.new{
  Histogram.new(-2, 2, 4){|pvalue_1| pvalue_1 }
}

PerfectosAPE::ResultShort.each_in_file(mutated_site_infos_cancer_filename).each do |site|
  motif = site.motif_name
  context = cancer_snv_contexts[ site ]
  histograms.add_element([motif, context], site.pvalue_1)
end

fitters = histograms.multiply(fitting_fold).fitter(raise_on_missing: false)

$stderr.puts "Loaded #{fitters.goal_total} original sites"

PerfectosAPE::ResultShort.each_in_file(mutated_site_infos_random_filename).each_with_index do |site, index|
  motif = site.motif_name
  context = random_snv_contexts[ site ]
  fitters.fit_element([motif, context], site.pvalue_1) do
    puts site.line
  end
end
$stderr.puts "Loaded #{fitters.current_total} fitted sites"

fitters.print_discrepancies(output_stream: $stderr)
