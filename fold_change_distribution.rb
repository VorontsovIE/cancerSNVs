$:.unshift File.absolute_path('lib', __dir__)
require 'perfectosape/results_short'
require 'optparse'
require 'fileutils'
require 'fitting/histogram'
require 'fitting/multi_histogram'
require_relative 'experiment_configuration'

def obtain_fold_change_frequencies(sites_filename, only_actual_sites:, pvalue_cutoff:, only_substitutions_in_core:, histogram_creator:)
  histograms = MultiHistogram.new(&histogram_creator)

  PerfectosAPE::ResultShort.each_in_file(sites_filename) do |site_info|
    next  if only_actual_sites && !site_info.site_before_substitution?(pvalue_cutoff: pvalue_cutoff)
    next  if only_substitutions_in_core && !site_info.substitution_in_core?
    
    motif_name = site_info.motif_name
    histograms.add_element([motif_name], site_info.fold_change)
  end
  histograms.frequencies
end

only_actual_sites = false
pvalue_cutoff = 0.0005
only_substitutions_in_core = false

OptionParser.new do |opts|
  opts.banner = "Usage: #{opts.program_name} <sites infos> <motif names> [options]\n" +
                "Calculate motif statistics for a list of site infos affected by SNVs."

  opts.on('--only-sites [CUTOFF]', 'Consider only actual sites (with P-value less than cutoff) in distribution.',
                                   'When CUTOFF not specified, cutoff of 0.0005 is taken.') {|value|
    only_actual_sites = true
    pvalue_cutoff = Float(value)  if value
  }
  opts.on('--only-core', 'Consider only substitutions in core.') {
    only_substitutions_in_core = true
  }
end.parse!(ARGV)

raise 'Specify output folder'  unless output_folder = ARGV.shift
raise 'Specify input file with mutation infos'  if (site_filenames = ARGV).empty? # ./results/intermediate/site_subsets/sites_cancer_cpg.txt


histogram_creator = ->{
  Histogram.new(1.0E-2, 1.0E+2, 1.0/4){|fold_change| Math.log(fold_change) / Math.log(2) }
}

histogram_sample = histogram_creator.call
bins = histogram_sample.each_bin(ignore_flanks: false).map{|bin,count| range_formatting(bin, rate: 3) }.to_a

frequencies_by_fn = site_filenames.map{|site_filename|
  freq = obtain_fold_change_frequencies(
    site_filename,
    histogram_creator: histogram_creator,
    only_actual_sites: only_actual_sites,
    pvalue_cutoff: pvalue_cutoff,
    only_substitutions_in_core: only_substitutions_in_core
  )
  [site_filename, freq]
}.to_h

FileUtils.mkdir_p(output_folder)
motif_names = File.readlines(LocalPaths::Secondary::MotifNames).map(&:strip).map(&:to_sym)

motif_names.each do |motif_name|
  File.open(File.join(output_folder, "#{motif_name}.csv"), 'w') do |fw|
    fw.puts ['Filename', *bins].join("\t")
    frequencies_by_fn.each do |filename, frequencies|
      motif_frequencies = frequencies[motif_name]
      freq_rounded = motif_frequencies ? motif_frequencies.map{|el| el.round(5)} : Array.new(bins.size, 0)
      fw.puts [filename, *freq_rounded].join("\t")
    end
  end
end
