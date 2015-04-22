$:.unshift File.absolute_path('lib', __dir__)
require 'perfectosape/results_short'
require 'optparse'
require 'set'

def count_each_element(values)
  result = {}
  values.each do |value|
    if result.has_key?(value)
        result[value] += 1
    else
        result[value] = 1
    end
  end
  result
end

mode = :all # by default consider not only disrupted/created sites, but all sites
fold_change_cutoff = 5

site_before_substitution = false # by default don't consider whether site was before substitution
pvalue_before_cutoff = 0.0005

site_after_substitution = false # by default don't consider whether site is after substitution
pvalue_after_cutoff = 0.0005
substitution_should_be_in_core = false
substitution_should_be_in_flank = false

OptionParser.new do |opts|
  opts.banner = "Usage: #{opts.program_name} <sites infos> <motif names> [options]\n" +
                "Calculate motif statistics for a list of site infos affected by SNVs.\n" +
                "By default it considers all sites in a list."

  opts.on('--site-before [PVALUE_CUTOFF]', 'Consider only sites with P-value before substitution less than cutoff.',
                                          '(default cutoff = 0.0005)') {|value|
    site_before_substitution = true
    pvalue_before_cutoff = value.to_f  if value
  }

  opts.on('--site-after [PVALUE_CUTOFF]', 'Consider only sites with P-value after substitution less than cutoff.',
                                          '(default cutoff = 0.0005') {|value|
    site_after_substitution = true
    pvalue_after_cutoff = value.to_f  if value
  }

  opts.on('--disrupted [FOLD_CHANGE]', 'Consider only sites which were disrupted by an SNV',
                                      '(default fold change cutoff = 5)') {|value|
    mode = :disrupted
    fold_change_cutoff = value.to_f  if value
  }

  opts.on('--emerged [FOLD_CHANGE]', 'Consider only sites which were emerged by an SNV',
                                    '(default fold change cutoff = 5)') {|value|
    mode = :emerged
    fold_change_cutoff = value.to_f  if value
  }

  opts.on('--substitution-in-core', 'Consider only sites where substitution was in core of the site') {
    substitution_should_be_in_core = true
  }
  opts.on('--substitution-in-flank', 'Consider only sites where substitution was in flanks of the site') {
    substitution_should_be_in_flank = true
  }
end.parse!(ARGV)

raise 'Specify input file with mutation infos'  unless mutated_site_infos_filename = ARGV[0] # ./results/intermediate/site_subsets/sites_cancer_cpg.txt
raise 'Specify input file with motif names'  unless motifs_filename = ARGV[1] # './source_data/motif_names.txt'

raise '--substitution-in-core and --substitution-in-flank should not be used together' if substitution_should_be_in_core && substitution_should_be_in_flank

if site_before_substitution
  site_before_checker = ->(site) { site.site_before_substitution?(pvalue_cutoff: pvalue_before_cutoff) }
else
  site_before_checker = ->(site) { true }
end

if site_after_substitution
  site_after_checker = ->(site) { site.site_after_substitution?(pvalue_cutoff: pvalue_after_cutoff) }
else
  site_after_checker = ->(site) { true }
end

case mode
when :disrupted
  fold_change_checker = ->(site) { site.disrupted?(fold_change_cutoff: fold_change_cutoff) }
when :emerged
  fold_change_checker = ->(site) { site.emerged?(fold_change_cutoff: fold_change_cutoff) }
when :all
  fold_change_checker = ->(site) { true }
else
  raise LogicError, 'Should never be here'
end

if substitution_should_be_in_core
  core_flank_checker = ->(site){ site.substitution_in_core? }
elsif substitution_should_be_in_flank
  core_flank_checker = ->(site){ site.substitution_in_flank? }
else
  core_flank_checker = ->(site) { true }
end

mutated_sites = PerfectosAPE::ResultShort.each_in_file(mutated_site_infos_filename)
                                .select(&site_before_checker)
                                .select(&site_after_checker)
                                .select(&fold_change_checker)
                                .select(&core_flank_checker)
                                .map(&:motif_name)

motif_names = File.readlines(motifs_filename).map(&:strip).map(&:to_sym)

number_of_motifs = count_each_element(mutated_sites)
motif_names.each do |motif_name|
  puts "#{motif_name}\t#{number_of_motifs.fetch(motif_name, 0)}"
end
