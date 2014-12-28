$:.unshift File.absolute_path('lib', __dir__)
require 'site_info'
require 'support'
require 'optparse'

mode = :all # by default consider not only disrupted/created sites, but all sites
fold_change_cutoff = 5

site_before_substitution = false # by default don't consider whether site was before substition
pvalue_before_cutoff = 0.0005

site_after_substitution = false # by default don't consider whether site is after substition
pvalue_after_cutoff = 0.0005

OptionParser.new do |opts|
  opts.banner = "Usage: #{opts.program_name} <sites infos> [options]\n" +
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
end.parse!(ARGV)

mutated_site_infos_filename = ARGV[0] # source_data/sites_cancer_cpg.txt

raise 'Specify input file with mutation infos'  unless mutated_site_infos_filename

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

mutated_sites = MutatatedSiteInfo.each_site(mutated_site_infos_filename)
                                .select(&site_before_checker)
                                .select(&site_after_checker)
                                .select(&fold_change_checker)
                                .map(&:motif_name)

motif_names = File.readlines('source_data/motif_names.txt').map(&:strip)

number_of_motifs = count_each_element(mutated_sites)
motif_names.each do |motif_name|
  puts "#{motif_name}\t#{number_of_motifs.fetch(motif_name, 0)}"
end
