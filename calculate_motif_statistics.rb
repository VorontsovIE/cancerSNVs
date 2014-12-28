$:.unshift File.absolute_path('lib', __dir__)
require 'site_info'
require 'support'
require 'optparse'

fold_change_cutoff = 5
mode = :all
OptionParser.new do |opts|
  opts.banner = "Usage: #{opts.program_name} <sites infos> [options]"
  opts.on('--all', 'Calculate motif statistics for all sites') {
    mode = :all
  }
  opts.on('--disrupted [FOLD_CHANGE]', 'Calculate motif statistics for disrupted sites (default fold change cutoff = 5)') {|value|
    fold_change_cutoff = value.to_f
    mode = :disrupted
  }
  opts.on('--emerged [FOLD_CHANGE]', 'Calculate motif statistics for emerged sites (default fold change cutoff = 5)') {|value|
    fold_change_cutoff = value.to_f
    mode = :emerged
  }
end.parse!(ARGV)

mutated_site_infos_filename = ARGV[0] # source_data/sites_cancer_cpg.txt

raise 'Specify input file with mutation infos'  unless mutated_site_infos_filename

if mode == :disrupted
  selector = ->(site) { site.site_before_substitution? && site.disrupted?(fold_change_cutoff: fold_change_cutoff) }
elsif mode == :emerged
  selector = ->(site) { site.site_after_substitution? && site.emerged?(fold_change_cutoff: fold_change_cutoff) }
elsif mode == :all
  selector = ->(site) { true }
else
  raise LogicError, 'Should never be here'
end

mutated_sites = MutatatedSiteInfo.each_site(mutated_site_infos_filename).select(&selector).map(&:motif_name)

motif_names = File.readlines('source_data/motif_names.txt').map(&:strip)

number_of_motifs = count_each_element(mutated_sites)
motif_names.each do |motif_name|
  puts "#{motif_name}\t#{number_of_motifs.fetch(motif_name, 0)}"
end
