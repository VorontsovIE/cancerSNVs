$:.unshift File.absolute_path('lib', __dir__)
require 'perfectosape/results_short'
require 'optparse'
require 'fileutils'
require_relative 'experiment_configuration'

def numbers_of_motifs_text(counters_by_motifs, motif_names)
  motif_names.map{|motif|
    "#{motif}\t#{counters_by_motifs.fetch(motif, 0)}"
  }.join("\n")
end

def empty_motif_counter
  Hash.new{|hsh, motif| hsh[motif] = Hash.new(0) }
end

fold_change_cutoff = 5
pvalue_cutoff = 0.0005

OptionParser.new do |opts|
  opts.banner = "Usage: #{opts.program_name} <sites infos> <motif names> [options]\n" +
                "Calculate motif statistics for a list of site infos affected by SNVs."

  opts.on('--pvalue-cutoff CUTOFF', 'Consider P-value less than cutoff as an actual site (default cutoff = 0.0005).') {|value|
    pvalue_cutoff = Float(value)
  }

  opts.on('--fold-change-cutoff CUTOFF', 'Consider only substitutions with P-value ratio greater than cutoff as disruption/emergence',
                                         '(default fold change cutoff = 5)') {|value|
    fold_change_cutoff = Float(value)
  }
end.parse!(ARGV)

raise 'Specify input file with mutation infos'  unless mutated_site_infos_filename = ARGV[0] # ./results/intermediate/site_subsets/sites_cancer_cpg.txt
raise 'Specify output folder'  unless output_folder = ARGV[1] # ./results/motif_statistics/slices/cpg/cancer/

FileUtils.mkdir_p(output_folder)  unless Dir.exist?(output_folder)

num_sites_before_substitution = empty_motif_counter
num_sites_after_substitution = empty_motif_counter
num_sites_disrupted = empty_motif_counter         # Was a site and disrupted
num_sites_emerged = empty_motif_counter           # Emerged and became a site
num_substitutions_in_core = empty_motif_counter   # Emerged and became a site
num_substitutions_in_flank = empty_motif_counter  # Emerged and became a site

PerfectosAPE::ResultShort.each_in_file(mutated_site_infos_filename) do |site_info|
  motif_name = site_info.motif_name

  if site_info.site_before_substitution?(pvalue_cutoff: pvalue_before_cutoff)
    num_sites_before_substitution[motif_name] += 1
    if site_info.disrupted?(fold_change_cutoff: fold_change_cutoff)
      num_sites_disrupted[motif_name] += 1
    end
  end

  if site_info.site_after_substitution?(pvalue_cutoff: pvalue_before_cutoff)
    num_sites_after_substitution[motif_name] += 1
    if site_info.emerged?(fold_change_cutoff: fold_change_cutoff)
      num_sites_emerged[motif_name] += 1
    end
  end

  if site.substitution_in_core?
    num_substitutions_in_core[motif_name] += 1
  else
    num_substitutions_in_flank[motif_name] += 1
  end
end

motif_names = File.readlines(LocalPaths::Secondary::MotifNames).map(&:strip).map(&:to_sym)

File.write  File.join(output_folder, 'sites_before.txt'), numbers_of_motifs_text(num_sites_before_substitution, motif_names)
File.write  File.join(output_folder, 'sites_after.txt'),  numbers_of_motifs_text(num_sites_after_substitution, motif_names)

File.write  File.join(output_folder, 'sites_disrupted.txt'), numbers_of_motifs_text(num_sites_disrupted, motif_names)
File.write  File.join(output_folder, 'sites_emerged.txt'),   numbers_of_motifs_text(num_sites_emerged, motif_names)

File.write  File.join(output_folder, 'substitutions_in_core.txt'),  numbers_of_motifs_text(num_substitutions_in_core, motif_names)
File.write  File.join(output_folder, 'substitutions_in_flank.txt'), numbers_of_motifs_text(num_substitutions_in_flank, motif_names)
