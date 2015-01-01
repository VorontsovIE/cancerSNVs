require 'optparse'
require 'fileutils'

pvalue = 0.0005
fold_change_cutoff = 5

OptionParser.new do |opts|
  opts.banner = "Usage: #{opts.program_name} <file with sites> <output_prefix> <motif names> [options]\n" +
                "Generate motif statistics for all sites and for disrupted/emerged sites.\n" +
                "`output_prefix` will be expanded into real filenames.\n" +
                "`output_prefix` should look like ./results/motif_statistics/cpg/cancer.txt\n" +
                "so that cancer_all.txt, cancer_sites_disrupted.txt etc will be created\n" +
                "in folder ./results/motif_statistics/cpg/"

  opts.on('--pvalue PVALUE', 'P-value to be treated as site') {|value|
    pvalue = value.to_f
  }
  opts.on('--fold-change FOLD_CHANGE', 'Fold change to be treated as disrupted/created') {|value|
    fold_change_cutoff = value.to_f
  }
end.parse!(ARGV)

sites_filename = ARGV[0] # ./results/intermediate/site_subsets/cancer_cpg.txt
output_prefix = ARGV[1] # ./results/motif_statistics/cpg/cancer.txt
motifs_filename = ARGV[2] # './source_data/motif_names.txt'

raise 'Specify sites file and prefix of output files (see --help)'  unless sites_filename && output_prefix

output_dir = File.dirname(output_prefix)
FileUtils.mkdir_p(output_dir)

extname = File.extname(output_prefix)
filename_prefix = File.basename(output_prefix, extname)

script_run = "ruby calculate_motif_statistics.rb #{sites_filename} #{motifs_filename}"

system("#{script_run} > " + File.join(output_dir, "#{filename_prefix}_all#{extname}"))
system("#{script_run} --site-before #{pvalue} > " + File.join(output_dir, "#{filename_prefix}_sites_before#{extname}"))
system("#{script_run} --site-after #{pvalue} > " + File.join(output_dir, "#{filename_prefix}_sites_after#{extname}"))
system("#{script_run} --site-before #{pvalue} --disrupted #{fold_change_cutoff} > " + File.join(output_dir, "#{filename_prefix}_sites_disrupted#{extname}"))
system("#{script_run} --site-after #{pvalue} --emerged #{fold_change_cutoff} > " + File.join(output_dir, "#{filename_prefix}_sites_emerged#{extname}"))
