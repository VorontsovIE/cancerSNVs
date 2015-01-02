$:.unshift File.absolute_path('../../lib', __dir__)
require 'optparse'
require 'fileutils'

flush_size = 1000
pvalue_cutoff = 0.0005
OptionParser.new do |opts|
  opts.banner = "Usage: #{opts.program_name} <sites> <output folder>\n" +
                "Extract each motif's site sequences (before substitution)\n" +
                "for all sites in file (which're sites before substitution)"
  opts.on('--flush-size SIZE', 'Flush motif sequences each SIZE sequences (default: #{flush_size}') {|value|
    flush_size = value.to_i
  }
  opts.on('--pvalue-cutoff PVALUE', 'P-value to treat motif as a site') {|value|
    pvalue_cutoff = value.to_f
  }
end.parse!(ARGV)

sites_filename = ARGV[0] # './results/intermediate/site_subsets/sites_cancer_any.txt'
fasta_output_folder = ARGV[1] # './results/motif_sites_fitted/'
raise "Specify file with sites' infos and output folder"  unless sites_filename && fasta_output_folder

FileUtils.mkdir_p(fasta_output_folder)  unless Dir.exist?(fasta_output_folder)
Dir.glob(File.join(fasta_output_folder, '*')).select{|fn| File.file?(fn) }.each{|fn| File.delete(fn) }

motifs = Hash.new{|h,k| h[k] = [] }


MutatatedSiteInfo.each_site(sites_filename).select{|info|
  info.site_before_substitution?(pvalue_cutoff: pvalue_cutoff)
}.each do |info|
  motifs[info.motif_name] << info.seq_1
  if motifs[info.motif_name].size == flush_size
    motif_filename = File.join(fasta_output_folder, "#{info.motif_name}.fasta")
    File.open(motif_filename, 'a') do |fw|
      fw.puts motifs[info.motif_name]
    end
    motifs[info.motif_name] = []
  end
end

motifs.each do |motif_name, seqs|
  motif_filename = File.join(fasta_output_folder, "#{motif_name}.fasta")
  File.open(motif_filename, 'a') do |fw|
    fw.puts seqs
  end
end
