$:.unshift File.absolute_path('../../lib', __dir__)
require 'support'
require 'optparse'
require 'breast_cancer_snv'
require 'site_info'
require 'support'

mode = :words
flank_length = nil
genome_folder = nil

OptionParser.new do |opts|
  opts.banner = "Usage: #{opts.program_name} <site infos file> <motif name> <-log2(pvalue) from> <-log2(pvalue) to> [options]"
  opts.on('--mode MODE', 'Specify output mode: `words` or `variant_ids`') {|value|
    mode = value.to_sym
    raise "Unknown mode `#{mode}`"  unless [:words, :variant_ids].include?(mode)
  }
  opts.on('--expand-flank LENGTH', 'Extracted words should be expanded (from a genome) by a specified length') {|value| flank_length = value.to_i }
  opts.on('--genome-folder FOLDER', 'Folder with a genome to load flanks from'){|value| genome_folder = value }
end.parse!(ARGV)


sites_filename, motif_name, from, to = ARGV.first(4)
raise "Specify file with sites, motif name, and bin boundaries (-log2 Pvalue) and (optional) mode: words or variant_ids" unless sites_filename && motif_name && from && to

from = from.to_f
to = to.to_f
range = from..to
sites = MutatatedSiteInfo.each_site(sites_filename).select(&:site_before_substitution?).select{|info|
  info.motif_name == motif_name && range.include?(-Math.log2(info.pvalue_1))
}

snvs = BreastCancerSNV.each_substitution_in_file('./source_data/SNV_infos.txt').map{|snv| [snv.variant_id, snv] }.to_h


case mode
when :words
  if flank_length
    raise 'Specify genome folder if you want to load flanks of a site'  unless genome_folder
    data = sites.map do |site|
      snv = snvs[site.variant_id]
      seq = snv.load_site_sequence(genome_folder, site, flank_length)
      site.orientation_1 == :direct ? seq : revcomp(seq)
    end
  else
    data = sites.map(&:seq_1)
  end
when :variant_ids
  data = sites.map(&:variant_id)
else
  raise ArgumentError, 'Unknown mode'
end

data.each{|x| puts x }

# puts words
