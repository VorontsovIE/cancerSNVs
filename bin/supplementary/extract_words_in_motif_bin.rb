# Caution! Script can work only with full perfectosape format, not with a short one
$:.unshift File.absolute_path('../../lib', __dir__)
require 'optparse'
require 'snv_info'
require 'perfectosape/results'
require_relative '../../experiment_configuration'

mode = :words
flank_length = nil

OptionParser.new do |opts|
  opts.banner = "Usage: #{opts.program_name} <site infos file> <motif> <-log2(pvalue) from> <-log2(pvalue) to> <SNV infos> [options]\n" +
                "Extract word sequences from `motif` sites which have a log2(P-value) between `from` and `to`."
  opts.on('--mode MODE', 'Specify output mode: `words` or `variant_ids`') {|value|
    mode = value.to_sym
    raise "Unknown mode `#{mode}`"  unless [:words, :variant_ids].include?(mode)
  }
  opts.on('--expand-flank LENGTH', 'Extracted words should be expanded (from a genome) by a specified length') {|value| flank_length = value.to_i }
end.parse!(ARGV)


sites_filename = ARGV[0] # sites_cancer.txt
motif_name = ARGV[1].to_sym # AHR_si
from = ARGV[2] # 15
to = ARGV[3] # 16
snv_infos_filename = ARGV[4] # './source_data/SNV_infos.txt'
raise "Specify file with sites, motif name, bin boundaries (-log2 Pvalue) and SNV infos file" unless sites_filename && motif_name && from && to && snv_infos_filename

from = from.to_f
to = to.to_f
range = from..to
sites = PerfectosAPE::Result.each_in_file(sites_filename).select(&:site_before_substitution?).select{|info|
  info.motif_name == motif_name && range.include?(-Math.log2(info.pvalue_1))
}

snvs = SNVInfo.each_in_file(snv_infos_filename).map{|snv| [snv.variant_id, snv] }.to_h

case mode
when :words
  if flank_length
    data = sites.map do |site|
      snv = snvs[site.variant_id]
      seq = snv.load_site_sequence(GENOME_READER, site, flank_length)
      site.orientation_1 == :direct ? seq : Sequence.revcomp(seq)
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
