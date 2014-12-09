$:.unshift File.absolute_path('../../lib', __dir__)
require 'support'

def each_site(mutation_filename, &block)
  return enum_for(:each_site, mutation_filename).lazy  unless block_given?
  each_mutated_site_info(mutation_filename).select{|mutated_site_info|
    mutated_site_info.pvalue_1 <= 0.0005
  }.each(&block)
end

sites_filename, motif_name, from, to = ARGV.first(4)
mode = (ARGV[4] || :words).to_sym
raise "Specify file with sites, motif name, and bin boundaries (-log2 Pvalue) and (optional) mode: words or variant_ids" unless sites_filename && motif_name && from && to

from = from.to_f
to = to.to_f
range = from..to
motifs = each_site(sites_filename).select{|info| 
  info.motif_name == motif_name && range.include?(-Math.log2(info.pvalue_1)) 
}

case mode
when :words
  data = motifs.map(&:seq_1)
when :variant_ids
  data = motifs.map(&:variant_id)
else
  raise ArgumentError, 'Inknown mode'
end

data.each{|x| puts x }

# puts words
