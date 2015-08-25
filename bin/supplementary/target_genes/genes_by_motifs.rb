$:.unshift File.absolute_path('../../../lib', __dir__)
require_relative '../../../experiment_configuration'
require 'uniprot_info'

uniprot_infos_by_motif = read_uniprot_infos_by_motif(LocalPaths::Secondary::HocomocoUniprots, LocalPaths::Secondary::UniprotDump)

puts ARGF.read.split.flat_map{|motif| uniprot_infos_by_motif[motif].primary_gene_name }
