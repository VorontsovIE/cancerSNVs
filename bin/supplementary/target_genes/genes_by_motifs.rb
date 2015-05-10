$:.unshift File.absolute_path('../../../lib', __dir__)
require_relative '../../../experiment_configuration'

genes_by_motif = File.readlines(LocalPaths::Secondary::GeneInfos).drop(1)
                      .map{|line| line.chomp.split("\t").first(2) }
                      .map{|motif,genes| [motif, genes.split.first] }.to_h

puts ARGF.read.split.flat_map{|motif| genes_by_motif[motif] }
