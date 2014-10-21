genes_by_motif = File.readlines('source_data/hocomoco_genes_more_info.csv').drop(1)
                      .map{|line| line.chomp.split("\t").first(2) }
                      .map{|motif,genes| [motif, genes.split] }.to_h

puts ARGF.read.split.flat_map{|motif| genes_by_motif[motif] }
