genes_by_motif = File.readlines('source_data/hocomoco_genes.csv').drop(1)
                      .map{|line| line.chomp.split("\t").first(2) }
                      .map{|motif,genes| [motif.sub(/\bHOCOMOCO\/(.+)\.pat\b/, '\1'), genes.split] }.to_h

puts ARGF.read.split.flat_map{|motif| genes_by_motif[motif] }
