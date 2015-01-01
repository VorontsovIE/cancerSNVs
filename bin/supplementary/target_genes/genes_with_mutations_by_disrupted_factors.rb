require 'optparse'
can_use_ensg = false
fold_change_cutoff = 1
OptionParser.new do |opts|
  opts.on('--min-foldchange FC', "minimal fold change for disruption. FC = 5 will be treated as 0.2") do |fold_change|
    fold_change = fold_change.to_f
    fold_change_cutoff = (fold_change > 1)  ?  (1.0 / fold_change) : fold_change
  end
  opts.on('--can-use-ensg', "Use ENSG identifiers when HGNC not available. Ignore such genes if option not enabled") do
    can_use_ensg = true
  end
end.parse!(ARGV)

motifs = ARGV

mutations = File.readlines('./source_data/MutatedGenes.csv').map do |line|
  line.chomp.split("\t")
end

genes = mutations.select{|motif, mut_name, chr, pos, mut_type, context, ensg, hgnc, fold_change, pvalue_1, pvalue_2|
  fold_change = fold_change.to_f
  motifs.include?(motif) && fold_change <= fold_change_cutoff # fold_change for disruption is not greater than 1.0
}.map{|motif, mut_name, chr, pos, mut_type, context, ensg, hgnc, fold_change, pvalue_1, pvalue_2|
  if can_use_ensg
    hgnc.empty? ? ensg : hgnc
  else
    hgnc
  end
}.compact.uniq.reject(&:empty?)

puts genes.join("\n")
