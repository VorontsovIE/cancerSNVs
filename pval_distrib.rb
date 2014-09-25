# провести сэмплирование случайных сайтов из бинов так, чтобы распределение pvalues было аналогично геномному

pvals = []

File.open('shuffled_mutations/all_res_SNPs_shuffle.txt') do |f|
# File.open('all_SNPs.txt') do |f|
  f.each_line do |l|
    row = l.chomp.split("\t")
    mot = row[1]
    pval = row[9].to_f
    pvals << pval  if mot == 'HIF1A_si'
  end
end

puts pvals.sort.map.with_index{|x,i| "#{i.to_f/pvals.size}\t#{x}" }.join("\n")
