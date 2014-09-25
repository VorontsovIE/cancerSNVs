# tsses should be sorted and have indices [[123,0],[199,1],[350,2],...]
def find_nearest_tss(tsses, pos)
  val, ind = tsses.bsearch{|val,ind| val >= pos }
  return []  unless ind
  # select tss-es in upstream and in downstream (because promoter region can be assymetric)
  tsses[ [ind - 1, 0].max, 2 ].map{|tss,ind| tss}
end

def pos_in_promoter?(chr, pos, gene_tss_sorted_by_chromosome)
  tsses_pos_plus = find_nearest_tss(gene_tss_sorted_by_chromosome["#{chr},+"], pos)
  tsses_pos_minus = find_nearest_tss(gene_tss_sorted_by_chromosome["#{chr},-"], pos)
  # p chr, pos, tsses_pos_plus, tsses_pos_minus
  tsses_pos_plus.any? {|tss_pos_plus| tss_pos_plus && pos >= tss_pos_plus - 2000 && pos <= tss_pos_plus + 500 } ||
    tsses_pos_minus.any? {|tss_pos_minus| tss_pos_minus && pos >= tss_pos_minus - 500 && pos <= tss_pos_minus + 2000 }
end

gene_tss_by_chromosome = Hash.new{|h,k| h[k] = [] }

File.open('gene_tss.txt') do |f|
  f.each_line do |line|
    chr, strand, tss = line.chomp.split("\t").last(3)
    if strand == '1'
      gene_tss_by_chromosome["#{chr},+"] << tss.to_i
    elsif strand == '-1'
      gene_tss_by_chromosome["#{chr},-"] << tss.to_i
    end
  end
end

gene_tss_sorted_by_chromosome = gene_tss_by_chromosome.map{|chr, tsses| [chr, tsses.uniq.sort.each_with_index.to_a] }.to_h
# p gene_tss_sorted_by_chromosome.keys

File.open('/home/ilya/iogen/BioSchool-summer2014/SNP/SUBSTITUTIONS_13Apr2012_snz.txt') do |f|
  puts f.readline
  f.each_line do |line|
    infos = line.chomp.split("\t")
    chr, pos = infos[2], infos[3].to_i

    if pos_in_promoter?(chr, pos, gene_tss_sorted_by_chromosome)
      mut_type = infos[17]
      infos[17] = [mut_type.strip.empty? ? nil : mut_type, 'Promoter'].compact.join(',')
    end
    puts infos.join("\t")
  end
end
