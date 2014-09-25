require 'set'

# snp = File.readlines('/home/ilya/iogen/BioSchool-summer2014/SNP/SNPs.txt')
# snp_splitted = snp.map{|el| el.split("\t")}
# # 3.times do |snp_shuffle|
#  snp_shuffle = snp_splitted.map do |el|
#     name = el[0]
#     seq = el[1]
#     snp_choose = seq[0..23].each_char.to_a.shuffle.join + seq[24..30] + seq[31..-1].chomp.each_char.to_a.shuffle.join
#     snp_new = [name + '_' + 5.to_s + "\t" + snp_choose]
#   end
#   puts snp_shuffle
# end

# snp = File.readlines('/home/ilya/iogen/BioSchool-summer2014/SNP/SNPs.txt')
# snp_splitted = snp.map{|el| el.split("\t")}

# snp_shuffle = snp_splitted.map do |el|
  
#   (1..8).map do |step_num|

#     name = el[0] + '_' + step_num.to_s
#     # name = "#{el[0]}_#{step_num}"
#     seq = el[1]
    
#     snp_choose = seq[0..23].each_char.to_a.shuffle.join + seq[24..30] + seq[31..-1].chomp.each_char.to_a.shuffle.join
#     name + "\t" + snp_choose
#   end

# end
# puts snp_shuffle.flatten

# sequences.select do |sequence|
#   sequence ...CpT

# end


# sequence.select{}

def count_each_element(more_5_motif)
  result = {}
  more_5_motif.each do |motif_name|
    if result.has_key?(motif_name)
        result[motif_name] += 1
    else
        result[motif_name] = 1
    end
  end
  return result
end  


filename = ARGV[0]   # real_mutations.txt
filename_result = ARGV[1] # results_real_mutations.txt --> cpg_promoter_results_real_mutations.txt, tpc_intronic_results_real_mutations.txt

motif_names = File.readlines('/home/ilya/iogen/BioSchool-summer2014/SNP/motif_names.txt').map(&:strip)


snp = File.readlines('/home/ilya/iogen/BioSchool-summer2014/SNP/SNPs.txt'); # ['name seq', 'name seq',...]
snp_splitted = snp.map{|el| el.split("\t")};  # [['name', 'seq'], ['name', 'seq'], ...]

# [['name','seq'],['name','seq'],['name', 'seq']... ] ---(select)---> [['name','seq'],..]  ---(map)---> ['name',...]

cpg_seq = snp_splitted.select do |name, sequence|
  ((sequence[26] == 'C') && (sequence[30] == 'G')) || ((sequence[26] == 'G') && (sequence[24] == 'C'))
end
cpg_names = Set.new( cpg_seq.map{|name, sequence| name } ) 

#   ((sequence[26] == 'C') && (sequence[24] == 'T')) || ((sequence[26] == 'G') && (sequence[30] == 'A')) 
tpc_seq = snp_splitted.select do |name, sequence|
  ((sequence[26] == 'C') && (sequence[24] == 'T')) || ((sequence[26] == 'G') && (sequence[30] == 'A'))
end
tpc_names = Set.new( tpc_seq.map{|name, sequence| name } )

#p cpg_names.size
#p tpc_names.size
#p cpg_names

mut = File.readlines('/home/ilya/iogen/BioSchool-summer2014/SNP/SUBSTITUTIONS_13Apr2012_snz_promoter_markup2.txt').drop(1);
mut_splitted = mut.map{|el| el.split("\t")};
#mut_types = mut_splitted[0], 
mut_types_intronic = mut_splitted.select do |el|
  el[17] == 'Intronic' || el[17] == 'Intronic,Promoter'  
end
mut_types_promoter = mut_splitted.select do |el|
  el[17] == 'Promoter' || el[17] == 'Intronic,Promoter' 
end;  
# mut_types_intronic_promoter = mut_splitted.select do |el|
#   el[17] == 'Intronic,Promoter'
# end;  
mut_types_name_intronic = Set.new(mut_types_intronic.map{|el| el[0]});
mut_types_name_promoter = Set.new(mut_types_promoter.map{|el| el[0]});


# $stderr.puts "CpG: #{cpg_names.size}\nTpC: #{tpc_names.size}"
# exit


fold_change_cutoff = 1.0 / 5.0
File.open('cpg_intronic_' + filename_result, 'w') do |fw|
  File.open(filename) do |f|
    cpg_more_5 = f.each_line.drop(1).select do |line| # "27610826_3 MAZ_f1  -7  direct  cggctgaGgaggaggag -7  direct  cggctgaCgaggaggag G/C 1.1218764110455249E-4 9.602413003842941E-4  0.11683275970285215"
      line_splitted = line.split("\t") # ["27610826_3", "MAZ_f1",  "-7",  "direct",  "cggctgaGgaggaggag", "-7",  "direct",  "cggctgaCgaggaggag", "G/C", "1.1218764110455249E-4", "9.602413003842941E-4",  "0.11683275970285215"]
      name_snp = line_splitted[0].split("_")[0]
      (line_splitted[-1].to_f < fold_change_cutoff) && (cpg_names.include? name_snp) && (mut_types_name_intronic.include? name_snp)
    end
    cpg_more_5_motif_intronic = cpg_more_5.map{|el| el.split("\t")[1]}
    counts = count_each_element(cpg_more_5_motif_intronic)
    # counts.each do |name,num|
    #   fw.puts name + "\t" + num.to_s
    # end
    motif_names.each do |motif_name|
      fw.puts "#{motif_name}\t#{counts[motif_name] || 0}"
    end
  end
end

File.open('cpg_promoter_' + filename_result, 'w') do |fw|
  File.open(filename) do |f|
    cpg_more_5_promoter = f.each_line.drop(1).select do |line| # "27610826 MAZ_f1  -7  direct  cggctgaGgaggaggag -7  direct  cggctgaCgaggaggag G/C 1.1218764110455249E-4 9.602413003842941E-4  0.11683275970285215"
      line_splitted = line.split("\t") # ["27610826", "MAZ_f1",  "-7",  "direct",  "cggctgaGgaggaggag", "-7",  "direct",  "cggctgaCgaggaggag", "G/C", "1.1218764110455249E-4", "9.602413003842941E-4",  "0.11683275970285215"]
      name_snp = line_splitted[0].split("_")[0]
      (line_splitted[-1].to_f < fold_change_cutoff) && (cpg_names.include? name_snp) && (mut_types_name_promoter.include? name_snp)
    end
    cpg_more_5_motif_promoter = cpg_more_5_promoter.map{|el| el.split("\t")[1]}
    counts = count_each_element(cpg_more_5_motif_promoter)
    # counts.each do |name,num|
    #   fw.puts name + "\t" + num.to_s
    # end
    motif_names.each do |motif_name|
      fw.puts "#{motif_name}\t#{counts[motif_name] || 0}"
    end
  end 
end

File.open('tpc_intronic_' + filename_result, 'w') do |fw|
  File.open(filename) do |f|
    tpc_more_5 = f.each_line.drop(1).select do |line| # "27610826 MAZ_f1  -7  direct  cggctgaGgaggaggag -7  direct  cggctgaCgaggaggag G/C 1.1218764110455249E-4 9.602413003842941E-4  0.11683275970285215"
      line_splitted = line.split("\t") # ["27610826", "MAZ_f1",  "-7",  "direct",  "cggctgaGgaggaggag", "-7",  "direct",  "cggctgaCgaggaggag", "G/C", "1.1218764110455249E-4", "9.602413003842941E-4",  "0.11683275970285215"]
      name_snp = line_splitted[0].split("_")[0]
      (line_splitted[-1].to_f < fold_change_cutoff) && (tpc_names.include? name_snp) && (mut_types_name_intronic.include? name_snp)
    end
    tpc_more_5_motif_intronic = tpc_more_5.map{|el| el.split("\t")[1]} 
    counts = count_each_element(tpc_more_5_motif_intronic)
    # counts.each do |name,num|
    #   fw.puts name + "\t" + num.to_s
    # end
    motif_names.each do |motif_name|
      fw.puts "#{motif_name}\t#{counts[motif_name] || 0}"
    end
  end
end

File.open('tpc_promoter_' + filename_result, 'w') do |fw|
  File.open(filename) do |f|
    tpc_more_5_promoter = f.each_line.drop(1).select do |line| # "27610826 MAZ_f1  -7  direct  cggctgaGgaggaggag -7  direct  cggctgaCgaggaggag G/C 1.1218764110455249E-4 9.602413003842941E-4  0.11683275970285215"
      line_splitted = line.split("\t") # ["27610826", "MAZ_f1",  "-7",  "direct",  "cggctgaGgaggaggag", "-7",  "direct",  "cggctgaCgaggaggag", "G/C", "1.1218764110455249E-4", "9.602413003842941E-4",  "0.11683275970285215"]
      name_snp = line_splitted[0].split("_")[0]
      (line_splitted[-1].to_f < fold_change_cutoff) && (tpc_names.include? name_snp) && (mut_types_name_promoter.include? name_snp)
    end
    tpc_more_5_motif_promoter = tpc_more_5_promoter.map{|el| el.split("\t")[1]}
    counts = count_each_element(tpc_more_5_motif_promoter) 
    # counts.each do |name,num|
    #   fw.puts name + "\t" + num.to_s
    # end
    motif_names.each do |motif_name|
      fw.puts "#{motif_name}\t#{counts[motif_name] || 0}"
    end
  end 
end
