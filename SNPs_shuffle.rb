require 'set'

def shuffle_string(str)
  str.each_char.to_a.shuffle.join
end
snp = File.readlines('source_data/SNPs.txt').map{|el| el.split("\t")}

(9..10).map do |step_num|
  snp_shuffle = snp.map do |snp_name, seq|
    shuffled_snp_name = "#{snp_name}_#{step_num}"
    shuffled_seq = shuffle_string(seq[0..23]) + seq[24..30] + shuffle_string(seq[31..-1].chomp)
    "#{shuffled_snp_name}\t#{shuffled_seq}"
  end
  File.write "shuffled_mutations/SNPs_shuffle_#{step_num}.txt", snp_shuffle.join("\n")
end
