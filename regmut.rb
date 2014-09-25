def seq_by_pos(pos, chr, flank_size)
  File.open ("/home/ilya/iogen/genome/hg19/chr#{chr}.plain") do |f|
    f.seek(pos - flank_size - 1)
    f.read(flank_size*2 + 1)
  end
end

def complement_letter(letter)
  letter = letter.upcase
  if letter == 'A' 
    return 'T'
  elsif letter == 'C'
    return 'G'
  elsif letter == 'G' 
    return 'C'
  elsif letter == 'T' 
    return 'A'
  else
    raise "#{letter} is not a nucleotide"  
  end
end

def complement(seq)
  seq_arr = seq.each_char.to_a
  seq_arr.map{|letter| complement_letter(letter)}.join
end

 hi = File.readlines('/home/ilya/iogen/BioSchool-summer2014/SNP/SUBSTITUTIONS_13Apr2012_snz.txt') [1..-1];
 # hi.map{|line| line.rstrip}
 stripped_lines = hi.map{|line| line.rstrip};
 splitted_lines = stripped_lines.map{|line| line.split("\t")};
 #chr_num = splitted_lines.map{|line| line[2]}
#chr = chr_num[1..-1]
#pos = splitted_lines.map{|line| line[3]}
#posit = pos.map{|el| el.to_f}

mut_coords = splitted_lines.map{|el| [el[2], el[3].to_i, el[6], el[11], el[0], el[7], el[10], el[5]]};

#p mut_coords





result = mut_coords.map do |mut|
  seq = seq_by_pos(mut[1], mut[0], 25).upcase
  seq_5 = mut[5].upcase
  seq_3 = mut[6].upcase
  seq_centre = mut[7].upcase
  if mut[3] == '+'
    raise 'Error'  if (seq_5 != seq[15..24]) || (seq_3 != seq[26..35]) || (seq_centre != seq[25])
  else
    seq_rev_5 = complement(seq_5).reverse
    seq_rev_3 = complement(seq_3).reverse
    if (seq_rev_3 != seq[15..24]) || (seq_rev_5 != seq[26..35]) || (seq_centre != seq[25])
      # raise "Error:\n#{mut[4]} #{seq}\n#{seq_rev_5}, #{seq_rev_3}, #{seq_centre}\n#{seq[26..35]}, #{seq[15..24]}, #{seq[25]}\n#{mut[3]}"
      raise "Error\n#{seq_rev_3} != #{seq[15..24]}\n or\n#{seq_rev_5} != #{seq[26..35]}\nor\n#{seq_centre} != #{seq[25]}"
    end
  end

  mut[4] + "\t" + seq[0..24] + '[' + seq[25] + '/' + mut[2] + ']' + seq[26..-1]
end;
puts result
