require 'set'

MutatatedSiteInfo = Struct.new( :line,
                                :variant_id, :motif_name,
                                :fold_change, :pvalue_1, :pvalue_2,
                                :pos_1, :orientation_1, :seq_1,
                                :pos_2, :orientation_2, :seq_2,
                                :variants ) do

  # "27610826_3 MAZ_f1  -7  direct  cggctgaGgaggaggag -7  direct  cggctgaCgaggaggag G/C 1.1218764110455249E-4 9.602413003842941E-4  0.11683275970285215"
  def self.from_string(line)
    variant_id, motif_name,
              pos_1, orientation_1, seq_1,
              pos_2, orientation_2, seq_2,
              variants,
              pvalue_1, pvalue_2, fold_change = line.split("\t")
    MutatatedSiteInfo.new(line,
                          variant_id, motif_name,
                          fold_change.to_f, pvalue_1.to_f, pvalue_2.to_f,
                          pos_1.to_i, orientation_1.to_sym, seq_1,
                          pos_2.to_i, orientation_2.to_sym, seq_2,
                          variants)
  end

  def normalized_snp_name
    variant_id.split("_")[0]
  end

  def length
    seq_1.length
  end

  def seq_1_direct_strand
    orientation_1 == :direct ? seq_1 : revcomp(seq_1)
  end

  def seq_1_five_flank_length
    -pos_1
  end

  def seq_1_three_flank_length
    length - 1 + pos_1
  end
end

def each_mutated_site_info_in_stream(stream, &block)
  stream.each_line.lazy.reject{|line|
    line.start_with?('#')
  }.map{|line|
    MutatatedSiteInfo.from_string(line)
  }.each(&block)
end

def each_mutated_site_info(all_mutations_filename, &block)
  return enum_for(:each_mutated_site_info, all_mutations_filename).lazy  unless block_given?
  File.open(all_mutations_filename) do |f|
    each_mutated_site_info_in_stream(f, &block)
  end
end

def each_site(mutation_filename, pvalue_cutoff: 0.0005, &block)
  return enum_for(:each_site, mutation_filename, pvalue_cutoff: pvalue_cutoff).lazy  unless block_given?
  each_mutated_site_info(mutation_filename).select{|mutated_site_info|
    mutated_site_info.pvalue_1 <= pvalue_cutoff
  }.each(&block)
end
