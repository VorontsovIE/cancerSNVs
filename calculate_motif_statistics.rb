require 'set'

def count_each_element(more_5_motif)
  result = {}
  more_5_motif.each do |motif_name|
    if result.has_key?(motif_name)
        result[motif_name] += 1
    else
        result[motif_name] = 1
    end
  end
  result
end  

def cpg_mutation?(sequence)
  (sequence[26] == 'C' && sequence[30] == 'G') || (sequence[26] == 'G' && sequence[24] == 'C')
end

def tpc_mutation?(sequence)
  (sequence[26] == 'C' && sequence[24] == 'T') || (sequence[26] == 'G' && sequence[30] == 'A')
end

def count_each_motif_mutations(all_mutations_filename, mutations_subset, fold_change_cutoff)
  File.open(all_mutations_filename) do |f|
    # "27610826_3 MAZ_f1  -7  direct  cggctgaGgaggaggag -7  direct  cggctgaCgaggaggag G/C 1.1218764110455249E-4 9.602413003842941E-4  0.11683275970285215"
    mutated_sites = f.each_line.drop(1).select do |line|
      line_splitted = line.split("\t")
      name_snp = line_splitted[0].split("_")[0]
      (line_splitted[-1].to_f < fold_change_cutoff) && mutations_subset.include?(name_snp)
    end
    mutated_motifs = mutated_sites.map{|el| el.split("\t")[1] }
    count_each_element(mutated_motifs)
  end
end

def output_counts_for_each_motif(output_filename, motif_names, motif_counts)
  File.open(output_filename, 'w') do |fw|
    motif_names.each do |motif_name|
      fw.puts "#{motif_name}\t#{motif_counts[motif_name] || 0}"
    end
  end
end

def filename_w_prefix(filename, prefix)
  File.join(File.dirname(filename), prefix + File.basename(filename))
end

def motif_statistics(filename_result, file_prefix, motif_names, input_file, mutations_subset, fold_change_cutoff)
  output_counts_for_each_motif(filename_w_prefix(filename_result, file_prefix), motif_names, 
                            count_each_motif_mutations(input_file, mutations_subset, fold_change_cutoff)) 
end

filename = ARGV[0]   # real_mutations.txt
filename_result = ARGV[1] # results_real_mutations.txt --> cpg_promoter_results_real_mutations.txt, tpc_intronic_results_real_mutations.txt

Dir.mkdir File.dirname(filename_result) unless Dir.exist?(File.dirname(filename_result))

raise 'Specify input file with mutations and output filename for results(actually 4 files will be created)'  unless filename && filename_result

motif_names = File.readlines('./motif_names.txt').map(&:strip)

snp_splitted = File.readlines('./SNPs.txt').map{|el| el.split("\t")}

cpg_names = Set.new( snp_splitted.select {|name, sequence| cpg_mutation?(sequence) }.map{|name, sequence| name } ) 
tpc_names = Set.new( snp_splitted.select {|name, sequence| tpc_mutation?(sequence) }.map{|name, sequence| name } )

$stderr.puts "CpG: #{cpg_names.size}\nTpC: #{tpc_names.size}"

mut_types = File.readlines('./SUBSTITUTIONS_13Apr2012_snz_promoter_markup2.txt').drop(1).map{|el| data = el.split("\t"); [data[0], data[17]] };
intronic_mutations = mut_types.select {|mut_name, mut_type| mut_type.split(',').map(&:strip).include? 'Intronic' }
promoter_mutations = mut_types.select {|mut_name, mut_type| mut_type.split(',').map(&:strip).include? 'Promoter' }

regulatory_mutations = mut_types.select {|mut_name, mut_type| types = mut_type.split(',').map(&:strip); types.include?('Promoter') || types.include?('Intronic') }

# mut_types_intronic_promoter = mut_types.select {|mut_name, mut_type| mut_type == 'Intronic,Promoter' }

intronic_mutation_names = Set.new(intronic_mutations.map{|mut_name, mut_type| mut_name} )
promoter_mutation_names = Set.new(promoter_mutations.map{|mut_name, mut_type| mut_name} )
regulatory_mutation_names = Set.new(regulatory_mutations.map{|mut_name, mut_type| mut_name} )

cpg_intronic_names = cpg_names & intronic_mutation_names
tpc_intronic_names = tpc_names & intronic_mutation_names
cpg_promoter_names = cpg_names & promoter_mutation_names
tpc_promoter_names = tpc_names & promoter_mutation_names
cpg_regulatory_names = cpg_names & regulatory_mutation_names
tpc_regulatory_names = tpc_names & regulatory_mutation_names


$stderr.puts "CpG intronic: #{cpg_intronic_names.size}\nCpG promoter: #{cpg_promoter_names.size}\nTpC intronic: #{tpc_intronic_names.size}\nTpC promoter: #{tpc_promoter_names.size}"
$stderr.puts "Regulatory: #{regulatory_mutation_names.size}"
$stderr.puts "CpG regulatory: #{cpg_regulatory_names.size}\nTpC regulatory: #{tpc_regulatory_names.size}"

fold_change_cutoff = 1.0 / 5.0


motif_statistics(filename_result, 'cpg_regulatory_', motif_names, filename, cpg_regulatory_names, fold_change_cutoff)
motif_statistics(filename_result, 'tpc_regulatory_', motif_names, filename, tpc_regulatory_names, fold_change_cutoff)

motif_statistics(filename_result, 'cpg_intronic_', motif_names, filename, cpg_intronic_names, fold_change_cutoff)
motif_statistics(filename_result, 'tpc_intronic_', motif_names, filename, tpc_intronic_names, fold_change_cutoff)

motif_statistics(filename_result, 'cpg_promoter_', motif_names, filename, cpg_promoter_names, fold_change_cutoff)
motif_statistics(filename_result, 'tpc_promoter_', motif_names, filename, tpc_promoter_names, fold_change_cutoff)
