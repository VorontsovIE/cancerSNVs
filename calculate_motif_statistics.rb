$:.unshift File.absolute_path('lib', __dir__)
require 'set'
require 'statistical_significance'
require 'mutation_features'
require 'support'

def count_each_motif_mutations(all_mutations_filename)
  File.open(all_mutations_filename) do |f|
    # "27610826_3 MAZ_f1  -7  direct  cggctgaGgaggaggag -7  direct  cggctgaCgaggaggag G/C 1.1218764110455249E-4 9.602413003842941E-4  0.11683275970285215"
    mutated_sites = f.each_line.drop(1).select do |line|
      name_snp, motif_name,
                pos_1,orientation_1,seq_1,
                pos_2,orientation_2,seq_2,
                variants,
                pvalue_1, pvalue_2, fold_change = line.split("\t")

      pos_1 = pos_1.to_i
      pos_2 = pos_2.to_i
      orientation_1 = orientation_1.to_sym
      orientation_2 = orientation_2.to_sym
      pvalue_1 = pvalue_1.to_f
      pvalue_2 = pvalue_2.to_f
      fold_change = fold_change.to_f

      yield name_snp, motif_name, fold_change, pvalue_1, pvalue_2
    end
    mutated_motifs = mutated_sites.map{|el| el.split("\t")[1] }
    count_each_element(mutated_motifs)
  end
end

mutation_infos_filename = ARGV[0]   # real_mutations.txt
filename_result = ARGV[1] # results_real_mutations.txt --> cpg_promoter_results_real_mutations.txt, tpc_intronic_results_real_mutations.txt

raise 'Specify input file with mutation infos and output filename for results(actually 4 files will be created)'  unless mutation_infos_filename && filename_result
# Dir.mkdir File.dirname(filename_result) unless Dir.exist?(File.dirname(filename_result))

motif_names = File.readlines('./motif_names.txt').map(&:strip)

output_configurator = BreastCancerSNPs::OutputConfigurator.new(filename_result, motif_names)


snp_splitted = File.readlines('./SNPs.txt').map{|el| el.split("\t")}

cpg_names = snp_splitted.select {|mut_name, sequence| cpg_mutation?(sequence) }.map{|mut_name, sequence| mut_name }.to_set
tpc_names = snp_splitted.select {|mut_name, sequence| tpc_mutation?(sequence) }.map{|mut_name, sequence| mut_name }.to_set
all_names = snp_splitted.map{|mut_name, sequence| mut_name }.to_set

# $stderr.puts "CpG: #{cpg_names.size}\nTpC: #{tpc_names.size}"

mut_types = File.readlines('./SUBSTITUTIONS_13Apr2012_snz_promoter_markup2.txt').drop(1).map{|el| data = el.split("\t"); [data[0], data[17]] };

intronic_mutation_names = mut_types.select{|mut_name, mut_type| intronic_mutation?(mut_type) }
                              .map{|mut_name, mut_type| mut_name}.to_set
promoter_mutation_names = mut_types.select{|mut_name, mut_type| promoter_mutation?(mut_type) }
                              .map{|mut_name, mut_type| mut_name}.to_set
regulatory_mutation_names = mut_types.select{|mut_name, mut_type| intronic_mutation?(mut_type) || promoter_mutation?(mut_type) }
                                .map{|mut_name, mut_type| mut_name}.to_set

# mut_types_intronic_promoter = mut_types.select {|mut_name, mut_type| mut_type == 'Intronic,Promoter' }

cpg_intronic_names = cpg_names & intronic_mutation_names
tpc_intronic_names = tpc_names & intronic_mutation_names
cpg_promoter_names = cpg_names & promoter_mutation_names
tpc_promoter_names = tpc_names & promoter_mutation_names
cpg_regulatory_names = cpg_names & regulatory_mutation_names
tpc_regulatory_names = tpc_names & regulatory_mutation_names


# $stderr.puts "CpG intronic: #{cpg_intronic_names.size}\nCpG promoter: #{cpg_promoter_names.size}\nTpC intronic: #{tpc_intronic_names.size}\nTpC promoter: #{tpc_promoter_names.size}"
# $stderr.puts "Regulatory: #{regulatory_mutation_names.size}"
# $stderr.puts "CpG regulatory: #{cpg_regulatory_names.size}\nTpC regulatory: #{tpc_regulatory_names.size}"

disrupted_and_in_set_checker = ->(mutations_subset) do
  combine_conditions(disrupted_site_checker(5), mutation_in_set_checker(mutations_subset))
end

disrupted_motifs_in_set = ->(mutations_subset) do
  count_each_motif_mutations(mutation_infos_filename, &disrupted_and_in_set_checker.call(mutations_subset))
end

output_configurator.motif_statistics('cpg_regulatory_', disrupted_motifs_in_set.call(cpg_regulatory_names))
output_configurator.motif_statistics('tpc_regulatory_', disrupted_motifs_in_set.call(tpc_regulatory_names))

output_configurator.motif_statistics('cpg_intronic_', disrupted_motifs_in_set.call(cpg_intronic_names))
output_configurator.motif_statistics('tpc_intronic_', disrupted_motifs_in_set.call(tpc_intronic_names))

output_configurator.motif_statistics('cpg_promoter_', disrupted_motifs_in_set.call(cpg_promoter_names))
output_configurator.motif_statistics('tpc_promoter_', disrupted_motifs_in_set.call(tpc_promoter_names))
