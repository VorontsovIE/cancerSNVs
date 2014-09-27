$:.unshift File.absolute_path('lib', __dir__)
require 'set'
require 'statistical_significance'
require 'mutation_features'
require 'support'

def count_each_motif_mutations(all_mutations_filename, mutations_subset, fold_change_cutoff)
  File.open(all_mutations_filename) do |f|
    # "27610826_3 MAZ_f1  -7  direct  cggctgaGgaggaggag -7  direct  cggctgaCgaggaggag G/C 1.1218764110455249E-4 9.602413003842941E-4  0.11683275970285215"
    mutated_sites = f.each_line.drop(1).select do |line|
      line_splitted = line.split("\t")
      name_snp = line_splitted[0].split("_")[0]
      fold_change = line_splitted[-1].to_f
      fold_change < fold_change_cutoff && mutations_subset.include?(name_snp)
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

intronic_mutations = mut_types.select{|mut_name, mut_type| intronic_mutation?(mut_type) }
                              .map{|mut_name, mut_type| mut_name}.to_set
promoter_mutations = mut_types.select{|mut_name, mut_type| promoter_mutation?(mut_type) }
                              .map{|mut_name, mut_type| mut_name}.to_set
regulatory_mutations = mut_types.select{|mut_name, mut_type| intronic_mutation?(mut_type) || promoter_mutation?(mut_type) }
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

fold_change_cutoff = 1.0 / 5.0

output_configurator.motif_statistics('cpg_regulatory_', count_each_motif_mutations(mutation_infos_filename, cpg_regulatory_names, fold_change_cutoff))
output_configurator.motif_statistics('tpc_regulatory_', count_each_motif_mutations(mutation_infos_filename, tpc_regulatory_names, fold_change_cutoff))

output_configurator.motif_statistics('cpg_intronic_', count_each_motif_mutations(mutation_infos_filename, cpg_intronic_names, fold_change_cutoff))
output_configurator.motif_statistics('tpc_intronic_', count_each_motif_mutations(mutation_infos_filename, tpc_intronic_names, fold_change_cutoff))

output_configurator.motif_statistics('cpg_promoter_', count_each_motif_mutations(mutation_infos_filename, cpg_promoter_names, fold_change_cutoff))
output_configurator.motif_statistics('tpc_promoter_', count_each_motif_mutations(mutation_infos_filename, tpc_promoter_names, fold_change_cutoff))
