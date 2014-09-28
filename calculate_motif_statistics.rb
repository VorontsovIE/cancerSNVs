$:.unshift File.absolute_path('lib', __dir__)
require 'set'
require 'statistical_significance'
require 'mutation_features'
require 'support'

mutation_infos_filename = ARGV[0]   # cancer.txt
filename_result = ARGV[1]           # results/cancer.txt --> results/cpg_promoter_cancer.txt,
                                    #                        results/tpc_intronic_cancer.txt, ...

raise 'Specify input file with mutation infos and output filename for results(actually 4 files will be created)'  unless mutation_infos_filename && filename_result

motif_names = File.readlines('./motif_names.txt').map(&:strip)

output_configurator = BreastCancerSNPs::OutputConfigurator.new(filename_result, motif_names)


snps_splitted = File.readlines('./SNPs.txt').map{|el| el.chomp.split("\t")}
cpg_names = mutation_names_by_mutation_context(snps_splitted){|mut_name, sequence| cpg_mutation?(sequence) }
tpc_names = mutation_names_by_mutation_context(snps_splitted){|mut_name, sequence| tpc_mutation?(sequence) }
all_names = mutation_names_by_mutation_context(snps_splitted){ true }

mut_types = File.readlines('./SUBSTITUTIONS_13Apr2012_snz_promoter_markup2.txt').drop(1).map{|el| data = el.split("\t"); [data[0], data[17]] };
intronic_mutation_names = mutation_names_by_mutation_type(mut_types){|mut_name, mut_type| intronic_mutation?(mut_type) }
promoter_mutation_names = mutation_names_by_mutation_type(mut_types){|mut_name, mut_type| promoter_mutation?(mut_type) }
regulatory_mutation_names = mutation_names_by_mutation_type(mut_types){|mut_name, mut_type| intronic_mutation?(mut_type) || promoter_mutation?(mut_type) }

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
