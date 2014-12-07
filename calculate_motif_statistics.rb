$:.unshift File.absolute_path('lib', __dir__)
require 'set'
require 'statistical_significance'
require 'mutation_features'
require 'support'
require 'output_configurator'

mutated_site_infos_filename = ARGV[0]   # cancer.txt
filename_result = ARGV[1]           # results/cancer.txt --> results/cpg_promoter_cancer.txt,
                                    #                        results/tpc_intronic_cancer.txt, ...

raise 'Specify input file with mutation infos and output filename for results(actually 4 files will be created)'  unless mutated_site_infos_filename && filename_result

motif_names = File.readlines('source_data/motif_names.txt').map(&:strip)

output_configurator = BreastCancerSNPs::OutputConfigurator.new(filename_result, motif_names)

# puts File.readlines('source_data/SNPs.txt').flat_map{|el| seq = el.chomp.split("\t")[1]; seq = seq.gsub(/\[(.)\/.\]/,'\1'); seq.each_char.each_cons(3).map{|x|x.to_a.join} }.each_with_object(Hash.new{|hsh,k| hsh[k]=0}){|tri,hsh| hsh[tri] += 1}.select{|k,v| k.start_with?('CG') || k.end_with?('CG')}.map{|k,v|v}.inject(0,&:+)
# exit 1
snps_splitted = File.readlines('source_data/SNPs.txt').map{|el| el.chomp.split("\t")}

cpg_names = mutation_names_by_mutation_context(snps_splitted){|mut_name, sequence| cpg_mutation?(sequence) }
tpc_names = mutation_names_by_mutation_context(snps_splitted){|mut_name, sequence| tpc_mutation?(sequence) }
not_cpg_tpc_names = mutation_names_by_mutation_context(snps_splitted){|mut_name, sequence| !tpc_mutation?(sequence) && !cpg_mutation?(sequence) }
any_context_names = mutation_names_by_mutation_context(snps_splitted){ true }
########
mut_types = File.readlines('source_data/SUBSTITUTIONS_13Apr2012_snz_promoter_markup2.txt').drop(1).map{|el| data = el.split("\t"); [data[0], data[17]] }
intronic_mutation_names = mutation_names_by_mutation_type(mut_types){|mut_name, mut_type| intronic_mutation?(mut_type) }
promoter_mutation_names = mutation_names_by_mutation_type(mut_types){|mut_name, mut_type| promoter_mutation?(mut_type) }
regulatory_mutation_names = mutation_names_by_mutation_type(mut_types){|mut_name, mut_type| intronic_mutation?(mut_type) || promoter_mutation?(mut_type) }
############

$stderr.puts "CpG: #{cpg_names.size}; TpC: #{tpc_names.size}; non-CpG-non-TpC: #{not_cpg_tpc_names.size}"
$stderr.puts  "regulatory: #{regulatory_mutation_names.size}; " +
              "regulatory CpG: #{(cpg_names & regulatory_mutation_names).size}; " +
              "regulatory TpC: #{(tpc_names & regulatory_mutation_names).size}; " +
              "non-CpG-non-TpC: #{(not_cpg_tpc_names & regulatory_mutation_names).size}"
# exit 1
# $stderr.puts "CpG intronic: #{cpg_intronic_names.size}\nCpG promoter: #{cpg_promoter_names.size}\nTpC intronic: #{tpc_intronic_names.size}\nTpC promoter: #{tpc_promoter_names.size}"
# $stderr.puts "Regulatory: #{regulatory_mutation_names.size}"
# $stderr.puts "CpG regulatory: #{cpg_regulatory_names.size}\nTpC regulatory: #{tpc_regulatory_names.size}"


# mutation_types = {regulatory: regulatory_mutation_names, intronic: intronic_mutation_names, promoter: promoter_mutation_names}
mutation_types = {regulatory: regulatory_mutation_names}
context_types = {cpg: cpg_names, tpc: tpc_names, not_cpg_tpc: not_cpg_tpc_names, any_context: any_context_names}

mutation_types.each do |mutation_type, mutation_type_nameset|
  context_types.each do |context_type, context_type_nameset|
    mutations_nameset = mutation_type_nameset & context_type_nameset
    output_configurator.motif_statistics( "#{mutation_type}_#{context_type}_disrupted_",
                                          created_motifs_in_set(mutated_site_infos_filename, mutations_nameset) )
    $stderr.puts "#{mutation_type}, #{context_type} -- disrupted"
    output_configurator.motif_statistics( "#{mutation_type}_#{context_type}_total_",
                                          improved_motifs_in_set(mutated_site_infos_filename, mutations_nameset) )
    $stderr.puts "#{mutation_type}, #{context_type} -- total"
  end
end
