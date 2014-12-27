$:.unshift File.absolute_path('../../lib', __dir__)
require 'set'
# require 'mutation_features'
require 'sequence_with_snp'
require 'support'
require 'optparse'
require 'mutation_context'
require 'import_information'

requested_mutation_types = nil
requested_mutation_contexts = nil
show_possible_mutation_types_and_exit = false
OptionParser.new do |opts|
  opts.banner = "Usage:\n" +
                "   #{opts.program_name} <site fold-changes file> [options]\n" + 
                "   <site fold-changes> | #{opts.program_name} [options]"

  opts.separator 'Options:'
  opts.on('--mutation-types POSS_TYPES', "Specify possible mutation types", "(e.g. `Intronic,Promoter`)") {|value|
    requested_mutation_types = value.split(',')
  }

  opts.on('--contexts POSS_CONTEXTS',  "Specify possible mutation contexts (e.g. `TCN,NCG`)",
                                        "`N` represents any-nucleotide placeholder. `N`-s are independent",
                                        "Context is a 3nt string, with center representing a mutation" ) {|value|
    requested_mutation_contexts = value.upcase.split(',').map{|context| MutationContext.from_string(context) }
  }

  opts.on('--possible-mutation-types', 'Show possible mutation types and exit') {
    show_possible_mutation_types_and_exit = true
  }

end.parse!(ARGV)

if requested_mutation_contexts
  requested_mutation_contexts = (requested_mutation_contexts + requested_mutation_contexts.map(&:revcomp)).uniq
else
  requested_mutation_contexts = [MutationContext.new('N','N','N')]
end

requested_mutation_types = requested_mutation_types.map(&:downcase)

mutations_markup_filename = 'source_data/SUBSTITUTIONS_13Apr2012_snz_promoter_markup2.txt'
if show_possible_mutation_types_and_exit
  puts possible_mutation_types(mutations_markup_filename)
  exit
end

snp_sequences_filename = 'source_data/SNPs.txt'

matching_type_mut_names = matching_type_mutation_names(mutations_markup_filename, requested_mutation_types)
matching_context_mut_names = matching_context_mutation_names(snp_sequences_filename, requested_mutation_contexts)
requested_mutation_names = matching_type_mut_names & matching_context_mut_names

output_filtered = ->(mutated_site_info) do
  next  unless requested_mutation_names.include?(mutated_site_info.normalized_snp_name)
  puts mutated_site_info.line
end

site_fold_changes_filename = ARGV[0] # 'source_data/cancer_SNPs.txt'

if site_fold_changes_filename
  MutatatedSiteInfo.each_site(site_fold_changes_filename, &output_filtered)
else
  MutatatedSiteInfo.each_site_in_stream($stdin, &output_filtered)
end
