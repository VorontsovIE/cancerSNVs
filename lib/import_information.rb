require 'sequence_with_snp'
require 'mutation_context'

def load_mutation_contexts(snp_sequences_filename)
  snps_splitted = File.readlines(snp_sequences_filename).map{|el|
    name, seq = el.chomp.split("\t")
    snv = SequenceWithSNP.from_string(seq, name: name)
    [snv.name, snv.context(before: 1, after: 1, allele_variant_number: 0)]
  }
end

def context_matches?(mutation_context, requested_mutation_contexts)
  requested_mutation_contexts.any?{|req_mut_context|
    req_mut_context.match?(mutation_context)
  }
end


## filter mutations by context
def matching_context_mutation_names(snp_sequences_filename, requested_mutation_contexts)
  mutation_contexts = load_mutation_contexts(snp_sequences_filename)
  mutation_contexts.select{|variant_id, mutation_context|
    context_matches?(mutation_context, requested_mutation_contexts)
  }.map{|variant_id, mutation_context| normalized_snp_name(variant_id) }.to_set
end

def normalized_snp_name(variant_id)
  variant_id.split('_')[0]
end
