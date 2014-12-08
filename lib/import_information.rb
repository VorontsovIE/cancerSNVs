require 'sequence_with_snp'
require 'mutation_context'

def load_mutation_types(mutations_markup_filename)
  File.open(mutations_markup_filename) do |f|
    header = f.readline.chomp.split("\t")
    snp_name_column = header.index('variant_id')
    mutation_type_column = header.index('mut_type')

    f.each_line.map do |line|
      data = line.chomp.split("\t")
      [data[snp_name_column], data[mutation_type_column].split(',')]
    end
  end
end

def load_mutation_contexts(snp_sequences_filename)
  snps_splitted = File.readlines(snp_sequences_filename).map{|el|
    name, seq = el.chomp.split("\t")
    snv = SequenceWithSNP.from_string(seq, name: name)
    [snv.name, snv.context(before: 1, after: 1, allele_variant_number: 0)]
  }
end

# requested_mutation_types should be downcased
def type_matches?(mutation_types, requested_mutation_types)
  return true  if requested_mutation_types.nil?
  mutation_types = mutation_types.map(&:downcase)
  requested_mutation_types.any?{|req_mut_type|
    mutation_types.include?(req_mut_type)
  }
end

def context_matches?(mutation_context, requested_mutation_contexts)
  requested_mutation_contexts.any?{|req_mut_context|
    req_mut_context.match?(mutation_context)
  }
end

## filter mutations by type
# requested_mutation_types should be downcased
def matching_type_mutation_names(mutations_markup_filename, requested_mutation_types)
  mutation_types = load_mutation_types(mutations_markup_filename)
  mutation_types.select{|variant_id, mutation_types|
    type_matches?(mutation_types, requested_mutation_types)
  }.map{|variant_id, mutation_types| normalized_snp_name(variant_id) }.to_set
end

def possible_mutation_types(mutations_markup_filename)
  mutation_types = load_mutation_types(mutations_markup_filename)
  mutation_types.flat_map{|variant_id, mutation_types|
    mutation_types
  }.uniq
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
